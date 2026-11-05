#!/usr/bin/env python3
"""Dispatch headless opencode subagents as slurm jobs.

Each subagent is a slurm job (built from agent.slurm.tmpl) that gets a whole
GPU node, starts its own llama-server, runs `opencode run` on a task prompt
inside a dedicated git worktree, and writes the agent's final message to
result.md. State for each run lives in <repo>/.slurm-agents/<name>/.

Typical manager loop:
    dispatch.py submit cice-foo --time 6:00:00 --prompt-file task.md
    dispatch.py status                 # all runs
    dispatch.py wait cice-foo          # blocks until the job leaves the queue
    dispatch.py result cice-foo        # agent's final report
    ... review + merge wt/cice-foo ...
    dispatch.py cleanup cice-foo       # remove worktree (+ branch with --delete-branch)
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parent
TEMPLATE = TOOL_DIR / "agent.slurm.tmpl"

# Slurm resource profiles. The llama.cpp build in
# /projects/illinois/eng/cee/ctessum/ctessum/llm is compiled for CUDA arch 75
# (Quadro RTX 6000) only, and launch-server.sh's v4 layout assumes 8 x 24 GB
# GPUs. Add a profile here for another partition once llama.cpp is rebuilt
# for its GPUs. A comma-separated --profile list submits to several partitions
# at once (slurm starts the job on whichever has room first); this only works
# when the profiles share the same gres/mem/ntasks.
#
# qos: on this cluster preemption is QOS-based (PreemptType=preempt/qos).
# `scav_high` (the ctessum account's default) preempts the `normal` QOS that
# every scavenger-partition job runs under, so an agent job outranks scavenger
# work no matter what its nice value is. It is set explicitly so that this does
# not depend on the account default.
PROFILES = {
    "ctessum": dict(partition="ctessum", account="ctessum", qos="scav_high",
                    gres="gpu:QuadroRTX6000:8", ntasks=40, mem=176992,
                    max_seconds=3 * 86400),
}
DEFAULT_PROFILE = "ctessum"

# Default slurm nice value for agent jobs. Among ctessum-partition jobs the
# priority is PriorityWeightPartition (2000) + age (up to 1000 after 7 days,
# about 6 per hour) - nice; the QOS priority weight is 0. A nice of 100 puts an
# agent job just below a group member's normal job of the same age (it takes
# about 17 h of extra waiting to catch up with a freshly submitted one). It
# does not affect preemption of scavenger jobs (see qos above). Users may only
# raise nice (lower priority); a negative value needs a slurm operator.
DEFAULT_NICE = 100

LLM_MODELS = {"v4": "local/deepseek-v4-flash", "v41": "local/deepseek-v4.1-flash"}

PROMPT_FOOTER = """

---
## Execution context (added by the dispatcher)

You are running unattended as a subagent in a slurm batch job; nobody will
answer questions. Make reasonable decisions and record them.

- Working directory: `{workdir}`{branch_line}
- Time budget: about {time_human} of wall time, after which you will be stopped.
  Commit work-in-progress regularly so nothing is lost if time runs out.
- Do not push to any remote, open pull requests, or modify the main checkout
  at `{repo}`.
- Do not edit `PROGRESS.md` or `jobs.json`; the manager maintains them from
  your report. If you find a bug or error in the CESM Fortran/C++ code, add an
  entry to `bugs.md` (following its template) and mention it in your report.
- If the task turns out to be impossible as specified, or much larger than it
  looks, stop early and explain why in your report rather than improvising a
  different task.
- When you are finished, end with a final message that is a concise report for
  the manager: what you did, the commits you made, test results (pass/fail
  counts), anything left unfinished, and any problems or open questions.
"""


def sh(*args, check=True, cwd=None):
    return subprocess.run(args, check=check, cwd=cwd, text=True,
                          capture_output=True).stdout.strip()


def repo_root():
    """Main checkout (not a worktree) of the repo containing the cwd."""
    common = Path(sh("git", "rev-parse", "--path-format=absolute",
                     "--git-common-dir"))
    return common.parent


def runs_dir():
    d = repo_root() / ".slurm-agents"
    d.mkdir(exist_ok=True)
    return d


def run_dir(name):
    d = runs_dir() / name
    if not d.is_dir():
        sys.exit(f"no run named {name!r} in {runs_dir()}")
    return d


def load_meta(d):
    return json.loads((d / "meta.json").read_text())


def parse_time(s):
    """Slurm-style time ('90' minutes, 'H:MM:SS', 'D-HH:MM:SS', or '2h'/'30m')."""
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([hm])", s)
    if m:
        return int(float(m[1]) * (3600 if m[2] == "h" else 60))
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        days = int(d)
    parts = [int(p) for p in s.split(":")]
    if len(parts) == 1:
        secs = parts[0] * 60
    elif len(parts) == 2:
        secs = parts[0] * 3600 + parts[1] * 60 if days else parts[0] * 60 + parts[1]
    else:
        secs = parts[0] * 3600 + parts[1] * 60 + parts[2]
    return days * 86400 + secs


def fmt_time(secs):
    d, r = divmod(secs, 86400)
    h, r = divmod(r, 3600)
    m, s = divmod(r, 60)
    return (f"{d}-" if d else "") + f"{h:02d}:{m:02d}:{s:02d}"


def slurm_state(jobid):
    out = sh("squeue", "-h", "-j", str(jobid), "-o", "%T|%M|%N|%R", check=False)
    if out:
        st, elapsed, node, reason = out.split("|")
        return st, elapsed, node if st == "RUNNING" else reason
    out = sh("sacct", "-n", "-X", "-P", "-j", str(jobid),
             "-o", "State,Elapsed,NodeList", check=False)
    if out:
        st, elapsed, node = out.splitlines()[0].split("|")
        return st, elapsed, node
    return "UNKNOWN", "", ""


# ---------------------------------------------------------------- commands

def cmd_submit(a):
    name = a.name
    if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        sys.exit("name may only contain letters, digits, '.', '_' and '-'")
    repo = repo_root()
    d = runs_dir() / name
    if d.exists():
        sys.exit(f"{d} already exists; pick another name or `cleanup {name}` first")

    if a.prompt_file:
        prompt = Path(a.prompt_file).read_text()
    elif a.prompt:
        prompt = a.prompt
    elif not sys.stdin.isatty():
        prompt = sys.stdin.read()
    else:
        sys.exit("give the task with --prompt-file, --prompt, or on stdin")
    if not prompt.strip():
        sys.exit("empty prompt")

    profiles = [PROFILES[p] for p in a.profile.split(",")]
    prof = dict(profiles[0])
    prof["partition"] = ",".join(p["partition"] for p in profiles)
    for k in ("partition", "account", "qos", "gres", "ntasks", "mem"):
        if getattr(a, k) is not None:
            prof[k] = getattr(a, k)

    secs = parse_time(a.time)
    if secs > prof["max_seconds"]:
        sys.exit(f"--time {a.time} exceeds the partition limit {fmt_time(prof['max_seconds'])}")
    if a.nice < 0:
        sys.exit("--nice must be >= 0 (a negative nice needs a slurm operator)")
    if secs < a.server_timeout + a.grace + 600:
        sys.exit(f"--time {a.time} is too short: model loading alone can take "
                 f"{a.server_timeout // 60} min")

    # Worktree
    branch = None
    if a.workdir:
        workdir = Path(a.workdir).resolve()
    else:
        branch = a.branch or f"wt/{name}"
        workdir = repo.parent / f"{repo.name}-wt-{name}"
        if workdir.exists():
            sys.exit(f"{workdir} already exists")
        if not a.dry_run:
            sh("git", "worktree", "add", "-b", branch, str(workdir), a.base, cwd=repo)

    agent_secs = secs - a.grace - 300  # rough, for the prompt only
    footer = PROMPT_FOOTER.format(
        workdir=workdir, repo=repo, time_human=fmt_time(max(agent_secs, 0)),
        branch_line=(f"\n- Git branch: `{branch}` (a worktree of `{repo}`); "
                     f"commit your work to this branch." if branch else ""))

    d.mkdir(parents=True)
    (d / "prompt.md").write_text(prompt.rstrip() + "\n" + ("" if a.no_footer else footer))

    subs = dict(
        PARTITION=prof["partition"], ACCOUNT=prof["account"], QOS=prof["qos"],
        NICE=a.nice, YIELD=0 if a.no_yield else 1, YIELD_AFTER=int(a.yield_after * 60),
        GRES=prof["gres"],
        NTASKS=prof["ntasks"], MEM=prof["mem"], TIME=fmt_time(secs),
        TIME_SECONDS=secs, GRACE_SECONDS=a.grace, SERVER_TIMEOUT=a.server_timeout,
        JOB_NAME=f"agent-{name}", RUN_DIR=d, WORKDIR=workdir,
        LLM_VARIANT=a.llm, OC_MODEL=a.model or LLM_MODELS[a.llm], TOOL_DIR=TOOL_DIR,
    )
    script = TEMPLATE.read_text()
    for k, v in subs.items():
        script = script.replace("{{" + k + "}}", str(v))
    left = re.findall(r"\{\{\w+\}\}", script)
    if left:
        sys.exit(f"unfilled template placeholders: {left}")
    (d / "job.slurm").write_text(script)

    meta = dict(name=name, workdir=str(workdir), branch=branch, base=a.base,
                time=fmt_time(secs), profile=a.profile, llm=a.llm, nice=a.nice, yield_on=not a.no_yield,
                submitted=time.strftime("%Y-%m-%d %H:%M:%S"), jobid=None)
    if a.dry_run:
        print(f"dry run: job script at {d / 'job.slurm'} (no worktree created, nothing submitted)")
        (d / "meta.json").write_text(json.dumps(meta, indent=2))
        return
    (d / "status").write_text("queued\n")
    out = sh("sbatch", "--parsable", str(d / "job.slurm"))
    meta["jobid"] = int(out.split(";")[0])
    (d / "meta.json").write_text(json.dumps(meta, indent=2))
    jobs = load_jobs()
    jobs.append(dict(
        name=name, jobid=meta["jobid"], purpose=a.purpose, milestone=a.milestone,
        parent=a.parent, branch=branch, workdir=str(workdir),
        base_commit=sh("git", "rev-parse", "--short", a.base, cwd=repo) if branch else None,
        time_limit=fmt_time(secs), nice=a.nice, llm=a.llm, prompt_chars=len(prompt),
        submitted=meta["submitted"], started=None, ended=None,
        wait_minutes=None, run_minutes=None, slurm_state=None, agent_status=None,
        outcome="pending", merge_commit=None, notes=""))
    save_jobs(jobs)
    print(f"submitted {name}: job {meta['jobid']}  (recorded in {jobs_path()})")
    print(f"  worktree: {workdir}" + (f"  (branch {branch})" if branch else ""))
    print(f"  run dir:  {d}")


def cmd_status(a):
    names = a.names or sorted(p.name for p in runs_dir().iterdir() if (p / "meta.json").exists())
    rows = []
    for n in names:
        d = run_dir(n)
        m = load_meta(d)
        st, el, where = slurm_state(m["jobid"]) if m.get("jobid") else ("DRYRUN", "", "")
        status = (d / "status").read_text().strip() if (d / "status").exists() else ""
        rows.append((n, str(m.get("jobid") or "-"), st, el, where, status))
    if not rows:
        print("no runs")
        return
    hdr = ("NAME", "JOBID", "SLURM", "ELAPSED", "NODE/REASON", "AGENT STATUS")
    w = [max(len(r[i]) for r in rows + [hdr]) for i in range(5)]
    for r in [hdr] + rows:
        print("  ".join(r[i].ljust(w[i]) for i in range(5)) + "  " + r[5])


def is_finished(d):
    m = load_meta(d)
    st, _, _ = slurm_state(m["jobid"])
    return st not in ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED", "SUSPENDED")


def cmd_wait(a):
    """Block until all named runs leave the queue (or, with --any, the first one)."""
    dirs = {n: run_dir(n) for n in a.names}
    pending = set(dirs)
    while pending:
        for n in sorted(pending):
            if is_finished(dirs[n]):
                status = (dirs[n] / "status").read_text().strip() if (dirs[n] / "status").exists() else "?"
                print(f"{n}: {status}", flush=True)
                pending.discard(n)
                if a.any:
                    return
        if pending:
            time.sleep(a.interval)


def cmd_result(a):
    d = run_dir(a.name)
    r = d / "result.md"
    if r.exists() and r.stat().st_size:
        print(r.read_text())
    else:
        print(extract_text(d) or f"(no result yet; status: {(d / 'status').read_text().strip()})")


def cmd_tail(a):
    """Human-readable view of the agent's event stream so far."""
    d = run_dir(a.name)
    lines = []
    for ev in iter_events(d / "opencode.jsonl"):
        for kind, text in describe_event(ev):
            lines.append(f"[{kind}] {text}")
    print("\n".join(lines[-a.n:]) if lines else "(no agent output yet)")


def cmd_cancel(a):
    for n in a.names:
        m = load_meta(run_dir(n))
        sh("scancel", str(m["jobid"]))
        print(f"cancelled {n} (job {m['jobid']})")


def _squeue_rows(partition):
    fields = ("id", "name", "state", "reason", "submit", "start", "nodes", "nice", "left", "nodelist")
    out = sh("squeue", "-h", "-p", partition, "-o", "%i|%j|%T|%r|%V|%S|%D|%y|%L|%N", check=False)
    rows = []
    for line in out.splitlines():
        r = dict(zip(fields, line.split("|")))
        try:
            r["nice"], r["nodes"] = int(r["nice"]), int(r["nodes"])
        except (KeyError, ValueError):
            continue
        rows.append(r)
    return rows


def _seconds_since(ts):
    try:
        return time.time() - time.mktime(time.strptime(ts, "%Y-%m-%dT%H:%M:%S"))
    except ValueError:
        return 0


def _seconds_left(s):
    try:
        return parse_time(s)
    except ValueError:
        return float("inf")  # UNLIMITED / NOT_SET


def _hostnames(spec):
    if not spec or spec == "(null)":
        return None
    return set(sh("scontrol", "show", "hostnames", spec, check=False).split())


def _node_limits(jobid):
    """(required nodes or None, excluded nodes) of a pending job."""
    out = sh("scontrol", "show", "job", "-o", str(jobid), check=False)
    f = dict(kv.split("=", 1) for kv in out.split() if "=" in kv)
    return _hostnames(f.get("ReqNodeList")), _hostnames(f.get("ExcNodeList")) or set()


def _agent_yield_settings(runs):
    """{jobid: yield_after_seconds or None (= never yields)} from each run's job script."""
    out = {}
    for m in runs.glob("*/meta.json"):
        try:
            jobid = str(json.loads(m.read_text()).get("jobid"))
            script = (m.parent / "job.slurm").read_text()
        except (OSError, ValueError):
            continue
        y = re.search(r"^YIELD=(\d)", script, re.M)
        after = re.search(r"^YIELD_AFTER=(\d+)", script, re.M)
        out[jobid] = int(after[1]) if y and y[1] == "1" and after else None
    return out


def cmd_should_yield(a):
    """Called once a minute by a running agent job. Prints the id of the waiting
    job this agent should make room for, or nothing.

    A job "outranks" an agent when its nice value is lower, whoever submitted
    it. Waiting jobs (PENDING on Resources/Priority), lowest nice and oldest
    first, are each given as many agent nodes as they request: running agent
    jobs (job name agent-*) with a higher nice, on a node the waiting job may
    use (its ReqNodeList/ExcNodeList), whose own yield-after time the job has
    waited (read from each run's job.slurm, so all watchers agree; agents
    submitted with --no-yield never qualify), highest nice and most recently
    started first, skipping agents with less than --min-left seconds left.
    Every watcher sees the same squeue snapshot, so they agree on who yields. A shared log in the runs dir records each
    yield, so a waiting job that still cannot start (e.g. it needs more than
    one node can give) never gets more yields than the nodes it asked for.
    """
    import fcntl
    rows = _squeue_rows(a.partition)
    me = next((r for r in rows if r["id"] == str(a.jobid)), None)
    if not me or me["state"] != "RUNNING":
        return
    runs = Path(a.runs_dir) if a.runs_dir else runs_dir()
    settings = _agent_yield_settings(runs)
    settings.setdefault(me["id"], a.after)
    finite = [v for v in settings.values() if v is not None]
    shortest = min(finite) if finite else a.after
    waiting = sorted((r for r in rows if r["state"] == "PENDING"
                      and r["reason"] in ("Resources", "Priority")
                      and _seconds_since(r["submit"]) >= shortest),
                     key=lambda r: (r["nice"], r["submit"]))
    if not waiting:
        return
    logp = runs / "yield-log.json"
    with open(runs / "yield-log.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        log = json.loads(logp.read_text()) if logp.exists() else {}
        pending_ids = {r["id"] for r in rows if r["state"] == "PENDING"}
        log = {w: v for w, v in log.items() if w in pending_ids}  # forget started/gone jobs
        already = {j for v in log.values() for j in v}
        agents = [r for r in rows if r["state"] == "RUNNING"
                  and r["name"].startswith("agent-") and r["id"] not in already
                  and _seconds_left(r["left"]) >= a.min_left]
        # highest nice first; within a nice, the most recently started first
        agents.sort(key=lambda r: r["start"], reverse=True)
        agents.sort(key=lambda r: r["nice"], reverse=True)
        taken = set()
        target = None
        for w in waiting:
            need = w["nodes"] - len(log.get(w["id"], []))
            if need <= 0:
                continue
            waited = _seconds_since(w["submit"])
            req, exc = _node_limits(w["id"])
            for ag in agents:
                if need <= 0:
                    break
                after = settings.get(ag["id"])
                if (ag["id"] in taken or ag["nice"] <= w["nice"] or after is None
                        or waited < after or ag["nodelist"] in exc
                        or (req is not None and ag["nodelist"] not in req)):
                    continue
                taken.add(ag["id"])
                need -= 1
                if ag["id"] == me["id"]:
                    target = w["id"]
            if target:
                break
        if target and not a.dry_run:
            log.setdefault(target, []).append(me["id"])
        if not a.dry_run:
            logp.write_text(json.dumps(log, indent=1) + "\n")
    if target:
        print(target)


def cmd_renice(a):
    """Change the nice value of queued/running runs (scontrol update)."""
    if a.nice < 0:
        sys.exit("nice must be >= 0 (a negative nice needs a slurm operator)")
    for n in a.names:
        d = run_dir(n)
        m = load_meta(d)
        r = subprocess.run(["scontrol", "update", f"jobid={m['jobid']}", f"nice={a.nice}"],
                           text=True, capture_output=True)
        if r.returncode:
            print(f"{n} (job {m['jobid']}): {(r.stderr or r.stdout).strip()}", file=sys.stderr)
            continue
        m["nice"] = a.nice
        (d / "meta.json").write_text(json.dumps(m, indent=2))
        jobs = load_jobs()
        for j in jobs:
            if j.get("jobid") == m["jobid"]:
                j["nice"] = a.nice
        save_jobs(jobs)
        print(f"{n} (job {m['jobid']}): nice={a.nice}")


OUTCOMES = ("merged", "partial", "rejected", "timeout", "failed", "cancelled", "superseded")


def jobs_path():
    return repo_root() / "jobs.json"


def load_jobs():
    p = jobs_path()
    return json.loads(p.read_text()) if p.exists() else []


def save_jobs(jobs):
    jobs_path().write_text(json.dumps(jobs, indent=2) + "\n")


def find_job(jobs, name):
    """Most recent jobs.json entry with this name."""
    for j in reversed(jobs):
        if j["name"] == name:
            return j
    sys.exit(f"no jobs.json entry named {name!r}")


def cmd_record(a):
    """Fill a jobs.json entry's timing from sacct and set its outcome/notes."""
    jobs = load_jobs()
    j = find_job(jobs, a.name)
    out = sh("sacct", "-n", "-X", "-P", "-j", str(j["jobid"]),
             "-o", "Submit,Start,End,State", check=False)
    if out:
        sub_t, start, end, state = out.splitlines()[0].split("|")
        parse = lambda s: time.mktime(time.strptime(s, "%Y-%m-%dT%H:%M:%S"))
        ok = lambda s: re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d", s)
        j["slurm_state"] = state.split()[0]
        if ok(start):
            j["started"] = start.replace("T", " ")
            j["wait_minutes"] = round((parse(start) - parse(sub_t)) / 60, 1)
            if ok(end):
                j["ended"] = end.replace("T", " ")
                j["run_minutes"] = round((parse(end) - parse(start)) / 60, 1)
    d = runs_dir() / a.name
    if (d / "status").exists():
        j["agent_status"] = (d / "status").read_text().strip()
    if (d / "yields.log").exists():
        # requeued to make room for other jobs; started/run_minutes cover the last attempt only
        j["yields"] = len((d / "yields.log").read_text().splitlines())
    if a.outcome:
        j["outcome"] = a.outcome
    if a.merge_commit:
        j["merge_commit"] = a.merge_commit
    if a.notes:
        j["notes"] = (j["notes"] + " " if j.get("notes") and a.append else "") + a.notes
    save_jobs(jobs)
    print(json.dumps(j, indent=2))


def cmd_cleanup(a):
    d = run_dir(a.name)
    m = load_meta(d)
    if m.get("jobid") and not is_finished(d):
        sys.exit(f"{a.name} is still in the queue; cancel it first")
    j = next((j for j in reversed(load_jobs()) if j["name"] == a.name), None)
    if j and j["outcome"] == "pending" and not a.force:
        sys.exit(f"jobs.json outcome for {a.name} is still 'pending'; "
                 f"run `record {a.name} --outcome ...` first (or --force)")
    repo = repo_root()
    if m.get("branch") and Path(m["workdir"]).exists():
        args = ["git", "worktree", "remove", m["workdir"]] + (["--force"] if a.force else [])
        r = subprocess.run(args, cwd=repo, text=True, capture_output=True)
        if r.returncode:
            sys.exit(f"git worktree remove failed (uncommitted changes? use --force):\n{r.stderr}")
        print(f"removed worktree {m['workdir']}")
    if a.delete_branch and m.get("branch"):
        r = subprocess.run(["git", "branch", "-D" if a.force else "-d", m["branch"]],
                           cwd=repo, text=True, capture_output=True)
        print(r.stdout.strip() or r.stderr.strip())
    if not a.keep_run_dir:
        archive = archive_run_logs(d, m, a.name)
        if archive:
            print(f"archived logs to {archive}")
        shutil.rmtree(d)
        print(f"removed {d}")


ARCHIVED_FILES = ("opencode.jsonl", "result.md", "prompt.md", "meta.json",
                  "slurm.out", "opencode.err", "status", "exit_code")


def archive_run_logs(d, m, name):
    """Copy a run's transcript and report into .work/agent-logs/<name>-<jobid>/
    (gitignored) before its run dir is deleted, so later log reviews can
    cover every job. opencode.jsonl is gzipped; the llama-server log is not
    kept (large, and only useful while debugging a live server)."""
    import gzip
    dest = repo_root() / ".work" / "agent-logs" / f"{name}-{m.get('jobid', 'nojob')}"
    try:
        dest.mkdir(parents=True, exist_ok=True)
        for fn in ARCHIVED_FILES:
            src = Path(d) / fn
            if not src.exists():
                continue
            if fn == "opencode.jsonl":
                with open(src, "rb") as fi, gzip.open(dest / (fn + ".gz"), "wb") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                shutil.copy2(src, dest / fn)
    except OSError as e:
        print(f"warning: could not archive logs for {name}: {e}", file=sys.stderr)
        return None
    return dest


# ------------------------------------------------- opencode event parsing

def iter_events(path):
    if not path.exists():
        return
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if line.startswith("{"):
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                pass


def event_part(ev):
    """The message part carried by one `opencode run --format json` event.

    Events look like {"type": "text"|"tool_use"|"step_start"|..., "part": {...}};
    text parts have part.type == "text", tool calls part.type == "tool".
    """
    p = ev.get("part")
    return p if isinstance(p, dict) else {}


def describe_event(ev):
    p = event_part(ev)
    if p.get("type") == "text" and p.get("text", "").strip():
        yield "text", p["text"].strip()
    elif p.get("type") == "tool":
        state = p.get("state") or {}
        inp = state.get("input") or {}
        title = inp.get("command") or inp.get("path") or inp.get("filePath") or state.get("title") or json.dumps(inp)
        yield f"tool:{p.get('tool')}", f"{state.get('status', '')} {str(title)[:200]}".strip()


def extract_text(d):
    """The agent's final assistant text: the text parts after its last tool call."""
    blocks = {}
    for ev in iter_events(d / "opencode.jsonl"):
        p = event_part(ev)
        if p.get("type") == "text" and p.get("text", "").strip():
            blocks[p.get("id", len(blocks))] = p["text"].strip()  # same id = streamed update
        elif p.get("type") == "tool":
            blocks = {}
    return "\n\n".join(blocks.values())


def cmd_extract(a):
    print(extract_text(Path(a.run_dir)))


# ------------------------------------------------------------------ CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("submit", help="create a worktree and submit a subagent job")
    s.add_argument("name", help="run name; also the worktree/branch suffix (wt/<name>)")
    s.add_argument("--time", "-t", required=True,
                   help="wall time: minutes, H:MM:SS, D-HH:MM:SS, or e.g. 6h / 90m")
    g = s.add_mutually_exclusive_group()
    g.add_argument("--prompt-file", "-f")
    g.add_argument("--prompt", "-p")
    s.add_argument("--base", default="HEAD", help="commit/branch the worktree starts from")
    s.add_argument("--branch", help="branch name (default wt/<name>)")
    s.add_argument("--workdir", help="run in this existing directory instead of a new worktree")
    s.add_argument("--profile", default=DEFAULT_PROFILE,
                   help=f"resource profile(s), comma-separated; known: {', '.join(PROFILES)}")
    s.add_argument("--partition"), s.add_argument("--account"), s.add_argument("--gres")
    s.add_argument("--qos", help="slurm QOS (default from the profile; scav_high for ctessum)")
    s.add_argument("--nice", type=int, default=DEFAULT_NICE,
                   help=f"slurm nice value, >= 0; higher = lower priority among ctessum jobs "
                        f"(default {DEFAULT_NICE}; 0 = same as a normal group job)")
    s.add_argument("--yield-after", type=float, default=10, metavar="MIN",
                   help="requeue to make room for a job with a lower nice that has waited "
                        "this many minutes for a node (default 10)")
    s.add_argument("--no-yield", action="store_true", help="never give up the node")
    s.add_argument("--ntasks", type=int), s.add_argument("--mem", type=int)
    s.add_argument("--llm", default="v4", choices=list(LLM_MODELS),
                   help="launch-server.sh model variant")
    s.add_argument("--model", help="opencode model id (default matches --llm)")
    s.add_argument("--server-timeout", type=int, default=1200,
                   help="seconds to wait for llama-server to become healthy")
    s.add_argument("--grace", type=int, default=300,
                   help="seconds reserved at the end of the job for cleanup")
    s.add_argument("--no-footer", action="store_true", help="don't append the execution-context footer")
    s.add_argument("--dry-run", action="store_true", help="write the job script but don't submit")
    s.add_argument("--purpose", required=True, help="one-line purpose, recorded in jobs.json")
    s.add_argument("--milestone", help="PLAN.md milestone / lane, recorded in jobs.json")
    s.add_argument("--parent", help="name of an earlier job this one retries or continues")
    s.set_defaults(func=cmd_submit)

    s = sub.add_parser("record", help="update a jobs.json entry (timing from sacct, outcome, notes)")
    s.add_argument("name")
    s.add_argument("--outcome", choices=OUTCOMES)
    s.add_argument("--merge-commit")
    s.add_argument("--notes")
    s.add_argument("--append", action="store_true", help="append to existing notes instead of replacing")
    s.set_defaults(func=cmd_record)

    s = sub.add_parser("status", help="show runs and their slurm/agent state")
    s.add_argument("names", nargs="*")
    s.set_defaults(func=cmd_status)

    s = sub.add_parser("wait", help="block until runs finish")
    s.add_argument("names", nargs="+")
    s.add_argument("--any", action="store_true", help="return when the first one finishes")
    s.add_argument("--interval", type=int, default=60)
    s.set_defaults(func=cmd_wait)

    s = sub.add_parser("result", help="print the agent's final report")
    s.add_argument("name")
    s.set_defaults(func=cmd_result)

    s = sub.add_parser("tail", help="readable summary of the agent's recent events")
    s.add_argument("name")
    s.add_argument("-n", type=int, default=40)
    s.set_defaults(func=cmd_tail)

    s = sub.add_parser("should-yield", help="(used by running jobs) print the waiting job to yield to")
    s.add_argument("--jobid", required=True)
    s.add_argument("--partition", required=True)
    s.add_argument("--after", type=int, default=600, help="seconds a job must have waited")
    s.add_argument("--min-left", type=int, default=1800,
                   help="don't yield with less than this many seconds of time limit left")
    s.add_argument("--runs-dir")
    s.add_argument("--dry-run", action="store_true", help="decide only; don't record the yield")
    s.set_defaults(func=cmd_should_yield)

    s = sub.add_parser("renice", help="change the slurm nice value of queued/running runs")
    s.add_argument("nice", type=int)
    s.add_argument("names", nargs="+")
    s.set_defaults(func=cmd_renice)

    s = sub.add_parser("cancel", help="scancel runs")
    s.add_argument("names", nargs="+")
    s.set_defaults(func=cmd_cancel)

    s = sub.add_parser("cleanup", help="remove a finished run's worktree and run dir")
    s.add_argument("name")
    s.add_argument("--delete-branch", action="store_true")
    s.add_argument("--force", action="store_true", help="discard uncommitted changes / unmerged branch")
    s.add_argument("--keep-run-dir", action="store_true")
    s.set_defaults(func=cmd_cleanup)

    s = sub.add_parser("extract", help=argparse.SUPPRESS)  # used by the job script
    s.add_argument("run_dir")
    s.set_defaults(func=cmd_extract)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
