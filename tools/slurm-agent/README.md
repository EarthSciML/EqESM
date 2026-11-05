# slurm-agent: opencode subagents as slurm jobs

A manager agent (Claude Code, or a person) can hand a task to a local
DeepSeek/opencode agent running on its own GPU node:

```
tools/slurm-agent/dispatch.py submit <name> --time 6h --prompt-file task.md --purpose "one line"
```

This:

1. creates a git worktree `../EqESM-wt-<name>` on a new branch `wt/<name>`
   (from `--base`, default `HEAD`). Use `--workdir DIR` to skip this.
2. writes `.slurm-agents/<name>/` (gitignored) with `prompt.md` (the task
   plus a footer about the worktree, time budget and the report format; turn
   the footer off with `--no-footer`) and `job.slurm` (filled in from
   `agent.slurm.tmpl`).
3. runs `sbatch` on it. The job takes a whole ctessum GPU node, starts
   `launch-server.sh v4` (model loading takes about 4 min), then runs
   `opencode run --standalone --auto` in the worktree with the prompt. The
   agent is stopped `--grace` seconds (default 300) before the slurm time
   limit.

## Priority

Agent jobs are submitted with `--nice 100` by default, so they rank just
below a normal job from someone in the ctessum group. Within the `ctessum`
partition, priority is 2000 (partition) + age (up to 1000 after 7 days, about
6 per hour) − nice. A job with nice 100 catches up with a group member's
freshly submitted job only after about 17 hours of extra waiting.
`--nice 0` gives normal group priority; a larger value lowers it further.
`dispatch.py renice <nice> <names…>` changes it for a queued job. Plain users
can only use values ≥ 0, and slurm may refuse to lower a value that has
already been raised.

### Yielding the node

Nice only orders waiting jobs, so a running agent job would otherwise keep its
node until it finishes. Agent jobs therefore give the node up themselves.
Once a minute a watcher in the job runs `dispatch.py should-yield`. It
checks the partition for a job with a **lower nice value** (any user, agent
or not) that has been pending on `Resources`/`Priority` for `--yield-after`
minutes (default 10). If this agent is the one chosen to make room, it:

1. stops opencode with SIGINT;
2. stops llama-server and copies the session DB to `opencode-data/`;
3. runs `scontrol requeue` on itself, so the status file reads `requeued`.

Agents with the highest nice yield first, and within the same nice the most
recently started. Agents with under 30 min left don't yield. A waiting job
gets at most as many yields as the nodes it requested; these are recorded in
`.slurm-agents/yield-log.json`. A job that still can't start therefore
doesn't cause repeated requeues.

When slurm restarts the job (`SLURM_RESTART_COUNT` > 0), the new attempt:

- restores the DB;
- continues the same opencode session with a "you were requeued, check git
  status" message;
- appends to `opencode.jsonl` and `slurm.out`.

The worktree is on shared disk, so its files are as the agent left them.
Each yield is a line in `yields.log` and is counted as `yields` in
`jobs.json`. `started`/`run_minutes` in `jobs.json` cover the last attempt
only. A restart costs the model load (~4 min) plus re-reading the session
context. Use `submit --no-yield` for a job that must not be interrupted.

Nice does not affect scavenger jobs. Preemption on this cluster is QOS-based,
and agent jobs run under QOS `scav_high` (set explicitly in the profile).
That QOS preempts the `normal` QOS every `scavenger`-partition job uses,
whatever its nice value.

## Files in the run dir

| file | contents |
|---|---|
| `status` | `queued` → `starting` → `loading-model` → `running` → `done` / `timeout …` / `failed: …` / `cancelled` |
| `result.md` | the agent's final message (its report to the manager) |
| `opencode.jsonl` | the full event stream (`--format json`) |
| `opencode.err`, `server.log`, `slurm.out` | opencode stderr, llama-server log, job log |
| `opencode-data/` | the job's opencode DB and logs, copied back when the job ends |
| `meta.json` | job id, worktree, branch, and so on |

## Manager commands

```
dispatch.py status [names…]          # table of runs: slurm state + agent status
dispatch.py wait <names…> [--any]    # block until they leave the queue
dispatch.py tail <name> [-n 40]      # readable view of recent tool calls/text
dispatch.py result <name>            # final report
dispatch.py record <name> --outcome merged|partial|rejected|timeout|failed|cancelled|superseded \
    [--merge-commit SHA] [--notes …]   # fill in jobs.json timing (from sacct) and outcome
dispatch.py renice <nice> <names…>    # change a queued run's slurm nice value
dispatch.py cancel <names…>
dispatch.py cleanup <name> [--delete-branch] [--force] [--keep-run-dir]
```

From Claude Code, run `dispatch.py wait <name>` with `run_in_background` to
get a notification when a subagent finishes. Then review the branch, merge
it, and run `cleanup`.

Before `cleanup` deletes a run dir, it archives the transcript
(`opencode.jsonl`, gzipped), `result.md`, `prompt.md`, `meta.json`,
`slurm.out`, `opencode.err`, `status`, and `exit_code` into the gitignored
`.work/agent-logs/<name>-<jobid>/`. The llama-server log is not kept. Read an
archived transcript with `zcat … | jq`, or use it for cross-job log reviews.

## jobs.json

`submit` appends an entry to `<repo>/jobs.json` (tracked in git). The entry
holds the name, jobid, purpose, milestone, parent, branch, base commit, time
limit and submission time. `record` adds `started`/`ended`,
`wait_minutes`/`run_minutes` from `sacct`, the slurm state, the agent's final
status, the outcome, the merge commit and notes. `cleanup` refuses to remove
a run whose outcome is still `pending`.

## Why opencode can run on several nodes at once here

By default opencode keeps its sqlite DB, background-service file and locks
in `~/.local/share/opencode` and `~/.local/state/opencode`, on the
NFS-shared home directory. Instances on different nodes fight over these
files. Each job therefore sets `XDG_DATA_HOME`, `XDG_STATE_HOME` and
`XDG_CACHE_HOME` to `/scratch.local/$USER/slurm-agents/$SLURM_JOB_ID`
(node-local), and uses `--standalone` so no shared background service is
involved. Config is still read from `~/.config/opencode`, so the provider,
model and permission settings are shared. To look at a finished session:

```
XDG_DATA_HOME=$(mktemp -d -p /scratch.local/$USER) && cp -r .slurm-agents/<name>/opencode-data $XDG_DATA_HOME/opencode
XDG_DATA_HOME=$XDG_DATA_HOME opencode session list
```

## Limits

- Only the `ctessum` profile works today. The llama.cpp build is compiled
  for CUDA arch 75 (Quadro RTX 6000) only, and the v4 layout assumes
  8 × 24 GB GPUs. To add another partition, rebuild llama.cpp for its
  GPUs and add an entry to `PROFILES` in `dispatch.py`. `--profile a,b`
  submits to several partitions at once, which requires matching gres.
  `--partition/--account/--gres/--mem/--ntasks` override individual fields.
- llama-server always uses port 8000, because `~/.config/opencode/opencode.json`
  points the `local` provider there. That is fine as long as each job has a
  whole node, and the job fails fast if the port is already taken.
- On ctessum nodes, scavenger jobs from other users are preempted when a
  job is submitted, so the job may wait in `PENDING (Resources)` for a few
  minutes.
