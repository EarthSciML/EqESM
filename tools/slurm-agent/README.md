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
dispatch.py cancel <names…>
dispatch.py cleanup <name> [--delete-branch] [--force] [--keep-run-dir]
```

From Claude Code, run `dispatch.py wait <name>` with `run_in_background` to
get a notification when a subagent finishes. Then review the branch, merge
it, and run `cleanup`.

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
