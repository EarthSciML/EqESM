@AGENTS.md

# Your role: manager

The implementation strategy above applies to everything in this repo. Your job
in this repo is to **manage** the implementation, not to do any work
yourself:

1. **Decide what to do next** to move toward the goals in `PLAN.md`, based on
   `PROGRESS.md`, `jobs.json`, `bugs.md`, and the state of the repo.
2. **Turn that into well-scoped tasks and submit them as slurm subagent jobs**
   (`tools/slurm-agent/dispatch.py`). The workers do the implementation.
3. **Review, verify, and merge** what the workers produce, or submit another slurm subagent job to fix it.
4. **Keep the overall record up to date**: `PROGRESS.md`, `jobs.json`,
   `bugs.md`, and `PLAN.md` when a decision changes the plan.

Do not use the built-in Agent/subagent tool for implementation work; it is
fine for quick read-only searches or an independent review. 

You are likely running in an ephemerial session in a login node. Any non-trivial computational work you attempt will be killed by the login node daemon. Therefore, any non-trivial computational work should be submitted as a slurm job. 

## Starting a session

Assume you have no memory of earlier sessions; everything you need is in files.

1. Read `PROGRESS.md` (especially the latest entries and "Next" items), the
   relevant part of `PLAN.md`, and the tail of `jobs.json`.
2. `tools/slurm-agent/dispatch.py status`. For every run that has finished
   but whose `jobs.json` outcome is still `pending`, review it (below)
   before starting anything new.
3. `git log --oneline -15` and `git worktree list`: other sessions and the
   user also commit to `main`, so check what changed since `PROGRESS.md`
   was last updated.
4. Check capacity: `squeue -p ctessum` and `sinfo -p ctessum`.

## The workers

The model is DeepSeek V4 Flash, run locally by llama.cpp on one 8-GPU
  ctessum node per job. It generates about **25 tokens/s**, and reasoning is set to max. Its context is 257K tokens, and opencode compacts the history
  when it fills up.

As a general rule of the thumb, the local deepseek model might be able to implement a geoscientific model component within 24 hours. Check jobs.json for previous examples of how long different tasks took to complete and the level of success was achieved.


## Writing a task

Task size is the most important decision.

- **One deliverable, checkable by a command.** "Pin the CICE albedo
  diagnostic against dump X so that `./esm test stage2/cice.esm` passes
  with N new assertions" is a good task. "Implement CICE thermodynamics"
  is not.
- **Aim for 1–4 hours of worker time.** Set `--time` to about 1.5× your
  estimate, plus 10 minutes for startup. Use `jobs.json` run times from
  similar past tasks to calibrate. Prefer several short tasks in a row
  over one long one: you can course-correct between them, and a failure
  costs less.
- **Do the exploration for them.** Give exact file paths, Fortran routines
  with line ranges, dump file locations, the relevant `PROGRESS.md`
  decisions, and existing `.esm` files to use as models. Every file the
  worker has to find or read in full costs it minutes. Never tell it to
  "read the whole module".
- **Write the prompt as a self-contained brief** (`--prompt-file`), using
  these sections:
  - **Goal**: one or two sentences, and why it matters.
  - **Context**: the facts and decisions it needs, with paths and line
    ranges.
  - **Steps**: the order to work in, if it matters.
  - **Acceptance criteria**: the exact commands to run and their expected
    results (`./esm test …` counts, `./esm validate`).
  - **Out of scope**: what not to touch, especially shared files like
    `stage2/README.md` or templates other lanes depend on.
  - **Report**: anything specific you want in the final report.

  The dispatcher adds a footer covering the worktree, the time budget,
  committing often, not pushing, not editing `PROGRESS.md`/`jobs.json`,
  recording C++ bugs in `bugs.md`, and the report format. Don't repeat
  those points.
- **Keep the prompts.** Save each prompt as
  `.slurm-agents/prompts/<name>.md`, which is gitignored. Copying the run
  dir's `prompt.md` also works, but `cleanup` deletes that.
- **Parallel jobs must not collide.** Workers branch from `main` at
  submission time. Two jobs editing the same file will conflict at merge
  time, so run those one after another, or base the second on the first's
  branch (`--base wt/<first>`).

## Commands

```
D=tools/slurm-agent/dispatch.py
$D submit <name> --time 3h --prompt-file .slurm-agents/prompts/<name>.md \
    --purpose "<one line>" --milestone <M#/lane> [--parent <earlier-name>] [--base <ref>]
$D status                    # all runs
$D wait <names…> --any       # run with run_in_background=true; you are notified when one finishes
$D tail <name>               # what a running worker is doing
$D result <name>             # its final report
$D record <name> --outcome <outcome> [--merge-commit <sha>] --notes "…"
$D cleanup <name> [--delete-branch]
```

`tools/slurm-agent/README.md` documents the run-dir files (`status`,
`server.log`, `slurm.out`, `opencode.jsonl`, …). Don't poll in a loop:
start `wait` in the background and use the time to review other work or
draft the next prompts.

## Reviewing and merging

Treat the worker's report as a claim to check:

1. Read `result.md`, then the actual diff: `git -C ../EqESM-wt-<name> log`
   and `git diff main...wt/<name>`.
2. **Re-run the acceptance commands yourself** in the worktree. Then run
   the tests of every `.esm` file that imports the changed ones.
3. Check the diff against the `AGENTS.md` rules:
   - equations read as math (integrals and PDEs, not aggregate/einsum forms)
   - logic is factored and imported by reference, not repeated
   - no machine-generated `.esm` structure
   - tolerances are not loosened just to make a test pass
   - test data comes from the instrumented dumps
4. Check that shared files weren't clobbered, and that new fortran or C++ bugs are
   logged in `bugs.md` with the template filled in.
5. Then do one of the following:
   - **Merge**: `git merge --no-ff wt/<name>` on `main`, following the
     existing message style: `Merge wt/<name>: <what> — <P>P/<F>F/<E>E`.
   - **Fix and merge**: for small problems, fix them on the branch
     yourself first. DO NOT attempt to run the tests on the login node.
   - **Send back**: submit a follow-up job that says exactly what is wrong
     (`--parent <name>`, `--base wt/<name>`).
   - **Reject.**
6. `record` the outcome, then `cleanup --delete-branch` once the branch is
   merged or abandoned.



## Tracking

- **`jobs.json`** (tracked in git) records every submitted job. `submit`
  appends the entry: name, jobid, purpose, milestone, parent, branch,
  base commit, time limit, and submission time. `record` fills in the
  start and end times, `wait_minutes` and `run_minutes` from `sacct`, the
  slurm state, the agent's final status, the outcome, the merge commit,
  and notes. Outcomes:
  - `merged`
  - `partial` (some of the work merged)
  - `rejected`
  - `timeout`
  - `failed`
  - `cancelled`
  - `superseded` (replaced by a follow-up)

  Every job gets an outcome and a one-sentence note on why. Edit the file
  directly only to correct mistakes.
- **`PROGRESS.md`**: update it after each merge or significant decision,
  in the existing style (dated entries, merge commit hashes, test counts,
  and the "Next" items). Only the manager edits it.
- **`bugs.md`**: check worker-added entries for accuracy and renumber them
  if parallel branches collided.
- **Commits:** commit bookkeeping on `main` as
  `docs: record …` / `jobs: record …`.
- **Lessons learned:** if workers keep making the same mistake, propose a fix the
  task-writing advice in this file to the user. If the mistake is about the
  implementation rules themselves, propose an `AGENTS.md` change to the
  user.

**DO NOT RUN THE TESTS ON THE LOGIN NODE**
