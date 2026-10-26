# PROGRESS — eqesm implementation

Tracking progress on PLAN.md (Stage 1 — instrument the CESM codes; Stage 2 —
stubs → physics → EarthSciModels PRs; Stage 3 — top-level `eqesm.esm`).
Last updated: 2026-10-25.

## Milestone 0 — Scaffolding (in progress)

- **Plan** ✅ — `PLAN.md` written (2026-10-25): CESM inventory, stage 1–3
  playbooks, milestones M0–M10. Target: CESM `cesm3.0-alphabranch` (submodule
  tags pinned at M1), first component **MOSART**, top-level model
  **`eqesm.esm`** (AGENTS.md naming fixed).
- **Tracking files** ✅ — `PROGRESS.md` (this file), `bugs.md` (template),
  `.gitignore` (esm binary, stage-1 artifacts).
- **Rust CLI** ⏳ — copy of the EarthSciAST CLI binary (`esm`, untracked) in
  the repo root, built from the sibling `../EarthSciAST` checkout
  (`cargo build --release -p earthsci-ast`) or copied from an existing
  binary (e.g. `../esast-perf-bins/`). Verify `./esm validate` / `test` /
  `units --check` on a trivial document.
- **Workflow check** ✅ — manager-worker arrangement with worktree
  subagents established (AGENTS.md); worktrees from `main`, up to 3 workers.

## Milestone 1 — CESM source access + build; reference runs; boundary catalog

Not started. See PLAN.md §4.1 for the playbook: pin `cesm3.0-alphabranch` +
submodule tags into `stage1/cesm` (gitignored); toolchain probe (local
netCDF/PIO build vs apptainer); CIME machine config; smoke runs (F-class,
X-class w/ MOSART, I-class, C-class, SCAM / Icepack / clm1pt harnesses);
minimal inputdata; `stage1/boundaries.md`.

## Milestone 2 — Instrumentation scaffolding + MOSART end to end

Not started. See PLAN.md §4.2–4.3 and §10.
