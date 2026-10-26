# PROGRESS — eqesm implementation

Tracking progress on PLAN.md (Stage 1 — instrument the CESM codes; Stage 2 —
stubs → physics → EarthSciModels PRs; Stage 3 — top-level `eqesm.esm`).
Last updated: 2026-10-25.

## Milestone 0 — Scaffolding (in progress)

- **Plan** ✅ — `PLAN.md` written (2026-10-25): CESM inventory, stage 1–3
  playbooks, milestones M0–M10. Target: CESM `cesm3.0-alphabranch` (submodule
  tags pinned at M1), first component **MOSART**, top-level model
  **`eqesm.esm`** (AGENTS.md naming fixed).
- **Format stance** ✅ — updated to **esm 1.3.0** (2026-10-25): all
  documents declare `"esm": "1.3.0"`. The local `../EarthSciAST` checkout
  was fast-forwarded to **v0.4.0** (released 2026-10-25; spec 1.3.0, schema
  `esm/1.3.0`). 1.3.0 adds `solver.external` (§2.2.6 — run-config
  reference, the EarthSciIntegrator hook) and `solver.time_span` (§2.2.7).
- **Tracking files** ✅ — `PROGRESS.md` (this file), `bugs.md` (template),
  `.gitignore` (esm binary, stage-1 artifacts).
- **Rust CLI** ⏳ — copy of the EarthSciAST CLI binary (`esm`, untracked) in
  the repo root, built from `../EarthSciAST` at v0.4.0
  (`cargo build --release -p earthsci-ast`). Verify `./esm validate` /
  `test` / `units --check` on a trivial document.
- **Workflow check** ✅ — manager-worker arrangement with worktree
  subagents established (AGENTS.md); worktrees from `main`, up to 3 workers.

## Milestone 1 — CESM source access + build; reference runs; boundary catalog

In progress. See PLAN.md §4.1 for the playbook.

- **Source access + pins** ✅ (2026-10-25) — CESM `cesm3.0-alphabranch` @
  `e59f251` checked out at `stage1/cesm` (gitignored, 4.1 GB) with
  `bin/git-fleximod update`; all submodule tags pinned and recorded in
  `stage1/runs/README.md` (cam6_4_212, ctsm5.4.063, cesm3_cice6_6_3_22,
  mi_260909, cismwrap_3_00_002, mosart1.1.13, cmeps1.1.73, cdeps1.0.112,
  cime6.7.1, …). Finding: CAM 6.4 radiation is **RTE-RRTMGP** (not RRTMG);
  photolysis is tuv-x v0.14.3.
- **Rust CLI** ✅ (2026-10-25) — `esm` 0.4.0 binary built from the
  fast-forwarded `../EarthSciAST` checkout, copied to the repo root
  (untracked); validates an esm 1.3.0 document (`./esm validate` clean).
- **Toolchain probe** ✅ (2026-10-25) — merged (`ee29816`):
  `stage1/runs/TOOLCHAIN.md`. Environment: SLURM, account/partition
  `ctessum` (5-node, 3-day limit), gcc 13.3.0 + openmpi 5.0.1 (sibling
  projects' pair), cmake 3.31.8/perl/python all present, apptainer 1.5.4
  (unprivileged). **No netCDF anywhere** — hand-built. `data.cesm.ucar.edu`
  does not resolve; no svn client (use `--protocol=https`); NCAR SVN +
  LCRC E3SM mirror + GitHub reachable.
- **netCDF stack** ✅ (2026-10-25) — hdf5 1.14.6 (serial) + netcdf-c 4.9.3 +
  netcdf-fortran 4.6.1 built into `stage1/build/local` via tracked
  `stage1/build/build-netcdf.sh` (gitignored prefix); Fortran write +
  ncdump round-trip verified. Note: `support.hdfgroup.org` 404s the classic
  release path — use the GitHub release asset (script updated).
  `stage1/runs/env.sh` provides the toolchain env (`NETCDF_PATH`).
- **Compset survey** ✅ (2026-10-25/26) — merged (`33243cb`):
  `stage1/runs/COMPSET_SURVEY.md` (368 lines). Smoke-run table: F-class
  `F2000climo`/`f19_f19_mg17` (CAM6+MOSART), CAM7 `FHISTC_LTso`/`ne3pg3`,
  I-class `I2000Clm60Sp` (CLM+MOSART coupled), D-class `D` (CICE),
  C-class `C_JRA`+`mom/tx10deg` (MOM6), B-class `B1850C_LTso`/`ne30pg3_t232`,
  SCAM `FSCAMARM97`, coupler testbed `X`, single-point `I1PtClm60SpRsGs`.
  Key findings: compsets split per-component; I-class → `D` compsets;
  no `ne4_ne4`/`T62_g17` grids (smallest SE = `ne3_ne3_mg37`); B-class is
  MOM6-only; machine config now CMake-based (`cmake_macros/*.cmake`,
  `$HOME/.cime/config_machines.xml` supported, `DIN_LOC_ROOT` required);
  inputdata via OSDF endpoint first. **Bug:** `R2000MOSART` broken as
  shipped (`DLND%LCPL` invalid) → `bugs.md` #1; PLAN.md grid note updated.
  Corrections (2026-10-26, `4c26cc8`, verified in `cam/bld/configure`):
  **SCAM requires the SE dycore** → FSCAM grid is `ne3_ne3_mg37`, not
  `f19_f19`; **radiation**: cam6 → RRTMG, cam7 → RTE-RRTMGP (my earlier
  "CAM 6.4 = RRTMGP" note was wrong — the F-class CAM6 reference runs use
  classic RRTMG). atmos_phys `0_31_000` pin re-verified via
  `git submodule status` (the vendored ChangeLog top tag `0_07_000` is just
  its first entry).
- **Boundary catalog** ✅ (2026-10-26) — merged (`6677181`):
  `stage1/boundaries.md` (398 lines): §0 CMEPS mediator phases/runseq/
  coupling intervals + esmFldsExchange field tables (~370 fields); §1
  per-component import/export caps (CAM, CTSM, CICE, MOM, MOSART, CISM,
  CDEPS) with D/I classification; §2 MOSART internals in detail
  (`mosart_run` orchestration, `Euler` hillslope→subnetwork→channel
  cascade, pointwise laws `CRVRMAN_nosqrt`/`CREHT_nosqrt`/`GR*` = first
  test targets, ID/dnID network handling; **no lakes/reservoirs/
  groundwater module in mosart1.1.13** — only `qgwl`; PLAN.md §3 row
  corrected); §3 grids/precision (real64 everywhere; CICE single-precision
  TBD); §4 instrumentation hooks. Caveat: cmeps is UFS-style
  (field-by-field API; coupling times from cime buildnml, no MED_NCPL
  symbol).
- **CIME machine config** ✅ (2026-10-26) — merged (`6ecaa1e`):
  `stage1/cime/machine/{config_machines.xml,config_batch.xml,README.md}` —
  machine `illinois` (gnu/openmpi, slurm, `ctessum` partition/account,
  QOS `scav_high`, 40 cores/node, NETCDF_PATH→stage1/build/local,
  ESMF_RUNTIME_PROFILE=ON, NODENAME_REGEX `^ccc02.*`); validated via
  CIME's Machines/Batch classes; `stage1/inputdata/` mirror dir created.
  TBDs: /scratch.local is node-local (CIME_OUTPUT_ROOT moved to
  `stage1/cime/cases` on /projects), login-node hostnames unverified,
  OMPI-5 mpirun slot discovery untested, first build compiles ESMF from
  source (long).
  `stage1/boundaries.md` (398 lines). §0 CMEPS mediator (run phases,
  runseq ordering, coupling intervals, esmFldsExchange tables ~370 fields);
  §1 per-component import/export caps (CAM/CTSM/CICE/MOM/MOSART/CISM/CDEPS)
  with D/I classification; §2 MOSART internals in detail (`mosart_run`
  orchestration incl. irrigation/flood/direct transfer/subcycling; `Euler`
  hillslope→subnetwork→channel cascade with CFL substeps; pointwise laws
  `CRVRMAN_nosqrt`/`CREHT_nosqrt`/`GR*` = preferred first test targets;
  ID/dnID network topology; **finding: no lakes/reservoirs/groundwater
  module in mosart1.1.13** — PLAN.md §3 row corrected); §3 grids/precision
  (real64 everywhere; CICE single-precision TBD); §4 instrumentation
  hooks. Caveat: this cmeps is UFS-style (field-by-field API; coupling
  times from cime buildnml, no MED_NCPL symbol).
- **CIME machine config** ✅ (2026-10-26) — merged (`6ecaa1e`):
  `stage1/cime/machine/{config_machines.xml,config_batch.xml,README.md}`
  — machine `illinois` (gnu/openmpi, slurm `ctessum` partition, account
  `ctessum`, QOS `scav_high`, 40 cores/node, NETCDF_PATH→stage1/build/local,
  ESMF_RUNTIME_PROFILE=ON, NODENAME_REGEX `^ccc02.*`), validated via
  CIME's own Machines/Batch classes; `stage1/inputdata/` mirror dir created.
  TBDs: /scratch.local is node-local (CIME_OUTPUT_ROOT may need /projects),
  login-node hostnames unverified, OMPI5 mpirun slot discovery in sbatch
  untested, first build compiles ESMF from source (long).
- **First case created** ✅ (2026-10-26) — `F2000climo_f19`
  (F2000climo/f19_f19_mg17) created + set up on machine `illinois`
  (36 tasks, 1 node, STOP nhours/6, wallclock 2 h). Operational note:
  `--extra-machines-dir` only applies at create time — `case.setup` etc.
  need the config in `$HOME/.cime/` (installed; README updated).
- **ESMF 8.9.1 built** ✅ (2026-10-26) — CMEPS needs ESMF
  (`find_package(ESMF REQUIRED)`; cmeps buildnml reads `$ESMFMKFILE`); no
  ESMF module on the cluster and not bundled in the checkout. Built from
  source via `stage1/build/build-esmf.sh` into `stage1/build/esmf`
  (gcc 13.3 + openmpi 5.0.1, netcdf split, no pnetcdf — must be *unset*,
  not "OFF": `ifdef ESMF_PNETCDF` in common.mk). Machine config sets
  `ESMFMKFILE` (note: real path is `lib/libO/.../esmf.mk`) + `ESMF_DIR`.
- **LAPACK** ✅ (2026-10-26) — link failed on LAPACK driver routines
  (`dgesv_` etc. from CLM/CICE); cluster `/sw/apps/lapack/3.12.1/lib64`
  has `liblapack.a` + `librefblas.a` but no `libblas` → symlinked
  `libblas.a`/`liblapack.a` into `stage1/build/local/lib` and set
  `LAPACK_LIBDIR` in the machine config.
- **First build** ✅ (2026-10-26) — `case.build` finished successfully
  (33 s incremental; total with ESMF ~1 h). `bld/cesm.exe` (50 MB ELF)
  built: CAM6 + CLM5 + CICE6 + DOCN + MOSART + CMEPS + CDEPS on gcc
  13.3/openmpi 5.0.1.
- **Inputdata** ✅ (2026-10-26) — all F2000climo/f19 files downloaded via
  wget (`check_input_data --download`): ESMF meshes (fv1.9x2.5, gx1v7),
  MOSART r05 routing, CLM surfdata/initdata/ndep, CAM physprops/ozone/
  solar/sst, etc.
- **First smoke run — SUCCESS** ✅ (2026-10-26) — `F2000climo_f19`
  (F2000climo/f19_f19_mg17, 36 tasks, 6 sim. hours) ran in 8 min wall
  on ccc0232: **"SUCCESSFUL TERMINATION OF CESM"**, exit 0. All components
  wrote restarts + history (cam.r 520 MB, clm2.r 198 MB, cpl.r, cice.r,
  **mosart.r + rh0a/rh0i**, docn.r). Run 1 failed at exec time
  (`libnetcdff.so.7` not found — CIME env vars *overwrite*
  LD_LIBRARY_PATH, so the machine config now lists every lib dir:
  netCDF/ESMF/gcc/openmpi; commit `1f952ac`). ESMF_Profile.summary shows
  expected coupling (12 LND run phases = 6 h @ ATM_NCPL 48; ctsm_run
  dominated by restart I/O).
- **I-class run — SUCCESS** ✅ (2026-10-26) — `I2000Clm60Sp_f19`
  (DATM%GSWP3v1 + CLM6 + MOSART coupled, 1 day, 36 tasks): GSWP3 forcing
  (~150 GB, 215 files) downloaded; **MOSART ran 8 routing steps**
  (ROF_NCPL=8), restarts + rh0a/rh0i archived. st_archive initially
  failed (`ncdump not found` — fixed by adding PATH to machine config).
- **D-class run — SUCCESS** ✅ (2026-10-26) — `D_f19g17` (active CICE,
  data atm/ocn/rof, 1 day, 18 tasks): CICE advanced 24 steps, restart at
  0001-01-02. First attempt failed at `ice_pio_init`: cold-start `ice_ic`
  file not in input_data_list (CIME quirk — see TOOLCHAIN.md §2.6 #11);
  fetched manually from svn-ccsm (240 MB).
- **SCAM run — SUCCESS** ✅ (2026-10-26) — `FSCAMARM97_ne3`
  (CAM60 SCAMARM97 IOP, ne3_ne3_mg37, **1 PE required** — first attempt
  failed with "single column mode must be run on 1 pe" after NTASKS was
  overridden to 18; fixed with NTASKS=1). 6 sim. hours, clean finalize.
- **X-class run** 🏃 (2026-10-26) — `X_f19g17` (all-stub coupler testbed)
  built + submitted (11827277).
- **C-class run** ⏳ (2026-10-26) — `C_JRA_tx10` (MOM6 active + tx10deg
  testmod, T62_t232) built successfully (mom 116 s); input data download
  (large MOM6 ICs/forcing) still in progress.
- **Next:** finish C-class + X-class; then M2 (instrumentation scaffolding
  + MOSART boundaries).

## Milestone 2 — Instrumentation scaffolding + MOSART end to end

Not started. See PLAN.md §4.2–4.3 and §10.
