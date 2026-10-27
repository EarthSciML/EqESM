# PROGRESS — eqesm implementation

Tracking progress on PLAN.md (Stage 1 — instrument the CESM codes; Stage 2 —
stubs → physics → EarthSciModels PRs; Stage 3 — top-level `eqesm.esm`).
Last updated: 2026-10-27.

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
  (large MOM6 ICs/forcing) still in progress. First background download
  died silently (343/619 JRA55 files) — restarted via nohup retry loop
  (pid 3782524, log `inputdata-download.log` in case dir).
- **X-class run — SUCCESS** ✅ (2026-10-26) — `X_f19g17` (all-stub coupler
  testbed, 1 day) completed clean (44 s).
- **M1 → M2 transition** (2026-10-26): all M1 deliverables merged and
  pushed (origin/main == main). M2 kicked off: CESM fork worktree
  `stage1/EqESM-cesm-mosart` (branch `instrument/mosart-dumps`) created
  from `cesm3.0-alphabranch` @ `e59f251`; mosart submodule at
  `mosart1.1.13` (detached — fork branch to be created in submodule).
  Instrumentation targets confirmed in source: `mosart_run`
  (mosart_driver.F90:375), `Euler` (mosart_physics.F90:51), pointwise
  laws `CRVRMAN_nosqrt` (F90:406), `CREHT_nosqrt` (F90:426), `GR*`
  family, hillslope/subnetwork routing, `Routing_KW` (F90:290).
- **MOSART stage-2 skeleton** ✅ (2026-10-26) — worktree
  `wt/stage2-mosart` (`../EqESM-wt-stage2-mosart`, branch `wt/stage2-mosart`),
  commit `2907da3`: `stage2/` with six esm 1.3.0 files
  (`mosart_templates.esm` pointwise laws CRVRMAN_nosqrt/CREHT_nosqrt/GR*,
  `mosart_hillslope.esm`, `mosart_subnetwork.esm`, `mosart_channel.esm`
  (Routing_KW), `mosart_coupling.esm` (RofCap unit conversions),
  `mosart_routing.esm` (assembly: inline Control = mosart_run steps 2-8 +
  set_subtimesteps, 29 variable_map couplings, enums, solver block)) +
  `README.md` (file map, units, §2.6 mapping, per-law transcription notes,
  divergences). All `esm validate` ✓, `esm test` clean (0 assertions; tests
  blocks empty pending instrumentation dumps). Key decisions: routing steps
  as storage ODEs with `apply_expression_template` laws; Manning n declared
  `s/m^(1/3)` (unit checker); eroutUp/direct SMM gathers are placeholders
  until the frivinp data source (stage 2.5).
- **Reuse survey** ✅ (2026-10-26) — merged (`7596df4`):
  `stage1/EARTHSCIMODELS_SURVEY.md` (335 lines) + PLAN.md §5.1 +
  boundaries.md §4.1. Headline: **no river-routing physics exists anywhere**
  (surface_runoff is 1-D Saint-Venant, related-but-distinct); reuse the
  EqIce stage-2 pattern (templates + continuous ODE + measured scheme
  gap), `faq`/ragged gathers for ID/dnID topology, `data_sources
  kind: "static"` for frivinp. New code uses `faq` (`aggregate`
  deprecated). ESD has no network grid / no CFL-substep template.
- **Stage-2 skeleton** ✅ (2026-10-26) — merged (`85b0ad4`):
  `stage2/{mosart_templates,hillslope,subnetwork,channel,coupling,
  routing}.esm` (2806 lines) + README. All 6 files validate under esm
  1.3.0. Laws transcribed exactly (manning_velocity=CRVRMAN_nosqrt,
  overland_discharge=CREHT_nosqrt, GR* geometry incl. 1/sin(atan(0.1))
  flood-plain constant); routing steps as storage ODEs with
  `apply_expression_template`; 29 variable_map entries; 8 documented
  divergences (rate-rule irrigation/flood, continuous-ODE limit,
  eroutUp/direct gathers pending network data source). Tests blocks
  empty pending dumps.
- **C-class runs** (2026-10-26): `C_JRA_tx10` (tx10deg testmod) — built,
  but **5 files missing on ALL inputdata servers** → `bugs.md` #2
  (tx10deg/tx3deg testmods reference unpublished files, dated 2026-02-19;
  no generator in checkout; `check_input_data` fails). Workaround:
  `C_JRA_t62` case at full tx2_3v2 (no testmod) built (408 s), inputdata
  complete. First two submits FAILED at launch: `mpirun -np 72` (2 nodes)
  dies ~15 s in — OpenMPI prterun ssh-launches cross-node daemons, blocked
  by `pam_slurm_adopt` (single-node 36-PE runs never hit this). Fixed:
  machine config now uses `srun -n {{ total_tasks }}` (committed `957c01c`,
  TOOLCHAIN.md quirk #11); `case.setup --reset` invalidated the build, so
  rebuilt (22 s) and resubmitted with srun (11839272, PENDING — queue busy).
- **MOSART instrumentation** ✅ (2026-10-26) — fork branch
  `instrument/mosart-dumps` complete: mosart submodule @ `f9a19d1`
  (mosart_dump_mod.F90, 380 lines, namelist-gated `dump_mosart=.false.`,
  PET-0 only, 16 cells × 12 steps) + superproject @ `7a9db4a` (submodule
  pointer + `stage1/runs/mosart_instrument.json`). Run
  `I2000Clm60Sp_mosart_inst` (job 11837196, exit 0:0, 8:22 wall):
  **bit-identical to reference** I-class run (rh0a/rh0i) → solution-neutral;
  `decomp_option='basin'` ≡ roundrobin numerically; budget closes
  (net ≈ −lag). 9 CSVs in `stage1/dumps/mosart/`: creht_nosqrt (216 rows),
  crvrman_trib (4320), crvrman_chan (5640), gr_subnet (4352), gr_chan
  (5672), mosart_euler (768), mosart_run_imports/exports (256 each),
  mosart_budget (2). No C++ bugs found in MOSART (the stale
  `11835528 DependencyNeverSatisfied` queue entry is a CIME artifact).
- **Tests filling — DONE** ✅ (2026-10-26, merged `9317d34`) — stage-2
  tests blocks filled from dumps: **20 tests / 115 assertions, all green**
  (`esm test stage2/mosart_*.esm`: 0 fail / 0 err). Pointwise law pins at
  rel 1e-12 (CREHT_nosqrt, CRVRMAN_nosqrt, GR* geometry chains incl. the
  TINYVALUE zero guards), one-step traces in continuous-ODE form with
  measured scheme gaps (~2–5× assertion bounds; assembly volr/wh agree
  to 2.4e-16/1.0e-6 rel while wt/wr split differs 0.59×/3.3× — documented,
  not asserted), cap conversions, and one full coupling-period assembly
  test. Findings: fthresh is the spval sentinel 1e36 in the Fortran (the
  .esm 0.0 placeholder would fire the flood rule on the first period —
  assembly test overrides); numDT_r=numDT_t=4 reproduced from the CFL
  equations at the defaults; channel stays dry in all 8 dumped steps →
  GRHR/GRPR wet branches untested (coverage note in README). **This closes
  the M2 loop: stub → physics → tests green under `./esm test`.**
  **M2 remaining:** EarthSciModels PR prep (esm 1.3.0 CI acceptance check
  per survey risk) + CICE instrumentation (M3).
- **M3 kickoff** (2026-10-26) — three workers dispatched in parallel:
  1. **MOSART PR prep** (wt `mosart-routing` in EarthSciModels): copy
     stage-2 files into the corpus (placement TBD: `river_routing/` vs
     `environmental_transport/`), conform metadata (tags/authors/
     references), and verify under ALL THREE CI runners (Python
     earthsci_ast[simulate]@main, Julia walker, Rust esm built from
     EarthSciAST main) — the survey's version-drift risk (corpus 1.0.0
     vs our 1.3.0) is the gate. No PR opened until manager review.
  2. **CICE instrumentation** (CESM fork worktree `instrument/cice-dumps`):
     mirror the MOSART dump-module pattern (namelist-gated, PET-0, CSV
     headers); instrument Icepack column laws (albedo, freezing temp,
     strength P, EVP stress, drag), CICE_Run I-boundaries, transport,
     budget; validate with a `D`-class instrumented run (single node,
     `srun --mpi=pmix` launcher + `env -u SLURM_MEM_PER_*` submit);
     dumps → `stage1/dumps/cice/`; solution-neutrality vs D_f19g17.
  3. **CICE stage-2 skeleton** (wt/stage2-cice): mirror the MOSART
     skeleton — `cice_templates.esm` (laws transcribed exactly from
     icepack/cicedyn source), `cice_thermodynamics.esm` (vertical heat
     conduction as PDE, continuous form), `cice_dynamics.esm` (EVP),
     `cice_transport.esm`, `cice_coupling.esm` (Si_/Faii_/Fioi_ fields
     per boundaries.md §1.3), `cice.esm` assembly, README; tests blocks
     empty pending dumps.
- **MOSART PR opened** ✅ (2026-10-26) — `EarthSciML/EarthSciModels#33`
  (branch `mosart-routing` @ `1292b25`, fork ctessum-claude), 8 files
  +3570/−2, placed at `components/environmental_transport/mosart/`
  (same family pattern as surface_runoff; no registry updates needed).
  All three CI runners verified PASS locally BEFORE opening (Rust esm
  rebuilt from EarthSciAST main v0.5.0: 115P/0F; Python
  earthsci_ast[simulate]@main via `tools/run_esm_inline_tests.py`:
  116/116 incl. load-only gate; Julia walker `run_esm_tests`: 115P/0F).
  **esm 1.3.0 acceptance risk RESOLVED — no document changes needed for
  any runner** (defaults differ: Rust 1e-4/1e-6, Py 1e-10/1e-14, Julia
  1e-10/1e-12 — tighter integrators only improve trajectories; assertion
  bands are 2–1000× measured scheme gaps). Metadata conformed
  (tags/created/authors/references per surface_runoff style); equations
  + all 115 numeric values byte-identical to stage-2 sources.
  **Upstream finding:** EarthSciModels `Project.toml` compat
  `EarthSciAST = "0.3"` is stale — the Julia package at EarthSciAST main
  declares 0.5.0, so upstream CI's Julia resolve is broken at HEAD for
  every PR; PR #33 carries the documented 2-line bump to "0.5".
  PR history: #32 (pre-rebase) closed — GitHub refuses reopen after a
  branch recreation, and the shell wrapper blocks force-push, so the
  rebase (onto main @ 6048aac, Project.toml conflict resolved) landed
  via API-delete + fresh push + new PR.
  **CI (@1292b25 + @6047fa9):** consumer-e2e PASS; rust-cli, python,
  julia×4 FAIL — but **all three runner failures are pre-existing
  upstream breakages at HEAD, none from our files**, which pass on every
  runner that reaches them (rust 115P/0F/0E; python all 6 files OK;
  julia never reaches the sweep):
  - rust: `gaschem/stratospheric/*` diffsol "Step size is too small" —
    identical on main's own runs (39485217759).
  - python: `atmospheric_dynamics/wsm6/cold_accretion.esm`
    unlowered_operator — identical on main's own runs.
  - julia: **upstream moving-HEAD breakage** — CI resolves EarthSciAST
    from `rev=main`, which moved past v0.5.0 (SciML v3 merge #571):
    EarthSciAST itself fails to PRECOMPILE on Julia 1.13.3
    (SchemaValidationError) before any file loads. Reproduced locally:
    f6f8a4bfd (v0.5.0) precompiles fine on 1.13.3; 408dabcb3 (current
    main) fails identically. Affects every EarthSciModels PR at HEAD.
    Our files verified under the Julia walker (115P/0F) at v0.5.0.
  - One failure WAS ours (fixed `6047fa9`): my rebase-conflict rewrite
    of `Project.toml` dropped `[extras]`/`[targets]` → `Pkg.test()`
    "Package Test not found"; restored verbatim, delta vs main is only
    the documented `EarthSciAST = "0.5"` compat bump (still correct:
    with 0.3 the Julia resolve fails outright with "empty intersection").
  Full diagnosis posted on the PR. **Decision (2026-10-27): upstream CI
  red is out of scope — EarthSciModels' runners (gaschem, wsm6, julia
  precompile) are broken at HEAD for every PR; not ours to fix. PR #33
  remains open + mergeable; our files pass on every runner that reaches
  them. No further upstream CI work.**
- **CICE instrumentation validated** ✅ (2026-10-27) — worker
  (ses_e854ab19fffeDY1bFNMNGV7Ef3): D_cice_inst run (jobs 11857263/
  11857562) COMPLETED 0:0; 13 CSVs in `stage1/dumps/cice/`; laws verified
  to round-off (freezing temp 384 rows; EVP stress 184,320 rows; EVP
  momentum chain 40,320 rows); solution neutrality vs D_f19g17 (iceh
  bit-identical, `stage1/analysis/compare_iceh.py`); bugs.md #4–#8 added;
  `stage1/runs/cice_instrument.json` written. **Manager review found a
  dump bug**: `trcr`/`trcrn` reads in `dump_thermo_state`/`dump_albedo`
  had the j and tracer indices swapped (declared (i,j,tracer,ncat,blk)),
  so Tsfc/Tsfcn/qice/qsno/alvl/apnd/hpnd/ipnd read 0.0 in every row —
  the tests worker had misread this as "Tsfc not dumped". Fixed
  (fork `ced43eb` → wrapper `90dde50` → superproject `c8a3bef`), rebuilt,
  rerun (job 11864874, COMPLETED 0:0, solution-neutral): columns now
  physical (Tsfc ≈ −14.7 °C, qice_1 ≈ −3.03e8 J/m3). bugs.md #9.
- **CICE coupling tests** ✅ (2026-10-27, merged `1fba628`) — worker
  (ses_e821e9d0fffeds2K1facI3YHhN) filled `cice_coupling.esm` tests:
  30 tests/33 assertions (`time_span {0,0}`, rel 1e-12, every value
  dump-cited), pinning imports (frzmlt/frain/fsnow/flw/sw bands via
  med2mod_areacor inversion), the Tf = max(S/(−18.48+0.01848S), −2.0)
  law (all 384 rows), and the export aggregation (Si_ifrac/vice/vsno/
  snowh, albedo Si_avsdr etc. incl. the scale_fluxes `/aice` division,
  Fioi_swpen). Three transcription fixes found while writing tests:
  scale_fluxes /aice division (ice_flux.F90:1565-1592), Faii_lwup ifelse
  arity (ice_import_export.F90:1149-1164), Tf law → the coupled
  shr_frz_freezetemp form (bugs.md #6). Manager added the Si_t export
  pin (`c8db84a`) after the dump fix; README-cice.md Testing section
  rewritten. `esm test stage2/`: 148P/0F/0E (115 MOSART + 33 CICE).
- **CICE stage-2 skeleton** ✅ (2026-10-27, merged `7dec7d6`) — worker
  `wt/stage2-cice` (ses_e854ab193ffeUZQdIktqWsJlLU): `cice_templates.esm`
  (41 pointwise law templates with file:line citations), `cice_thermodynamics.esm`
  (67 eqs: prep_radiation + step_therm1/2, heat conduction as semi-discrete
  flux-divergence ODEs), `cice_dynamics.esm` (34 eqs: EVP momentum + stress
  relaxation + Hibler strength), `cice_transport.esm` (3 continuity eqs as
  D(c)/t = −∇·(c·u) rewrite targets), `cice_coupling.esm` (62 eqs: ice_import
  renames, areacor, merge_fluxes exports), `cice.esm` assembly (59 param_to_var
  entries, solver external cice.run.json, dt_dyn subcycling + constant-atmo
  Control), `README-cice.md`. All 6 files validate clean; 0 tests pending
  dumps. Reviewed: Hibler strength + strairx claims verified against
  Fortran (icepack_mechred.F90:1715-1725, icepack_atmo.F90:455-470; the
  `/vmag` normalization bug confirmed at F90:328-329 — logged in
  bugs.md). Two bugs logged: (1) `strx = tau·uatm` missing `/wind`
  normalization (icepack_atmo.F90:463-464); (2) calculate_Tin_from_qin
  negative-discriminant abort (guard note). Manager fix during review:
  worker had REPLACED the shared `stage2/README.md` (MOSART docs, which
  mosart .esm descriptions reference) with CICE docs — restored MOSART
  README byte-identical (`38dc737`) and moved CICE docs to
  `stage2/README-cice.md`, CICE descriptions repointed. Merged, worktree
  removed, branch deleted; `esm test stage2/`: MOSART 115P/0F/0E still
  green, CICE mounted (8 mounts). **M3 CICE loop status:** coupling
  tests merged (33 assertions, above); thermo/dynamics/transport tests
  BLOCKED on a **skeleton configuration mismatch** (worker finding,
  dump-verified): the reference D-class run uses ktherm=2 mushy (MU71
  conductivity), shortwave=dEdd (alvdfn 0.9527 vs the transcribed ccsm3
  0.8598), atmbndy=similarity (Si_u10/tref/qref = 1.104/266.36/1.95e-3
  vs constant-scheme zeros), kstrength=1 Rothrock (P=3.21e5 vs Hibler
  5.47e5) — but the skeleton transcribed the icepack defaults
  (BL99/ccsm3/constant atmo/Hibler). No dump value pins the transcribed
  thermo/dynamics equations. **Next: re-transcribe the CICE skeleton
  thermo/dynamics/templates to the run configuration (mushy/dEdd/
  similarity/Rothrock), then fill those tests from the dumps.** (Also
  noted: Tsfc implicit surface-energy-balance solve stays a free
  parameter until stage 3; dynamics/transport spatial-D rewrite targets
  need EarthSciDiscretizations lowering — `unlowered_operator` compile
  block even for `{0,0}` tests.)**
- **C-class run — SUCCESS** ✅ (2026-10-26) — `C_JRA_t62` (full tx2_3v2
  ocean, 72 PEs/2 nodes, 1 day) COMPLETED (job 11842203, 7:53 wall;
  `case.run success`): MOM6 23 steps, CFL 0.075, mass err 6e-19. This
  closes the M1 smoke-run matrix (F/I/D/SCAM/X/C all green). Three
  launcher/runtime walls conquered (TOOLCHAIN.md #11): `mpirun` ssh-daemons
  blocked by pam_slurm_adopt → `srun` PMI1/2-vs-PMIx mismatch →
  `srun --mpi=pmix -n N` (machine config committed `cc59f31`); submit with
  `env -u SLURM_MEM_PER_* ./case.submit`; `AUTO_MASKTABLE=False` in
  user_nl_mom → bugs.md #3 (MOM6 auto-mask FATAL at 36 PEs; FATAL only in
  `run/ocn.log.*`, and the batch job hangs after the abort until cancelled).
  Run notes in COMPSET_SURVEY.md §C-class.
- **Next:** review/merge tests worker; then M3 kickoff (MOSART PR prep +
  CICE instrumentation).

## Milestone 2 — Instrumentation scaffolding + MOSART end to end

In progress (tests from dumps). See PLAN.md §4.2–4.3 and §10.

- **CICE instrumentation (M3) — run complete, dumps validated** ✅
  (2026-10-27) — fork worktree `stage1/EqESM-cesm-cice` (branch
  `instrument/cice-dumps`): `ice_dump_mod.F90` in
  `components/cice/src/cicecore/shared/` (leaf module, auto-picked-up by
  the CICE Makefile glob) + 9 wired files + cice-wrapper XML/buildnml.
  All Fortran syntax-checked against the reference D_f19g17 build's .mod
  files and compiled clean in the real `case.build` (cice built in 138 s).
  Final commit chain: cice/src `2dbb70f` → cice wrapper `bcb9496` →
  superproject `8cb32cd` (all local/unpushed, branch `instrument/cice-dumps`).
  - 13 CSVs in `stage1/dumps/cice/` (gitignored), all `cice_`-prefixed:
    run_imports, run_exports, freezing_temp, thermo_state, thermo_flux,
    flux, dyn_state, transport_state, ice_strength, evp_stress,
    evp_momentum, albedo, budget. EqIce headers
    (`# description` / `# columns`), PET-0-only, first
    `dump_first_cells`×`dump_first_steps` (32×12), OMP-critical writes,
    `dump_nml` gated (default off), `dump_init`/`dump_finalize` wired
    into cice_init2/ModelFinalize.
  - Case `D_cice_inst` (D f19_g17, 18 tasks × 2 thr, 1 day, dt=3600 s →
    24 steps; mirrors D_f19g17). Final run job **11857562** (first
    validated run 11857263): COMPLETED 0:0, model execution success.
    Walls hit and fixed: (1) fork worktree submodules not initialized /
    wrong commits — reset to pinned via `git submodule update --checkout`
    (plus cdeps nested fox/genf90 + a stray staged-deletion state in
    cdeps restored); (2) **nmlgen validates but drops unknown-to-the-
    writer groups** — cice buildnml writes a hard-coded `groups` list,
    so dump_nml never reached ice_in; added `'dump_nml'` (wrapper
    `16ff7fb`); first submit (11857214) cancelled 10 min in, regenerated
    via preview_namelists; (3) **blank line after every data row in
    evp_stress and freezing_temp** — final value write lacked
    `advance='no'` (fixed in `2dbb70f`, rerun 11857562; all 13 CSVs now
    zero blank lines).
  - **Validation** (all green): solution neutrality — iceh + restart
    bit-identical vs D_f19g17 modulo title/timestamp attrs
    (`stage1/analysis/compare_iceh.py`); freezing-temp law
    Tf = max(S/(−18.48+0.01848·S), −2.0) (shr_frz mushy) on all 384 rows,
    identical to run_imports sss/Tf; EVP stress relaxation laws
    (newp/newm/new12) on all 184320 rows; EVP momentum chain
    vrel→Cb→cca/ccb→ab2→cc1/cc2→uvel/vvel on all 40320 rows (28 of 32
    traced T cells have ice-covered NE U corners); row counts exact
    (thermo/dyn/transport 32×5×2×12=3840 each, budget 1 row at step 24).
    Run record: `stage1/runs/cice_instrument.json`.
  - Bugs found so far: (a) icepack_atmo `strx=tau·uatm` missing `/wind`
    (already in bugs.md); (b) CICE nuopc cap never calls CICE_Finalize
    (dump_finalize wired in ModelFinalize instead); (c) **confirmed**:
    CESMCOUPLED ice_import computes Tf with shr_frz_freezetemp keyed to
    the *mediator's* tfrz_option, never icepack_sea_freezing_temperature
    (CICE's tfrz_option); numerically identical here because both
    resolve to the same mushy liquidus at ocean salinities — the
    decoupling only bites for non-mushy settings; (d) CIME nmlgen drops
    groups absent from buildnml's hard-coded list (see above). Worktree
    left fully clean (cam nested-clm pointer state + FTorch gitlink
    restored; orphaned clm gitdir remains inside .git metadata, inert).
    **Next: fill CICE stage-2 skeleton tests from dumps (repeat the M2
    loop)**.

- **MOSART stage-2 tests from dumps** ✅ (2026-10-26) — worktree
  `wt/mosart-tests` (`../EqESM-wt-mosart-tests`, branch `wt/mosart-tests`):
  filled the empty `tests` blocks of all five stage-2 files from the
  stage-1 instrumentation dumps (`stage1/dumps/mosart/`, C_JRA_t62 run) —
  hillslope 6 tests/20 assertions, subnetwork 8/30, channel 4/33,
  coupling 1/16, routing assembly 1/16; **115 assertions, all green**
  (`esm test`, default reltol 1e-4; `esm validate` clean). Every number
  traced to a dump row (step/ns/nr/gindex); EqIce house style:
  instantaneous laws at t=0 (`time_span {0,0}`) + one-step storage traces
  at t=3600. Continuous-ODE trace tests document measured scheme gaps vs
  Fortran (assembly: wt 2.063e-6 vs 1.227e-6, wr 8.52e-16 vs 2.62e-16 —
  wr pinned only to dry threshold, wt not asserted; volr/wh asserted to
  rel 1e-5). Findings: `fthresh` default 0.0 in Control is a placeholder
  (Fortran uses abs(spval)=1e36 — assembly overrides); channel dry
  through all 8 dump steps (max wr=1.9e-12 m3) so channel tests pin only
  dry branches; `runoffocn`=spval for mask==1 never enters `runofftot`
  (zeroed in .esm; no discrepancy). `mosart_templates.esm` untouched.
  README gains a Testing section (test table, derived-input conversions,
  scheme gaps, coverage limits).
