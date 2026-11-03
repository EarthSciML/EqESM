# PROGRESS — eqesm implementation

Tracking progress on PLAN.md (Stage 1 — instrument the CESM codes; Stage 2 —
stubs → physics → EarthSciModels PRs; Stage 3 — top-level `eqesm.esm`).
Last updated: 2026-11-02.

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
  2026-11-01: manager instructions moved to `CLAUDE.md`; workers are now
  DeepSeek V4 Flash opencode slurm jobs (`tools/slurm-agent/`), tracked in
  `jobs.json`.
- **esm CLI 0.7.0 + repo cleanup** (2026-11-01) — `./esm` upgraded
  0.6.0 → 0.7.0 (release binary, sha256 verified; breaking unit/table
  checks, new `units_dimensionless_literal_mismatch` warning). Regression
  check `esm070-regress` (15 min): **no regressions** — `esm test stage2/`
  249P/0F/0E, `esm validate` clean on all 12 files (needs `ESD_ROOT`
  exported for the CICE B-grid rule imports). 0.7.0 warns about 66
  dimensionless literals (cice_coupling 23, cice_thermodynamics 26,
  cice_dynamics 3, cice Control 3, mosart channel 5 / subnetwork 5 /
  hillslope 1; 129 counting mounts). **Fixed** by `units-literals`
  (merged `06e80cf`, 4 h): EarthSciAST codemod + hand review, warnings
  66 → 0, `esm test stage2/` 249P/0F/0E unchanged; values verified
  unchanged by a const-unwrapped JSON comparison except two **transcription
  fixes the unit checker exposed** (both untested branches):
  `cice_thermodynamics` ustar was `(τ²/ρw)^(1/4)` instead of
  `sqrt(|τ|/ρw)` (icepack_therm_vertical.F90:665), and
  `meltpond_sealvl_dvn` mixed rates and per-step amounts
  (`(meltt·ρi + melts·ρs + frain·dt)/ρfresh` with the .esm's m/s rates
  and pre-divided frain; icepack_meltpond_sealvl.F90:157). Also: MOSART
  `hydraulic_radius` template param `m` renamed `area` (a param named
  like a unit symbol shadows unit strings, esm-spec §9.6.3); manager
  unfolded two folded ρ·Lvap constants. **Follow-ups found:** (a) MOSART
  `channel_depth` trapezoid factor `0.1*(2.0/1.0)` = 0.2 vs Fortran
  SLOPE1 0.05 (mosart_physics.F90:501-532) — an untested .esm
  deviation, also in EarthSciModels#33 (**fixed 2026-11-02**, `e013d2c`); (b) salinity is declared `ppt`,
  which the esm unit registry reads as parts-per-trillion (CICE means
  parts per thousand, g/kg) — harmless while nothing converts it, but
  wrong; (c) `atmo_constant_tau` coefficient has velocity units while
  Control declares the stress in N/m². Cleanup: merged worktrees
  `wt/cice-dedd` / `wt/stage2-clm` removed; the uncommitted CLM skeletons
  from the latter rescued to `wt/clm-skeletons` (`485b260`);
  `units/dimensionless-literals` (codemod on a stale base, conflicts in
  CICE files) dropped locally — origin copy kept; redo the codemod on
  current main after the 0.7.0 check. CTSM W1b (`clm-w1b`, base
  `wt/clm-skeletons`) and W2 dump extension (`clm-dumps-w2`) dispatched.

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
- **CICE similarity atmo + Rothrock strength** ✅ (2026-10-28, merged
  `20d305f`) — worktree `wt/cice-atmo` (branch `wt/cice-atmo`): the
  Monin-Obukhov similarity boundary layer transcribed into
  `cice_templates.esm` as 18 atmo templates (neutral drag, relative wind,
  stability parameter, stable/unstable psi, exchange coefficients, 2-m
  psix2, Tref/Qref/Uref) with the 5-iteration unroll for the D-class
  natmiter = 5, plus the Rothrock strength chain (cumulative distribution,
  exponential participation, hrmin/hrexp/krdg, strength term). Coupling
  tests extended to 36 tests/47 assertions (Si_tref/Si_qref at rel 1e-9,
  Si_u10 with a documented rel 1e-4 gap from the uvelT_icep
  reconstruction, strength at rel 1e-9); emissivity fixed to the D-class
  0.95. `esm test stage2/`: 162P/0F/0E.
- **CICE thermodynamics re-transcribed to the D-class config** ✅
  (2026-10-28, merged `dfe2671`) — worktree `wt/cice-thermo` (branch
  `wt/cice-thermo`): the thermodynamics skeleton transcribed the icepack
  namelist defaults (BL99 7/1, bubbly, ccsm3, cesm pond), which the
  D-class run does not exercise; re-transcribed `cice_thermodynamics.esm`
  + `cice_templates.esm` to the run config: ktherm = 2 mushy
  (nilyr = 8/nslyr = 3), conduct = MU71, shortwave = dEdd,
  frzpnd = hlid sea-level ponds, snwredist = 'none', hs0 = 0,
  phi_i_mushy = 0.45, tfrz_option = mushy, congel_freeze = one-step.
  Thermodynamics model: mushy enthalpy states qicen/qsnon as unknowns
  with defaults (MOSART wr precedent), vertical conduction as
  semi-discrete flux-divergence ODEs, dEdd surface laws (snow fraction,
  effective pond area, band-weighted ice/snow albedos, absorbed/
  penetrating shortwave conservation), sea-level melt-pond volume change.
  Templates 63 → 76 (mushy liquidus/temperature-from-enthalpy/
  conductivity/enthalpy-of-melting/brine/liquid-fraction, dEdd laws).
  9 tests/10 assertions pin the dEdd surface laws against
  `cice_albedo.csv` to round-off (rel 1e-12): snow fraction on bare ice
  and snow, apeffn through the infiltration/hypsometry/hpmin/taper
  branches, ice and snow broadband albedos, shortwave conservation
  fswint = swabs − fswsfc − fswthru. Findings: (a) assertion syntax —
  the runner materializes shaped unknowns only under `coords`
  (element names `u[1]` are not found); (b) Fortran unit bug at
  icepack_shortwave.F90:3587-3588 — `dhs = hsnlvl − fsnow·dt` subtracts
  the snowfall mass flux (kg/m²/s) from a snow depth (m); the .esm
  reproduces the numerics via a dedicated m/s parameter `fsnow_dEdd`
  (bugs.md #10); (c) cice_transport.esm/cice.esm metaparameters moved
  to nilyr = 8/nslyr = 3 so the assembly's index-set merge stays
  deep-equal. Worktree had been left mid-write by the worker (compaction
  failure): file truncated to invalid JSON — repaired by the manager
  (missing metadata close brace + surplus trailing braces), then the
  fsnow_dEdd binding and coords conversions. `esm test stage2/`:
  172P/0F/0E (115 MOSART + 57 CICE).
- **CICE dynamics pointwise laws pinned in the coupling cap** ✅
  (2026-10-28) — worktree `wt/cice-dynpins` (branch `wt/cice-dynpins`):
  the compile block was confirmed empirically first (a trivial `{0,0}`
  `zetax2` test in `cice_dynamics.esm` fails whole-file compilation with
  `unlowered_operator: rewrite-target operator 'D' reached evaluation
  without being lowered`, esm-spec §4.2/§9.6.8 — reverted), so the
  dynamics pointwise laws were pinned in `IceCap` as algebraic probes
  (12 probe equations, 166 total): the classic-EVP stress maps
  (newp/newm/new12, ice_dyn_evp.F90:1584-1609), the strain invariant
  Delta and viscosities/replacement pressure (zetax2/etax2/rep_prs,
  ice_dyn_shared.F90:2156-2159, 2467-2471), and the stepu ocean-ice
  stress (vrel/taux/tauy/Cb, F90:817-818, 931-937) via the existing
  templates. Constants verified from the dump fit, not assumed: e_factor =
  e_yieldcurve²/e_plasticpot⁴ = 0.25 (all 184,320 stress rows),
  epp2i = 0.25, capping = 1, Ktens = 0 — no mismatch. Findings: the
  dumped taux/tauy are the stepu partials vrel·waterx/y (not the
  dyn_finish full form of `cice_dynamics.esm` eqs 26-27); the corner
  strain rates are area-weighted (m²/s) so the corner viscosities carry
  kg/(m²·s) (nominal kg/s in the Fortran); TbU = 0 in every traced row
  (Cb pins are wiring checks). 4 tests/22 assertions at rel 1e-12, all
  inputs bound via `parameter_overrides`, numbers from
  `extract_dynamics_tests.py`. `esm test stage2/`: 200P/0F/0E (115 MOSART
  + 85 CICE).
- **CICE thermodynamics flux-law pins** ✅ (2026-10-29, worktree
  `wt/cice-thermoflux`, branch `wt/cice-thermoflux`) — pinned the
  pinnable thermodynamics flux laws of `cice_thermo_flux.csv` into
  `cice_thermodynamics.esm` (3 new tests / 6 assertions):
  `thermo_surface_fluxes_step1_cell3_5` pins the cell-level aggregates
  flwout/fsens (step 1, cell (3,5)) at rel 1e-9 — the similarity-derived
  shcoef/lhcoef at the ENTRY Tsfcn (the atmo boundary layer runs before
  the Tsfc solve), the fluxes at the EXIT Tsfcn, the exit-aicen
  aggregation of icepack_flux.F90:204/209-211; the 1e-9 tolerance covers
  the ITD-weighting residual (~1e-10) of the merged dump values.
  `thermo_radiation_scale_step2_cell3_5`/`..._cell4_5` pin fswabs/fswthru
  (step 2) at rel 1e-12 through the prep_radiation rescale: the step-2
  thermo consumes the stored step-1 dEdd outputs × netsw(2)/netsw(1)
  (icepack_shortwave.F90:3893-3922), netsw(1) = the step-1 cice_flux exit
  scale_factor ("Save net shortwave", CICE_RunMod.F90:618-624), netsw(2)
  from the step-2 bands and the step-1 albedo state with the quirky
  per-unit-grid-area Σ aicen·al aggregation (NOT /aice); reproduced to
  1.5e-16. Model fixes: emissivity 0.985 → 0.95 in the
  `surface_heat_flux` template body/description and the model's flwoutn
  (ε = 0.95 is the D-class forcing_nml value the cap's 0.05·flw
  correction already assumed; ε = 0.985 misses the dump by 3.3e-3);
  freshn/fsaltn rewritten to the corrected rate laws — the F90:510-513
  evapn and new-snowfall terms cancel exactly, leaving
  freshn = −(ρi·(congeln+meltbn+melttn) + ρs·meltsn) and
  fsaltn = −ρi·dhi·ice_ref_salinity·0.001 with
  dhi = congeln+meltbn+melttn+evapin/ρi (no separate /dt); the rain
  drainage fresh = fresh + frain·aice (icepack_therm_itd.F90:2043) and
  the per-category dsnown = fsnow/ρs − meltsn + evapsn/ρs added, the
  merged fresh verified to 3e-17 and fsalt exactly at step-1 (3,5).
  Per-category meltt/meltb/congel/melts/fswthru/fswint renamed to
  melttn/meltbn/congeln/meltsn/fswthrun/fswintn (the cice.esm
  `Thermodynamics.fswthru → IceCap.fswthrun` variable_map updated), and
  18 cell-level merged aggregates added (aice, fsurf, fcondtop, fhocn,
  fresh, fsalt, flwout, fsens, flat, evap, fswabs, fswthru, fswint,
  congel, meltb, meltt, melts, dsnow — faq sums WITHOUT /aice, matching
  the F90 merges). Unpinnable laws documented in README-cice.md with
  residuals (flat ~6e-8 step 1; the step-2 surface fluxes and the
  growth/melt closures ~1e-6..1e-3; fhocn/fresh/fsalt/evap not pinnable).
  Numbers from `extract_thermoflux_tests.py` (hand-authored structure,
  in-script verification of every expected value). `esm test stage2/`:
  206P/0F/0E (115 MOSART + 91 CICE).
- **CICE thermodynamics flux-law pins, part 2 — latent heat +
  sublimation/condensation** ✅ (2026-10-30, worktree `wt/cice-fluxchains`,
  branch `wt/cice-fluxchains`) — extended the `cice_thermo_flux.csv` pins
  to flat/evap/dsnow (2 new tests / 4 assertions):
  `thermo_surface_fluxes_step1_cell2_4` pins flwout (rel 1e-8 observed
  2.6e-9) and flat (rel 1e-6 observed 1.0e-7) at cell (2,4) step 1;
  `thermo_sublimation_condensation_step1_cell2_4` pins evap and dsnow
  (rel 1e-6, observed 1.2e-7/1.1e-7; dsnow expected = dump/dt = 7.88e-10
  m/s since the dump volume tendency is [m/step]). `esm test stage2/`:
  210P/0F/0E (115 MOSART + 95 CICE). Findings (all dump-proven, README
  transcription fixes #9-#12): (a) **evapn condensation branch was
  missing** — the Fortran adds the condensation `dhs = econ/qsub`
  (icepack_therm_vertical.F90:1355-1399) to the same evapn accumulator as
  the sublimation `dhs = esub/qsub` (F90:1502-1521), so the closure is
  `evapn = (esub + econ)/qsub = −flatn·dt/qsub`; the esub-only model gave
  −5.6e-8 vs the dump +2.04e-7 kg/m²/s at (2,4) (net condensation; the
  `econ` law was dead code). Fixed in both branches. (b) **The merged
  flux aggregates are ENTRY-aicen weighted** — merge_fluxes sums
  `aicen_init·xn` (F90:3169, save_init ice_step_mod.F90:80), correcting
  the earlier exit-aicen assumption of the (3,5) pin (harmless there:
  remap ≤ 1.8e-5, exit weighting approximates the merge to ~1e-10;
  dump-proven at (2,4) where aicen remaps up to 3.7e-4: entry-weighting
  residuals 2.6e-9/1.0e-7 vs exit-weighting 4e-6/9.9e-5). (c) **The
  residual mechanism for all surface-flux columns**: merged values are
  evaluated at the converged-solve Tsfc, and the exit-dump Tsfc differs
  by ~1e-7..3e-6 K (linear_itd tracer remap) — this, not the weighting,
  sets the 1e-6 tolerances. (d) **dsnow dump units are [m/step]** — the
  rate-law expected values divide by dt (documented in Units). (e)
  **meltsn bare-ice guard added** — kh_top = 2·kslyr/hslyr diverges for
  the bare-ice cat (hslyr → 0) and the model's unconverged etop_mlt
  poisoned dsnown; the Fortran clamps `dhs = max(−dzs, etop_mlt/zqsn)`
  (F90:1530) so the .esm now guards with the hs > puny branch (0 for
  bare ice). The evap/dsnow closures evaluate qsub with the post-solve
  state, so the test carries the EXIT qicen/qsnon (initial_conditions)
  and EXIT zSin (parameter) as the solved-state proxy (entry state:
  evap rel 1.4e-4). Unpinnable and documented: fsurf/fcondtop
  (solve-internal zTsn(1)/Tsfc; reconstructed fswsfc fails fsurf by rel
  0.78), congel/meltb/meltt (every dumped zTin > Tbot; meltt = solver
  residual 3.35e-20), melts (no dump column), fresh/fsalt/fhocn
  (lateral-melt contamination at step 2), step-2 flat (history-file
  uvel vs actual, rel 17%), dsnow at (2,6)-(2,8) (solve did not
  converge there: dump fsurf − fcondtop ≈ +0.3 W/m²). Numbers from
  `extract_fluxchains_tests.py` (`/scratch.local/ctessum`, hand-authored,
  %.17g, entry/exit dual-state verification to stderr). The (3,5)
  description was corrected to the entry-weighting framing.
- **CICE Tsfc surface-energy-balance pin** ✅ (2026-10-31, worktree
  `wt/cice-tsfc`, branch `wt/cice-tsfc`) — landed milestone-2 part (a):
  the per-category surface energy balance residual of the implicit Tsfc
  solve, as the diagnostic equation `fsurf_bal = fsurfn − fcondtopn`
  (eq 85; 84 → 85 equations, 83 algebraic + 2 enthalpy PDEs) pinned at
  the step-2 exit state of cell (2,4) (1 new test / 4 assertions, cats
  2-5, `abs 2e-2`). Key finding: the balance is `fsurfn − fcondtopn = 0`
  (the converged solve has fsurfn = fcondtopn, SAME sign — the
  downward-positive conductive flux enters the surface heat budget as
  −fcondtopn; the two-stage consistency checks `fcondtop − fsurfn <
  ferrmax`, icepack_therm_mushy.F90:535/588/834/892, and the etop_mlt
  closure F90:1337-1338). The earlier worker's ~−29 W/m² "residual" was
  the WRONG sum fsurfn + fcondtopn ≈ 2×fcondtopn; the difference closes
  to 6.8e-4..6.7e-3 W/m² at the exit state (the ferrmax = 1e-3 solve
  tolerance plus the step-2 bottom-melt regrid shifts, meltb(2) =
  8.76e-4 m; the dumped merged solve residual is fsurf − fcondtop =
  −3.3e-10). The measured residuals sit at 6.8e-4 (cat 2), 6.7e-3 (cat
  3), 2.0e-3 (cat 4), 8.8e-4 (cat 5); abs tolerance 2e-2 ≈ 3× the max.
  Test inputs: EXIT Tsfcn/zTsn (via exit qicen/qsnon + zSin) with ENTRY
  aicen/vicen/vsnon, fswsfc(2) = albedo_dump(1) × netsw(2)/netsw(1) =
  ×0.8451509697 (the prep_radiation rescale chain, pinned to round-off
  by the thermo_radiation_scale_* tests), similarity shcoef/lhcoef at the
  ENTRY Tsfcn with the uarea-weighted ice_ic uvelT. Cat 1 excluded:
  vsnon = 0 exactly at entry (bare ice) where the .esm kh_top =
  2·kslyr/hslyr diverges (the Fortran switches to the bare-ice branch
  kcstar(1) = km(1), calc_intercell_conductivity F90:2038-2051).
  Transcription fix #13: the file description and the Tsfc parameter
  description stated the balance as `fsurfn + fcondtopn = 0` (the
  equations were always correct); fixed with the dump proof (merged
  fsurf = fcondtop to round-off at steps 1-3; the step-2 per-category
  sum would read 2×fcondtopn ≈ −29 W/m²). Part (b), the conduction ODE
  trajectory test, remains BLOCKED and documented as such: the per-layer
  Iswabs/Sswabs are never dumped, and all 384 traced rows are daytime
  (coszen 0.158..0.572, fswintn 2..106 W/m² comparable to the
  conduction fluxes), so no zero-shortwave conduction-only trajectory
  exists to
  compare; the step-1 balance is also unpinnable (the step-1 thermo
  consumed the init-time dEdd solve, init_shortwave CICE_InitMod.F90:
  192-201, not dumped; the albedo-dump shortwave leaves the residual at
  ~−15 W/m²). Unblock list: (i) extend the stage-1 dumps with the
  per-layer Iswabs/Sswabs (or transcribe the dEdd 3-band solve), and
  (ii) trace a night step or otherwise zero the shortwave, for the
  conduction-only trajectory; plus (iii) the esm-spec implicit-equation
  construct or an EarthSciIntegrator-side nonlinear solve for the Tsfc
  coupling at stage 3. Harness probe (milestone-1 item 3, scratch
  `time_span {0,3600}` runs): the default bdf solver integrates the
  55-state enthalpy ODE system over a full step when every category has
  a finite kh (the state evolves ~1.8e6 J/m³ in the top ice layer), but
  collapses to step size ~0 on the bare-category kh_top = 2·kslyr/hslyr
  ≈ 6e14 stiffness (vsnon ≈ 0; the Fortran's bare-ice branch
  kcstar(1) = km(1) avoids it) — transcribing the lsnow guard into the
  .esm kh_top/fcondtopn is a prerequisite for any ODE-trajectory test.
  `esm test stage2/`: 214P/0F/0E (115 MOSART + 99 CICE). Numbers from
  `extract_tsfc_tests.py` (`/scratch.local/ctessum`, hand-authored,
  %.17g; the balance-residual diagnosis, similarity coefficients, and
  the fswsfc chain verified in-script).
- **CICE B-grid spatial-D lowering (dynamics + transport)** ✅ (2026-11-01,
  worktrees `EarthSciDiscretizations-wt-lowering` / `EqESM-wt-cice-lowering`)
  — the milestone that unblocks the dynamics/transport/assembly tests:
  EarthSciDiscretizations gains `grids/cice_bgrid` (ESD commit `11fe246`,
  unpushed): the 2-D B-grid with the strain-rate rules
  (`strain_rates_divu/tension/shear` lowering `D(u,x)+D(v,y)`,
  `D(u,x)−D(v,y)`, `D(u,y)+D(v,x)` to the four-corner staggered combinations
  of strain_rates, ice_dyn_shared.F90:2094-2175), the stress-divergence rules
  (`stress_divergence_strintx/strinty` lowering `D(σ11,x)+D(σ12,y)` /
  `D(σ12,x)+D(σ22,y)` to the stepu `uarear·(str1−str2+str3−str4)` integrals,
  ice_dyn_evp.F90:1698-1765, 958-961), and the conservative flux-divergence
  rule (`transport_flux_divergence` lowering `D(c·u,x)+D(c·v,y)` to the mass
  update `mm − w1·tarear`, ice_transport_remap.F90:3586-3590, with the
  centered edgearea_e/n reconstruction F90:658-666), all transcribed with
  line citations and six ast conformance fixtures (N_I = N_J = 4, NCAT = 2)
  whose goldens all five bindings reproduce byte-identically (54/54).
  Consumer (EqESM worktree): `cice_dynamics.esm` rewritten to the 3×5 traced
  subdomain (local N_I/N_J/NCAT meta + i/j/corner/ncat index sets, 34
  equations, four-corner stresses on a corner axis — no corner collapse),
  importing the 5 rules by reference via `${ESD_ROOT}`; `cice_transport.esm`
  rewritten to 3×5 with the flux-divergence rule, 3 algebraic div_flux
  unknowns (per-tracer: div_flux 1/s for aicen, div_flux_v/div_flux_s m/s
  for vicen/vsnon — a shared div_flux failed unit validation) + 3 continuity
  D-t equations; the assembly description updated (the unlowered-D block is
  gone; the assembly tests stay deferred to the shape-mismatched edges). New
  tests: dynamics 1 test/19 assertions (`time_span {0,0}`, rel 1e-12 —
  corner strain rates at T (3,5), strintx/strinty at U(2,4), Hibler P =
  546804.85434316727, full dyn_finish vrel/taux/tauy/Cb with the dump's
  aiX/uocn/vocn), transport 2 tests/16 assertions (integrated [0,3600]
  continuity trace to the t_exit state at (3,7): aicen ≤ 3.04e-5 / vicen ≤
  3.33e-5 / vsnon ≤ 1.78e-4 relative, per-assertion tols 1e-4/1e-4/2e-4 —
  the centered-reconstruction-vs-incremental-remap scheme difference
  O(C·Δc); plus the uniform-field instantaneous div_flux pin = 1.1240e-7
  1/s, rel 1e-12, integrator-independent after the bdf-vs-Euler
  exp(−D·t) vs (1−D·dt) restructuring). **`esm test stage2/`: 249P/0F/0E**
  (115 MOSART + 134 CICE = 69 coupling + 30 thermo + 19 dynamics + 16
  transport). README-cice.md rewritten (testing table, ESD_ROOT prerequisite,
  units, divergences 4/6, coverage limits, Fortran-bug table rows 4-8).
  Bugs and findings from this milestone (see also bugs.md #13-#18 (plus #19, the CLM urban entry)):
  (a) **5 Fortran stencil asymmetries** transcribed exactly — strain_rates
  shear corners 2/3/4 (ice_dyn_shared.F90:2170-2175) carry sign-inconsistent
  v/u terms vs corner 1 (F90:2168-2169) and the stress() str_dy k7/k8 atoms
  (ice_dyn_evp.F90:1760-1765) weight the p/m pair sums west-dominant
  (`p333·ssigpw + p166·ssigpe`) while k5/k6 (F90:1748-1749) are
  east-dominant; a symmetric transcription fails the pins; (b) **esm CLI
  RangeSpec bug**: CLI 0.4.0 AND 0.6.0 fail with "JSON parse error: data did
  not match any variant of untagged enum RangeSpec" whenever a faq range
  contains a string metaparameter ref (`[2, "N_I"]`) AND any template-library
  import carries bindings — both `esm validate` and `esm test` fail;
  workaround: declare the metaparameters/index_sets locally and import ONLY
  the rules (no grid.esm, no bindings) via `${ESD_ROOT}` refs (env expansion
  works when exported); (c) **extractor bug (fixed)**: the transport
  uvel/vvel pins were read from the MOMENTUM dump (m1) instead of the
  transport dump's t_entr phase — fixed with `source="transport"` and the
  pins file regenerated (only the uvel/vvel sections changed); (d) builder
  index() bug in the scratch extractor tooling (fixed; the uniform-field test
  is hand-authored, not injected); (e) **retracted the earlier claim of rel
  6.3e-4 agreement between the Hibler transcription and the dumped P** — the
  dump's P is the kstrength = 1 Rothrock form (321248.866048), not Hibler
  (546804.85434316727); the Hibler pin is a formulation pin, the Rothrock
  strength is pinned in the coupling cap; (f) the ocean-ice stress consumes
  the dump's **aiX = the cap's aice_init-based concentration**
  (ice_dyn_shared.F90:928-930) — deriving aiX from the model aicen's 4-cell
  U-point average misses the cap value by up to 2.6e-3, so aiX stays a dump
  parameter; (g) **U(4,4)-vs-uold RESOLVED**: the momentum IC convention —
  the dump cell U(4,4) = (0,0) velocity, and the stress-dump uvel column
  equals the momentum uold at U(i,j) exactly at every overlapping point; no
  discrepancy. Native compiler note: hand-written mixed-rank faqs are refused
  ("combine: array operand boxes differ") — rule-fired faqs and plain
  equations broadcast mixed ranks fine, so tmpcalc & co. are plain equations.
- **CICE mushy/snow temperature-inversion pins + nuopc-cap constants fix**
  ✅ (2026-10-28) — dump extension: the instrumented fork now dumps the
  mushy layer salinities sice_1..8 and the inversion outputs zTin_1..8 /
  zTsn_1..3 (icepack `6ce32d6` exposes `temperature_snow` via
  icepack_intfc; cice `44d019e` adds the `nt_sice` query and the columns;
  wrapper `8ef0fce`, superproject `80b0eef`; rerun job 11883231:
  COMPLETED, solution-neutral — all 30 pre-existing columns bit-identical;
  the new columns are physical: sice ∈ [0.10, 12.27] ppt, zTin ∈
  [−7.0, −0.008] °C, zTsn ∈ [−10.45, 0] °C, with 31/3840 no-snow rows
  dumping zTsn = +157.6 °C = Lfresh/cp_ice at qsn = 0 — dump guard note,
  bugs.md #12). **Root cause found (bugs.md #11):** the CESM nuopc cap
  (ice_comp_nuopc.F90:415-450, CESMCOUPLED) overrides ~20 icepack
  constants with SHR_CONST values: cp_ice = 2117.27, cp_ocn = 3996,
  cp_air = 1004.64, zvir = 0.6077930728241565, Lsub = 2,834,700,
  Lfresh = 333,700, snowpatch = 0.005, Tocnfrz = −1.836 — the first
  zTin/zTsn transcriptions used the icepack defaults (2106/4218/1005/
  0.606/2835000/334000/0.02) and were off by 0.3-2.6% against the dump.
  Fixed: all 16 affected templates in `cice_templates.esm` (incl. the
  mushy quadratic constants with cprho = cp_ocn·rhow = 4,099,896, F1/G1/
  H1/F2/G2/H2, D_liq/E_liq, I_liq), `cice_coupling.esm` inline zvir/
  cp_air (eqs 68, 69, 75), `cice.esm` Control Lheat binding (Lsub);
  descriptions cite the cap override. Validation: a Python transcription
  of `icepack_mushy_temperature_mush` with the run constants reproduces
  all 8 dumped zTin layers at rel 0.00e+00 (vs ~4e-3 with defaults), and
  the patched template bodies reproduce zTin/zTsn at rel ≤ 5e-16. New
  pins (`thermo_mushy_temperature_step1_cell2_4_ncat1` asserting zTin at
  layers 1/4/8, `thermo_mushy_snow_temperature_step1_cell2_4_ncat2`
  asserting zTsn at all 3 snow layers; qicen/qsnon via
  `initial_conditions`, zSin/aicen/vicen/vsnon via `parameter_overrides`,
  rel 1e-12; numbers from `extract_inversion_tests.py`, structure
  hand-authored) added to `cice_thermodynamics.esm`; liquidus Tmlts not
  directly pinnable (no dump column — the liquidus high-T branch is
  pinned via Tf(sss) in the coupling cap and its constants are exercised
  inside the zTin quadratic). `esm test stage2/`: 178P/0F/0E (115 MOSART
  + 63 CICE).
- **CICE upstream PR branches prepared (not pushed)** ✅ (2026-11-02,
  `cice-upstream-prep`, merged `0f3fcd4`, notes `stage2/UPSTREAM-cice-pr.md`)
  — local EarthSciDiscretizations `cice-bgrid-rules` @ `a0c0421` (B-grid
  rules + 6 conformance fixtures; all five bindings 1521P/0F/0S) and
  EarthSciModels `cice` @ `f7e05c5` (`components/cryosphere/cice/`, 6 `.esm`
  + README, strip-compare identical to `stage2/`), each with an untracked
  `PR_DESCRIPTION.md`. Rust CLI: 134P/0F/0E with the ESD branch, 99P/0F/3E
  against ESD `main` → **ESD PR must merge first**. Runner findings (files
  are spec-correct): python runner silently drops array-valued
  `parameter_overrides` (9-line repro; 69E coupling + 7E thermo); Julia
  tree-walk does not finish `cice_coupling.esm` (166 observeds) within
  any CI budget and rejects dynamics `sigma12` shape / thermo `index` op;
  RangeSpec bug (#16) still present on 0.8.0, workaround kept;
  EarthSciModels `Project.toml` still pins EarthSciAST 0.3 (no bump in
  this PR to avoid conflicting with #33). **Opened 2026-11-02** from the
  `ctessum-claude` fork: EarthSciML/EarthSciDiscretizations#50 (merge
  first) and EarthSciML/EarthSciModels#37 (description corrected for the
  final python/julia runner results). **Next: shepherd #50 then #37 through
  review; file the python-runner array-override and Julia-walker issues
  upstream.**
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

## CTSM lane (parallel to CICE; M5) — see `stage1/runs/CTSM_LANE.md`

Started 2026-10-31 as a second, non-overlapping lane alongside the peer
agent's CICE work. Scope: **all of CTSM**, with the explicit stance that the
existing EarthSciModels land/urban corpus is **presumptively unverified**.

- **Recon** ✅ (2026-10-31) — reports in `/scratch.local/ctessum/`:
  `clm_source_map.md`, `land_reuse_audit.md`, `ctsm_harness_recon.md`.
  Findings: CLM driver is `src/main/clm_driver.F90` (`Biogeophysics1/2/3` no
  longer exist); ~20 dump targets at file:line; urban is **on by default for
  clm6_0** and grid-driven (the global `I2000Clm60Sp_f19` run already
  exercises urban); `1x1_mexicocityMEX`/`urbanc_alpha` single-point runs need
  missing DATM inputdata (`urbanc_alpha` also a known upstream setup fail).
- **Reuse/trust audit** ✅ (2026-10-31) — `urban_canopy/` = 44 files, esm
  1.0.0, **0/44 cite a CESM `*.F90` line**: migrated from `UrbanCanopy.jl`
  (a Julia reimplementation of Oleson et al. 2010). It has tests (392/1032,
  green) but **none are CESM-grounded**; 2/44 fail `./esm validate`
  (`urban_canopy_model.esm` 63 eq / 0 unknowns; `urban_radiation.esm` θ_0
  units). Green tests ≠ CLM-correct.
- **CLM instrumentation (W1a)** ✅ (2026-10-31) — fork `stage1/EqESM-cesm-ctsm`
  (branch `instrument/clm-dumps`), `clm_dump_mod.F90` (559 lines) + namelist
  gating (`dump_clm`) + 7 wiring sites (additive dump calls only →
  solution-neutral). superproject `4181331`, clm submodule `ac5192b7d`.
- **Instrumented run** ✅ (2026-10-31) — case `I2000Clm60Sp_clm_inst`
  (I2000Clm60Sp/f19_f19_mg17, 36 tasks, 1 day), job `11961132` COMPLETED
  8:24. 8 CSVs in `stage1/dumps/clm/` (urban-radiation/albedo + surfrad +
  bareground + soilfluxes + qsat), one row per line, urban rows present
  (ctypes 71–75, l=252). Two dump-code defects fixed en route (urbalb format
  arity split rows; surfalbedo header order). Run record
  `stage1/runs/clm_instrument.json`.
- **CLMU urban radiation verified against CLM dumps** ✅ (2026-11-01) —
  EarthSciModels branch `clmu-verify` @ `e21ee29` (+365/−20 in
  `urban_canopy/urban_radiation.esm`; validate clean, 59/59 tests). Verdict:
  the radiation/albedo physics **is a correct transcription of the CESM CLM
  laws** — view factors bit-identical, roof SW/LW rel ≤7e-9, canyon LW
  multi-re-emission abs ≤2.1e-4 W/m². Fixed the objective `θ_0` units bug
  (`asin` result multiplied by `one_rad` → angle²; numerically a no-op).
  **No physics transcription error found.** The one systematic deviation is
  CLM's own loose absolute convergence cutoff (bugs.md §13: canyon
  multiple-reflection fixed-point exits at 0.001 W/m² LW / 1e-5 SW, ~1e-5
  rel); the `.esm` uses the exact closed-form inverse and is *more* accurate,
  so the LW-canyon tests assert at the Fortran cutoff. Capability gaps
  documented (per-surface `albgr*`/`lwup_*` and separate pervious/impervious
  roads are not exposed by the `.esm`). Not pinnable: daytime canyon
  shortwave — the only daytime urban landunit (255) has only its roof column
  traced; **dump follow-up: trace a daytime landunit with all 5 ctypes.**
- **CLM W2 dumps + layer-index fix** ✅ (2026-11-02, `clm-dumps-w2`,
  merged `8f0c24d`, 8.7 h) — 3 new slices (`clm_canopyflux`,
  `clm_run_imports`, `clm_soiltemp` 13056 rows), urban slices now by
  landunit id (daytime l=255 has all 5 ctypes, coszen up to 0.92). Root
  cause of shifted layer indices: gfortran passes arrays to assumed-shape
  `(:,:)` dummies as 1-based positional views → fixed with an explicit dim2
  lower bound (fork clm `2d8e0c697` / superproject `7896e49`); soiltemp now
  uses the model layer j (file slot j+nlevsno) and soilfluxes t_soisno_1..3
  are the real top-soil temperatures. Entry state matches finidat 34/34
  layers; neutrality **verified bit-identical** vs the f19 reference
  restart. Snow-layer inter-step jumps explained as DivideSnowLayers /
  SnowCapping behavior. bugs.md **#20**: SoilWaterMovementMod prints
  negative h2osoi_liq (−20.9 kg/m², col 182) with `endrun` commented out.
- **CLM stage-2 skeleton (W1b)** ✅ (2026-11-02, `clm-w1b`, merged
  `91c3be3`, 4.9 h) — `stage2/clm_canopy.esm` (CanopyHydrology +
  CanopyFluxInterface, 47 equations, 14 tests / 47 assertions at rel
  1e-12), `clm_templates.esm` trimmed to the 15 exercised templates,
  `README-clm.md`. `esm test stage2/` **296P/0F/0E** (249 + 47). Fixes
  made while testing: `roverg` = Rw/g·1000 mm/K (not Rw), `qsat_qsdT` had
  an extra 0.622. bugs.md #21 (dump index shift — the W2 defect above,
  now marked fixed) and #22 (VanGenuchten1980 module is a verbatim
  Clapp-Hornberg copy), renumbered from the colliding #20/#21. Caveats:
  pins were taken from the pre-fix dumps (values valid, j labels in test
  descriptions refer to the old layout); canopy tests mostly hit the
  bare-patch / zero-precip branches — the vegetated interception branch,
  liquid QSat branch, and h2osfc branch are unpinned (README coverage gaps).
- **Next (W2 cont.)** — (a) re-pin `clm_canopy.esm` against the W2 dumps:
  vegetated canopyflux rows (frac_veg_nosno=1, nonzero storage) and
  soilfluxes t_soisno_1..3, update provenance to the new j convention;
  (b) verify the remaining CLMU leaves (`temps/`, `hydro/`,
  `clmu_atmosphere`, `offline_mode_forcing`) and daytime canyon shortwave
  now that l=255 is fully traced; (c) next CTSM slices (soil temperature
  solve from `clm_soiltemp`, bare-ground / canopy fluxes).

## MOM6 lane (M7 stage-1) — see `stage1/runs/MOM6_LANE.md`

Started 2026-11-02 as the MOM6 stage-1 instrumentation lane (this worktree,
`wt/mom6-w1`). Complete: fork + dumps + validation + run record.

- **Instrumented fork** ✅ — worktree `stage1/EqESM-cesm-mom6` (branch
  `instrument/mom6-dumps`, local/unpushed): new module
  `MOM6/src/framework/mom6_dump_mod.F90` (auto-picked-up via the MOM6
  Filepath list), namelist-gated (`DUMP_MOM6`/`DUMP_DIR`/`DUMP_FIRST_CELLS`/
  `DUMP_FIRST_STEPS` read via get_param from MOM_override — no CIME
  namelist-definition change needed, unlike CICE), OMP-critical writes,
  PET-0-only. Three boundaries: (B1) Wright_FULL EOS ρ(T,S,p) + analytic
  derivatives per (cell, layer, thermo step) from inside KPP_compute_BLD;
  (B2) CVMix KPP OBL computation inputs/outputs + full KPP_calculate
  result (Kd_heat/Kd_salt/Kv, NLT, final Kt/Ks/Kv, smoothed OBLdepth);
  (B3) surface forcing import → fluxes/forces mapping per coupling step.
  Chain: MOM6 submodule `9955b6b31` (dump fixes 84f9c9a66 + LangEnhK init)
  → components/mom `f8bd526` → superproject `a6a28ad`.
- **Instrumented run** ✅ (2026-11-02) — case `C_JRA_t62_mom6_inst`
  (mirror of the validated C-class reference `C_JRA_t62`, T62_t232
  tx2_3v2, 72 PEs/2 nodes, 1 day), job `11992759` COMPLETED 0:0 (7:43
  wall). 4 CSVs in `stage1/dumps/mom6/`: eos 11700 rows (12 cells × 75
  layers × 13 thermo steps), kpp_obl/kpp 156 rows each (12 × 13), sfc
  forcing 144 rows (12 × 12). Traced set: all 12 ocean cells of PET 0
  (SW corner block; the 77 southernmost grid rows are Antarctic land).
  Trace pattern: CESM startup ocean-lag (mom_cap.F90:1835-1856) → the
  forcing dump fired 12× (OS%Time tods 0, 7200..43200) and the KPP dumps
  on all 13 thermo steps covered (tods 3600..46800; cpl_step 1 twice).
  Two dump-code defects fixed en route: dump_open missing path separator
  (files landed as `mom6mom6_*.csv` one dir up) + per-cell step counter →
  per-interval counter (84f9c9a66); uninitialized LangEnhK denormal
  1.636e-311 in the kpp CSV (LT_K_ENHANCEMENT=False) → explicit else
  branch. Run record: `stage1/runs/mom6_instrument.json`.
- **Validation** ✅ (all green, `validate_mom6_dumps.py`): solution
  neutrality — all 6 restarts + ocean.stats.nc + Vertical_coordinate.nc
  bit-identical to the reference (`compare_mom6_restarts.py`, 0 diffs);
  EOS law bit-exact on 11700 rows (max rel err 0.0); hydrostatic
  dp/dh = RHO_0·g = 10143 constant over 11544 interfaces; KPP update rule
  (KPP_IS_ADDITIVE=False) on all 11856 interfaces; NLT structure (all
  traced columns stably stratified); OBL/kOBL consistency + z_inter
  cross-match on 156 columns; forcing latent/sw/evap/p/ice/liq/melt/tau/
  ustar/omega laws 144/144 and salt_flux semantics pinned
  (RESTORE_SALINITY=True: salt_flux_o = salt_flux_in + restoring flux;
  the mediator's import is 0 in this run).
- **Next:** stage-2 EOS component first (`mom6_eos.esm` — closed-form
  pointwise law, the dump pins it to 1e-12 relative), then KPP, then the
  surface forcing mapping (PLAN.md M7).

## Lane round 2026-11-02 (manager review of 6 workers)

- **MOSART channel_depth fix** ✅ (`mosart-channel-depth`, merged `e013d2c`,
  0.9 h) — GRHR trapezoid branch had two transcription defects: the
  `0.1·(2/1)` factor (4× too deep) **and** a dropped `rdepth` offset.
  SLOPE1 is now a named dimensionless template parameter; 5 new tests / 20
  assertions pin every nonzero GRHR/GRPR/GRRR branch against a compiled
  gfortran driver of the verbatim F90 routines (rel 1e-12). `stage2/`
  269P/0F/0E. Ported to EarthSciModels `mosart-routing` as local commit
  `705518c` (**not pushed**; PR #33 still needs the push).
- **CDEPS DATM** ✅ (`cdeps-datm`, merged `64b7b38`) — instrumented
  strdata + GSWP3/CLMNCEP export dumps (4800 + 192 + 192 rows),
  solution-neutral on a different PE layout, 4992 rows law-checked. Stage 2:
  `stage2/cdeps_datm_gswp3.esm` (Streams/Clmncep/GSWP3Files, 4 netcdf
  `data_sources`) + `cdeps_templates.esm` (10 templates), 58P/0F/0E. First
  working netcdf `data_sources` recipe (see README-cdeps.md). Small
  follow-ups: `Faxa_lwdn`'s `e2` repeats the `vapor_pressure` law; bare
  literals lack units (check `./esm units`).
- **CISM lane start** ✅ (`cism-w0`, merged `94c330f`, 11.9 h) — glissade
  instrumentation fork (`instrument/cism-dumps`), 7 validated traces,
  solution-neutral; `stage1/runs/CISM_LANE.md`, `cism_instrument.json`.
  bugs.md #23–26.
- **CAM6 SCAM lane start** ✅ (`cam-scam-w1`, merged `5f7f674`) — 26
  site-calls/step ptend/cam_in/cam_out/qsat dumps (FSCAMARM97, 48 steps),
  bit-neutral vs pristine; qsat law round-off on 3072 rows. Reconstructed
  the missing `domain.lnd.ne3np4_gx3v7.230718.nc` (ocean column). bugs.md
  #27 (CMEPS single-column driver unchecked netCDF errors). Open: the
  land-column flwds→0 abort is unexplained.
- **MOM6 lane start** ✅ (`mom6-w1`, merged `48db9bd`) — see the MOM6
  section above.
- **CMEPS aoflux** ◐ (`cmeps-aoflux`, partial merge `da3c97c`) — stage 1
  merged: 288-row dump of flux_atmOcn_large on C_JRA_t62_inst, all 46
  computed columns replayed to 4.4e-15, restarts md5-identical; bugs.md
  #28 (rainc uninitialized allocation) / #29 (gfortran es26.17 format).
  Stage 2 **sent back**: the stability iteration was unrolled by copying a
  16-equation block 5 times (136 equations) and the commit title "6P/312F"
  mislabeled 312 passes. Follow-up `cmeps-aoflux-refactor` (job 12016134)
  rewrites it as a recurrence over an iteration index (esm-spec §4.3.1.1).
- **Bookkeeping:** three parallel branches all claimed bugs.md #20+; they
  were renumbered on the branches before merge (CISM 23–26, CAM 27, CMEPS
  28–29). Three `st_archive.FSCAMARM97_ne3*` jobs (11987659, 11987710,
  12001102) from failed CAM runs are stuck at DependencyNeverSatisfied
  and need `scancel`.
- **Follow-ups (same day):** `cmeps-aoflux-refactor` failed silently (a
  reasoning turn hit the 32,768-token output cap; status "done", empty
  result). It found that coupled multi-variable recurrences are rejected as
  observed cycles (esm-spec §4.9.6). Resubmitted as `cmeps-aoflux-refactor2`
  (job 12017252: a one-step model chained 5× as subsystems). MOSART
  `705518c` **pushed** to EarthSciModels#33 (comment explains the GRHR fix;
  CI pending). **CISM dumps were lost:** `dump_dir` pointed inside the
  cism-w0 worktree, which `cleanup` deleted; `cism-basal` regenerates them
  into `$M/stage1/dumps/cism/` (1-minute rerun of the existing exe). All
  prompts now say that dumps live in the main checkout, and that bugs.md
  entries are written as `## TBD.` and numbered at merge.
- **Submitted:** `mom6-eos` (12017293), `cam-qsat` (12017294, includes the
  cdeps vp/units cleanup), `cism-basal` (12017295), `clm-w2-repin`
  (12017296), `clm-soiltemp` (12017297).
