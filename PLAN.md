# EqESM — CESM Implementation Plan (EarthSciAST `.esm`)

This document is the working plan for re-implementing the Community Earth
System Model (CESM, https://github.com/ESCOMP/CESM) as a composition of
EarthSciAST `.esm` components, following the staged process described in
`AGENTS.md`. It is a living document: update it as components are
instrumented, migrated, and merged. Progress is tracked in `PROGRESS.md`;
bugs found in the CESM Fortran/C++ code are tracked in `bugs.md`.

Sibling projects using the same workflow (see `AGENTS.md`):
[EqIce](https://github.com/EarthSciML/EqIce) (PISM) and
[EqWeFiC](https://github.com/EarthSciML/EqWeFiC) (WRF + fire + chemistry).
Their `PLAN.md`s, instrumentation patterns, and `.esm` factoring are the
prior art for everything below; this plan reuses their conventions rather
than re-inventing them.

## 1. Context and objectives

- Goal: a top-level `eqesm.esm` model, assembled by reference from
  individual `.esm` components, that reproduces CESM simulations.
- The original CESM implementation (hereafter "the Fortran codes"; CESM is
  Fortran, with some C++ in the infrastructure) is the reference for physics
  and numerics. Components that already exist in
  [EarthSciModels](https://github.com/EarthSciML/EarthSciModels) must be
  reused (and their tests extended from instrumented CESM outputs), not
  duplicated.
- All model logic lives in `.esm` files. Tests and examples live in the
  `tests` and `analysis` sections of `.esm` files.
- All local tests run through the EarthSciAST Rust CLI binary, kept untracked
  in this repo's root. Final CI gates run each file's inline `tests` through
  the EarthSciModels `tools/run_esm_inline_tests.py` / `run_esm_tests.jl`
  runners (the EarthSciModels single-pathway rule).
- **Time stepping for the fully coupled model is done with
  [EarthSciIntegrator](https://github.com/EarthSciML/EarthSciIntegrator)**
  (AGENTS.md). Its run-config concept — operators, macro-step control,
  partitioning — maps directly onto CESM's component/coupler architecture
  (§6).

## 2. Tooling and conventions

- **Format:** `.esm`, JSON, current version **1.3.0** (esm-spec §2;
  EarthSciAST **v0.4.0**, released 2026-10-25 — the local `../EarthSciAST`
  checkout is fast-forwarded to it). Every document we author declares
  `"esm": "1.3.0"`. Variables are `unknown` or `parameter`; everything else
  (ODE states, observed quantities, algebraic unknowns, parameter kinds) is
  derived (esm-spec §6.3). Version gates: a document uses a construct only
  at the version that introduces it — 1.1.0 brought the `solver` block
  (§2.2), 1.2.0 the `const` node `units` (§4.8.5), **1.3.0** the `solver`
  fields `external` (§2.2.6) and `time_span` (§2.2.7). The `solver` block
  declares stiffness/tolerances in the document itself (never in a gate),
  and `solver.external` references an external run-configuration document
  that EarthSciAST never opens — the hook for the EarthSciIntegrator run
  config at stage 3 (§6.2). `solver.time_span` declares the interval a
  document is meant to be simulated over.
- **Authoring stance (esm-spec §1.1):** AST first, closed-function registry
  second, factoring third. Use `expression_templates` liberally; import
  shared templates by reference (esm-spec §9.6/§9.7) instead of repeating
  logic. Components are mounted into subassemblies by reference with
  `variable_map` coupling edges (esm-spec §9.7.10), with import renaming /
  metaparameter rebinding (esm-spec §9.7.7) for grid/axis collisions.
- **PDEs:** write them as PDEs, not einsum/aggregate ops. Reuse (or submit
  PRs to)
  [EarthSciDiscretizations](https://github.com/EarthSciML/EarthSciDiscretizations)
  grid+stencil+rule templates; boundary conditions live inside the rule.
  Relevant grids today: `cartesian_uniform_{1d,2d}`, `cartesian_nonuniform_1d`,
  `column_nonuniform_1d` (atmospheric/soil columns), `latlon`, `latlon3d`
  (tracer transport), `duo` (gnomonic cubed sphere), MPAS TRiSK; plus
  conservative overlap regridding and reprojection templates (the coupler's
  core). **Known gaps to submit as ESD PRs when needed:** a tripole/
  displaced-pole grid (CICE, MOM), a GLL cubed-sphere / spectral-element
  grid (CAM-SE), and regridding onto/off of those grids.
- **Reuse:** check EarthSciModels (`components/<domain>/`), its `lib/`, and
  EarthSciDiscretizations first (see the inventory table in §3 for the
  per-component picture). The EqIce project's ice-sheet components (SIA,
  SSA, energy, hydrology, calving — `stage2/` of the sibling EqIce repo,
  pending merge into EarthSciModels) are the reuse path for CISM.
- **No scripts for equation authoring.** Scripts may only (a) extract numeric
  test tuples from CESM dumps into a hand-authored `tests` block, and (b)
  translate tabular data — KPP `.eqn` reaction lists (CAM chemistry),
  lookup tables (MOSART routing parameters, CICE albedo tables, CLM
  soil/veg tables), forcing datasets — into `reaction_systems`,
  `function_tables`, or `data_sources`.
- **Rust CLI:** a copy of the EarthSciAST CLI (`esm`) lives untracked in the
  repo root, built from the sibling `../EarthSciAST` checkout at **v0.4.0**
  (`cargo build --release -p earthsci-ast`; the released binary supports
  esm 1.3.0). Local testing is `./esm validate`, `./esm test`,
  `./esm units --check`. Julia/Python runners are used only to debug a
  suspected CLI fault or to match the EarthSciModels CI gate.
- **EarthSciIntegrator:** coupled-model time stepping (spec 0.1.0, Julia +
  Rust implementations with a shared conformance suite; local checkout
  `../EarthSciIntegrator`). Stage 3 consumes its run-config format
  (`esm_integrator` key, operators, macro-step control, partition,
  precision); the Rust implementation
  (`rust/earthsci-integrator*` crates) is the local path. The run config
  is referenced from `eqesm.esm`'s `solver.external` (esm 1.3.0, §2.2.6).
  Revisit the spec's status before stage 3 starts; gaps found are filed
  upstream, not worked around locally.
- **Workflow (AGENTS.md):** manager-worker arrangement. The manager
  dispatches up to three subagents at a time, each working in a git
  worktree; the manager checks, merges, and cleans up as they complete.
  Worktrees are created from `main`; every task closes with a merged commit
  and (for `.esm` work) passing inline tests via the Rust CLI.
- **Hygiene:** no large files in `/tmp` (RAM-backed); temporary files in
  `/scratch.local/ctessum` (may disappear); medium-horizon non-git files in
  a `.gitignored` subdirectory of this repo. Dumps and builds are gitignored
  (see §7). No "Co-Authored-By" or other LLM-authorship claims in commit
  messages.

## 3. CESM component inventory

Reference implementation: CESM **cesm3.0-alphabranch** (default branch),
with component submodules pinned at stage 1 (M1) to the tags listed below
(the `.gitmodules` values at the time of writing), the way EqIce pinned
PISM v2.3.2. Components map to `.esm` components/subassemblies; group by
the CESM source tree so each maps to a well-defined input/output boundary
for instrumentation.

| CESM area | Submodule (tag) | Main processes | `.esm` components (draft) | Existing impl? |
|---|---|---|---|---|
| Coupler | `components/cmeps` (cmeps1.1.73) | NUOPC mediator: import/export field mediation, merging/averaging over coupling intervals, conservative regridding between component grids, coupler flux computations (ocean/ice albedo, atm–ocn and sea-ice fluxes), mediator phases | `coupler_mediation.esm`, `coupler_fluxes.esm`, regridding via ESD templates | ESD regridding/reprojection yes; mediator no |
| Atmosphere | `components/cam` (cam6_4_212) | Dynamics (SE cubed-sphere spectral element, FV3 dycore submodule, vertical remap); physics (radiation per phys pkg: **cam6 → RRTMG, cam7 → RTE-RRTMGP** (rte-rrtmgp v1.9.2 submodule); CLUBB PBL; MG2/P3 microphysics; ZM deep + UW shallow convection; cloud macro/micro; gravity waves; dry deposition); chemistry (MOZART family via KPP `chem_proc`, MAM4 aerosol, photolysis via **tuv-x** v0.14.3; GEOS-Chem/HEMCO in the CAM-chem stack) | `cam_dynamics.esm`, `cam_radiation.esm`, `cam_pbl.esm`, `cam_microphysics.esm`, `cam_convection.esm`, `cam_cloud.esm`, `cam_chemistry.esm` (reaction_systems from KPP `.eqn`), `cam_aerosol.esm`, `cam_photolysis.esm` | RRTM LW (WRF) in `atmospheric_radiation/` — RRTMGP is a different scheme (k-distribution, GPU-capable), but the RRTM work's structure/constants transfer; `gaschem/fastjx` photolysis (CAM-chem path; CESM's tuv-x is the CAM photolysis library — verify which the reference runs use); `atmospheric_deposition/wesley_dry_gas.esm` (Wesely dry dep); `atmospheric_dynamics/holtslag_boville` (CAM4-era PBL family); `aerosol/dynamics` (check modal-aerosol overlap with MAM); no dycore, no CLUBB/MG/P3/ZM |
| Land | `components/clm` (ctsm5.4.063) | Biogeophysics (canopy, snow — incl. SNICAR, soil hydrology/thermal), biogeochemistry (CN, fire), urban (CLMU), irrigation, data assimilation | `clm_biogeophysics.esm`, `clm_snow.esm`, `clm_soil.esm`, `clm_canopy.esm`, `clm_urban.esm`, `clm_bgc.esm` | `urban_canopy/clmu_atmosphere.esm` (check provenance — candidate CLMU reuse); `land_surface/noah` (WRF's Noah — different model, same structural patterns); `vegetation/` (stand dynamics — check) |
| Sea ice | `components/cice` (cesm3_cice6_6_3_22) | Icepack column thermodynamics (Bitz–Lipscomb, mushy), EVP dynamics, transport (incremental remapping/advection), fluxes, coupling | `cice_thermodynamics.esm` (icepack), `cice_dynamics.esm` (evp), `cice_transport.esm`, `cice_fluxes.esm` | none (new `sea_ice` domain) |
| Ocean | `components/mom` (mi_260909) | MOM6: layered isopycnal/z\* dynamics on Arakawa C-grid (FMS framework), mixed-layer physics (KPP), isopycnal/neutral diffusion, mesoscale parameterizations, tides; coupling via FMS cap | `mom6_dynamics.esm`, `mom6_kpp.esm`, `mom6_mixing.esm`, `mom6_coupling.esm` | none substantial (check `environmental_transport/` for mixing pieces) |
| Ice sheets | `components/cism` (cismwrap_3_00_002) | SIA, SSA/higher-order stress balance, enthalpy/energy, calving, basal resistance, coupling | `cism_sia.esm`, `cism_ssa.esm`, `cism_energy.esm`, `cism_calving.esm`, `cism_coupling.esm` | **EqIce reuse path:** SIA/SSA/energy/hydrology/calving exist in the EqIce repo (`stage2/`, 3443+ tests) pending merge into EarthSciModels; add CISM-derived tests, don't re-derive |
| River routing | `components/mosart` (mosart1.1.13) | Network river routing: hillslope (overland) → subnetwork (tributary) → main channel, kinematic wave/variable storage with CFL substepping; irrigation + flood terms; direct transfer (glc/frozen/negative runoff to outlet); **no separate lakes/reservoirs/groundwater module in 1.1.13** (the `qgwl` stream from CTSM is routed directly per `qgwl_runoff_option`; see `stage1/boundaries.md` §2.5) | `mosart_routing.esm`, `mosart_coupling.esm` | `environmental_transport/surface_runoff` is Saint-Venant overland flow — related but NOT river routing; treat as new |
| Data components | `components/cdeps` (cdeps1.0.112) | datm/dlnd/docn/dice prescribed forcing | `data_sources` entries (JRA55, GSWP3, CPC, HadISST, …) | `earthsci_data/` has the pattern (era5, geosfp, ncep_ncar, …); CESM's own datasets are new `data_sources` |
| Legacy/optional | `components/rtm` (rtm1_0_89), `components/mizuRoute` (v3.1.1), `components/ww3` (main_0.1.1) | River transport (legacy), alternative routing, waves | defer; revisit if a reference compset needs them | — |
| Infra (build/test only, not modeled) | cime (cime6.7.1), ccs_config (ccs_config_cesm1.0.92), share (share1.1.22), libraries/parallelio (pio2_6_8), libraries/FMS (fi_260310), libraries/tuv-x (v0.14.3), libraries/mpi-serial, FTorch (optional), tools/* | Case creation, build, run, I/O, framework | not components; harness only | — |

Ordering rationale (also the §10 milestone order): **MOSART first** (smallest,
standalone-runnable, well-defined import/export boundaries — establishes the
full pipeline end to end), then **CICE** (column harness, bounded physics,
central to the coupled system), then **CISM** (reuses EqIce), then **CTSM**,
**CAM** (largest), **MOM** (largest + FMS), and finally the **CMEPS coupler**
and **CDEPS data** pieces, whose coupling tests depend on real component
imports/exports.

Each row becomes one or more component PRs into EarthSciModels (stage 2).
The top-level `eqesm.esm` (stage 3) composes the merged components.

## 4. Stage 1 — Instrument the CESM codes

Purpose: produce authoritative input/output traces to build tests from.
Status is tracked in `PROGRESS.md`; this section is the playbook.

### 4.1 CESM source access + build (M1)

1. Check out CESM `cesm3.0-alphabranch` with submodules (`git-fleximod`
   or `manage_externals`) into `stage1/cesm` (gitignored). **Pin** the
   submodule hashes/tags (see §3) and record them in
   `stage1/runs/README.md` for reproducibility.
2. **Toolchain.** The cluster has `gcc/13.3.0` + `openmpi/5.0.1` modules
   (the EqIce/ESD builds already use them). CESM needs netCDF-C/Fortran,
   PIO, and (for MOM) the FMS stack. Two routes, decided at M1 by probe:
   (a) build the netCDF/PIO stack with spack or by hand under
   `stage1/build/` (the EqIce local-deps pattern — expat/udunits2 were
   built this way); (b) an apptainer image (the EqWeFiC `wrf-build.def`
   pattern, which is proven in setuid mode on this cluster). **CIME machine
   config:** add a `config_machines.xml` entry for this cluster (or use
   the generic Linux target) — record it in `stage1/build/`.
3. **Smoke reference runs.** Create and run small cases with CIME
   (`create_newcase`/`create_test`): the compsets actually used at stage 1
   are chosen and verified here (exact names/options recorded in
   `stage1/runs/*.json`). Candidates:
   - F-class compset (e.g. `F2000climo`) at a coarse grid (e.g. `f19_f19`
     or the smallest SE alias `ne3_ne3_mg37`; there is **no `ne4_ne4`**
     alias on `cesm3.0-alphabranch` — see `stage1/runs/COMPSET_SURVEY.md`
     §4.1) for a few months — the atmosphere+land reference.
   - An all-data "X"-class compset with `rof = mosart` — the standalone
     MOSART reference (plus any CIME component-test route that drives
     MOSART alone). Note: the shipped `R2000MOSART` compset is broken
     (`DLND%LCPL` invalid; see `bugs.md` #1) — use the `DLND%RCPL` fix or
     an F/I-class coupled route.
   - I-class compset (active CICE, data atm/ocean) for sea ice.
   - C-class compset (active MOM, data ice/atm/land) for the ocean.
   - B/G-class fully coupled compsets for the coupler and stage 3.
   - **Single-column / single-point harnesses** (the EqWeFiC kernel-harness
     analog): CAM single-column mode (SCAM/FSCAM compset), CICE's Icepack
     column driver, CLM single-point (`clm1pt`), MOM6's own idealized test
     cases, and CISM's standalone mode. These give fast, small, branchable
     reference traces for stage 1/2 and are preferred over full 3-D runs
     for leaf components.
4. **Input data.** CESM needs `inputdata` (datasets, forcing, initial
   conditions — many GB). Curate the minimal set each reference run needs;
   record the source (NCAR servers / local mirror) and checksums in
   `stage1/runs/README.md`. This is a stage-1 risk item (§11).

### 4.2 Boundary catalog (`stage1/boundaries.md`)

Deliverable of stage-1 step 1 (the EqIce pattern): map the §3 inventory to
concrete CESM modules, types, and methods with signatures, inputs→outputs,
and **D** (instantaneous/derivative — preferred) vs **I** (integrated
trajectory) classification. For CESM the natural boundary layers are:

- **Component boundaries:** each component's NUOPC cap import/export field
  bundles as exchanged with CMEPS (`cam_import/export`, `lnd_import/export`,
  `ice_import/export`, `ocn_import/export`, `rof_import/export`, and the
  mediator's `esmFldsExchange` field list). These are the subassembly
  boundaries for the coupled model and for data-component coupling.
- **Subassembly boundaries inside each component:** physics parameterization
  interfaces (CAM physics buffer / `cam_in`/`cam_out`; CLM's land-unit /
  column / pft decomposition; CICE's icepack column API; MOM6's
  parameterization entry points; MOSART's routing/irrigation/flood/direct-transfer steps;
  CMEPS's mediator phases and flux computations).
- Record the exact run configuration (grid, timestep, coupling interval,
  namelist, precision — CESM components can run single- or double-precision)
  for reproducibility. Precision matters: single-precision reference values
  have a real32 floor, exactly as EqWeFiC found.

### 4.3 Instrumentation

1. **Fork strategy (the EqWeFiC pattern):** fork the CESM repo and the
   needed component repos under the `ctessum-claude` account with
   `earthsciml-instrumented` branches; keep the forks in
   `stage1/` checkouts (gitignored). Instrumentation lives in the fork,
   never in this repo.
2. **Dump mechanism:** a small ESM dump module (JSON or NetCDF per dumped
   call, the EqWeFiC `module_esm_dump.F` pattern) writing the boundary's
   inputs and outputs; env-var selection of what to dump. Where a boundary
   is already in history/restart output, prefer turning on history output
   over new hooks.
3. **Which dumps:** per component, a set of regimes covering every branch
   (the EqIce/EqWeFiC discipline: night/day, stable/unstable, land/water,
   ice/no-ice, wet/dry, cold/warm, snow/no-snow, …). Prefer
   **instantaneous derivatives**; integrated trajectories only where
   derivatives are unavailable. Also dump inputs/outputs of significant
   internal subassemblies (e.g. CAM's radiation vs. convection vs. PBL;
   CICE's thermodynamics vs. dynamics) to support later factoring.
4. Store curated dumps in `stage1/dumps/<component>/` (gitignored), with
   a `README.md` documenting layout and regeneration commands, and the run
   configs in `stage1/runs/`. Extract a select subset into numeric test
   tuples (script allowed, §2).

## 5. Stage 2 — Stubs → physics → EarthSciModels PRs

Bottom-up, leaf physics first, then subassemblies, in the §3 order
(MOSART → CICE → CISM → CTSM → CAM → MOM → CMEPS/CDEPS). For each
component:

1. **Hand-author a stub** `.esm` declaring the component's unknowns,
   parameters, and a `tests` block with the Stage-1 test tuples. The stub
   equations are placeholders (or absent) so tests initially fail/define
   the contract.
2. **Fill in the physics** so the tests pass:
   - Author equations as math (PDEs where the Fortran discretizes PDEs),
     importing discretization rules from EarthSciDiscretizations; submit
     ESD PRs for missing grids/operators (tripole, GLL cubed sphere,
     coupler regridding variants) as they are needed.
   - Factor repeated math into `expression_templates` and import them by
     reference. Shared CESM-wide constants/thermodynamics (e.g. the CESM
     `share/` module constants, the CAM/CLM/CICE shared physical constants)
     become a `lib/` family in the EqWeFiC `lib/wrf_constants.esm` style.
   - Translate KPP `.eqn` reaction lists (CAM chemistry), MOSART/CLM lookup
     tables, and forcing datasets with scripts into `reaction_systems`,
     `function_tables`, and `data_sources` (§2); equations and factoring
     are always hand-authored.
3. **Review + merge** each PR into EarthSciModels (`components/<domain>/`)
   after human review, following the existing layout and documentation
   style of sibling components. EarthSciModels CI has no skip/xfail: a PR
   opens only when its inline tests pass locally under the Rust CLI.
   Components that already exist in EarthSciModels get the Stage-1 test
   tuples added to the existing `.esm` instead of a new file.
4. **Cascade:** once leaf components are merged, create subassembly `.esm`
   files that import the leaves **by reference** and couple them
   (`variable_map` edges), adding subassembly tests. Repeat until the
   component's subassemblies pass.
5. **EqIce coordination:** the EqIce project's SIA/SSA/energy/hydrology/
   calving components are the CISM reuse path. Coordinate so that CISM
   work either (a) imports those components by reference once they are
   merged into EarthSciModels, or (b) adds CISM-derived tests to them in
   the same PR wave. Do not duplicate.

## 6. Stage 3 — Top-level `eqesm.esm`

### 6.1 Design: CESM-matching operator splitting, EarthSciIntegrator stepping

The coupled model is a set of PDE components exchanging fields through a
mediator at a coupling interval — exactly CESM's CIME/CMEPS architecture,
and exactly the shape EarthSciIntegrator's run config models (operators,
macro-step control, partition, precision). Mirror the CIME driver order:

1. **Components are operators.** Each merged component (or subassembly)
   is one leaf operator with its own state, its own timestep (component
   timesteps differ: e.g. CAM ~30 min, CICE ~1 h with subcycling, MOM ~
   minutes, MOSART ~hours), and its import/export fields.
2. **The mediator is the coupling layer.** Field exchange, merging and
   time averaging over coupling intervals, and regridding between
   component grids are expressed with ESD regridding templates and
   `variable_map` coupling edges in the CMEPS phase order (init → prep
   hydrology → run phases → post-run), reproducing the mediator's flux
   computations (ocean/ice albedo, atm–ocn and sea-ice fluxes) as
   components.
3. **Algebraic/equilibrium pieces** (where the harness refuses implicit
   algebraic equations — CONFORMANCE_SPEC §5.39, and CESM solves them
   iteratively anyway, e.g. the EVP ice dynamics) are realized as
   pseudo-transient relaxation to steady state with slow variables pinned
   (the EqIce stage-3 validated pattern: relaxed steady state = the
   Picard/Newton fixed point of the Fortran solver), or as EarthSciIntegrator
   partition entries — decided per piece at implementation time.
4. **Stiff pieces** keep the EqIce block pattern (implicit on small states,
   explicit on full grids; document stiffness via the `solver` block).
5. **The full-model C++/Fortran driver (reference data):** an
   `instrument_cesm.cc`-style driver is not needed — the instrumented
   component forks (§4.3) plus the coupled reference runs already dump
   every boundary per coupling step into `stage1/dumps/fullmodel/`
   (gitignored), with the run configuration in
   `stage1/runs/fullmodel.json`. These dumps are the stage-3 test
   references.

### 6.2 The `eqesm.esm` orchestration

- `stage3/eqesm.esm` mounts the merged components/subassemblies by
  reference and couples them with `variable_map` edges in the CIME/CMEPS
  phase order: atm ⇄ mediator ⇄ (ocn, ice, lnd, rof), with data
  components where a reference compset uses them.
- Time stepping of the coupled model runs through
  **EarthSciIntegrator** (AGENTS.md): a run config naming the operators,
  the macro-step (coupling interval), the partition, and the precision
  (real32/real64, matching the reference run). The run config is
  referenced from the model's `solver.external` (esm 1.3.0), and the
  simulation interval is declared in `solver.time_span`.
- Tests cover 3–5 orchestrated coupling steps against the per-step dumps
  (pinned inputs; documented tolerances — exact for pure functions,
  measured for scheme gaps and for implicit-vs-instantaneous gaps, the
  EqWeFiC lesson).
- Discrete events (CICE transport steps, calving, snow aging) stay step
  maps unless EarthSciIntegrator events mature; revisit per spec status.

### 6.3 Sequencing

1. Instrumented coupled reference runs + `dumps/fullmodel/` (first).
2. `eqesm.esm` orchestration skeleton + first coupled steps vs dumps.
3. Iterate on coupling order/timestep and subassembly discrepancies;
   extend to the full reference compset (e.g. a B/G-class compset) and
   compare integrated trajectories.

## 7. Repository layout (this repo)

```
EqESM/
  AGENTS.md            # process rules (authoritative)
  PLAN.md              # this document
  PROGRESS.md          # milestone tracking
  bugs.md              # centralized list of bugs/errors found in CESM code
  <rust CLI binary>    # untracked local EarthSciAST runner ("esm")
  stage1/
    boundaries.md      # subassembly boundary catalog (step 1 deliverable)
    dumps/             # instrumented traces (gitignored; referenced by tests)
    runs/              # run configurations + env
    instrument/        # instrumentation hooks/drivers for the forks
    build/             # build scripts + artifacts (gitignored)
    cesm/              # CESM + component checkouts, instrumented forks (gitignored)
  stage2/
    <component>.esm    # working stubs/physics before merging upstream
  stage3/
    eqesm.esm          # top-level coupled model
```

Components that pass review are removed from `stage2/` once merged into
EarthSciModels and imported by reference from there.

## 8. Centralized CESM bug list

Maintain `bugs.md` (format as in EqIce `bugs.md`). Record for each bug:
location (file:line/function), description, impact on results, whether the
`.esm` reproduces or corrects it, and status. CESM is mostly Fortran —
"the C++ code" in AGENTS.md means the reference implementation regardless
of language.

## 9. Test strategy summary

- Leaf physics: instantaneous-derivative tuples (preferred), from
  single-column/single-point harness traces where available.
- Subassemblies: composed-component tests with the same tuples, plus
  subassembly dumps (per-physics-driver traces, mediator field bundles).
- Full model: integrated trajectories vs. instrumented coupled CESM runs.
- All tests are inline `tests` blocks runnable by the Rust CLI locally and
  by the EarthSciModels CI gate on merge (all three bindings sweep the
  corpus; a red job names the binding — upstream fixes, never local
  workarounds).

## 10. Milestones / sequencing

Each milestone closes with merged PRs and passing inline tests.

1. **M0 — Scaffolding** ✅ (2026-10-25): PLAN.md, PROGRESS.md, bugs.md,
   .gitignore; AGENTS.md top-level-model naming fixed to `eqesm.esm`;
   next: Rust CLI binary in repo root (untracked) from `../EarthSciAST`.
2. **M1 — CESM source access + build; reference runs; boundary catalog.**
   Pin `cesm3.0-alphabranch` + submodule tags; toolchain (netCDF/PIO or
   apptainer) + CIME machine config; smoke runs (F-class, X-class with
   MOSART, I-class, C-class, single-column harnesses: SCAM, Icepack,
   clm1pt); minimal inputdata curated; `stage1/boundaries.md` published.
3. **M2 — Instrumentation scaffolding + MOSART end to end.** Dump module +
   fork branches; MOSART boundaries instrumented (routing/hillslope/
   subnetwork/channel one-step traces + pointwise laws; irrigation/flood
   rate rules; direct transfer); `stage1/dumps/mosart/`;
   `stage2/mosart_routing.esm` stub → physics, tests green under
   `./esm test`. This closes the loop for the whole pipeline.
4. **M3 — MOSART PR → EarthSciModels (merged); CICE instrumented.**
   Subassembly coupling of MOSART pieces; CICE boundaries
   (thermodynamics/icepack, EVP dynamics, transport, fluxes) instrumented
   via Icepack column + I-class runs; dumps + stubs.
5. **M4 — CICE physics; CISM via EqIce reuse.** CICE `.esm` components
   green and PR'd; CISM work starts from the EqIce components (coordinate
   merge timing), instrumented CISM traces added as tests.
6. **M5 — CTSM.** Biogeophysics (energy/water, snow, soil, canopy) →
   urban (CLMU reuse) → BGC; single-point + F-class traces; PRs.
7. **M6 — CAM.** Physics parameterizations (radiation/RRTMG, CLUBB,
   MG2/P3, convection, cloud) on SCAM traces → dynamics (FV lat-lon first;
   SE cubed sphere via ESD `duo`-family work) → chemistry
   (KPP `.eqn` → reaction_systems; photolysis verified vs tuv-x/fastjx)
   → aerosol (MAM4, reusing `aerosol/` pieces). Largest component; PRs
   land per parameterization, not all at once.
8. **M7 — MOM.** MOM6 dynamics + KPP/mixing on MOM6 test cases and C-class
   traces; FMS coupling boundary. PRs.
9. **M8 — Coupler + data.** CMEPS mediator (field mediation, flux
   computations) and CDEPS `data_sources`; coupling tests against
   instrumented mediator field bundles; ESD regridding PRs as needed.
10. **M9 — Subassemblies across components.** Cross-component subassemblies
    (atm+land, atm+ice+ocn flux exchange, land+river) import components by
    reference and couple; tests vs the coupled-run dumps.
11. **M10 — Stage 3: full model.** Instrumented coupled reference runs →
    `dumps/fullmodel/`; `stage3/eqesm.esm` orchestration +
    EarthSciIntegrator run config; 3–5-step coupled tests; then extended
    runs vs full CESM simulations; regression suite on the merged corpus.

## 11. Risks and open questions

- **Scale.** CAM alone is ~10⁶ lines; MOM6 and CTSM are similar. Mitigation:
  single-column harnesses for leaf physics, per-parameterization PRs, and
  heavy reuse of EarthSciModels/ESD/EqIce content. The full model will take
  years; milestones are open-ended per component and the plan is ordered so
  each stage delivers merged, tested content.
- **CESM 3.0 alpha is a moving target.** Pin tags at M1 and treat version
  bumps as deliberate updates (EqIce pinned PISM v2.3.2 for the same
  reason).
- **Build/infrastructure risk.** CIME case creation on a non-NCAR cluster
  needs a machine config; netCDF/PIO/FMS toolchain is the first hard
  dependency. The EqWeFiC apptainer route and the EqIce local-build route
  are both proven on this cluster; probe both at M1.
- **Input data.** CESM inputdata is large and versioned; curate the minimal
  set per reference run and record provenance. Check whether the cluster
  has a local CESM inputdata mirror.
- **Precision.** Single-precision components (CICE, CAM options) set a
  real32 floor on reference values; record precision per dump and expect
  the EqWeFiC tolerance story (real64 replays, dt → 0 limits, cancellation
  floors).
- **Grid coverage gaps.** Tripole (CICE/MOM), GLL cubed sphere (CAM-SE),
  and coupler regridding need ESD PRs; schedule them ahead of the
  components that need them (§5 step 2, §6).
- **EarthSciIntegrator maturity.** Spec is at 0.1.0 (M0 draft); confirm
  status and the Rust implementation's capabilities before stage 3, and
  file gaps upstream rather than working around them.
- **Open question:** which CAM6.4 chemistry/physics option set and which
  photolysis scheme the reference compsets use — pinned at M1 from the
  actual namelists of the reference runs.
