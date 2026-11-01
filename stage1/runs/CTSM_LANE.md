# CTSM parallel lane (land, milestone M5)

Parallel to the peer agent's CICE work (M3/M4). Scope: **all of CTSM**, with the
explicit stance that the existing EarthSciModels land/urban corpus is
**presumptively unverified** and must be checked against instrumented CESM CLM
dumps before we rely on it.

## Stance — do not assume the existing CLMU is correct

- `EarthSciModels/components/urban_canopy/` = 44 files, esm 1.0.0, migrated
  from **`UrbanCanopy.jl@9200edd4`** (a *Julia reimplementation* of the Oleson
  et al. 2010 CLMU tech note) — **0/44 cite a CESM `*.F90` line**. It is a
  reimplementation of a reimplementation.
- It has tests (392 tests / 1032 assertions, all green) but **zero are
  CESM-grounded**; they reproduce the upstream Julia test files. A green run
  says nothing about CLM correctness.
- 2/44 fail `./esm validate`: `urban_canopy_model.esm` (63 equations / **0
  unknowns**; subsystems inlined not imported) and `urban_radiation.esm`
  (`θ_0 [rad] ~ asin(...) [rad²]` unit bug). Both still pass `esm test`.
- No README; the top-level model inlines its subsystems instead of importing
  the 40 leaf files.
- Rule (AGENTS.md): existing implementations get **new tests from instrumented
  CESM outputs**; where a dump disagrees, the Fortran is the reference; if the
  Fortran itself is wrong, log it in `bugs.md`.

## Ownership / collision guards

- **Ours:** `stage2/clm_*.esm`, fork `stage1/EqESM-cesm-ctsm` (branch
  `instrument/clm-dumps`), `stage1/dumps/clm/`, `stage1/runs/clm_instrument.json`,
  `stage1/analysis/*clm*`, this file.
- **Peer's (CICE):** `stage2/cice_*.esm`, `stage1/EqESM-cesm-cice`, `wt/cice-*`.
- **Coordinate:** `PROGRESS.md` (contended — add only a short pointer),
  `lib/cesm_constants.esm`, ESD templates, the shared cluster queue and
  `stage1/build`.

## Recon (2026-10-31)

Reports: `/scratch.local/ctessum/{clm_source_map,land_reuse_audit,ctsm_harness_recon}.md`.

- **Call tree:** driver `src/main/clm_driver.F90` (`clm_drv` @108), entered from
  `lnd_comp_nuopc.F90:913`. `Biogeophysics1/2/3.F90` **no longer exist** (folded
  into `clm_drv` + `BiogeophysPreFluxCalcsMod` + `SoilFluxesMod` +
  `HydrologyNoDrainage/DrainageMod`). Index space is a flat 1-D
  gridcell→landunit→column→patch hierarchy per clump; verticals
  `-nlevsno+1:nlevmaxurbgrnd` (nlevsno 12, nlevgrnd 15, nlevsoi 10, mxpft 78).
- **Namelist gating is manual:** `controlMod.F90` (`clm_inparm` group) +
  `clm_varctl.F90` + `bld/namelist_files/namelist_definition_ctsm.xml` +
  `namelist_defaults_ctsm.xml`. New `dump_clm` needs all four; OMP-critical
  writes (clump loop).
- **Top dump targets:** SurfaceAlbedo, SurfaceRadiation, CanopyFluxes,
  BareGroundFluxes, SoilFluxes, SoilTemperature/SoilThermProp/Phasechange,
  HydrologyNoDrainage, Infiltration/TotalSurfaceRunoff/Drainage, SoilWater,
  soil_hk/soil_suction, SnowWater/SnowCompaction/Combine/Divide, SnowAge_grain,
  SNICAR_RT, UrbanAlbedo/UrbanRadiation/UrbanFluxes, CNDriver.
- **Harness:** urban is **on by default for clm6_0** (`urban_explicit_ac=.true.`,
  `urban_hac=ON_WASTEHEAT`) and grid-driven via surface data — the existing
  global `I2000Clm60Sp_f19` run already exercises urban. Clean urban
  single-point = **`I1PtClm60SpRsGs` / `1x1_mexicocityMEX`** (1 PE);
  `1x1_urbanc_alpha` is a known upstream SETUP failure. 1x1 surfdata + DATM 1PT
  forcing are **missing** from `stage1/inputdata`.

## Waves

- **W0 — recon** ✅ (2026-10-31): source map, reuse/trust audit, harness recon.
- **W1 — instrumentation + skeletons.**
  - W1a: fork `stage1/EqESM-cesm-ctsm` + `clm_dump_mod.F90` (namelist-gated,
    PET-0, CSV) + first slice (urban/radiation/energy) + instrumented global run.
  - W1b: `stage2/clm_{biogeophysics,snow,soil,canopy}.esm` skeletons + README.
  - W1c: uninstrumented single-point urban reference case + 1x1 inputdata.
- **W2 — dumps → independent tests.** Extract `stage1/dumps/clm/`; fill tests;
  **verify CLMU and the other reused pieces against the dumps**; `bugs.md`.
- **W3 — PRs** into EarthSciModels (new `clm_*`; CLMU as test/fix PRs onto the
  existing files, plus fixing the 2 invalid files).

## Risks

- Cluster queue is shared with the peer's CICE reruns; serialize submissions.
- `PROGRESS.md` is contended — keep edits minimal/delimited.
- Surface-dataset-driven subgrid: a dump is only meaningful with the subgrid
  census recorded (urban columns exist only where surfdata has urban landunits).
