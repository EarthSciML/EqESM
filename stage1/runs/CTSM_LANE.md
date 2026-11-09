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
  - W2a (2026-11-01): added 4 dump slices to `clm_dump_mod.F90` (fork clm
    `cef3e7256` / superproject `9e349ee`, branch `instrument/clm-dumps`):
    `clm_canopyflux.csv` (canopy interception/throughfall + CanopyFluxes
    boundary, stages ent|hyd|flx), `clm_run_imports.csv` (pre-flux column
    state, stages preflux|surfhum), `clm_soiltemp.csv` (snow/soil temperature
    solve, stages entry|exit, one row per column layer), and switched the
    urban slices (`clm_urbalb`, `clm_urban_rad`) from first-16-columns to
    landunit-by-id selection (all 5 column types of the first 16 urban
    landunits, so the daytime landunit is fully traced; `clm_urban_rad` gained
    a coszen column). Instrumented rerun: job **11971155** (in-node
    `case.submit --no-batch` fallback because the 5-node ctessum partition was
    fully occupied by agent sessions; backfill slot was 00:31).
  - Law spot-checks pass: fpiliq/fpisnow effective coefficients recovered
    from the traced fluxes match the CLM formula
    `fpiliq = 0.25*(1-exp(-0.5*(elai+esai)))` (`use_clm5_fpi=.false.`
    default), `fpisnow = 1-exp(-0.5*(elai+esai))`; `liqcan/snocan` update is
    `max(0, storage + dtime*interception)`; `qg` composition
    `frac_sno*qg_snow + (1-frac_sno-frac_h2osfc)*qg_soil + frac_h2osfc*qg_h2osfc`;
    soiltemp entry/exit `t_soisno` continuity across the solve.
  - W2b (2026-11-02): root-caused + fixed the soiltemp/soilfluxes layer
    indexing bug (fork clm `2d8e0c697` / superproject `7896e49`):
    **gfortran rebases assumed-shape dummy args to 1-based positional
    views**, so `(:,:)` dummies read every layer shifted by nlevsno slots
    (snow appeared at j=4..12; j=-8..0 read OOB garbage; W1 soilfluxes
    t_soisno_1..3 read unused snow slots = 0). Fix: explicit dim2 lower
    bound `-nlevsno+1` on all 2-D dummies; qflx_snomelt_lyr/snofrz_lyr
    (snow-only arrays) reads guarded `j <= 0` (NaN for soil). Also fixed
    the clm_soiltemp `# columns` header order (h2osoi before qflx, to match
    the data emission) and documented the j-layer convention (file slot =
    j+nlevsno). Rebuilt + rerun twice with `--no-batch` in allocation
    11967609 (final CSVs LID 261102-004337). Verified: entry(step 1)
    matches finidat exactly (34/34 layers, c=1 and c=14); liq+ice conserved
    entry→exit to 0.0; soil-layer cross-step continuity 0 bad; snow-layer
    differences (≤4e-3 K) fully explained by designed inter-step snow
    hydrology (DivideSnowLayers dzmax_u cascade + snow capping at
    h2osno_max=10 m WE + snowfall; reproduced analytically to ~3%);
    neutrality re-verified bit-identical vs the f19 reference restart
    (290/292 vars identical, 2 differ only in embedded path strings).
    bugs.md: SoilWaterMovementMod negative-h2osoi_liq print (endrun
    commented out) candidate entry.
- **W3 — soiltemp dump extension (2026-11-03, lane wt/clm-soiltemp2).**
  `dump_soiltemp` extended (clm `1c0d375dd` / superproject `687355c`):
  per-layer tkmg/tkdry/watsat/csol, sabg_lyr_col, fn; per-column nbedrock,
  frac_sno_fluxes, hs_top_snow, hs_soil, dhsdT, eflx_bot; the `zi` column
  split into correctly-labelled `zi_top` = zi(j−1) and `zi_bot` = zi(j)
  (the old column carried the interface *above* layer j — see the
  clm_instrument.json `w3.zi_labelling_defect`). Rerun in-node `--no-batch`
  in allocation 12027055 (LID 261103-111900); neutrality re-verified
  bit-identical vs the I2000Clm60Sp_f19 reference restart (290/292 vars,
  2 path-string-only diffs). NOTE: the fork checkout is shared with the
  concurrent clm-surfalbedo lane; their clm `3d7038304` / superproject
  `527ab5d` commits (surfalbedo albsnd/albsni + soil color + h2osoi_vol_1)
  are the parents of the soiltemp W3 commits and their clm_dump_mod.F90
  diff contains both lanes' dump changes.
  - W2c (2026-11-03): surfalbedo dump extension (fork clm `3d7038304` /
    superproject `527ab5d`): `clm_surfalbedo.csv` gained the SNICAR snow
    albedos albsnd/albsni (direct/diffuse, vis/nir), the top-layer
    volumetric soil water h2osoi_vol_1 and the soil color class isoicol
    per column (urban columns write NaN in the albsnd/albsni slots — they
    do not run SNICAR). Instrumented rerun inside allocation 12027802
    (in-node `--no-batch`), LID 261103-133002. The first attempt (LID
    261103-131613) aborted in `dump_surfalbedo`'s urban branch: the module
    passed the shr_infnan derived-type NaN parameter (`nan => shr_infnan_nan`)
    directly in a formatted es26.17 write list, which gfortran rejects at
    runtime ("Expected REAL ..., got LOGICAL"); fixed by materializing a
    real local through the module's `assignment(=)` interface (bugs.md #31
    entry). Verified: albsnd == albgrd and albsni == albgri bit-identical
    on all nourban rows (frac_sno_albedo = 1 endpoint of the ground-albedo
    blend), h2osoi_vol_1 = 0.1505, isoicol = 10, and the shared columns
    bit-identical to the W3 CSVs (the extension is dump-only; the model
    trajectory is unchanged).
  - W2d (2026-11-08): vegetated-patch dumps (fork clm `afbae49dd` on top of
    `65f12f0d8`, branch `instrument/clm-dumps`; worktree lane
    `wt/clm-veg-instrument`): new namelist `dump_veg_patches` (default 0),
    dynamic per-step top-N selection of the most-vegetated patches
    (qualifying: frac_veg_nosno=1, elai+esai>0, liqcan+snocan>0; ranked by
    elai+esai then liqcan+snocan then patch index), traced in 4 new CSVs:
    `clm_canopyflux_veg.csv` (CanopyFluxes vegetated-branch exit: leaf-temp
    iteration num_iter/det/dele + end-of-routine canopy storage update),
    `clm_twostream.csv` (TwoStream radiation law; coszen>0 patches;
    sun/shade vis only), `clm_interception.csv` (interception entry/hyd:
    fpiliq/fpisnow, throughfall, liqcan/snocan update, fwet/fcansno/fdry),
    `clm_soilwater.csv` (SoilWater solve entry/exit per selected column).
    Writer PET chosen once at step 1 by two integer allreduces (max
    qualifying count, tie-break lowest PET) — on this cold start all PETs
    tie at nveg=0, so PET 0 writes. Bring-up fixes: (1) the five dump_*
    namelist vars were never mpi_bcast in control_spmd, so only PET 0 read
    them and PETs 1..35 returned early from dump_select_veg, deadlocking
    PET 0's allreduce; (2) shr_mpi_min aliasing fix; (3) dump_twostream was
    inside `if (use_SSRE)` and never fired; (4) gfortran rebase fix (veg
    dummies take begp/begc/begg, arrays declared (begp:)/(begc:)/(begg:)).
    Final rerun in-node (`--no-batch`), LID 261108-072738 (run log
    /scratch.local/ctessum/clmveg/run_veg_final.log; diagnostic LIDs
    261108-064127/261108-071237 with markers, reverted). Verified:
    validator `stage1/analysis/verify_clm_veg_dumps.py` — ALL CHECKS
    PASSED (selection consistency incl. step-1 cold-start skip (traced
    steps 2..12), interception laws, canopyflux laws incl. the
    liqcan(hyd)+(tran-evap)·dtime identity, full TwoStream recomputation
    vs ctsm60_params xl, soilwater structure; `--existing` 11 old CSVs
    bit-identical; `--restart` neutrality 292/292 vars, 0 differ).
    `soilwater_movement_method=1` (moisture_form) in this case: the
    zengdecker2009 hk loop never runs, so clm_soilwater.csv has entry/exit
    rows only and the hk-law columns are NaN (documented in
    `stage2/README-clm-veg.md`; the pre-existing `layer, h2osoi_liq: 182 2
    -20.878` debug print at SoilWaterMovementMod:1422 and the step-1
    balance warnings are baseline behavior, identical in both runs).
- **W4 — PRs** into EarthSciModels (new `clm_*`; CLMU as test/fix PRs onto the
  existing files, plus fixing the 2 invalid files).

## Risks

- Cluster queue is shared with the peer's CICE reruns; serialize submissions.
- `PROGRESS.md` is contended — keep edits minimal/delimited.
- Surface-dataset-driven subgrid: a dump is only meaningful with the subgrid
  census recorded (urban columns exist only where surfdata has urban landunits).
