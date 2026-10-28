# Bugs and issues found in the CESM Fortran/C++ codes

Centralized list per PLAN.md §8. Each entry records: location, description,
impact, whether the `.esm` reproduces or corrects it, and status. (CESM is
mostly Fortran; "the C++ code" in AGENTS.md means the reference
implementation regardless of language.)

Format template:

## N. <one-line summary>

- **Location:** `<file>:<line>/<function>`
- **Description:** what the code does vs. what it should do.
- **Impact on results:** none / quantitative effect.
- **Reproduced/corrected in .esm:** which `.esm` component matches or
  corrects it (or n/a).
- **Status:** open / patched locally / reported upstream / fixed upstream.

## 1. R2000MOSART standalone compset is non-functional as shipped

- **Location:** `components/mosart/cime_config/config_compsets.xml`
  (`R2000MOSART` alias, lname
  `2000_SATM_DLND%LCPL_SICE_SOCN_MOSART_SGLC_SWAV`), found via the M1
  compset survey (2026-10-25, CESM `cesm3.0-alphabranch` @ `e59f251`,
  mosart 1.1.13).
- **Description:** the compset is commented "NOT yet functional" and
  references `DLND%LCPL`, which is not a valid dlnd mode on this branch
  (valid modes per `components/cdeps/dlnd/cime_config/config_component.xml`:
  `SCPL`, `GCPL`, `RCPL`). `RCPL` ("drive mosart with data from a previous
  run's coupler history") is tested in the dlnd testlist, so the intended
  fix is `DLND%LCPL` → `DLND%RCPL` (or use an F/I-class compset, which run
  MOSART coupled to CLM).
- **Impact on results:** cannot create a standalone-MOSART reference case
  as shipped; workaround: `usermods_dirs`/local compset override with
  `DLND%RCPL`, or couple via F-class/I-class.
- **Reproduced/corrected in .esm:** n/a (run-configuration bug, not a
  physics bug).
- **Status:** open (recorded; workaround documented in
  `stage1/runs/COMPSET_SURVEY.md` §3).

## 2. MOM `tx10deg`/`tx3deg` testmods reference inputdata files that are not published

- **Location:** `components/mom/cime_config/testdefs/testmods_dirs/mom/{tx10deg,tx3deg}/`
  (added `f0f3f07`, 2026-02-21; grid-independence fix `e5edf77`, 2026-07-09).
- **Description:** the testmods' `user_nl_mom`/`shell_commands`/`user_nl_cice`
  reference six files per grid under
  `$DIN_LOC_ROOT/ocn/mom/tx2_3v2/tx10deg/` (e.g. `ocean_hgrid_tx10deg_20260219.nc`,
  `ESMF_mesh_tx10deg_20260219.nc`, `cice_grid_tx10deg_20260219.nc`,
  `ocean_topog_tx10deg_20260219.nc`, `seawifs-clim-1997-2010-tx10deg_20260219.nc`,
  `state_restore_tx10deg_20260219.nc`; same for `tx3deg`), but **no `tx10deg`/`tx3deg`
  directory exists on any inputdata server** — checked 2026-10-26:
  `svn-ccsm-inputdata.cgd.ucar.edu` (directory listing of `ocn/mom/tx2_3v2/` shows no
  `tx10deg`), the OSDF endpoint (`osdf-data.gdex.ucar.edu`, cert-expired anyway), the
  LCRC E3SM mirror (`web.lcrc.anl.gov`), and `ftp.cgd.ucar.edu` (all 404). The
  `check_input_data --download` run for a `C_JRA`/`T62_t232` case with the `tx10deg`
  testmod ends with "Could not find all inputdata on any server" (5 missing refs, all
  the same `ESMF_mesh_tx10deg_20260219.nc` plus the MOM files). No grid-generation
  tool exists in the checkout (`MOM6/` has no coarse-grid generator; the testmod
  commit message says nothing about generating them).
- **Impact on results:** cannot create the ultra-coarse MOM smoke cases
  (`C_JRA`+`tx10deg`, `mom/tx3deg` variants) as shipped; workaround: run the
  C-class case at full `tx2_3v2` resolution (inputdata present) or a coarser
  atmosphere grid without the testmod.
- **Reproduced/corrected in .esm:** n/a (run-configuration/inputdata gap).
- **Status:** open (recorded; may be worth reporting upstream — the files
  were presumably generated once and never uploaded).

## 3. MOM6 `AUTO_MASKTABLE=True` FATAL on 36 OCN PEs at tx2_3v2

- **Location:** `components/mom/MOM6/src/framework/MOM_domains.F90`
  `gen_auto_mask_table` (FATAL: "Couldn't auto-eliminate any land blocks.
  Try to increase the number of MOM6 PEs or set AUTO_MASKTABLE to False.").
- **Description:** a `C_JRA` case at `T62_t232` (full `tx2_3v2` ocean) with
  the default `AUTO_MASKTABLE=True` (MOM_input.yaml default for
  `TOPO_CONFIG="file"`) and 36 OCN PEs crashes in initialization: MOM6's
  mask-table layout search cannot eliminate any land blocks at this PE
  count/grid, calls `MOM_error(FATAL)`, then the error path segfaults
  (recurse into `mpp_error_basic`; backtrace frames repeat), and the batch
  job HANGS (srun never returns — the model aborted but the batch script
  waits) until manually cancelled. The FATAL is only visible in
  `run/ocn.log.*`, not `cesm.log`. The `mom/tx10deg`/`tx3deg` testmods set
  `AUTO_MASKTABLE = False` in `user_nl_mom` — the known workaround.
- **Impact on results:** C-class (active MOM6) cases at default PE counts
  fail to initialize; workaround: `AUTO_MASKTABLE = False` in `user_nl_mom`
  (no mask_table file needed — MOM6 then includes all-land blocks).
  Possibly related: `TARGET_IO_PES` layout search interacting with the
  36-PE/2-node decomposition.
- **Reproduced/corrected in .esm:** n/a (run-configuration issue; the .esm
  MOSART work is unaffected). Verified fix: `C_JRA_t62` runs to completion
  (23 MOM6 steps, CFL 0.075, mass err 6e-19).
- **Status:** open (recorded; workaround in `user_nl_mom`).

## 4. CICE `strx`/`stry` wind stress missing the `/wind` (unit-vector) normalization

- **Location:** `components/cice/src/icepack/columnphysics/icepack_atmo.F90:328-329` (high-frequency branch) and `:463-464` (standard branch, `icepack_atm`).
- **Description:** the commented intent is `strx = tau * uatm / vmag` with `tau = rhoa * ustar^2` (a stress magnitude, N/m^2). The standard branch computes `tau = rhoa*ustar*rd` (still proportional to wind speed) and then `strx = tau*uatm`, `stry = tau*vatm` — the division by the wind magnitude is missing, so the surface stress scales with wind speed instead of being a unit-vector-projected stress. The high-frequency branch is worse: `strx = tau*sqrt((uatm-uvel)^2+(vatm-vvel)^2)*(uatm-uvel)` with `tau = rhoa*0.0012*wind` — an extra wind-magnitude factor on top of the missing normalization.
- **Impact on results:** standalone CICE atmosphere surface stress (and hence ice drift forcing) is wrong by a factor of the wind speed when `calc_strair` uses these branches; the CESM-coupled D-class runs receive the atmosphere stress from the coupler (`Faii_taux`), so this path is not exercised there.
- **Reproduced/corrected in .esm:** the stage-2 CICE skeleton (`stage2/cice.esm` Control `strairx`/`strairy` parameters, cited F90:463-464) transcribes the active (buggy) law `strairx = tau*uatm` with its honest `kg/s^3` dimension and logs the missing-`/wind` quirk in `stage2/README-cice.md`; the correct `strx = tau*uatm/wind` unit-vector law is documented there as the divergence.
- **Status:** open (recorded; reported to the EqESM bug list per PLAN.md §8).

## 5. CICE `calculate_Tin_from_qin` aborts on negative quadratic discriminant

- **Location:** `components/cice/src/icepack/columnphysics/icepack_therm_shared.F90:89-94`.
- **Description:** for brine ice (`l_brine`), `Tin` comes from the quadratic `aa1*T^2 + bb1*T + cc1 = 0`; when the discriminant `csqrt = bb1^2 - 4*aa1*cc1` is negative the code calls `icepack_warnings_setabort(.true.,...)` and returns, aborting the whole run, instead of clamping or handling the unphysical (qin, Tmltk) combination. The guard comment in the stage-2 skeleton notes the risk.
- **Impact on results:** an unphysical enthalpy state anywhere aborts the model run rather than producing a local, bounded answer.
- **Reproduced/corrected in .esm:** n/a (guard noted in the CICE stage-2 skeleton docs).
- **Status:** open (recorded).

## 6. CESMCOUPLED CICE freezing temperature decoupled from CICE's `tfrz_option`

- **Location:** `components/cice/src/cicecore/drivers/nuopc/cmeps/ice_import_export.F90:817-832` (`ice_import`, `#ifdef CESMCOUPLED` branch).
- **Description:** under CESMCOUPLED the imported ocean freezing temperature is computed with `shr_frz_freezetemp(sss)` (share/shr_frz_mod.F90.in, Assur-1958 liquidus, clamped to >= -2 C), keyed to the **mediator's** `tfreeze_option`; `icepack_sea_freezing_temperature` (which honors CICE's own `tfrz_option`) is compiled only in the `#else` branch. CICE's `tfrz_option` namelist setting therefore has no effect on the Tf used in the coupled import path.
- **Impact on results:** with the default `tfrz_option='mushy'` on both sides the two functions coincide numerically for ocean salinities (verified against the stage-1 dump: all 384 `cice_freezing_temp.csv` rows match `Tf = max(S/(-18.48+0.01848*S), -2.0)` to round-off, and `icepack_liquidus_temperature_mush`'s high-temperature branch is the same law), so no difference in this run; a non-mushy CICE `tfrz_option` (or mediator mismatch) would silently use the wrong law.
- **Reproduced/corrected in .esm:** the stage-2 CICE coupling skeleton documents the import-side Tf law as actually used (shr_frz mushy); flagged so the .esm can reproduce either law explicitly.
- **Status:** open (recorded; confirmed numerically via the M3 dump).

## 7. CICE nuopc cap never calls `CICE_Finalize`

- **Location:** `components/cice/src/cicecore/drivers/nuopc/cmeps/ice_comp_nuopc.F90` (cap `Finalize`); the model-level `CICE_Finalize` exists in the other drivers (e.g. `cicecore/drivers/direct/hadgem3/CICE_FinalMod.F90:30`).
- **Description:** the nuopc cap's finalize routine does not call the model-level `CICE_Finalize`, so finalize-side work (diagnostics, and in our instrumentation `dump_finalize`) never runs in the CESM context. The stage-1 dump module's `dump_finalize` was therefore wired into the cap's `ModelFinalize` instead.
- **Impact on results:** none for the science; finalize-side bookkeeping (e.g. final output flushes) is skipped.
- **Reproduced/corrected in .esm:** n/a (harness-level; documented in the instrumentation notes).
- **Status:** open (recorded).

## 8. CIME nmlgen silently drops namelist groups absent from the writer's hard-coded list

- **Location:** `components/cice/cime_config/buildnml` (the `groups` list) + cime `nmlgen`/`namelist_definition_cice.xml`.
- **Description:** CIME's namelist generator validates `namelist_definition_cice.xml` (including unknown-to-the-writer groups) but only *writes* groups that appear in buildnml's hard-coded `groups` array; a new group (e.g. `dump_nml`) that is perfectly valid in the XML silently disappears from `ice_in` with no error. First D_cice_inst submit ran without `&dump_nml` (job 11857214 cancelled); fixed by adding `'dump_nml'` to the list (cice wrapper `16ff7fb`).
- **Impact on results:** namelist settings silently ignored — a footgun for any new CICE namelist group.
- **Reproduced/corrected in .esm:** n/a (tooling).
- **Status:** open (recorded; workaround: keep buildnml's `groups` list in sync).

## 10. CICE dEdd pond code subtracts the snowfall mass flux (kg/m2/s) from a snow depth (m)

- **Location:** `components/cice/src/icepack/columnphysics/icepack_shortwave.F90:3587-3588` (`shortwave_dEdd_set_eff`).
- **Description:** `fsnow` is documented "snowfall rate (kg/m^2 s)" (F90:3511) and is passed straight through from the coupler import (icedrv_flux `fsnow`, kg/m2/s). But the pond-lid snow-cover initialization computes `dhs = hsnlvl - fsnow*dt` (F90:3588) and tests `fsnow*dt > hs_min` (F90:3587), treating `fsnow*dt` as a snow **depth** (m). `hsnlvl` is a snow depth in m, so `fsnow*dt` has dimension kg/m2 — mixing units; the conversion `fsnow*dt/rhos` (as thickness_changes does at icepack_therm_vertical.F90:1632-1660) is missing. Numerically the code uses the kg/m2 value as if it were meters, overestimating the equivalent snowfall depth by the factor rhos (~330).
- **Impact on results:** the branch fires when a refrozen pond lid exists (`ipn > puny`) with `dhs` (the dhsn tracer) below puny and `fsnow*dt > hs_min` — in the D-class run 92/1920 dumped (cell, category) rows meet the condition (fsnow ~ 0-1.05e-5 kg/m2/s, dt = 3600 s, so fsnow*dt up to ~0.038 against hs_min = 1e-4). Where it fires, `dhs` (and hence the snow-on-pond taper `asnow = min(spn/hs1, 1)` and the retained-meltwater remap) is computed from an over-large snowfall depth. The effect is usually masked because the trailing guard `if (ipn*spn < puny) dhs = c0` (F90:3590) zeroes `dhs` for small ponds, but for larger ponds the taper uses the inflated value.
- **Reproduced/corrected in .esm:** the stage-2 dEdd template (`stage2/cice_templates.esm` `dEdd_effective_pond`, used by `stage2/cice_thermodynamics.esm` eq 55) reproduces the Fortran numerics exactly: the coupled mass flux is carried into the template through a dedicated `fsnow_dEdd` parameter declared in m/s, so `hsnlvl - fsnow*dt` types as length minus length and evaluates to the same number the Fortran produces (bugs.md #10 documented in the parameter description). The dEdd test pins (`thermo_dedd_apeffn_*`) pass to round-off against the dump.
- **Status:** open (recorded; transcribed as-run in the .esm with the quirk documented).

## 9. CICE dump module: `trcr`/`trcrn` tracer reads swapped the j and tracer indices

- **Location:** fork-only instrumentation `components/cice/src/cicecore/shared/ice_dump_mod.F90` (instrument/cice-dumps, `ced43eb`), `dump_thermo_state` (5 reads) and `dump_albedo` (5 reads, minus the corrected ones).
- **Description:** the tracer arrays are declared `trcr(nx,ny,ntrcr,nblk)` / `trcrn(nx,ny,ntrcr,ncat,nblk)` — index order (i, j, tracer, ncat, block) — but the dump module read `trcrn(i, nt, j, n, blk)` (i, tracer, j, n, block): the j and tracer indices swapped. Every tracer-sourced column therefore read the wrong array slot: Tsfc/Tsfcn/qice_1..8/qsno_1..3 in `cice_thermo_state.csv` and Tsfcn/alvl/apnd/hpnd/ipnd in `cice_albedo.csv` came out exactly 0.0 for all 3840/1920 rows. State-array reads (aicen/vicen/vsnon, albedo outputs) were unaffected. This is why the first tests-filling pass concluded "Tsfc not dumped" — it was a dump bug, not a model state.
- **Impact on results:** none on the model (dump-only code); the affected CSV columns were unusable until the fix. The stage-2 coupling Si_t pin had to wait for the fix.
- **Reproduced/corrected in .esm:** n/a (instrumentation). Fixed in the fork (`ced43eb` → wrapper `90dde50` → superproject `c8a3bef`); rerun job 11864874 verified solution-neutral and the columns now carry physical values (Tsfc ≈ −14.7 °C, qice_1 ≈ −3.03e8 J/m3). The Si_t export pin (`icecap_export_sit_step1_cell2_4`) validates the restored values to round-off.
- **Status:** fixed (2026-10-27; lesson: validate that every dumped column is non-degenerate against known physics before writing tests).
