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

## 9. CICE dump module: `trcr`/`trcrn` tracer reads swapped the j and tracer indices

- **Location:** fork-only instrumentation `components/cice/src/cicecore/shared/ice_dump_mod.F90` (instrument/cice-dumps, `ced43eb`), `dump_thermo_state` (5 reads) and `dump_albedo` (5 reads, minus the corrected ones).
- **Description:** the tracer arrays are declared `trcr(nx,ny,ntrcr,nblk)` / `trcrn(nx,ny,ntrcr,ncat,nblk)` — index order (i, j, tracer, ncat, block) — but the dump module read `trcrn(i, nt, j, n, blk)` (i, tracer, j, n, block): the j and tracer indices swapped. Every tracer-sourced column therefore read the wrong array slot: Tsfc/Tsfcn/qice_1..8/qsno_1..3 in `cice_thermo_state.csv` and Tsfcn/alvl/apnd/hpnd/ipnd in `cice_albedo.csv` came out exactly 0.0 for all 3840/1920 rows. State-array reads (aicen/vicen/vsnon, albedo outputs) were unaffected. This is why the first tests-filling pass concluded "Tsfc not dumped" — it was a dump bug, not a model state.
- **Impact on results:** none on the model (dump-only code); the affected CSV columns were unusable until the fix. The stage-2 coupling Si_t pin had to wait for the fix.
- **Reproduced/corrected in .esm:** n/a (instrumentation). Fixed in the fork (`ced43eb` → wrapper `90dde50` → superproject `c8a3bef`); rerun job 11864874 verified solution-neutral and the columns now carry physical values (Tsfc ≈ −14.7 °C, qice_1 ≈ −3.03e8 J/m3). The Si_t export pin (`icecap_export_sit_step1_cell2_4`) validates the restored values to round-off.
- **Status:** fixed (2026-10-27; lesson: validate that every dumped column is non-degenerate against known physics before writing tests).

## 10. CICE dEdd pond code subtracts the snowfall mass flux (kg/m2/s) from a snow depth (m)

- **Location:** `components/cice/src/icepack/columnphysics/icepack_shortwave.F90:3587-3588` (`shortwave_dEdd_set_eff`).
- **Description:** `fsnow` is documented "snowfall rate (kg/m^2 s)" (F90:3511) and is passed straight through from the coupler import (icedrv_flux `fsnow`, kg/m2/s). But the pond-lid snow-cover initialization computes `dhs = hsnlvl - fsnow*dt` (F90:3588) and tests `fsnow*dt > hs_min` (F90:3587), treating `fsnow*dt` as a snow **depth** (m). `hsnlvl` is a snow depth in m, so `fsnow*dt` has dimension kg/m2 — mixing units; the conversion `fsnow*dt/rhos` (as thickness_changes does at icepack_therm_vertical.F90:1632-1660) is missing. Numerically the code uses the kg/m2 value as if it were meters, overestimating the equivalent snowfall depth by the factor rhos (~330).
- **Impact on results:** the branch fires when a refrozen pond lid exists (`ipn > puny`) with `dhs` (the dhsn tracer) below puny and `fsnow*dt > hs_min` — in the D-class run 92/1920 dumped (cell, category) rows meet the condition (fsnow ~ 0-1.05e-5 kg/m2/s, dt = 3600 s, so fsnow*dt up to ~0.038 against hs_min = 1e-4). Where it fires, `dhs` (and hence the snow-on-pond taper `asnow = min(spn/hs1, 1)` and the retained-meltwater remap) is computed from an over-large snowfall depth. The effect is usually masked because the trailing guard `if (ipn*spn < puny) dhs = c0` (F90:3590) zeroes `dhs` for small ponds, but for larger ponds the taper uses the inflated value.
- **Reproduced/corrected in .esm:** the stage-2 dEdd template (`stage2/cice_templates.esm` `dEdd_effective_pond`, used by `stage2/cice_thermodynamics.esm` eq 55) reproduces the Fortran numerics exactly: the coupled mass flux is carried into the template through a dedicated `fsnow_dEdd` parameter declared in m/s, so `hsnlvl - fsnow*dt` types as length minus length and evaluates to the same number the Fortran produces (bugs.md #10 documented in the parameter description). The dEdd test pins (`thermo_dedd_apeffn_*`) pass to round-off against the dump.
- **Status:** open (recorded; transcribed as-run in the .esm with the quirk documented).

## 11. CESM nuopc cap silently overrides ~20 icepack thermophysical parameters at init

- **Location:** `components/cice/src/cicecore/drivers/nuopc/cmeps/ice_comp_nuopc.F90:415-450` (`#ifdef CESMCOUPLED` block in the cap init; the values come from `share/src/shr_const_mod.F90.in` SHR_CONST_*).
- **Description:** at cap init the CESMCOUPLED branch calls `icepack_init_parameters` with CESM shared values, overriding the icepack namelist/defaults: cp_ice = 2117.27 (SHR_CONST_CPICE, vs icepack's 2106), cp_ocn = 3996 (SHR_CONST_CPSW, vs 4218), cp_air = 1004.64 (SHR_CONST_CPDAIR, vs 1005), zvir = 0.6077930728241565 (SHR_CONST_ZVIR, vs 0.606), Lsub = 2,834,700 (SHR_CONST_LATSUB, vs 2,835,000) and Lvap = 2,501,000 (unchanged) so Lfresh = Lsub − Lvap = 333,700 (vs 334,000), snowpatch = 0.005 (vs 0.02), Tocnfrz = −1.836 (vs −1.8). rhoi/rhow/rhos/gravit/vonkar/stefan_boltzmann/cp_wv/Tffresh/ice_ref_salinity are passed with the same values icepack already uses. The driver's own `icepack_init_parameters` call (icedrv_init.F90:1027) does not touch these (it only updates args present), so the cap's values are what the run uses.
- **Impact on results:** every enthalpy/temperature conversion shifts: e.g. `zTsn = qsn/(rhos·cp_ice) + Lfresh/cp_ice` is off by up to ~2.6% if transcribed with icepack defaults (cp_ice 2106 vs 2117.27, Lfresh 334000 vs 333700), and the mushy quadratic inversion constants (which embed cp_ocn·rhow − rhoi·cp_ice) are off by ~0.4%. The first zTin/zTsn transcription used the icepack defaults and failed the dump pins; with the run values all 8 mushy layers and 3 snow layers match the dumped zTin/zTsn to round-off.
- **Reproduced/corrected in .esm:** transcribed as-run: all 16 affected templates in `stage2/cice_templates.esm` (cp_ice/cp_ocn/cp_air/zvir/Lfresh/snowpatch, the mushy quadratic constants with cprho = cp_ocn·rhow = 4,099,896), the coupling model's inline zvir/cp_air (`cice_coupling.esm` eqs 68, 69, 75), and the `Control` Lheat binding (`cice.esm` eq 6, Lsub). Template descriptions cite the cap override. Verified to round-off by the `thermo_mushy_*` dump pins (2026-10-28).
- **Status:** open (recorded; a footgun for any CESM-coupled CICE transcription — always read the cap's overrides, not the icepack defaults).

## 12. CICE dump module: no-snow rows dump a garbage snow temperature (+157.6 °C)

- **Location:** fork-only instrumentation `components/cice/src/cicecore/shared/ice_dump_mod.F90` (`dump_thermo_state`, the zTsn_1..3 columns added in `44d019e`).
- **Description:** for the 31/3840 rows with no snow (qsno = 0), the dumped `zTsn` reads +157.6086... °C = Lfresh/cp_ice with the cap's constants — the `temperature_snow` result `qsn/(rhos·cp_ice) + Lfresh/cp_ice` evaluated at qsn = 0 — i.e. the temperature of "zero-enthalpy snow", not a physical snow temperature. The Fortran's snow temperature is only meaningful where snow exists; the dump module does not guard on vsnon.
- **Impact on results:** none on the model (dump-only); the affected rows must be excluded from any zTsn pin (the stage-2 test selects the ncat with vsnon > 0).
- **Reproduced/corrected in .esm:** n/a (instrumentation note; the extraction script picks a snow-present row).
- **Status:** open (recorded; harmless for tests, documented in the extraction script).

## 13. CICE strain_rates shear corner 2 sign asymmetry (ice_dyn_shared.F90:2170-2171)

- **Location:** `components/cice/src/cicecore/cicedyn/dynamics/ice_dyn_shared.F90:2170-2171` (`strain_rates`, the shear corner-2 combination).
- **Description:** the four corner shearing strain rates are not mutually consistent under a uniform B-grid pattern: corner 1 (F90:2168-2169) is all-negative `shear1 = -cym*v(i,j) - dyT*v(i-1,j) - cxm*u(i,j) - dxT*u(i,j-1)`, but corner 2 carries a sign-flipped v term `shear2 = -cyp*v(i-1,j) + dyT*v(i,j) - cxm*u(i-1,j) - dxT*u(i-1,j-1)` (the `+dyT*v(i,j)` would be `-dyT*v(i-1,j)`-pattern under a symmetric transcription).
- **Impact on results:** the computed shear (and hence Delta, the viscosities and the plastic stress) differs from the symmetric form; the pins cannot distinguish "bug" from "intended" because every dumped row is produced by the active code.
- **Reproduced/corrected in .esm:** reproduced exactly — the EarthSciDiscretizations `grids/cice_bgrid/stencils/strain_shear_corner2.esm` atom (ESD commit `11fe246`) and the `cice_dynamics.esm` pin (12 strain-rate assertions at rel 1e-12) match the dump; a symmetric guess fails the pins.
- **Status:** open (recorded; transcribed as-run, asymmetries 13-15 flagged as a cluster).

## 14. CICE strain_rates shear corners 3/4 sign asymmetries (ice_dyn_shared.F90:2172-2175)

- **Location:** `components/cice/src/cicecore/cicedyn/dynamics/ice_dyn_shared.F90:2172-2175` (`strain_rates`, the shear corner-3 and corner-4 combinations).
- **Description:** same cluster as #13: corner 3 `shear3 = -cyp*v(i-1,j-1) + dyT*v(i,j-1) - cxp*u(i-1,j-1) + dxT*u(i-1,j)` (F90:2172-2173) and corner 4 `shear4 = -cym*v(i,j-1) - dyT*v(i-1,j-1) - cxp*u(i,j-1) + dxT*u(i,j)` (F90:2174-2175) each carry sign flips (the `+dyT*v` and `+dxT*u` terms) relative to the corner-1 pattern.
- **Impact on results:** as #13 — the active code's corner shear values are what they are; a symmetric transcription misses the dump.
- **Reproduced/corrected in .esm:** reproduced exactly — the ESD `strain_shear_corner3/4.esm` atoms (commit `11fe246`), pinned at rel 1e-12.
- **Status:** open (recorded; transcribed as-run).

## 15. CICE stress() str_dy k7/k8 west-dominant p/m weighting (ice_dyn_evp.F90:1760-1765)

- **Location:** `components/cice/src/cicecore/cicedyn/dynamics/ice_dyn_evp.F90:1760-1765` (`stress`, the str(i,j,7)/str(i,j,8) y-divergence corner atoms).
- **Description:** the NW/SW corner atoms weight the p/m east-west pair sums west-dominant — `strp_tmp = p25*dxT*(p333*ssigpw + p166*ssigpe)` (F90:1754-1755) — while the NE/SE atoms k5/k6 use the east-dominant `p333*ssigpe + p166*ssigpw` (F90:1748-1749); the U-point stress interpolation is not symmetric under reflection.
- **Impact on results:** strinty differs from a symmetric reconstruction at the ~16.7%-weight level on the traced subdomain; the dump pins reproduce the active code.
- **Reproduced/corrected in .esm:** reproduced exactly — the ESD `str_dy_cell_k7/k8.esm` atoms (commit `11fe246`) and the `cice_dynamics.esm` strinty pin (rel 1e-12).
- **Status:** open (recorded; transcribed as-run).

## 16. esm CLI: faq-range metaparameter refs + bound template imports break RangeSpec parsing

- **Location:** tooling — the EarthSciAST `esm` CLI (local v0.4.0 and v0.6.0 binaries; `src/schema`/faq range handling).
- **Description:** `esm validate` and `esm test` both fail with `JSON parse error: data did not match any variant of untagged enum RangeSpec` whenever a faq range contains a string metaparameter reference (e.g. `[2, "N_I"]`) AND the model imports a template library with bindings. Either feature alone works; the combination breaks the RangeSpec round-trip parse.
- **Impact on results:** the intended design — a 3×5 consumer importing `grid.esm` plus the rules with metaparameter bindings — cannot compile; all stage-2 dynamics/transport workarounds depend on the local-metaparameter pattern.
- **Reproduced/corrected in .esm:** workaround only (2026-11-01): the consumer declares its own metaparameters (N_I/N_J/NCAT) and index_sets (i/j/corner/ncat) and imports ONLY the rule files (no grid.esm import, no bindings) via `${ESD_ROOT}` refs (env expansion works when exported). The workaround is documented in the `cice_dynamics.esm`/`cice_transport.esm` descriptions and stage2/README-cice.md.
- **Status:** open (tooling bug; probe files in /tmp/opencode/).

## 17. CICE dump extractor: transport uvel/vvel read from the momentum dump

- **Location:** fork-only scratch tooling `extract_bgrid_consumer.py` (/scratch.local/ctessum/cice_verify).
- **Description:** the transport test pins read uvel/vvel from the MOMENTUM dump (phase m1) instead of the transport dump's t_entr phase; the momentum uvel differs from the transport-side velocity the remap consumes (the transport dump's ncat-1 rows are the post-dynamics field the remap actually uses).
- **Impact on results:** none on the model (extraction-only); the transport pins generated before the fix carry the wrong velocity inputs.
- **Reproduced/corrected in .esm:** n/a (tooling). Fixed with a `source="transport"` parameter; the pins file was regenerated (only the uvel/vvel sections changed; the 249-assertion suite re-verified green).
- **Status:** fixed (2026-11-01).

## 18. CICE dump extractor: builder index() lookup bug (scratch tooling)

- **Location:** fork-only scratch tooling `build_consumer_tests.py` (/scratch.local/ctessum/cice_verify).
- **Description:** an index() lookup in the test-builder failed on the reordered/extended pin file (ValueError on a missing element after the uvel/vvel source fix), aborting the builder mid-injection.
- **Impact on results:** none on the model (tooling); the uniform-field test is hand-authored rather than builder-injected.
- **Reproduced/corrected in .esm:** n/a (tooling). Fixed in the script.
- **Status:** fixed (2026-11-01).

## 19. CLM urban canyon longwave/shortwave multiple-reflection solvers stop at a loose absolute tolerance, leaving ~1e-5 relative error

- **Location:** `components/clm/src/biogeophys/UrbanRadiationMod.F90:596-674` (`net_longwave`, `crit = max(road_a, sunwall_a, shadewall_a); if (crit < .001_r8) exit`); same pattern in `UrbanAlbedoMod.F90:1109-1250` (`net_solar`, `errcrit = .00001_r8`).
- **Description:** the CLM urban canyon multiple-reflection exchange is solved by fixed-point iteration with a per-iteration *absolute* convergence test on the absorbed flux (0.001 W/m² for longwave, 1e-5 per-unit-incident-flux for solar). Because the absorbed flux shrinks by a factor of only ~0.05-0.08 per iteration but the test is on the raw increment (not normalised), the iteration exits several orders of magnitude before machine convergence. The physically exact solution is a closed-form matrix inverse, so the dumped `lwnet_*`/`sabs_*` values carry a systematic truncation error relative to the exact law.
- **Impact on results:** measured against the exact linear-system solution of the same equations, the CLM dump `lwnet_sunwall`/`lwnet_shadewall`/`lwnet_canyon` for the traced urban rows are off by up to ~1.0e-4 W/m² (2.6e-5 W/m² on a 2.49 W/m² net, i.e. ~1e-5 relative); tight-converging the iteration to 1e-12 reproduces the exact values to <1e-14. This is a precision artefact, not a logic error — it is largely below history-file significance — but it must be accounted for when pinning `*.esm` against these dumps.
- **Reproduced/corrected in .esm:** `components/urban_canopy/urban_radiation.esm` (EarthSciModels `clmu-verify` branch) computes the exact closed-form solution (`det_lw`, `det_dir_*` and the `band_*` expression templates are the analytic inverse of `(I - T·A)`), so it is *more* accurate than the dumped values. The CESM-dump tests `cesm_dump_canyon_longwave_hwr1p4`/`_hwr4p8` therefore assert with an absolute tolerance of 0.001 W/m² (the Fortran cutoff) rather than machine precision; the `net_solar` comparison uses the exact system solve.
- **Status:** open (numerical-precision note; no code change proposed — tightening the tolerance would only cost runtime).
- **Note (2026-11-01):** renumbered from #13 to #19 when the CICE lowering bugs (#13-#18) landed from the same numbering space.


## 20. CLM soil-water negative h2osoi_liq is printed but not stopped (endrun commented out)

- **Location:** `components/clm/src/biogeophys/SoilWaterMovementMod.F90:1412-1416`
  (private `soilwater_moisture_form`, adaptive-substep Richards solve; the
  check `if(h2osoi_liq(c,j) < -1e-6_r8)` prints
  `'layer, h2osoi_liq: ', c, j, h2osoi_liq(c,j)` and then the
  `call endrun(subname // ':: negative soil moisture values found!')` is
  **commented out**).
- **Description:** when the Richards substep solve produces negative liquid
  soil moisture, the code writes a diagnostic line and continues with the
  unphysical negative value instead of aborting (or clamping). In the
  instrumented I2000Clm60Sp run (cesm.log.11967609.261101-200125 and
  .261102-004337, lines 14-19) column 182 prints layers 2-5:
  `-20.878`, `-1.957e-2`, `-2.258e-2`, `-7.776e-3` kg/m2 (deterministic,
  4 lines per run). The -20.9 kg/m2 layer-2 value is far below the
  physical range (a 0.02 m soil layer at porosity ~0.4 holds at most
  ~8 kg/m2), so the state is corrupt at that point.
- **Impact on results:** the traced soil-moisture-dependent outputs
  (infiltration, drainage, transpiration stress, soil heat capacity via
  h2osoi in SoilThermProp) inherit the negative water mass for the
  affected column/layers; the effect is confined to the flagged column in
  this run but the code path is reachable in any run with the adaptive
  substep scheme.
- **Reproduced/corrected in .esm:** n/a for stage 2 so far (no soil-water
  .esm component yet); the stage-2 clm soil skeleton should reproduce the
  Richards solve and clamp/flag negatives.
- **Status:** open (recorded; endrun deliberately disabled upstream — the
  surrounding code already moved the over-saturation excess upward just
  before this check, so the negative is a residual from the substep
  solve, not from that redistribution).

## 21. CLM dump: t_soisno value columns shifted 12 layers below the written j label (clm_dump_mod.F90 dump_soiltemp / dump_soilfluxes)

- **Location:** fork instrumentation `components/clm/src/main/clm_dump_mod.F90` — `dump_soiltemp` (loop and writes at :885-909) and `dump_soilfluxes` (`t_soisno(c,1..3)` at :588); artifacts `stage1/dumps/clm/clm_soiltemp.csv`, `clm_soilfluxes.csv`.
- **Description:** `dump_soiltemp` loops over the true layer indices `j = snl(c)+1 .. jbot` (documented in its header comment), but the t_soisno/dz/z/zi value columns at row `j` carry true layer `j−12`: the dummy arguments read an array whose slot s holds true layer s−12, so the top snow layer t_soisno(c,snl+1) appears at row j = 4 and the top soil layer at row j = 13 (verified: row j = 4 equals the preflux t_grnd to round-off on all 192 rows, row j = 13 equals the surfhum instant's t_soisno_top_soil, and rows j = 13..25 reproduce the true soil dz grid 0.02..0.44 m for true layers 1..13). Rows j = 1..3 read 0.0-filled unused slots below the snow layers; rows j ≤ 0 read out-of-bounds/stale memory (garbage dz, stale ~225.7 K temperatures; the dz column at j ≤ 0 shows denormals/NaN/1e3 values). The same shift makes `dump_soilfluxes`' t_soisno_1..3 identically 0.0 on all 192 rows.
- **Impact on results:** any consumer using the raw j labels misattributes every layer; layer-resolved pins (e.g. the surfhum top-layer temperatures) must decode j → j−12 and restrict to the window j = 4..25. The defect affects only the layer index columns of the two dumps; the t_soisno *values* are clean once decoded (all 13,056 soiltemp rows fit the decode exactly).
- **Reproduced/corrected in .esm:** stage2/clm_canopy.esm reconstructs t_soisno_top_snow (entry j = 4) and t_soisno_top_soil (entry j = 13) from `clm_soiltemp.csv`; the three surfhum tests reproduce the dumped qg_snow/dqgdT to round-off (rel 1e-12). Documented in stage2/README-clm.md (Dump defects).
- **Status:** fixed in the instrumentation (fork clm `2d8e0c697`, clm-dumps-w2): root cause was gfortran passing whole arrays to assumed-shape `(:,:)` dummies as 1-based positional views; the 2-D dummies now declare dim2 lower bound `-nlevsno+1`, and the regenerated dumps (LID 261102-004337) use the true model layer index j (file slot j+nlevsno) with valid soilfluxes t_soisno_1..3. The clm_canopy.esm pins were extracted from the pre-fix files; their values are unaffected (the decode was exact), only the j labels cited in descriptions refer to the old layout.
- **Note (2026-11-02):** renumbered from #20 to #21 (and VanGenuchten #21→#22) when clm-dumps-w2's #20 landed in parallel.

## 22. CLM VanGenuchten1980 soil retention option is a verbatim Clapp-Hornberg copy (SoilWaterRetentionCurveVanGenuchten1980Mod.F90)

- **Location:** `components/clm/src/biogeophys/SoilWaterRetentionCurveVanGenuchten1980Mod.F90:45-158` (`soil_hk`, `soil_suction`, `soil_suction_inverse`) vs `SoilWaterRetentionCurveClappHornberg1978Mod.F90:45-159`.
- **Description:** with `use_soil_water_retention_curve = 'VanGenuchten1980'`, the Fortran still evaluates the Clapp-Hornberg equations: the Van Genuchten module is a verbatim copy of the Clapp-Hornberg module — a file diff shows only the module/type names and comments differ — so `smp = -sucsat*s**(-bsw)` (VanGenuchten1980Mod.F90:115) is the Clapp-Hornberg power law, and `soil_hk` (F90:75) is `hk = imped*hksat*s^(2*bsw+3)`. The actual Van Genuchten 1980 curve is never evaluated in this code path.
- **Impact on results:** none on the run (the default ClappHornberg1978 and the VanGenuchten1980 option are numerically identical in this source), but any `.esm` implementing the Van Genuchten curve from the literature would fail dump pins; the template library must transcribe the active (Clapp-Hornberg) form regardless of the namelist option.
- **Reproduced/corrected in .esm:** stage2/clm_templates.esm's `soil_suction` (used by the hr law of clm_canopy.esm) transcribes the Clapp-Hornberg form with the provenance noted; the Van Genuchten-named duplicate templates were removed from the library (see stage2/README-clm.md, Removed templates and the Van Genuchten 1980 note).
- **Status:** open (recorded; transcribe as-run).

## 23. CISM pseudo-plastic basal law (HO_BABC_PSEUDO_PLASTIC_OLD) uses cell-center topg where stagtopg was intended

- **Location:** `components/cism/source_cism/libglissade/glissade_basal_traction.F90:262-277` (`glissade_calcbeta`, PSEUDO_PLASTIC_OLD branch).
- **Description:** the phi(bed) ramp and the resulting tau_c/beta are computed from `bed = topg(ew,ns) - eus` — the **cell-center** bed elevation — while beta lives on the **staggered** grid (vertices). The code's own comment admits it: "There is a minor bug in the loop below. beta(ew,ns) is computed based on topg(ew,ns); should use stagtopg(ew,ns) instead. Leaving the bug as is for back compatibility, given that this method will be deprecated." (The newer BPMP branch at :1191 does use `stagtopg`.) The stage-1 basal_traction trace pins the as-run behavior: `bed = topg - eus` with cell-center topg (validated by cross-referencing the velo trace's topg/eus columns).
- **Impact on results:** tau_c/beta are evaluated at a slightly different bed elevation than the vertex-consistent one; the effect is small for smooth bed topography but systematic. The stage-1 dump reproduces the active code exactly.
- **Reproduced/corrected in .esm:** reproduced as-run in the stage-2 CISM basal law test (bed from cell-center topg − eus).
- **Status:** open (recorded; transcribed as-run).

## 24. CISM cap import turns the mediator's 0 K missing-data sentinel into a +0 °C surface air temperature

- **Location:** `components/cism/drivers/cpl/nuopc/glc_import_export.F90:323-329` (`import_fields`).
- **Description:** after `tsfc = tsfc - tkfrz`, the cap applies the "Jer hack fix": `where (tsfc < -250.d0) tsfc = 0.d0`. CLM sends 0 K (i.e. -273.15 °C after conversion) for land points where CLM sees ocean and for all ocean points, so those cells get a surface air temperature of **0 °C** rather than being masked/excluded. The `-250 °C` threshold is a hard-coded sentinel test.
- **Impact on results:** on the gris20 active grid, boundary cells adjacent to the ocean receive artm = 0 °C instead of a physical value; the tsfc trace pins the post-processing values exactly as consumed by CISM.
- **Reproduced/corrected in .esm:** reproduced as-run in the cap-boundary test (tsfc column of the cism_cap trace).
- **Status:** open (recorded; legacy behavior kept for bit-reproducibility).

## 25. CISM beta_internal unit conversion at the solver boundary is undocumented (Pa/(m/yr) internal vs Pa/(m/s) stored)

- **Location:** `components/cism/source_cism/libglissade/glissade_velo_higher.F90:3989-3990` (`glissade_velo_higher_scale_output`); field docs `libglide/glide_types.F90:1457-1459,2345-2347`.
- **Description:** the HO velocity solver works internally in Pa/(m/yr) (the units in which the sliding laws are evaluated, matching the dumped beta column), but `scale_output` silently converts the stored field with `beta_internal = beta_internal * scyr` to Pa/(m/s). The `glide_types` field comments describe beta as "Pa yr/m" without mentioning the stored-field unit, so the two unit conventions coexist without documentation.
- **Impact on results:** none on the trajectory (self-consistent); matters for anyone reading the stored beta field or writing an `.esm` that must match the solver-internal value (the stage-1 beta trace deliberately reports the pre-scale_output Pa/(m/yr) value).
- **Reproduced/corrected in .esm:** the basal-law pin uses the Pa/(m/yr) solver-internal units (documented in the cism_basal_traction.csv header).
- **Status:** open (documentation/units note; recorded).

## 26. CISM geomderv%stagtopg is vestigial (never staggered in the Glissade path)

- **Location:** `components/cism/source_cism/libglissade/glissade_velo.F90:137,180` (passed to `glissade_velo_driver`), `libglissade/glissade_types.F90` (geomderv).
- **Description:** `model%geomderv%stagtopg` is passed through the Glissade velocity driver but is never computed by any stagger call; the higher-order solve uses its own local `stagtopg` (from `glissade_stagger(topg, stagtopg)` at `glissade_velo_higher.F90:1520`), and the basal-traction routines use their own local copies. The derived-type field is dead weight.
- **Impact on results:** none (dead field; the stage-1 velo trace dumps the solve-local stagtopg, which is the value actually used).
- **Reproduced/corrected in .esm:** n/a (no behavioral effect).
- **Status:** open (dead-code note; recorded).

## 27. CMEPS single-column driver: unchecked netCDF errors (SIGSEGV on missing file) and undetected mesh-match failure (garbage scol_* attributes, run continues)

- **Location:** `components/cmeps/cesm/driver/esm.F90` `esm_get_single_column_attributes` (lines ~1312-1455, cam6_4_212-era CMEPS; called once per component in the cap init for SCAM/single-column mode).
- **Description:** every `nf90_*` call checks its status only by logging (`call shr_log_error(...)`), which does not stop the run; the code then proceeds to use the (uninitialized) return buffers. Two concrete failures seen in this lane:
  (a) **Missing domain file:** when `single_column_lnd_domainfile` does not exist, `nf90_open` fails, `get_var xc` reads into an uninitialized buffer, and the run dies with a SIGSEGV inside ESMF/io (observed for `domain.lnd.ne3np4_gx3v7.230718.nc`; both instrumented and pristine runs crashed identically).
  (b) **Mesh-element search failure:** after the nearest-neighbor search in the domain file, the code searches the ESMF mesh for an element whose center matches the domain point within 1e-4; if no element matches, `scol_mesh_n` stays 0 and `scol_lon = lonMesh(0)` / `scol_lat = latMesh(0)` / `scol_area = dataptr(0)` read out of bounds (bounds-checking off), silently leaving the attributes as garbage (observed: 0.00000, 0.00000). The run then continues with `scol_lon=scol_lat=0`, `scol_lndmask` garbage, `scol_lndfrac=0`, and only the components' own mesh-based grids keep the solution sane.
- **Impact on results:** with a domain file that does not contain the model mesh element centers, the single-column attributes are garbage yet the run proceeds; CAM/CLM/CICE/DOCN still build their grids from the meshes, but the `scol_*` attributes (used e.g. by CLM's column fraction, by the IOP-column selection and by diagnostics) carry wrong values. In this lane it produced a run whose land fraction was 0 (all-ocean single column) while the archive shows the correct ARM97 column coordinates — i.e., the run "worked" only by accident.
- **Reproduced/corrected in .esm:** n/a (run-configuration/driver robustness bug, not a physics bug).
- **Status:** open (upstream candidate). Not patched: the lane works around it with a reconstructed domain file that matches the mesh element centers.
- **Note (2026-11-02):** the same unchecked-error pattern (`shr_log_error` + continue) recurs throughout this subroutine (dim/var inquires, `get_var` for xc/yc/area/mask/frac); each should be `ChkErr`/abort. A second observation from this lane: with a land single column (domain frac=1) in the FSCAMARM97 case, CAM's longwave-down export (`flwds`) collapses to 0.0 from step ~5 onward and CLM's `check_for_errors` (`components/clm/src/cpl/utils/lnd_import_export_utils.F90:107-130`) then aborts the run ("Longwave down sent from the atmosphere model is negative or zero"); with the all-ocean column (frac=0, as the original inputdata file produced) radiation stays healthy and the check never fires. Whether the flwds=0 collapse with a land column is a CAM bug or a SCAM-configuration quirk is unconfirmed; recorded here for the inputdata owner and CAM developers.

## 28. CMEPS: aoflux rainc fallback allocation is uninitialized when add_gusts=.false.

- **Location:** `components/cmeps/mediator/med_phases_aofluxes_mod.F90:1616-1622` (`set_aoflux_in_pointers`).
- **Description:** when `add_gusts` is false the code allocates the `aoflux_in%rainc` pointer with `allocate(aoflux_in%rainc(1))` and never initializes it, so `aoflux_in%rainc(1)` is an uninitialized allocation. In the pristine code this is benign (rainc is never read on the no-gusts path — `flux_atmOcn_large` only touches it inside `if (add_gusts)`), but the stage-1 dump traces rainc as a boundary input and read garbage: values seen across runs include a subnormal 1.27887920104815652e-315 and a 2.5e180 (both unparseable with gfortran es26.17, see #29).
- **Impact on results:** none on the coupled model (rainc unused without gusts); affects only the instrumentation dump's rainc column.
- **Reproduced/corrected in .esm:** n/a (fork-only). Patched in the fork (med_phases_aofluxes_mod.F90, commit 4750efc5): the fallback allocation is initialized to 0.0, so the dump writes rainc = 0 for the no-gusts C_JRA_t62 run.
- **Status:** patched locally (fork); benign in pristine code.

## 29. gfortran es26.17 omits the exponent letter for subnormals and 3-digit exponents

- **Location:** toolchain — gfortran (GNU Fortran) list-directed-free formatted output; triggered by the stage-1 CMEPS dump (`med_dump_mod.F90`).
- **Description:** gfortran's `es26.17` writes subnormal values (|x| < 2.2e-308) as e.g. `1.27887920104815652-315` and values with 3-digit exponents (|x| >= 1e100) as e.g. `2.51968629031625108+180`, both without the `E` — unparseable by standard CSV/float readers. `es26.17e3` always prints the exponent letter (reproduced with a 12-line standalone test).
- **Impact on results:** none on the model; breaks machine-readable dumps when such values occur (rainc garbage per #28 triggered the large-exponent case).
- **Reproduced/corrected in .esm:** n/a (toolchain). Fork patched (med_dump_mod.F90, commit 4750efc5): dump writes `es26.17e3` and flushes subnormals to 0.0.
- **Status:** patched locally (fork); open upstream.

## 30. CLM dump: clm_canopyflux.csv ent/hyd stages trace bare soil patches, so the vegetated interception flux split is not dump-covered

- **Location:** fork instrumentation `components/clm/src/main/clm_dump_mod.F90`
  `dump_canopyflux` — called with `num_soilp, filter_soilp` from
  CanopyHydrologyMod.F90:228,459 (ent/hyd stages) and with
  `num_exposedvegp, filter_exposedvegp` from CanopyFluxesMod.F90:1764 (flx
  stage); artifact `stage1/dumps/clm/clm_canopyflux.csv` (whose description
  reads "One row per traced vegetated patch per step").
- **Description:** the ent/hyd stages dump the first `dump_first_cols`
  patches of the soil-patch filter — in the traced run the 16 bare patches
  1-16 (frac_veg_nosno = 0, check_point = 0 everywhere) — while only the flx
  stage (CanopyFluxes exit, end of step) traces the vegetated patches 48-65
  (frac_veg_nosno = 1). The header's "vegetated patch" wording applies only
  to the flx rows. Because the interception/throughfall fields
  (fpiliq/fpisnow, qflx_intercepted_*, qflx_through_*) are written only at
  the hyd stage, they are NaN on every flx row and never dumped for a
  vegetated patch: no row anywhere carries both a vegetated patch and the
  interception fluxes. The flx-stage fwet/fdry/fcansno are also stale
  relative to the flx liqcan: they were computed by this step's hydrology
  from the pre-evap storage, while the flx liqcan is the post-evap state
  (CanopyFluxesMod.F90:1618-1631) — on the 10 flx rows with fwet < 0.05 the
  post-evap liqcan is 0, so those dumped fractions cannot be reproduced from
  any dumped state.
- **Impact on results:** none on the model (dump-only). The stage-2 pins
  cover the vegetated fraction laws only where the fwet clamp saturates
  (182/192 flx rows at exactly 0.05) and cannot pin the interception laws,
  the excess-fall fluxes, the snow-unloading branch, or the storage ODE with
  nonzero interception; the 10 fwet < 0.05 rows are unpinnable (pre-evap
  storage not dumped).
- **Reproduced/corrected in .esm:** stage2/clm_canopy.esm
  `clm_canopy_fractions_veg_step*` tests pin the fraction laws at the dumped
  post-evap state via initial_conditions (clamp-saturated, so the pin is
  robust to the pre/post-evap ambiguity); the interception/throughfall flux
  split and the excess-fall fluxes are documented as unpinned in
  stage2/README-clm.md (Coverage gaps, Dump extensions). Closing the gap
  needs the hyd-stage filter to include the vegetated patches (or the
  interception fields to be written at the flx stage).
- **Status:** open (instrumentation coverage note; recorded).

## 31. CLM fork dump: derived-type NaN passed directly in a formatted write list crashes the instrumented run

- **Location:** fork instrumentation `components/clm/src/main/clm_dump_mod.F90`
  `dump_surfalbedo` (urban branch, fork clm `3d7038304`): the write
  statement passed the module's `nan` (the renamed derived-type parameter
  `shr_infnan_nan` of `shr_infnan_mod`, a `type(shr_infnan_nan_type)` whose
  single component is LOGICAL) directly in an `es26.17` data list.
- **Description:** `shr_infnan_nan` is a derived-type parameter, not a
  real: it can only be converted to a real through the module's
  `assignment(=)` interface (the module header says "cannot be passed to
  functions that expect real arguments"). gfortran's formatted transfer
  therefore failed at runtime with "Expected REAL for item 14 in formatted
  transfer, got LOGICAL" (item 14 = the first albsnd/albsni slot) the first
  time the urban branch ran with the extended dump (LID
  12027802.261103-131613, step 2), aborting the whole run: the run's CSVs
  were truncated mid-step-2. (The earlier runs never hit it because they
  predated the extension's albsnd/albsni columns.)
- **Impact on results:** the aborted run wrote truncated dumps (the
  complete pre-extension dumps were preserved). Fix verified by a rerun
  (LID 12027802.261103-133002): with the NaN materialized to a real local
  (`nanval = nan` via the assignment interface) all 11 dumps complete and
  the urban albsnd/albsni slots read `NaN` as designed.
- **Reproduced/corrected in .esm:** n/a (instrumentation defect, dump-only;
  no model state affected).
- **Status:** patched locally (uncommitted in the fork working tree; to be
  committed with the surfalbedo dump work).

## 32. GLC cap: hflx export sign convention contradicts the code comments and field dictionary

- **Location:** `components/cism/source_cism/libglad/glad_output_fluxes.F90:134-139`
  (`accumulate_output_fluxes`, `hflx_sum = hflx_sum - ucondflx`) vs
  `components/cism/source_cism/libglad/glad_main.F90:625` (the `glad_gcm`
  `hflx` intent comment "output heat flux (W/m^2, positive down)") and the
  CMEPS field dictionary `components/cmeps/mediator/esmFldsExchange_cesm_mod.F90:816`
  ("Flgg_hflx: downward heat flux from glacier interior").
- **Description:** `ucondflx` is the upper-surface conductive flux, positive
  down (`libglissade/glissade_therm.F90:669`); the accumulation code flips the
  sign with the comment "Flip the sign so that hflx is positive up"
  (F90:134-136), so the exported `hflx` is **positive up**. The `glad_gcm`
  interface comment (glad_main.F90:625) and the CMEPS field dictionary both
  label `Flgg_hflx` as **positive down**. The two labels cannot both describe
  the exported value. (The value is currently unconsumed: the mediator errors
  out if `Flgg_hflx` is actually routed to the land model,
  `med_phases_post_glc_mod.F90:343-344`, so no downstream component pins the
  convention.)
- **Impact on results:** none in the reference run — every traced export row
  of `cism_cap.csv` has hflx = 0 (ucondflx = 0 in the 2-year cold-base
  window), so the sign is unobservable in the dumps. Would matter for a run
  with surface temperature evolution.
- **Reproduced/corrected in .esm:** `stage2/cism_cap.esm` (GlcCapExport)
  transcribes the code as written (`hflx = -ucondflx`, i.e. positive up per
  the F90:135-136 comment) and documents the label contradiction; the sign is
  unpinned by the dump.
- **Status:** open (comment/label inconsistency; needs an upstream decision
  on the intended convention).

## 33. CLM Sturm1997 snow thermal conductivity is evaluated with single-precision literals (missing `_r8` suffix)

- **Location:** `components/clm/src/biogeophys/SoilTemperatureMod.F90:774-778` (Sturm1997 branch of `SoilThermProp`, non-glacier land units; the identical glc branch at :754-758 and the same copy in `LakeFluxesMod.F90:497` / `LakeTemperatureMod.F90:1195` carry the same defect).
- **Description:** the five Sturm et al. (1997) fit literals `0.023`, `0.234`, `0.138`, `1.01`, `3.233` have no kind suffix. Under default `real(4)` (gfortran, CESM default) they are stored as single precision and promoted to `real(8)` at evaluation, so the executed coefficients are the float32-rounded values — 0.1379999965429306, 1.0099999904632568, 3.2330000400543213, 0.023000000044703484, 0.23399999737739563 — not the paper's coefficients. Verified in the compiled object (`SoilTemperatureMod.o` `.rodata.cst8`: the constants are the float32 values; the branch threshold 156 is an exact integer double): a double-precision evaluation of the written formula disagrees with the executed snow interface conductivities at rel ~2.6e-8 (layer thk rel ~1.2e-7 at the traced densities ~340-460 kg/m3), while the float32-constant formula reproduces the dump to rel 1.6e-16 on all 144 traced snow interfaces (16 columns x 9). The neighbouring Jordan1991 branch uses suffixed literals (`7.75e-5_r8`, `1.105e-6_r8`) — the Sturm literals simply lost their suffix when the Jordan formula was replaced (commit ab94ab484, "Replace thermal conductivity of snow of Jordan with Sturm").
- **Impact on results:** snow layer thermal conductivity is ~1e-7 relative above the intended fit everywhere, silently changing the effective Sturm coefficients; tiny for climate but a real precision defect, and it makes the "as-written" formula fail exact pins.
- **Reproduced/corrected in .esm:** reproduced exactly — `stage2/clm_templates.esm` `snow_thermal_cond_sturm1997` carries the float32-rounded coefficients with this provenance, and the `clm_soiltemp.esm` snow interface-tk pins pass at rel 1e-12 to the dump's full double precision.
- **Status:** open (upstream candidate; the fix is appending `_r8` to the five literals).

## 34. CLM snow hydrology leaves the top snow layer's (z, dz) inconsistent with its zi interfaces

- **Location:** SUSPECTED attribution — `components/clm/src/biogeophys/SnowHydrologyMod.F90` (`DivideSnowLayers` / `AddNewSnow`; called from `SnowWater` before `SoilTemperature`; clm_driver.F90:921/:971 order rules out in-step snowfall as the cause). Not traced to the specific call that leaves the state inconsistent; the observed defect is in the state arrays `col%z`, `col%dz`, `col%zi` (SoilStateType) at the SoilTemperature entry.
- **Description:** after snow-layer division/addition the layer centers and thicknesses (z, dz) of the top snow layer are updated but the interface depths (zi) are not re-derived from them, so `z(c,snl) - 0.5*dz(c,snl) != zi(c,snl)`. In the W3 traced January columns the top-layer mismatch is ~3.02e-5 m (column 1: `zi(-9) = -22.13156338` vs `z(-8) - dz(-8)/2 = -22.13159361`; the layer thickness from dz is 0.01557835 m vs the interface spacing zi(-8) - zi(-9) = 0.01551790 m, i.e. a 6.05e-5 m discrepancy) and up to ~1.6e-4 m (columns 11-15); all other layers satisfy the identity to machine precision.
- **Impact on results:** small but real: the top snow row's `fact` (ComputeHeatDiffFluxAndFactor, SoilTemperatureMod.F90:1896) uses `z(c,j) - zi(c,j-1)` directly, so the capr-weighted surface factor carries the inconsistency; it does not cancel elsewhere because the matrix/rhs rows use zi only through the top-row factor. The W3 dump now writes `zi_top`/`zi_bot` directly, so the harmonic-mean interface pins feed the true interface depths and reproduce the dump to round-off (rel 1e-12) despite the inconsistency.
- **Reproduced/corrected in .esm:** the harmonic-mean pins in `stage2/clm_soiltemp.esm` (model `SoilTemperature`) use the dumped `zi_bot` (true interface below each layer), which absorbs the inconsistency exactly; the top-snow factor is transcribed verbatim for the PDE step (the local capr-weighted top-row rule `stage2/clm_soiltemp_discretization.esm`). The inconsistency itself is not corrected (out of scope; it reproduces the F90 state).
- **Status:** open (upstream candidate; recorded).

## 35. CMEPS ocnalb init: else branch sets the local use_min_ocnalb instead of the module use_min_albedo

- **Location:** `components/cmeps/mediator/med_phases_ocnalb_mod.F90:220-229`
  (`med_phases_ocnalb_init`, the `ocean_albedo_limit` attribute read).
  Line numbers below refer to the pristine `stage1/cesm` tree; the
  instrumented fork `stage1/EqESM-cesm-cmeps` is offset by +2 (and by +24
  in the albedo loop).
- **Description:** when the `ocean_albedo_limit` attribute is absent, the
  else branch executes `min_albedo = 0.0_R8; use_min_ocnalb = .false.`
  (F90:227-228). `use_min_ocnalb` is a LOCAL of the init subroutine
  (declared F90:107) that is never read anywhere; the albedo loops
  (F90:250 in the init pass and F90:483) gate the `max(anidr, min_albedo)` clamp on the MODULE
  variable `use_min_albedo` (declared F90:64). The intent of the dead
  assignment was clearly to switch the clamp off in this branch. The bug is
  latent but undefined-behaviour: the module variable has no initializer,
  so in this branch its value is formally undefined (zero, i.e. `.false.`,
  in practice with the CESM compilers), and the else branch behaves as
  intended only by accident; the `isPresent .and. isSet` branch sets
  `use_min_albedo = .true.` correctly (F90:225).
- **Impact on results:** none in the traced run (`ocean_albedo_limit` not
  set, so `use_min_albedo` stays `.false.` either way; the min_albedo clamp
  is inactive). Would matter on a compiler/flags combination that does not
  zero-initialize module variables.
- **Reproduced/corrected in .esm:** `stage2/cmeps_ocnalb.esm` transcribes
  the clamp as `ifelse(use_min_albedo, max(anidr, min_albedo), anidr)`
  driven by the dumped `use_min_albedo` flag, i.e. it pins the intended
  semantics, not the dead assignment; the traced rows all have
  `use_min_albedo = 0` so the clamp branch is unpinned.
- **Status:** open (dead assignment; benign in all CESM configurations on
  this branch).

## 36. DATM JRA55 datamode: unsuffixed real literals 0.608/0.069/0.011 silently round to single precision

- **Location:** `components/cdeps/datm/datm_datamode_jra_mod.F90:292` (qcoef 0.608), `:324` (LY08 albedo 0.069, 0.011); every other literal in the datamode is `_R8`-suffixed.
- **Description:** the three coefficients are written as unsuffixed Fortran literals (default real kind = single precision), so the compiled laws use the r4-rounded values r4(0.608) = 0.607999980449676513671875, r4(0.069) = 0.068999998271465301513671875, r4(0.011) = 0.010999999940395355224609375 rather than the intended decimal values. Found 2026-11-05 while matching the C_JRA_t62_datm_inst stage-1 dumps: the Sa_dens and avg_alb values are ~2e-8 relative off from the exact-r8 laws (density 1.6e-11, albedo 2.2e-8). qcoef = 0.608 is itself an intentional rounding of eps1/eps = 0.607717... (documented in the code comment), but the r4 rounding is an additional 8th-decimal effect that survives into the forcing.
- **Impact on results:** none practically (forcing-level rounding), but the laws differ subtly from the CLMNCEP/GSWP3 equivalents and from any naive transcription, so tests must use the r4-rounded constants.
- **Reproduced/corrected in .esm:** reproduced exactly by `density_jra` and `jra_albedo` in `stage2/cdeps_templates.esm` (constants carry the r4-rounded values with provenance notes); `stage2/cdeps_datm_jra.esm` tests pin the dumps to 1e-12 relative.
- **Status:** open (cosmetic; suggest suffixing `_R8` upstream).

## 37. DATM JRA55 datamode: dead code cosfactor/phs_c0/dLWarc in advance

- **Location:** `components/cdeps/datm/datm_datamode_jra_mod.F90:64-65,274`.
- **Description:** `cosfactor = cos((2π·rday)/365 − phs_c0)` is computed every advance (F90:274) but never used; the constants `phs_c0 = 0.298` and `dLWarc = −5.0` are declared but never referenced. The rday computation at F90:271-273 is therefore also dead weight (rday is used only by cosfactor).
- **Impact on results:** none (dead code).
- **Reproduced/corrected in .esm:** not reproduced (omitted from the Jra model in `stage2/cdeps_datm_jra.esm`).
- **Status:** open (cleanup; suggest removing upstream).

## 38. ZM: zm_conv_evap_run ignores ke_lnd and landfrac (kemask = ke unconditionally)

- **Location:** `components/cam/src/atmos_phys/schemes/zhang_mcfarlane/zm_conv_evap.F90:160` (`kemask = ke`); `ke_lnd` and `landfrac` are declared `intent(in)` (F90:23, 25, 50, 55) but never referenced in the body. The caller (`physics/cam/zm_conv_intr.F90:671-679`) passes `zmconv_ke` and `zmconv_ke_lnd` and the land fraction, so the namelist option `zmconv_ke_lnd` is wired through but has no effect.
- **Description:** upstream CAM6.4's ZM evaporation applies one evaporation coefficient to all columns (`kemask = ke`), despite the code carrying the `ke_lnd`/`landfrac` arguments. (The CESM2-era intent was presumably `kemask = ke_lnd` over land, the analog of the cldprp `c0_lnd`/`c0_ocn` split, but the branch was never completed.)
- **Impact on results:** none in the reference run — all dumped columns of the SCAM window are ocean g = 1 columns, so the distinction is unobservable in the dumps; the esm model transcribes the code as written (one coefficient `ke`, dumped 5e-6).
- **Reproduced/corrected in .esm:** `stage2/cam_zm_evap.esm` (ZMevap) reproduces the code as written; the `ke` variable documents the dead arguments. Verified identical in the vanilla CESM checkout (fork-vs-vanilla diff empty — not an instrumentation artifact).
- **Status:** open (dead arguments; the zmconv_ke_lnd namelist option has no effect on the simulation).

## 39. CLM CanopyFluxes below-freezing storage transfer lacks the max(0,.) clamp of its above-freezing mirror

- **Location:** `components/clm/src/biogeophys/CanopyFluxesMod.F90:1630` (the
  all-snow-evaporated branch of the end-of-CanopyFluxes storage update;
  above-freezing mirror with the clamp at `:1622`).
- **Description:** in the above-freezing branch, when the evaporation demand
  exceeds liqcan the remainder is taken from snocan with a positivity clamp —
  `snocan = max(0, snocan + liqcan + (tran-evap)*dtime)` (`:1622`) — and
  liqcan itself is clamped by `max(0, ...)` (`:1623`). The below-freezing
  mirror (`:1628-1634`) instead transfers to liqcan **without any clamp**:
  `liqcan = liqcan + snocan + (tran-evap)*dtime` (`:1630`), and only the final
  snocan line is clamped (`:1634`). Since the branch condition only
  guarantees `(tran-evap)*dtime < -snocan`, the transfer is negative whenever
  the evaporation demand exceeds the total canopy water liqcan + snocan, and
  the negative liqcan persists (the `truncate_small_values` call at
  `:1640-1643` touches only snocan). The stage-1 validator
  (`stage1/analysis/verify_clm_veg_dumps.py` `check_canopyflux`) silently
  assumes the clamp (`max(0, liqcan_h + snocan_h + net)`), i.e. it encodes
  the intended semantics, not the Fortran as written.
- **Impact on results:** none on the traced data — all 176 vegetated rows
  have t_veg > tfrz, so the below-freezing branch never runs; potentially
  negative liqcan on below-freezing patches with sublimation demand
  exceeding the canopy storage.
- **Reproduced/corrected in .esm:** `stage2/clm_canopy_storage.esm`
  transcribes the Fortran as written (the below-freezing transfer at
  `:1630` without the clamp); the branch is unpinned (see
  `stage2/README-clm-canopy-veg.md`).
- **Status:** open (upstream asymmetry; unexercised by this lane's dumps).
