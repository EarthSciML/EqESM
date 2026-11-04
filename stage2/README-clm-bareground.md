# CLM stage 2 — bare-ground turbulent fluxes

Stage-2 implementation of the CLM bare-ground surface-layer fluxes
(ctsm5.4.063, cesm3.0-alphabranch; checkout `stage1/cesm/components/clm`) as
esm 1.3.0 files: the pointwise laws of `BareGroundFluxesMod.F90` and
`FrictionVelocityMod.F90` (MoninObukIni, FrictionVelocity, StabilityFunc1/2)
for the non-exposed-vegetation patches, with the shared pointwise laws
factored into the expression-template library `clm_templates.esm`. The
`tests` blocks pin the stage-1 instrumentation dump
`stage1/dumps/clm/clm_bareground.csv` (192 rows: 12 steps × 16 bare patches,
dtime = 1800 s) with its feeding dumps `clm_run_imports.csv` (surfhum:
qg_snow/qg_soil/dqgdT; preflux: thv) and `clm_soiltemp.csv` (entry stage,
j = −8 snow-top / j = 1 soil-top temperatures) — see
[Testing](#testing).

Sources: `stage1/cesm/components/clm/src/biogeophys/BareGroundFluxesMod.F90`,
`FrictionVelocityMod.F90`, `clm_varcon.F90` (constants), the parameter file
`ctsm60_params.c260518.nc` (wind_min, zlnd), the traced case namelist
`stage1/cime/cases/I2000Clm60Sp_clm_inst/run/lnd_in`, and the QSat laws
already in `clm_templates.esm` (QSatMod.F90).

## Traced configuration

The run is the single-clump column experiment of README-clm.md (16 patches on
16 columns, bare, `frac_veg_nosno = 0`, deep snowpacks). The flux-relevant
namelist/parameter choices transcribed:

| Setting | Value | Source |
|---|---|---|
| `z0param_method` | `'Meier2022'` | lnd_in (note: README-clm.md's `'ZengWang2007'` assumption applies to the vegetation roughness of that component, not to the ground roughness here) |
| `use_z0m_snowmelt` | `.true.` | lnd_in — z0mg comes from the snowmelt branch of `SetRoughnessLengthsAndForcHeightsNonLake` (FrictionVelocityMod.F90:620-623) |
| `zetamaxstable` | `2.0` | lnd_in `friction_velocity` — the stable-branch clamp of the iteration (F90:370) and of MoninObukIni (F90:1201) |
| `soil_resis_method` | `1` (sl14) | lnd_in — the Swenson & Lawrence 2014 branch of raiw (F90:434-435), inactive on the pinned rows |
| `do_soilevap_beta` | inactive | the Lee & Pielke 1992 branch of raiw (F90:427-433) is not reached |
| forcing heights | 30 m (GSWP3) for u/t/q | `forc_hgt_u = forc_hgt_t = forc_hgt_q = 30` |
| `wind_min` | `1.0` m/s | param file |
| `zlnd` | `8.5e-4` m | param file (bare-soil z0mg, untested branch) |
| Meier2022 constants | `beta_param = 7.2`, `nu_param = 1.5e-5` m²/s, `b1_param = 1.4`, `b4_param = −0.31`, `meier_param1 = 0.23`, `meier_param2 = 0.08`, `meier_param3 = 70` | clm_varcon.F90:164-171 |
| `zetam` / `zetat` | `1.574` / `0.465` | FrictionVelocityMod.F90:795-796 |
| `zii` / `beta` | 1000 m / 1 | BiogeophysPreFluxCalcsMod.F90:362-363 |
| `niters` | 3 | BareGroundFluxesMod.F90:342 (the do-while is unrolled as three mounts) |

## File map

| File | Contents | Fortran source |
|---|---|---|
| `clm_templates.esm` | 10 new templates added to the library (25 total): `stability_func_momentum` (StabilityFunc1, exact SHR_CONST_PI·0.5 constant), `stability_func_scalar` (StabilityFunc2, **reused by reference** from the cmeps `psixhu_unstable`, esm-spec §9.7.3), `monin_obukhov_initial_wind`, `bulk_richardson`, `monin_obukhov_initial_zeta` (the [0.01, zetamaxstable] / [−100, −0.01] clamps), `friction_velocity_momentum` (four-branch ustar), `friction_velocity_scalar` (four-branch temp1, shared by the humidity and 2-m profiles), `z0hg_meier2022`, `z0mg_meier2022_snowmelt` (the log10/atan snowmelt branch), `dewpoint_magnus` (Lawrence 2005 eq. 7) | `FrictionVelocityMod.F90:1131-1134,1151-1152,1189-1205,857-871,956-970,620-623`; `BareGroundFluxesMod.F90:365,548-578` |
| `clm_bareground_step.esm` | Model `BareGroundStep` — **one** pass of the stability do-while, authored once and mounted 3× (esm-spec §4.7): the three profile heights zldis_u/zldis_t/zldis_2m, ustar/temp1/temp12m from the four-branch profile templates, temp2 = temp1 and temp22m = temp12m (equality branches), tstar/qstar, the Meier2022 z0h, thvstar, zeta_raw, the sign-preserving stability clamp, the wc/um branch update and obu; 17 equations | `BareGroundFluxesMod.F90:342-388` |
| `clm_bareground.esm` | Model `BareGroundFluxes` — the derived inputs (ur, dth, dqh, dthv, z0mg, zldis), the MoninObukIni initialization (um_ini, rib, zeta_ini, obu_ini), the exit state (= step3 outputs), and the flux laws: ram/rah/raw, raih, forc_esat/forc_e/forc_dewpoint, the raiw selection, ram1, cgrnds/cgrndl/cgrnd, taux/tauy, eflx_sh_grnd/soil/snow/h2osfc, qflx_evap_soi/ev_snow/ev_soil/ev_h2osfc, t_ref2m/q_ref2m/rh_ref2m, kbm1; 45 equations + top-level mounts step1..step3 and 72 `variable_map` couplings (20 update-invariant fields × 3 steps feed every mount; um/obu/z0h/wc thread the state with 12 chained edges) + 5 tests | `BareGroundFluxesMod.F90:325-329,331,342-388,396-499` |

The parent imports `./clm_templates.esm` and `./cmeps_templates.esm` by
reference; no law is duplicated. The 3-mount composition mirrors the
cmeps_aoflux precedent (README-cmeps.md): top-level mounts are required
because the pinned CLI drops couplings in the single-top-level-model build
path.

## Units

| Quantity | Symbols | Units |
|---|---|---|
| Roughness lengths | `z0mg`, `z0h`, `zlnd`, `zldis`, `zldis_u/t/2m`, `obu*` | m |
| Wind speeds | `ur`, `um*`, `wc*`, `ustar`, `forc_u/v` | m/s |
| Temperatures | `forc_t`, `forc_th`, `thm`, `t_grnd`, `t_soisno_*`, `t_h2osfc`, `thv`, `dth`, `dthv`, `tstar`, `thvstar`, `t_ref2m`, `forc_dewpoint`, `tfrz` | K |
| Humidities | `forc_q`, `qg*`, `dqh`, `qstar`, `q_ref2m`, `qsat_ref2m` | kg/kg |
| Pressures | `forc_pbot`, `forc_esat`, `forc_e`, `e_ref2m` | Pa |
| Resistances | `ram1`, `ram`, `rah`, `raw`, `soilresis` | s/m |
| Conductances | `raih`, `cgrnds`, `cgrnd` | W/m²/K |
| Latent derivative | `cgrndl` | kg/m²/s/K |
| Stresses | `taux`, `tauy` | kg/m/s² |
| Heat fluxes | `eflx_sh_grnd/soil/snow/h2osfc` | W/m² |
| Evaporation | `qflx_evap_soi`, `qflx_ev_soil/snow/h2osfc`, `raiw` | kg/m²/s |
| Dimensionless | `zeta*`, `rib`, `zetam`, `zetat`, `zetamaxstable`, `vkc`, `beta*`, `temp1/2`, `temp12m/22m`, `kbm1`, `rh_ref2m`, `frac_sno_albedo`, `snomelt_accum` (m) | 1 (%, for rh_ref2m) |
| Meier2022 constants | `beta_param` | s^0.5/m^0.5/K^0.25 |
| | `nu_param` | m²/s |
| | `meier_param1/2/3`, `b1_param`, `b4_param` | 1 |

The fractional-exponent units of `beta_param` are carried through the unit
checker (verified with a probe). One unit-algebra note: the convective
velocity law `wc = beta·(−grav·ustar·thvstar·zii/thv)^0.333` (F90:372) is
written with the exponent 1/3 in the esm because the unit checker only
reduces exact rational exponents on dimensional bases; 0.333 differs from 1/3
by ~1e-4 relative, but wc is only read on the unstable branch (never pinned —
all traced rows are stable), so no pinned value is affected.

## Testing

The tests block of `clm_bareground.esm` pins 5 rows of the 192-row dump, one
per traced configuration regime; each test asserts all 14 dumped output
columns (ram1, cgrnds, cgrndl, taux, tauy, eflx_sh_grnd/soil/snow/h2osfc,
qflx_evap_soi, qflx_ev_soil, qflx_ev_snow, t_ref2m, q_ref2m) plus z0mg
(dumped input recomputed from the snowmelt template), 15 assertions per row,
75 total, tolerance rel 1e-12 at time_span {0,0}. Inputs per row: the
clm_bareground.csv row values (forc_t/forc_q/forc_u/forc_v/forc_pbot/
forc_rho/t_grnd/qg/thm/beta), thv from clm_run_imports.csv (preflux) with
forc_th = thv/(1 + 0.61·forc_q) (forc_th is not dumped; the inversion is
exact for the dumped dthv), qg_snow/qg_soil/dqgdT from clm_run_imports.csv
(surfhum), t_soisno_top = clm_soiltemp.csv j = −8 and t_soisno_soil = j = 1,
t_h2osfc = t_soisno_soil (equality verified on all 192 rows), and the fixed
configuration (snomelt_accum = 0, frac_sno_albedo = 1, zetamaxstable = 2.0,
soil_resis_method = 1, do_soilevap_beta = .false., GSWP3 heights, wind_min).

| Test id | Dump row | Regime |
|---|---|---|
| `row_1_1` | step 1, p 1, c 1 | MoninObukIni initial zeta **unclamped** (the only 5 dump rows with zeta0 < zetamaxstable; all other rows clamp at 2.0) |
| `row_2_1` | step 2, p 1, c 1 | steady stable regime |
| `row_6_8` | step 6, p 8, c 8 | mid-run, warm column |
| `row_9_4` | step 9, p 4, c 4 | mid-run, low wind (ur clamped at wind_min on forc_u = forc_v = 1.26 m/s) |
| `row_12_16` | step 12, p 16, c 16 | last traced step, coldest column |

Branch census of the full 192-row dump (verified by the Python replica
`scratch/bareground_replica.py`, which reproduces all 192 rows to rel 0.0):
every row is stable (dthv ≥ 0, so um = max(ur, 0.1) and the zeta clamp at
2.0), every row is in the dew branch (dqh > 0 with t_grnd < forc_dewpoint, so
raiw = forc_rho/raw), the iteration momentum/heat profiles are in the
zeta > 1 branch and the 2-m profiles in the zeta ≤ 1 branch, z0mg is in the
snomelt_accum < 1e-5 branch, z0hg is the Meier2022 branch, the dewpoint is in
the liquid branch and the temp2 = temp1 / temp22m = temp12m equality branches
hold. A stage-1 Python replica of the full component (kept untracked in
scratch/) reproduces all 192 dumped rows exactly and was used to verify the
transcription before writing the .esm.

## Coverage gaps

- **Unstable branches** (ζ < 0): the MoninObukIni unstable clamp, the wc/um
  unstable update (and with it the 1/3-vs-0.333 exponent note above), the
  StabilityFunc1/2 unstable forms (via the profiles' ζ < 0 and ζ < −ζₘ/ζₜ
  branches with their 1.14/0.8 free-convection terms) — no traced row is
  unstable, so none of these are pinned.
- **Dewpoint value**: only the branch selection (the `t_grnd > forc_dewpoint`
  gate, which pins the sign of t_grnd − forc_dewpoint, margin ≥ 1.41 K) is
  pinned; the ice branch (t < tfrz selects the ice coefficients) is untested.
- **raiw dqh ≤ 0 side**: the sl14 (`forc_rho/(raw + soilresis)`) and beta
  (`soilbeta·forc_rho/raw`) sub-branches are transcribed but untested
  (traced rows are all dew rows); soilresis/soilbeta are bound only via the
  dump's beta column.
- **z0mg snowmelt branch** (snomelt_accum ≥ 1e-5, the log10/atan form) and
  the bare-soil `zlnd` branch of the Meier2022 select case: untested.
- **ZengWang2007 z0hg branch** (BareGroundFluxesMod.F90:364): not transcribed
  — the traced run is Meier2022.
- **Unpinned outputs**: qflx_ev_h2osfc (qg_h2osfc is not dumped; bound to qg
  in the tests), cgrnd (htvp is not dumped; bound to hvap), rh_ref2m, kbm1 —
  all transcribed with the Fortran-exact formulas.
- **Not transcribed** (dead or undumped in this configuration): the
  human-stress indices (F90:502-528), the www soil-wetness variable (dead —
  only used by the ZengWang2007 a_coef z0hg branch), the rssun/rssha
  photosynthesis resistances, the u10 diagnostics of FrictionVelocity, the
  `use_lch4` grnd_ch4_cond output, and the qflx_evap_tot/qflx_tran_veg/
  qflx_evap_veg sums (zero for bare ground).
- **Documented simplifications**: forc_th recovered from thv (exact);
  t_h2osfc = t_soisno(c,1) (verified on all rows); temp2 = temp1 and
  temp22m = temp12m by the equality branches (hold because
  forc_hgt_q_patch == forc_hgt_t_patch and z0qg == z0hg in the traced run).

## Suite

Baseline (before this component): 1171P/0F/0E. After:
`SUITE: 1246P/0F/0E  validate-failures=0  units-issues=0` — the 75 new
assertions are the bare-ground pins; nothing else changed. No CESM defects
found during this transcription (no bugs.md entries).
