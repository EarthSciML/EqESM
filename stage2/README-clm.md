# CLM stage 2 — `.esm` skeleton (canopy hydrology + flux interface)

Stage-2 implementation of the CTSM/CLM land component (ctsm5.4.063,
cesm3.0-alphabranch; checkout `stage1/cesm/components/clm`) as esm 1.3.0
files: the canopy interception/throughfall hydrology of CanopyHydrologyMod.F90
and the pre-flux surface interface of BiogeophysPreFluxCalcsMod.F90 plus
SurfaceHumidityMod.F90, with the shared pointwise laws factored into an
expression-template library. The `tests` blocks pin the stage-1
instrumentation dumps (`stage1/dumps/clm/` — the clm-dumps-w2 regeneration,
fork clm `2d8e0c697`/LID 261102-004337, with layer-exact
`clm_soiltemp.csv` j labels and real `clm_soilfluxes.csv` t_soisno_1..3;
see [Dump defects](#dump-defects)) — see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/clm/src/biogeophys/`
(`CanopyHydrologyMod.F90`, `BiogeophysPreFluxCalcsMod.F90`,
`SurfaceHumidityMod.F90`, `QSatMod.F90`, `FrictionVelocityMod.F90`,
`SoilWaterRetentionCurveClappHornberg1978Mod.F90`), the parameter file
`ctsm60_params.c260518.nc` and `clm_varcon.F90`/`shr_const_mod.F90` for the
constants, and the instrumentation in
`stage1/EqESM-cesm-ctsm/components/clm/src/main/clm_dump_mod.F90`.

Default configuration transcribed (the D-class reference run's namelist):
`use_clm5_fpi = .false.` (the fpiliq = 0.25·(1−e^(−0.5·(elai+esai))) law,
CanopyHydrologyMod.F90:538), `z0param_method = 'ZengWang2007'`
(z0m = z0mr·htop, displa = displar·htop, BiogeophysPreFluxCalcsMod.F90:175-176),
`use_soil_water_retention_curve = 'ClappHornberg1978'` (see the
[Van Genuchten note](#van-genuchten-1980-note)), the vegetated-landunit
branch of CalculateSurfaceHumidity (SurfaceHumidityMod.F90:41), and the
non-lake, non-urban patch filter of the traced run (see
[Coverage gaps](#coverage-gaps)). The traced run is a single-clump column
experiment: 16 patches on 16 columns (mostly bare, `frac_veg_nosno = 0`,
`elai = esai = htop = 0`, deep snowpacks `snl = −9` on the traced columns),
12 steps of dtime = 1800 s. The canopyflux dump's ent/hyd stages trace the
bare patches 1-16; its flx stage (CanopyFluxes exit, end of step) traces the
vegetated patches 48-65 (`frac_veg_nosno = 1`, nonzero elai/esai,
forc_snow = 0, liqcan > 0 post-evap).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `clm_templates.esm` | Expression-template library (15 templates): the QSat saturation laws — `poly8` (the Flatau et al. 1992 power-sum form of the Fortran's Horner fits), `qsat_esat`/`qsat_esat_dT` (water/ice polynomials a0..a8, b0..b8, c0..c8, d0..d8), `qsat_qs` (0.622·es/(p−0.378·es)) and `qsat_qsdT` (0.622·esdT·p/(p−0.378·es)²) | `QSatMod.F90`: 18-56 (fits), 90-103, 108-120 |
| | soil matric potential `soil_suction` = −sucsat·s^(−bsw) | `SoilWaterRetentionCurveClappHornberg1978Mod.F90:115` |
| | canopy interception/storage/excess laws `canopy_intercept_liq`/`_snow`, `canopy_storage_max`, `canopy_excess_flux` | `CanopyHydrologyMod.F90:538,541,755-758` |
| | roughness/displacement `pft_z0m_zengwang2007`/`pft_displa_zengwang2007`, emissivities `vegetation_emissivity`/`ground_emissivity`, latent heat `latent_heat_of_vapor` | `BiogeophysPreFluxCalcsMod.F90:175-176,344-358,396-397` |
| `clm_canopy.esm` | Model `CanopyHydrology` — the canopy interception window and throughfall/interception fluxes (check_point = frac_veg_nosno==1 ∧ (forc_snow+forc_rain+irrigation) > 0, fpiliq/fpisnow, throughfall, interception), the canopy storage ODEs d(liqcan)/dt = qflx_intercepted_liq − qflx_liqcanfall and d(snocan)/dt = qflx_intercepted_snow − qflx_snocanfall − qflx_snow_unload (the continuous limit of the Fortran's discrete AddInterceptionToCanopy → RemoveCanfallFromCanopy → RemoveSnowUnloading chain), the excess-storage fall fluxes gated on check_point, the Roesch et al. (2001) snow unloading, the wet/dry/snow-covered fractions fwet/fcansno/fdry, and the ground fluxes qflx_snow_grnd/qflx_liq_grnd including the drip-irrigation term; 23 equations | `CanopyHydrologyMod.F90`: CanopyInterceptionAndThroughfall 167-435, interception law 479-566, storage updates 670-762, snow unloading 873-927, ground fluxes 1023-1116, wet/dry fractions 1120-1195 |
| | Model `CanopyFluxInterface` — the pre-flux state consumed by CanopyFluxes/BareGroundFluxes/SoilFluxes: z0m/displa, the ground temperature t_grnd (snow/soil/h2osfc weighted mean), ground/vegetation emissivity, latent heat htvp, thv/thm, the beta/zii constants, and the surface humidity block (s_top, hr with the smpmin clamp, soilalpha, qsatg_soil/qsatgdT_soil with the supersaturated-window clamp, qg_soil, the snow/h2osfc humidities, qg, dqgdT); 23 equations | `BiogeophysPreFluxCalcsMod.F90`: SetZ0mDisp 122-221, CalcInitialTemperatureAndEnergyVars 225-405; `SurfaceHumidityMod.F90`: CalculateSurfaceHumidity 41, 132-216; `FrictionVelocityMod.F90:632,664` (forc_hgt_t_patch) |

The two models import all pointwise laws by reference from
`./clm_templates.esm`; no law is duplicated. The canopy storage ODEs are the
continuous limit of the Fortran's discrete updates — with check_point = 0 and
zero storage the Fortran's `max(0, .)` clamps are inert, so the ODE trace pins
the wiring exactly (see [Testing](#testing)).

## Units

| Quantity | Symbols | Units |
|---|---|---|
| Water fluxes | `forc_rain`, `forc_snow`, `qflx_through_*`, `qflx_intercepted_*`, `qflx_*_grnd`, `qflx_liqcanfall`, `qflx_snocanfall`, `qflx_snow_unload`, `qflx_irrig_*` | mm/s |
| Canopy storage | `liqcan`, `snocan`, `h2ocan`, `liqcanmx`, `snocanmx` | mm |
| Temperatures | `t_soisno_top_*`, `t_h2osfc`, `t_grnd`, `thm`, `thv`, `forc_t`, `forc_th`, `TKFRZ` | K |
| Humidities | `qg`, `qg_soil`, `qg_snow`, `qg_h2osfc`, `qsatg_*`, `forc_q` | kg/kg; `dqgdT`, `qsatgdT_*` | 1/K |
| Pressure | `forc_pbot` | Pa |
| Soil | `s_top`, `watsat_top`, `bsw_top`, `sucsat_top`, `frac_sno_*`, `frac_h2osfc`, `snl`, `soilalpha`, `hr` | 1 (sucsat_top in mm) |
| Emissivities | `emg`, `emv`, `emg_soil` = 0.96, `emg_snow` = 0.97 | 1 |
| Energy | `htvp`, `hvap` = 2.501e6, `hsub` = 2.8347e6 | J/kg |
| Time | `dtime` = 1800 s | s |

**Unit notes**

- `roverg` (the hr law's constant) is `mm/K`, not J/kg/K: clm_varcon.F90:56
  defines `roverg = SHR_CONST_RWV/SHR_CONST_G*1000` = 47062.72789962227 mm/K
  (Rw/g·1000 converts the soil water potential in mm of head to temperature
  units), so `hr = exp(psit/(roverg·T))` is dimensionless.
- `snowcan_unload_wind_fact` is dimensionless (a scaling of the fixed
  1.56e5 s/m divisor of Roesch et al. 2001, CanopyHydrologyMod.F90:918),
  not a rate.
- `sucsat_top` follows the CLM convention of mm (head), hence `smpmin` =
  −1.e8 mm (SoilStateInitTimeConstMod.F90:246).
- `qflx_evap_soi`/`eflx_lh_grnd` in clm_soilfluxes.csv carry htvp = 2,834,700
  J/kg exactly (SoilFluxesMod.F90:394), pinning `hsub = LATICE + LATVAP` =
  3.337e5 + 2.501e6 (SHR_CONST_LATSUB), not the 2.835e6 rounded value.

## Dump defects

The stage-1 CLM dumps were produced by the fork instrumentation
(`stage1/EqESM-cesm-ctsm/.../clm_dump_mod.F90`). The clm-dumps-w2
regeneration (fork clm `2d8e0c697`, LID 261102-004337) fixed the layer
indexing defect of the original dumps; one defect remains (both recorded in
`bugs.md`):

1. **t_soisno lower-bound shift — fixed in W2 (bugs.md #21).** The original
   `dump_soiltemp`/`dump_soilfluxes` wrote every t_soisno/dz/z/zi value
   column 12 layers below its j label (the top snow layer appeared at row
   j = 4, the top soil layer at j = 13, and `dump_soilfluxes`' t_soisno_1..3
   were identically 0.0). W2 declares the 2-D dummy arguments with the dim2
   lower bound `-nlevsno+1`, so the regenerated files carry the true model
   layer index j = snl(c)+1 .. jbot: the top snow layer is j = snl+1 = −8
   (for snl = −9; equal to the preflux t_grnd to round-off on all 192 rows)
   and the top soil layer is j = 1, and the soilfluxes t_soisno_1..3 are
   real (t_soisno_1 equals the soiltemp exit j = 1 row to round-off on all
   192 rows). The clm_canopy.esm pins were re-decoded to the new convention
   (entry j = −8 / j = 1); the values themselves are unchanged — the old
   decode was exact, so only the cited j labels moved (see
   [Testing](#testing)).
2. **qsat derivative columns zero.** `clm_qsat.csv`'s qsdT/esdT are
   identically 0: the traced QSat callers do not request the optional
   derivatives, and the dump writes the unset optionals. The derivative law
   is pinned instead through dqgdT by the surfhum tests (which trace
   callers that do request qsdT).

The surfalbedo dump limitation of the W2 files (no albsnd/albsni columns)
was closed by the W3 surfalbedo extension (fork clm `3d7038304` /
superproject `527ab5d`): `clm_surfalbedo.csv` now also carries the SNICAR
snow albedos, the top-layer volumetric soil water and the soil color class
per column (see `README-clm-albedo.md`, the surface-albedo stage-2
component).

## Coverage gaps

The 18 tests pin exactly what the dumps can support; the following branches
are transcribed but not dump-pinned (all inputs/behaviors are the active
Fortran's, and no tolerance has been loosened or quantity invented to force a
pin):

- **Liquid qsat branch.** All 192 clm_qsat.csv rows are at T < 273.15 K, so
  the td ≥ 0 branches of qsat_esat/qsat_esat_dT (water polynomials) are
  untested.
- **hr/qg_soil/soilalpha interior.** The traced surfhum rows have
  frac_sno_fluxes = 1, so qg = qg_snow and the soil branch of
  CalculateSurfaceHumidity never fires; h2osoi_liq/h2osoi_ice, sucsat, bsw,
  watsat and dz are not dumped, so hr ≈ 0.991 (inferred from the
  qg_soil/qsat(T_soil) ratio of the dump's internal values) cannot be pinned
  directly. The soil-branch QSat *wiring* is now pinned by
  `clm_soilfluxes_soil_qsat_step1_c1` (qsatg_soil_raw/qsatg_soil/qsatgdT_soil
  at the post-solve top-soil temperature from clm_soilfluxes.csv t_soisno_1,
  clamp inactive since qsatg_soil = 4.9e-5 < forc_q = 6.5e-4), but the hr
  equation, its smpmin clamp and the supersaturated-window clamp
  (qsatg_soil → forc_q when qsatg > forc_q > hr·qsatg,
  SurfaceHumidityMod.F90:185-190) are transcribed with the verified constants
  but carry no assertion.
- **qg_h2osfc and the snl = 0 branches.** frac_h2osfc = 0 and snl < 0 on all
  traced rows, so qg_h2osfc/qsatg_h2osfc and the qg_snow = qg_soil fallback
  are untested.
- **Canopy interception/throughfall flux split.** The ent/hyd stages trace
  only the bare patches 1-16 (check_point = 0 everywhere there), and the flx
  stage — which traces the vegetated patches 48-65 — writes NaN for
  fpiliq/fpisnow/qflx_intercepted_*/qflx_through_* (fields defined only at
  the hyd stage), so the interception fractions themselves, the
  throughfall/interception split, the excess-fall fluxes (qflx_liqcanfall/
  qflx_snocanfall) and the nonzero snow unloading (snocan > 0) are untested;
  qflx_irrig_drip is never nonzero in the dump. The vegetated *fraction
  laws* (fwet/fdry/fcansno) are pinned at nonzero canopy storage by the
  flx-stage tests `clm_canopy_fractions_veg_step*` (fwet clamp-saturated at
  0.05 on 182/192 flx rows; the 10 rows with fwet < 0.05 are unpinnable —
  their fwet was evaluated from the pre-evap storage, which is not dumped —
  see Testing). qflx_liq_grnd for the vegetated columns is a wtcol-weighted
  column aggregate (ratio to forc_rain constant per column), not a
  patch-level value, so it is not pinned either. The ODE integration test
  still pins the zero branch only.
- **thv.** forc_th is not dumped, so thv = forc_th·(1+0.61·forc_q) carries no
  assertion.
- **t_grnd mixing.** frac_sno_fluxes = 1 and frac_h2osfc = 0, so only the
  t_grnd = t_soisno_top_snow endpoint is pinned; the weighted mean and the
  snl = 0 branches are untested.
- **htvp/emg/emv/z0m/displa alternatives.** Only the hsub branch (frozen top
  layer with no liquid; pinned with the real top-snow h2osoi_ice =
  5.39439117719703098 kg/m2), the emg_snow endpoint (frac_sno_albedo = 1),
  and the zero elai/esai/htop endpoints are exercised; the hvap branch
  (liquid top layer) is untested.
- **Landunit/column-type branches.** The icol_road_perv branch of
  CalculateSurfaceHumidity and the urban road/roof columns of
  CanopyHydrology are outside the traced filters.
- **Constants not dumped.** beta = 1 and zii = 1000 m
  (BiogeophysPreFluxCalcsMod.F90:362-363) are transcribed from the code.

## Van Genuchten 1980 note

`use_soil_water_retention_curve = 'VanGenuchten1980'` is a namelist option,
but `SoilWaterRetentionCurveVanGenuchten1980Mod.F90` is a verbatim copy of
`SoilWaterRetentionCurveClappHornberg1978Mod.F90` (only the type/module names
differ): its soil_hk (F90:45-84), soil_suction (F90:87-124) and
soil_suction_inverse (F90:127-158) are the Clapp-Hornberg equations, so the
two options are numerically identical in the current source (bugs.md #22).
The `soil_suction` template transcribes the Clapp-Hornberg form (the active
code); a true Van Genuchten implementation would be a stage-3 concern once
the source actually differs.

## Testing

The stage-2 `tests` blocks pin the stage-1 instrumentation dumps
(`stage1/dumps/clm/`, the clm-dumps-w2 regeneration): `clm_qsat.csv`,
`clm_run_imports.csv` (preflux + surfhum stages), `clm_soiltemp.csv`
(entry/exit stages, layer-exact j labels: top snow j = −8, top soil j = 1),
`clm_bareground.csv`, `clm_soilfluxes.csv` (real t_soisno_1..3) and
`clm_canopyflux.csv` (ent/hyd stages for the bare patches 1-16, flx stage
for the vegetated patches 48-65). Run `./esm test stage2/` for the full
suite (404 assertions: 58 DATM/GSWP3 + 135 MOSART + 134 CICE + 77 CLM).
The CLM files themselves are pointwise (no ESD rules), so they validate and
test without ESD_ROOT; the suite-level command still needs
`export ESD_ROOT=/path/to/EarthSciDiscretizations` for the lowered cice/
mosart files.

| file | tests | assertions | what is pinned |
|---|---|---|---|
| `clm_canopy.esm` — CanopyHydrology | 6 | 44 | the no-interception branch against `clm_canopyflux.csv` (hyd stage): throughfall carrying the full snowfall (qflx_through_snow = qflx_snow_grnd = forc_snow to round-off on all 192 hyd rows), zero interception/excess/canopy-fraction rows for a zero-precipitation patch, `clm_canopy_ode_zeros_step6_p3` — the storage ODEs d(liqcan)/dt = d(snocan)/dt = 0 integrated over the full 1800 s step from the ent-stage state (liqcan = snocan = 0), pinning the ODE wiring and the check_point gate (time_span {0, 1800}, `initial_conditions`, MOSART precedent; no scheme gap: the Fortran's discrete max(0, .) clamps are inert at zero storage) — and the vegetated fraction laws at nonzero canopy storage against the flx-stage rows (1,48)/(1,60)/(12,52) (`clm_canopy_fractions_veg_step*`): fwet = 0.05 exactly (the raw law (h2ocan/(vegt·0.1))^0.666666666666 evaluates to 0.436/0.116/0.982, so the clamp saturates — the pin is robust to the pre/post-evap storage ambiguity; 182/192 flx rows carry fwet = 0.05, and the 10 rows with fwet < 0.05 are unpinnable since their pre-evap storage is not dumped), fdry = (1−fwet)·elai/(elai+esai) to round-off, fcansno = 0 (snocan = 0), and the h2ocan/vegt/liqcanmx/snocanmx storage equations (expected values are f64 evaluations at the dumped elai/esai/liqcan/snocan), with the state set via `initial_conditions` at time_span {0, 0} (the flx liqcan is post-evap; the evaporation is outside CanopyHydrology) |
| `clm_canopy.esm` — CanopyFluxInterface | 12 | 33 | the ice-branch QSat law against `clm_qsat.csv` rows 0/48/191(Tmin)/7(Tmax) (qsat_qs = qs to round-off; qg_snow = qsatg_snow for snl < 0); the surface-humidity boundary against `clm_run_imports.csv` surfhum rows (1,1)/(6,8)/(12,16) with t_soisno_top_snow from `clm_soiltemp.csv` entry j = −8 and forc_pbot from `clm_bareground.csv` (qg_snow = qsat_qs(T_snow, p), dqgdT = qsat_qsdT(T_snow, p), qg = frac_sno_fluxes·qg_snow, soilalpha = 1, t_grnd = t_soisno_top_snow — max rel err 3.1e-16 over all 192 surfhum rows); the soil-branch QSat wiring `clm_soilfluxes_soil_qsat_step1_c1` against `clm_soilfluxes.csv` t_soisno_1 = 225.739712678586642 K (= soiltemp exit j = 1 to round-off on all 192 rows): qsatg_soil_raw = qsatg_soil = qsat_qs(T,p) and qsatgdT_soil = qsat_qsdT(T,p) with the supersaturated-window clamp inactive (qsatg = 4.9e-5 < forc_q = 6.5e-4); the thm law thm = forc_t + 0.0098·forc_hgt_t with forc_hgt_t = 30 + z0mg (FrictionVelocityMod.F90:632,664, z0mg = 8.13402050576016129e-05 m) against `clm_bareground.csv` rows (1,1)/(7,9); htvp = hsub = 2,834,700 J/kg against the preflux row and `clm_soilfluxes.csv` (eflx_lh_grnd = qflx_evap_soi·hsub exactly on all 192 rows), with the real top-snow h2osoi_ice = 5.39439117719703098 kg/m2 (soiltemp entry j = −8); and the emg/emv/z0m/displa endpoints emg = 0.97, emv = z0m = displa = 0 against the preflux row (dumped exactly) |

All assertions use `tolerance {rel 1e-12, abs 0}` (the ODE test adds
`abs 1e-18`); every expected value traces to a specific dump row (step,
patch/column) cited in the test description, and the instantaneous tests use
`time_span {0, 0}` (house style). The surfhum tests are reconstructions, not
direct dumps: the input temperatures are taken from the entry-stage
`clm_soiltemp.csv` of the same step/column (the surfhum instant's t_soisno —
the entry dump is the previous step's post-solve array, and the preflux t_grnd
equals the entry j = −8 value to round-off on all 192 rows), and forc_pbot
comes from the same step's `clm_bareground.csv` row. The soil-QSat test is
also a reconstruction: t_soisno_top_soil is the post-solve value from
`clm_soilfluxes.csv` t_soisno_1 (the soilfluxes dump point, one step-solve
increment after the surfhum instant's entry value), with forc_pbot from the
same step's `clm_bareground.csv` row. The flx-stage fraction tests set the
liqcan/snocan ODE states through `initial_conditions` (flx liqcan is
post-evap, so the storage cannot be integrated by CanopyHydrology alone);
the fraction laws are pinned at that state with time_span {0, 0}.

### Transcription fixes found while writing the tests

1. **qsat_qsdT extra 0.622 factor.** The first draft of the derivative
   template squared the whole 0.622/(p−0.378·es) factor, giving
   esdT·0.622²·p/(p−0.378·es)² — the surfhum dqgdT pins failed by exactly
   M_wv/M_da = 0.622. The Fortran is qsdT = esdT·vp1·vp·p with
   vp1 = 0.622·vp (QSatMod.F90:118-120), i.e. 0.622·esdT·p/(p−0.378·es)².
   Fixed in the template; all three surfhum dqgdT pins now pass at rel 1e-12.
2. **roverg units/constant.** The hr law's roverg was initially carried as
   461.5 J/kg/K (the water-vapor gas constant); the unit checker rejected the
   exp argument as length⁻¹·time². clm_varcon.F90:56 defines
   roverg = SHR_CONST_RWV/SHR_CONST_G·1000 = 47062.72789962227 mm/K; with the
   constant and its mm/K units the hr equation closes dimensionally.

### What the ODE test does and does not cover

`clm_canopy_ode_zeros_step6_p3` integrates d(liqcan)/dt and d(snocan)/dt over
[0, 1800] s from zero storage with check_point = 0, asserting the end-of-step
state and ground fluxes against the hyd-stage dump row. With the interception
branch off, all source terms are identically 0, so the integration is exact —
it pins the ODE wiring, the check_point gate and the initial_conditions
plumbing, not the interception physics (the dumps contain no vegetated patch
with precipitation at a stage that carries the interception fluxes; see
Coverage gaps). The Fortran's discrete storage update chain
(AddInterceptionToCanopy → RemoveCanfallFromCanopy → RemoveSnowUnloading,
CanopyHydrologyMod.F90:686-703, 835, 978) becomes the continuous ODEs; at zero
storage the max(0, .) positivity clamps are inert, so the continuous and
discrete forms agree exactly on the pinned trace. The flx-stage fraction
tests (see the table above) complement the ODE test by pinning the
vegetated branch of the fraction laws at nonzero storage — but with the fwet
clamp saturated, so the (h2ocan/(vegt·0.1))^0.666666666666 exponent itself is
only pinned to the inequality raw ≥ 0.05, not to round-off.

## Removed templates

The template library originally carried the snow/soil physics laws
transcribed during stage-1 analysis (soil_hk/soil_dhkds/soil_dsmpds/
soil_suction_inverse, soil/snow thermal conductivity and heat capacity,
frac_sno/fsca laws, snow grain/compaction/layer laws, the delta-Eddington
albedo set, soil/ground albedo laws, topmodel baseflow/aquifer
yield/impedance, soil_water_flux, frozen_water_potential, supercool_water,
snow_percolation_flux, monin_obukhov_zeta, dewpoint_magnus, aerodynamic
ram/rah, snow_albedo_two_stream_approx). None of these is exercised by the
canopy hydrology/flux-interface component, so they were removed to keep the
library to what is exercised. Four of them have since returned with the
surface-albedo stage-2 component (`soil_albedo_wetness`, `soil_albedo`,
`snow_soil_albedo`, `absorbed_shortwave`, `reflected_shortwave` — see README-clm-albedo.md);
the rest will return with their stage-2 components (soil hydrology and
temperature, snow hydrology), where they will be tested against the
corresponding dumps. Two of the removed laws are documented in bugs.md: the
Van Genuchten 1980 duplicates (#22) and the soil_hk/soil_suction
Clapp-Hornberg provenance (kept in the template descriptions).

## Dump extensions that would close the gaps

The clm-dumps-w2 regeneration (bug #21 fixed) already made the soiltemp dump
directly usable layer-by-layer and the soilfluxes t_soisno_1..3 real — the
clm_canopy.esm pins are re-decoded and extended accordingly (vegetated
canopy fractions from the flx stage, soil-branch QSat wiring from
t_soisno_1). To pin the remaining branches, the stage-1 instrumentation
should additionally dump (per step, per traced patch/column): the
interception/throughfall fluxes (or fpiliq/fpisnow) at a stage that traces
the vegetated patches — the hyd stage currently traces only the bare patches
1-16, while the flx stage writes NaN for them (the vegetated fraction pins
are exact only because the fwet clamp saturates; a vegetated hyd-stage row
would pin the interception laws, the excess-fall fluxes and the storage ODE
with evap off, and the pre-evap storage that the 10 rows with fwet < 0.05
were evaluated from); the preflux h2osoi_liq/h2osoi_ice, sucsat, bsw, watsat
and dz of the top layer (for hr/soilalpha/qg_soil and the supersaturated-
window clamp); a T > 273.15 QSat call (for the liquid branch); a
frac_h2osfc > 0 state (for the h2osfc humidity branch); a snow-unloading
state (snocan > 0 with forc_t/forc_wind); the patch-level qflx_liq_grnd for
the vegetated columns (the dumped value is the wtcol-weighted column
aggregate); and forc_th (for thv). The frac_sno_albedo + albsnd/albsni
extension (ground-albedo mixing) is done: the W3 surfalbedo dump extension
(see README-clm-albedo.md) adds the SNICAR snow albedos, h2osoi_vol_1 and
isoicol to clm_surfalbedo.csv; what is still missing for the *interior* of
the mixing law is a traced state with 0 < frac_sno_albedo < 1.

## Next steps (stage 2.5 / stage 3)

The canopy hydrology/flux-interface skeleton and the surface albedo
component (README-clm-albedo.md) are done; the remaining CLM pieces (snow
hydrology and soil temperature/moisture, canopy fluxes/biogeophysics, and
the top-level clm assembly) follow the same pattern — transcribe the Fortran
laws into templates, mount them in the component models, and pin them against
the remaining `stage1/dumps/clm/` files (`clm_soiltemp.csv` exit rows,
`clm_canopyflux.csv` flx stage (t_veg/eflx_sh_veg/qflx_evap_veg for the
canopy-flux energy balance), `clm_bareground.csv`, `clm_soilfluxes.csv`).
The W2 soiltemp/soilfluxes files are now layer-exact, so the soiltemp
component can pin its layer temperatures directly (entry/exit j = −8..13+).
Stage 3 couples the merged components into the top-level `eqesm.esm` model.
