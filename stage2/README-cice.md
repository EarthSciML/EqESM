# CICE stage 2 — `.esm` skeleton

Stage-2 implementation of the CESM CICE sea-ice component (cesm3_cice6_6_3_22,
checkout `stage1/cesm/components/cice`) as esm 1.3.0 files: component structure,
imports/exports, solver block, the pointwise column/dynamics/transport laws of
the default configuration transcribed exactly, subcomponents, and top-level
assembly. The `tests` blocks are empty pending the stage-1 instrumentation
dumps (`stage1/dumps/cice/`) — see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/cice/src/cicecore/drivers/nuopc/cmeps/`
(`ice_import_export.F90`, `CICE_RunMod.F90`),
`stage1/cesm/components/cice/src/cicecore/cicedyn/general/`
(`ice_step_mod.F90`, `ice_forcing.F90`, `ice_flux.F90`),
`stage1/cesm/components/cice/src/cicecore/cicedyn/dynamics/`
(`ice_dyn_evp.F90`, `ice_dyn_shared.F90`),
`stage1/cesm/components/cice/src/cicecore/cicedyn/analysis/`
(`ice_transport_driver.F90`, `ice_transport_remap.F90`),
`stage1/cesm/components/cice/src/icepack/columnphysics/`
(`icepack_shortwave.F90`, `icepack_therm_shared.F90`, `icepack_therm_bl99.F90`,
`icepack_therm_vertical.F90`, `icepack_meltpond_lvl.F90`, `icepack_atmo.F90`,
`icepack_mechred.F90`, `icepack_mushy_physics.F90`, `icepack_flux.F90`,
`icepack_parameters.F90`), and the boundary plan in `stage1/boundaries.md` (§1.3).

Default configuration transcribed (the D-class reference run's namelist,
`stage1/cime/cases/D_cice_inst/CaseDocs/ice_in`, unless noted): `shortwave = 'dEdd'`
(albedo/shortwave laws — the delta-Eddington 3-band solver itself is deferred,
its outputs carried as parameters), `kdyn = 1` (EVP), `grid_ice = 'B'`, `ndtd = 1`, `ndte = 120`,
`revised_evp = .false.`, `elasticDamp = 0.36`, `Ktens = 0`,
`e_yieldcurve = e_plasticpot = 2`, `capping = 1`, `u0 = 5e-5`, `cosw = 1`,
`sinw = 0`, `kstrength = 1` (the Rothrock potential-energy strength is
transcribed and dump-pinned in the coupling cap; the ridge redistribution
itself is a stage-3 dynamics concern), `atmbndy = 'similarity'` (the
Monin-Obukhov iteration is transcribed and dump-pinned in the coupling cap;
the constant scheme remains for the assembly Control), `ktherm = 2` (mushy), `conduct = 'MU71'`,
`tfrz_option = 'mushy'`, `saltflux_option = 'constant'`, `frzpnd = 'hlid'`
(sea-level ponds, `tr_pond_sealvl`), `snwredist = 'none'`, `hs0 = 0.0`,
`phi_i_mushy = 0.45`, `congel_freeze = 'one-step'`, `emissivity = 0.95`
(forcing_nml), `natmiter = 5`, `ncat = 5`, `nilyr = 8`, `nslyr = 3`,
`l_brine = .true.`, `dt = 3600 s`.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `cice_templates.esm` | Expression-template library (76 templates): snow cover / albedo / absorbed-penetrating shortwave, the dEdd snow fraction and effective pond area | `icepack_shortwave.F90`: compute_albedos 466-597, absorbed_solar 701-896, shortwave_dEdd_set_snow 3448-3454, shortwave_dEdd_set_eff 3495-3637 |
| | saturation humidity, surface heat flux, liquidus (linear + mushy), enthalpy↔temperature (BL99 brine quadratic, mushy quadratic, fresh, snow) | `icepack_therm_shared.F90`: surface_heat_flux 147-165, liquidus 383-418; `icepack_mushy_physics.F90`: liquidus_temperature_mush 227-261, icepack_mushy_temperature_mush 348-452, enthalpy_of_melting 264-286, liquid_fraction_mush 288-310; calculate_Tin_from_qin 68-101; snow temperature 470-484 |
| | conductivity (bubbly Pringle 2007 / MU71 / mushy) and interface conductances, congelation enthalpy, ocean heat flux, melt-pond volume change (cesm and sea-level schemes) | `icepack_therm_bl99.F90`: conductivity 874-961 (MU71 911-916); `icepack_mushy_physics.F90`: conductivity_mush 455-485; `icepack_therm_vertical.F90`: qbot 1432-1443, fbot 662-678; `icepack_meltpond_lvl.F90`: 154-174; `icepack_meltpond_sealvl.F90`: 76, 157 |
| | constant atmo boundary layer (tau, shcoef, lhcoef) and the similarity scheme (neutral drag, relative wind, stability parameter, stable/unstable psi, exchange coefficients, 2-m psix2, Tref/Qref/Uref) | `icepack_atmo.F90`: atmo_boundary_const 398-496 (460-464, 493-494); atmo_boundary_layer 52-389 (215, 176-210, 985-989, 1015-1075, 289-291, 364-378) |
| | Hibler and Rothrock strength (cumulative distribution, exponential participation, ridge redistribution hrmin/hrexp/krdg, aksum, strength term), Delta invariant, viscosities/replacement pressure, classic-EVP stress maps, ocean-ice stress, seabed stress | `icepack_mechred.F90`: 1720, ridge_itd 848-1036, icepack_ice_strength 1660-1712; `ice_dyn_shared.F90`: strain_rates 2081-2161, visc_replpress 2444-2473, set_evp_parameters 451-484; stepu 923-964 |
| `cice_thermodynamics.esm` | Model `Thermodynamics` — the D-class column physics: prep_radiation dEdd shortwave/albedo, the mushy (ktherm = 2) vertical heat conduction and surface energy balance, growth/melt/ponds (frzpnd = hlid sea-level) and ocean heat/fresh/salt fluxes, plus the cell-level merged aggregates of icepack_flux.F90:196-257 (aice, fsurf, fcondtop, fhocn, fresh, fsalt, flwout, fsens, flat, evap, fswabs, fswthru, fswint, congel, meltb, meltt, melts, dsnow); 84 equations over ncat × nilyr/nslyr (82 algebraic pointwise laws + 2 enthalpy PDEs) | `icepack_shortwave.F90` (step_radiation, prep_radiation 3893-3922, dEdd set_eff/set_snow), `icepack_therm_mushy.F90` (two-stage mushy tridiagonal 1800-1985, semi-discrete flux-divergence ODEs), `icepack_mushy_physics.F90` (temperature/enthalpy inversions), `icepack_therm_vertical.F90` (thickness changes, freshn/fsaltn 510-522, fhocnn 1625-1630, sublimation 1355-1403, 1502-1559), `icepack_flux.F90` (merge_fluxes 180-260), `icepack_meltpond_sealvl.F90` (compute_ponds_sealvl, dvn) |
| `cice_dynamics.esm` | Model `Dynamics` — B-grid momentum balance (continuous limit of stepu) and classic-EVP stress relaxation ODEs, Hibler strength, viscosities; 34 equations over one cell | `ice_dyn_shared.F90`: stepu 923-964, dyn_prep2 815-833, dyn_finish 1345-1348, Tmass 540; `ice_dyn_evp.F90`: stress 1584-1609 |
| `cice_transport.esm` | Model `Transport` — per-category conservative advection d(a)/dt + div(a·u) = 0 (and volume forms) as spatial-D rewrite targets; enthalpy-tracer statement documented, deferred to stage 3 | `ice_transport_driver.F90`: 251-700; `ice_transport_remap.F90`: 317-831 |
| `cice_coupling.esm` | Model `IceCap` — med↔ice cap renames, areacor scaling, Sa_pbot potT/rhoa branch, export aggregation (ailohi, Tsrf, merge_fluxes sums), Faii_lwup fix; the similarity atmo boundary layer (5-iteration unroll) with the aicen_init-weighted reference-height exports, the Rothrock strength law, and the dynamics pointwise-law probes (EVP stress maps, strain invariant, viscosities/replacement pressure, stepu ocean-ice stress) as algebraic probes of the dynamics physics; 166 equations | `ice_import_export.F90`: ice_import 414-700, ice_export 1030-1317; `icepack_flux.F90`: merge_fluxes 180-260; `icepack_atmo.F90`: atmo_boundary_layer 52-389; `icepack_mechred.F90`: ridge_itd 848-1036, icepack_ice_strength 1660-1712 |
| `cice.esm` | Top-level assembly: mounts the four subcomponents by reference + inline `Control` (dt_dyn subcycling, constant atmo boundary layer), 59 `variable_map` coupling entries, solver block | `CICE_RunMod.F90`: ice_step 109-369; `ice_step_mod.F90`: step_dyn_horiz 990-1074; `ice_forcing.F90`: prepare_forcing 1606-1825; `icepack_atmo.F90`: 398-496 |

The PDEs are written as spatial-D rewrite targets
(`div(c·u) = D(c·uvel, x) + D(c·vvel, y)`) to be lowered by discretization
templates from `github.com/earthsciml/EarthSciDiscretizations` in stage 3,
exactly as the incremental-remapping scheme replaces the continuous flux form
with the departure-cell integrals of the Fortran. The momentum equation is the
continuous limit of the Fortran's implicit stepu solve (the cca/ccb/ab2 algebra
is solver machinery, documented rather than transcribed); the classic-EVP
stress maps become relaxation ODEs d(σ)/dt = (σ_target − σ)/Td with
Td = 2·elasticDamp·dt.

## Units

All fields are SI; temperatures are °C internally (`degC` ≡ K in the unit
checker's algebra) with `Tffresh = 273.15` added where Kelvin is needed;
salinity is ppt; the CMEPS flux imports are W/m² and kg/m²/s (the latter equals
mm/s of water).

| Quantity | Symbols | Units |
|---|---|---|
| Thicknesses / volumes / depths | `hi`, `hs`, `hilyr`, `hslyr`, `vicen`, `vsnon`, `hpnd`, `vice`, `vsno`, `Si_snowh` | m (volumes per unit grid-cell area) |
| Enthalpies | `qicen`, `qsnon`, `qbot` | J/m³ |
| Temperatures | `zTin`, `zTsn`, `Tmlts`, `Tbot`, `Tsfc`, `Tf`, `sst` | degC; `potT`, `Tair`, `Si_t` in K |
| Conductivities / conductances | `kilyr`, `kslyr`; `kh_top`, `kh_snowice`, `kh_bot` | W/m/K; W/m²/K |
| Fluxes | `fswsfc`, `fswpen`, `fswthru`, `fsurfn`, `fsensn`, `flatn`, `flwoutn`, `fbot`, `fhocnn`, `frzmlt`, sw bands, Faii_*/Fioi_* | W/m² |
| Water fluxes | `evapn`, `freshn`, `fsaltn`, `fsnow` | kg/m²/s; `frain` in m/s (cap conversion by rhofresh) |
| Stresses | `strength`, `stressp/m/12`, `sigma*`, `rep_prs` | N/m; `tau`, `strairx/y`, `taux/y`, `strintx/y`, `strtltx/y`, `Cb` | N/m² |
| Strain rates | `divu`, `tension`, `shear`, `Delta`, `deltaminEVP` | 1/s |
| Mass | `m` (Tmass) | kg/m² |
| Viscosities | `zetax2`, `etax2`, `tmpcalc` | kg/s |
| Dynamics | `uvel`, `vvel`, `uocn`, `vocn`, `uatm`, `vatm`, `wind`, `u0` | m/s |
| Air | `rhoa` | kg/m³; `Qa`, `Qsfc` kg/kg; `Sa_z`/`zlvl` m; `Sa_pbot` Pa |
| Time | `dt`, `dt_dyn`, `Td` | s |
| Constants | `puny = 1e-11`, `Tffresh = 273.15`, `rhoi = 917`, `rhos = 330`, `rhow = 1026`, `cp_ice = 2106`, `cp_ocn = 4218`, `cp_air = 1005`, `Lfresh = 3.34e5`, `Lsub = 2.835e6`, `Lvap = 2.501e6`, `depressT = 0.054`, `emissivity = 0.95` (D-class forcing_nml), `stefan_boltzmann = 5.67e-8`, `gravit = 9.80616`, `vonkar = 0.4`, `zref = 10`, `iceruf = 5e-4`, `zvir = 0.606`, `cp_wv = 1810`, `Cf = 17`, `Gstar = 0.15`, `astar = 0.05`, `maxraft = 1`, `mu_rdg = 4` (D-class namelist) | (icepack_parameters.F90:36, 82, 62, 61, 65, 74-75, 281, 83-84, 77, 159, 56, 190, 284-285, 283, 282, 261, 269; icepack_mechred.F90:70-72, 74; icepack_parameters.F90:2814) |

**Unit notes**

- `"C"` parses as Coulomb in the esm unit checker, so all Celsius fields use
  `degC` (identical to K in the checker's algebra).
- `lhcoef` is W/m² (per specific-humidity difference, not per K): the Fortran
  `lhcoef = latncoef·Lheat·rhoa·wind` (icepack_atmo.F90:494) carries the latent
  heat inside the coefficient.
- `evapn` requires the ρ in the numerator: the Fortran
  `evapn += dhs·rhos` (icepack_therm_vertical.F90:1502-1521) is transcribed as
  the continuous rate `evapn = esub·rhos/((qsn − rhos·Lvap)·dt)`; `fsaltn` is
  the rate form `fsaltn = −rhoi·dhi·ice_ref_salinity·0.001` with the same net
  thickness change `dhi = congeln + meltbn + melttn + evapin/rhoi` as freshn
  (F90:513, 521; the per-step dt divides out of both rate forms).
- `freshn` in rate form: the Fortran's `freshn = evapn − (rhoi·dhi +
  rhos·dhs)/dt` (F90:510-513, per-step) cancels the sublimation/condensation
  and new-snowfall terms exactly — `evapn = evapin + evapsn` and
  `dhs = −meltsn − evapsn/rhos` with `hsn_new` added back to dzs(1) before
  hin/hsn are summed (F90:1641-1672) — leaving
  `freshn = −(rhoi·(congeln + meltbn + melttn) + rhos·meltsn)`; the rain is
  added at the cell level (`fresh = fresh + frain·aice`,
  icepack_therm_itd.F90:2043, with the model's m/s frain × rhofresh).
- All icepack module constants are inlined as bare literals inside the
  templates; the unit checker treats bare literals as dimension-Unknown and
  skips products/comparisons containing them, which is how the many
  Fortran-sloppy unit mixes (e.g. `meltt·rhoi + frain·dt` in the melt-pond
  law) validate.
- `sign()` preserves its argument's dimension, so the Fortran's
  `sign(1, fm)` of the ocean-stress rotation is carried as the `hemis`
  parameter (+1/−1).

## §1.3 mapping, as implemented

| Boundary (boundaries.md §1.3) | Proposed | Implemented |
|---|---|---|
| Import (med→ice): atm state, rad/precip, ocn state, Fioo_q | `cice_coupling.esm` | ✓ `IceCap`: renames + areacor + Sa_pbot potT/rhoa branch; uocn/vocn/ss_tltx/ss_tlty wired to Dynamics; potT/Qa/rhoa/sw*/flw/sst/sss/frzmlt/frain/fsnow wired to Thermodynamics; uatm/vatm/rhoa wired to the assembly `Control` atmo scheme |
| Export (ice→med): Si_* state, Faii_* atm/ice fluxes, Fioi_* ice/ocn fluxes | `cice_coupling.esm` | ✓ `IceCap`: ailohi, Tsrf, aggregate albedos, vice/vsno/snowh, merge_fluxes sums, areacor; Fioi_taux/tauy sign flip from the Dynamics full ocean stress |
| Wave exports (Si_thick, Si_floediam), aerosols (Fioi_bcpho/bcphi, Fioi_flxdst), isotopes, per-band Fioi_swpen_* | `cice_coupling.esm` | ◐ deferred — carried as zero parameters with file:line notes (wave/aerosol/isotope tracers off; per-band fswthru_* not in the stage-2 column skeleton) |
| EVP dynamics + icepack column thermodynamics (the run step) | `cice_dynamics.esm` + `cice_thermodynamics.esm` | ✓ mounted in `cice.esm` with the ice_step ordering documented (step_prep → prep_radiation → step_therm1 → step_therm2 → update_state → ndtd loop step_dyn_horiz(evp → transport_remap) → step_dyn_ridge → update_state) |
| Transport (advection) | `cice_transport.esm` | ✓ continuity system as spatial-D rewrite targets; the enthalpy-tracer equations documented, deferred to the stage-3 incremental-remap lowering |
| Ridging (kstrength = 1 Rothrock) | `cice_coupling.esm` (algebraic probe) | ◐ the Rothrock strength law (icepack_mechred.F90:1660-1712 with the ridge_itd state, krdg_partic = krdg_redist = 1, mu_rdg = 4) is transcribed and dump-pinned in `IceCap` (the dynamics skeleton cannot compile tests on the single-cell geometry); the ridge redistribution of area/volume/enthalpy (ridge_shift, F90:1068+) is a stage-3 dynamics concern |
| Tsfc implicit surface energy balance | `cice_thermodynamics.esm` | ◐ the balance fsurfn(Tsfc) + fcondtopn = 0 (icepack_therm_bl99.F90:1116-1128) is transcendental in Tsfc (saturation-humidity exponential); esm 1.3.0 cannot express implicit equations, so Tsfc is a parameter fed by the coupled solve in stage 3 |

## Pointwise laws — Fortran vs esm transcription notes

All templates are match-less, invoked by name (`apply_expression_template`
with bindings), and transcribed exactly with file:line citations in the
descriptions. Representative transcriptions:

| Template | Fortran | esm |
|---|---|---|
| `snow_cover_fraction` | `asnow = hs/(hs+snowpatch)` for `hs > puny` (icepack_shortwave.F90:574-579) | `ifelse(hs > 1e-11, hs/(hs+0.02), 0)` |
| `thickness_albedo_factor` | `fh = min(atan(hi·4)/atan(ahmax·4), 1)` (F90:541-542) | `min(atan(hi·4)/atan(0.3·4), 1)` |
| `layer_absorbed_shortwave` | `Iswabs = fswpen·(exp(−kappav·ztop) − exp(−kappav·zbot))` (F90:862-880) | `fswpen·(transmittance(ztop) − transmittance(zbot))`, calling `shortwave_transmittance` by reference |
| `surface_heat_flux` | `fsurfn = fswsfc + ε·flw − ε·σ·TsfK⁴ + shcoef·(potT−TsfK) + lhcoef·(Qa−Qsfc)` (icepack_therm_shared.F90:147-165) | same, with `Qsfc` from `saturation_specific_humidity_ice` |
| `ice_temperature_from_enthalpy` | brine quadratic `cp_ice·T² + ((cp_ocn−cp_ice)·Tmltk − qin/ρi − L)·T + L·Tmltk = 0` solved by the negative root (F90:68-101) | `ifelse(l_brine > 0.5, min((−bb1−√(bb1²−4·cp_ice·cc1))/(2·cp_ice), Tmltk), (L+qin/ρi)/cp_ice)`; the Fortran abort on a negative discriminant is a noted guard; the D-class run uses the mushy inversion instead (`icepack_mushy_temperature_mush`) |
| `ice_conductivity_bubbly` | `k = max((2.11 − 0.011·T + 0.09·S/min(−puny,T))·ρi/917, kimin)` (icepack_therm_bl99.F90:918-924) | same with `min(−1e-11, T)` |
| `ice_conductivity_mu71` | `k = max(rhoi·(1.88 + 0.0102·T) − 0.9e-3·T², kimin)` (icepack_therm_bl99.F90:911-916) | same; the D-class `conduct = 'MU71'` law |
| `ice_conductivity_mush` | mushy two-stage `k = max(k1, k2)` with the solid-fraction weighted law (icepack_mushy_physics.F90:455-485) | same; used by the ktherm = 2 solver |
| `congelation_enthalpy` | `qbot = −ρi·(cp_ice·(Tmlts−Tbot) + L·(1−Tmlts/Tbot) − cp_ocn·Tmlts)`, capped at `−0.5·ρi·L` (icepack_therm_vertical.F90:1432-1443) | same, `ifelse` brine/fresh, `min` cap |
| `ocean_heat_flux` | `fbot = max(−cp_ocn·rhow·0.006·max(sst−Tf,0)·ustar, frzmlt)` (F90:662-678) | same; ustar floored at 0.005 in the model |
| `meltpond_cesm_dvn` | `dvn = rfrac/rhofresh·(meltt·ρi + melts·ρs + frain·dt)·aicen − volpn·(1−exp(rexp·dTs/Tp))` (icepack_meltpond_lvl.F90:154-174) | same with `Tp = −2`, `rexp = 0.01`, `dTs = max(−2−Tsfcn, 0)`; the D-class run uses the sea-level scheme instead (`meltpond_sealvl_dvn`) |
| `atmo_constant_tau/shcoef/lhcoef` | `tau = ρa·0.0012·wind`; `shcoef = 0.0012·1005·ρa·wind`; `lhcoef = 0.0015·Lheat·ρa·wind` (icepack_atmo.F90:460, 493-494) | identical, `Lheat = Lsub = 2.835e6` for ice |
| `atmo_neutral_drag` | `rdn = vonkar/log(zref/iceruf)` (icepack_atmo.F90:215) | identical (`vonkar = 0.4`, `zref = 10`, `iceruf = 5e-4`) |
| `atmo_stability_parameter` | `hol = vonkar·gravit·zlvl·(tstar/thva + qstar/(1/zvir+Qa))/ustar²`, clamped to |hol| ≤ 10 (F90:985-989) | identical (sign-preserving clamp via min/max) |
| `atmo_psi_momentum_unstable` / `atmo_psi_scalar_unstable` | Paulson-1970 forms with `capital_X = sqrt(max(sqrt(|1−16·hol|), 1))` (F90:1040-1075) | identical, `pih = π/2`; the stable branch `psi_stable = −(0.7·hol + 0.75·(hol−14.3)·exp(−0.35·hol) + 10.7)` (F90:1017-1019) |
| `atmo_reference_temperature/humidity/wind` | `Tref = potT − delt·(rh/vonkar)·(alzs+al2−psixh+psix2) − 0.01·zTrf`; `Qref = Qa − delq·(re/vonkar)·(…)`; `Uref = vmag_rel·rd/rdn` (F90:366-378) | identical; psix2 = −5·hols·stable + (1−stable)·ψ_scalar_unstable(hols) at hols = hol·2/zlvl (F90:364-365) |
| `ice_strength_hibler` | `strength = Pstar·vice·exp(−Cstar·(1−aice))` (icepack_mechred.F90:1720) | identical (`Pstar = 2.75e4`, `Cstar = 20`) |
| `ridge_cumulative_area` | `Gsum(n) = Gsum(n−1) + aicen(n)` for `aicen(n) > puny`, `Gsum(0) = aice0` for `aice0 > puny` (F90:850-879) | nested prefix sum over m ≤ n with the aicen/aice0 masks |
| `ridge_participation_exponential` | `apartic = [exp(−G(n−1)/astar) − exp(−G(n)/astar)]/(1 − exp(−astari))`, `astari = 20` (F90:919-934) | identical (`astar = 0.05`) |
| `ridge_hrmin/hrexp/krdg` | `hrmin = min(2·hi, hi+maxraft)`; `hrexp = mu_rdg·√hi`; `krdg = (hrmin+hrexp)/hi` for `aicen > puny`, else `krdg = 1` (F90:1004-1010) | identical (`maxraft = 1`, `mu_rdg = 4`), `hi = max(vicen/aicen, puny)` |
| `ridge_strength_term` | `apartic·(h2rdg/krdg − hi²)` with `h2rdg = hrmin² + 2·hrmin·hrexp + 2·hrexp²` (F90:1697-1706) | identical; the cap sums it into `strength = Cf·Cp·Σ/aksum` (F90:1710) |
| `viscosity_scale` | `tmpcalc = capping·P/max(Δ,Dmin) + (1−capping)·P/(Δ+Dmin)` (ice_dyn_shared.F90:2467-2468) | identical; `Dmin = deltaminEVP = 1e-11` |
| `evp_stressp_update` | `stressp_new = (stressp·(1−arlx1i·revp) + arlx1i·(zetax2·divu − rep_prs))·denom1` (ice_dyn_evp.F90:1584-1591) | transcribed as the relaxation ODE `D(stressp)/t = (stressp_target − stressp)/Td` with `Td = 2·elasticDamp·dt` (see divergences) |

## Documented divergences (all annotated in file metadata)

1. **Heat conduction as semi-discrete ODEs.** The Fortran assembles and solves
   an implicit backward-Euler tridiagonal for Tsfc and the layer temperatures —
   the mushy two-stage solver (icepack_therm_mushy.F90:1800-1985) under the
   D-class ktherm = 2, with the kh conductances of the conductivity
   subroutines (icepack_mushy_physics.F90:455-485) and the per-layer absorbed
   shortwave of absorbed_solar (icepack_shortwave.F90:862-880). The .esm writes
   the equivalent energy-conserving method-of-lines form
   `d(q_k)/dt = (kh_up·(T_up−T_k) + kh_dn·(T_dn−T_k) + Iswabs_k)/hilyr`
   (volumetric-enthalpy balance, no ρ·cp factor). The stage-3 integrator
   reproduces the backward-Euler trajectory; the layer-energy and growth/melt
   rates are pinned at stage 3 once the Tsfc implicit solve is coupled.
2. **Tsfc is a parameter.** The implicit surface energy balance
   `fsurfn(Tsfc) + fcondtopn = 0` (icepack_therm_bl99.F90:1116-1128) is
   transcendental in Tsfc; esm 1.3.0 has no implicit equations, so Tsfc is a
   free parameter of the thermodynamics model fed by the coupled solve.
3. **Momentum equation as continuous limit.** stepu's implicit-solve algebra
   (cca/ccb/ab2/cc1/cc2, ice_dyn_shared.F90:949-959) is solver machinery; the
   skeleton writes `m·du/dt = τa + τtilt + τo + ∇·σ − Cb·u + m·f·k̂×u` with the
   full ocean stress of dyn_finish (F90:1345-1348) and the tilt stress
   `strtltx = −gravit·m·ss_tltx` (F90:826).
4. **EVP stress evolution as relaxation ODEs.** The classic-EVP under-relaxed
   subcycle maps (ice_dyn_evp.F90:1584-1609, `arlx = 2·elasticDamp·ndte = 86.4`,
   `brlx = ndte`) are written as `d(σ)/dt = (σ_target − σ)/Td` with
   `Td = 2·elasticDamp·dt`; the four-corner (NE/NW/SW/SE) staggered stresses
   are collapsed to one cell value (documented; the B-grid stencil is a stage-3
   discretization concern).
5. **Similarity iteration unrolled to 5 fixed updates.** The Monin-Obukhov
   do-while (icepack_atmo.F90:269) is unrolled for the D-class `natmiter = 5`
   and `atmiter_conv = 0` (unset), which always runs the full 5 updates; the
   always-5 unroll reproduces the dumped reference-height exports to machine
   precision for Tref/Qref. The assembly `Control` still uses the constant
   scheme for the stresses (strairxT/strairyT stay parameters; the similarity
   stress wiring is a stage-3 concern).
6. **Dynamics physics in the coupling cap.** The kstrength = 1 Rothrock
   strength law (icepack_mechred.F90:1660-1712) and the dynamics pointwise
   laws — the EVP stress maps (`evp_stressp/m/12_update`,
   ice_dyn_evp.F90:1584-1609), the strain invariant and
   viscosities/replacement pressure (`strain_invariant_delta`,
   `viscosity_scale`, `bulk_viscosity`, `shear_viscosity`,
   `replacement_pressure`, ice_dyn_shared.F90:2156-2159, 2467-2471), and the
   stepu ocean-ice stress (`ocean_ice_stress_magnitude`,
   `seabed_stress_coefficient` plus the waterx/watery projections, F90:817-818,
   931-937) — are transcribed and dump-pinned in `IceCap` because the dynamics
   skeleton cannot compile tests (unlowered spatial D); the `cice_evp_stress.csv`
   / `cice_evp_momentum.csv` / `cice_ice_strength.csv` rows reproduce to rel
   1e-12 / 1e-9. The Hibler kstrength ≠ 1 form stays available as
   `ice_strength_hibler`. The spatial-D-dependent ODEs and stencils of
   `cice_dynamics.esm` remain stage-3 deferrals (see the coverage limits).
7. **Enthalpy-tracer transport documented, not carried.** The conservative
   forms `d(vicen·qicen)/dt + div(vicen·qicen·u) = 0` are the Fortran's type-2
   tracer transport (ice_transport_driver.F90:98-102); the product-form LHS
   cannot define an unknown under the esm equation-balance rule, and the
   tracers are parameters fed by the thermodynamics states, so the equations
   are documented and asserted at stage 3 when the incremental-remap
   discretization lands.
8. **Grid and stagger collapse.** The single-cell skeleton collapses the B-grid
   staggering (T/U points), the ANGLET POP-grid rotations of the export
   stresses (ice_import_export.F90:962-968), the `aiU` division of the
   T-point ocean-stress average (ice_step_mod.F90:1051-1072), and the tmask
   ocean-mask guards; each collapse is annotated at the affected variable.
9. **`frain` in m/s.** The cap converts Faxa_rain kg/m²/s to m/s by rhofresh
   for the melt-pond term `frain·dt` (icepack_meltpond_lvl.F90:160), matching
   the skeleton's rate-form meltt/melts; the Fortran's per-step lengths vs
   the esm rate convention is documented in the melt-pond template.

## Bugs found in the Fortran code

| # | Location | Bug / quirk |
|---|---|---|
| 1 | `icepack_atmo.F90:463-464` (atmo_boundary_const) | `strx = tau·uatm` omits the `/wind` direction-cosine normalization — the vmag division is commented out at F90:328-329 and the active stress is dimensionally `kg/s³` (a 10 m/s wind gives a 10× stress). Transcribed as written (the stage-3 tests must pin the active code); the esm skeleton declares the honest `kg/s^3` dimension at Control.strairx/strairy and the dynamics consumes the same value with its N/m² momentum balance annotated. |
| 2 | `icepack_therm_shared.F90:90-94` (calculate_Tin_from_qin) | The Fortran aborts on a negative discriminant of the brine quadratic; the .esm skeleton cannot express the abort and carries it as a documented guard for the stage-2 tests (a negative discriminant in a dump would expose a thermodynamics inconsistency upstream). |
| 3 | `icepack_shortwave.F90:3587-3588` (shortwave_dEdd_set_eff) | `dhs = hsnlvl − fsnow·dt` subtracts the snowfall **mass flux** (kg/m²/s, F90:3511) from a snow **depth** (m) — the `/rhos` conversion is missing (bugs.md #10). The .esm reproduces the Fortran numerics by carrying the coupled flux in an m/s parameter (`fsnow_dEdd`); the branch fires in 92/1920 dumped rows but is usually masked by the `ipn·spn < puny` reset at F90:3590. |

## Testing

The stage-2 `tests` blocks pin the stage-1 instrumentation dumps
(`stage1/dumps/cice/`) with instantaneous laws (`time_span {0,0}`, MOSART
pattern). Run `./esm test stage2/` for the full suite (115 MOSART + CICE):

| file | tests | assertions | what is pinned |
|---|---|---|---|
| `cice_coupling.esm` | 40 | 69 | the import laws (areacor scaling, frain/rhofresh, potT/rhoa branch wiring), the freezing-temperature law, the export state fields, the per-category albedo aggregation, the scale_fluxes `/aice` flux aggregation, the air-ice stress exports, the similarity reference-height exports Si_tref/Si_qref/Si_u10 (4 cells × 2 steps), the Rothrock strength law (2 dump rows), the EVP stress pointwise laws Delta/zetax2/etax2/rep_prs and the newp/newm/new12 subcycle maps (2 dump rows, different cells/corners/ksub), and the stepu ocean-ice stress laws vrel/taux/tauy/Cb (2 dump rows) |
| `cice_thermodynamics.esm` | 14 | 22 | the dEdd shortwave laws against `cice_albedo.csv`: the snow fraction on bare ice and snow (vsnon = 0 / > 0), the effective pond area apeffn through the infiltration, hypsometry, hpmin and taper branches, the ice/snow broadband albedos, and the shortwave conservation fswintn = swabs − fswsfc − fswthrun; the mushy enthalpy↔temperature inversions zTin/zTsn (eqs 1-2) against the `cice_thermo_state.csv` zTin_1..8/zTsn_1..3 columns (all `rel 1e-12`); the step-1 surface energy-block aggregates flwout/fsens against `cice_thermo_flux.csv` (the similarity-driven shcoef/lhcoef at the entry Tsfcn, fluxes at the exit Tsfcn, ε = 0.95; `rel 1e-9`, ITD-weighting residual); and the step-2 radiation aggregates fswabs/fswthru testing the prep_radiation rescale `albedo_dump(1) × netsw(2)/netsw(1)` (icepack_shortwave.F90:3893-3922; `rel 1e-12`, cells (3,5) and (4,5)) |
| `cice_dynamics.esm` | 0 | 0 | none — unlowered spatial-D stencils block compilation of any test (see coverage limits) |
| `cice_transport.esm` | 0 | 0 | none — same spatial-D compile block |
| `cice.esm` | 0 | 0 | none — the assembly inherits the dynamics/transport compile block |

All 47 pre-probe coupling assertions reproduce the dumped values at the
tolerances documented per test: the imports/exports/aggregations at round-off
(`rel 1e-12`, every expected value tracing to a specific dump row (step, cell,
ncat) cited in the test description), the similarity exports (Si_tref/Si_qref
at `rel 1e-9`, machine precision) and the Rothrock strength (two dump rows at
`rel 1e-9`); Si_u10 carries a documented `rel 1e-4` gap: the uvelT_icep
reconstruction (ice_step_mod.F90:114-128, the uarea-weighted average of the
U-point uvel) reproduces the run's T-point velocity only to ~1e-4 m/s, and
Uref = |Uatm−U|·rd/rdn amplifies that difference by rd/rdn ≈ 8. The source of
the residual is unidentified (the run-side value is consumed inside
step_therm1, not dumped); the pin is kept at the observed precision rather
than tightening past it. The 22 dynamics-probe assertions (4 tests) reproduce
`cice_evp_stress.csv` and `cice_evp_momentum.csv` to round-off (`rel 1e-12`):
all probe inputs are bound from the dump rows, and a Python transcription of
the templates' laws reproduces every asserted column over all 184,320 stress
rows and 40,320 momentum rows (see "Transcription fixes found while writing
the tests" for the e_factor/epp2i/capping/Ktens fit results). The 22
thermodynamics assertions (14 tests) likewise reproduce `cice_albedo.csv`,
the `cice_thermo_state.csv` zTin/zTsn columns, and the `cice_thermo_flux.csv`
flwout/fsens (step 1, cell (3,5), `rel 1e-9`) and fswabs/fswthru (step 2,
cells (3,5) and (4,5), `rel 1e-12`) to the documented tolerances.

### Derived inputs (the raw CMEPS fields are not dumped)

The run dumps the fields *after* the import mapping (`cice_run_imports.csv`).
The raw CMEPS inputs are reconstructed by inverting the same law, and the .esm
must reproduce the dumped value:

- `Faxa_*` / `Fioo_q` = dumped import / `med2mod_areacor`. The mesh areacor ≈
  1.00033 (not 1), so the scaling genuinely bites (e.g. `flw = Faxa_lwdn·areacor`).
- `Faxa_rain` = dumped frain × 1000 / areacor (the kg/m²/s → m/s rhofresh
  conversion; pinned with the nonzero step-2 rain rows — step 1 is rain-free).
- `Sa_ptem` / `Sa_dens` = the dumped potT / rhoa. The dump rows have
  `potT = Tair` exactly, so the run used the Sa_ptem/Sa_dens branch; the
  Sa_pbot fallback is not exercised (inverting the dumped rhoa gives
  pbot = 99250.47 Pa, which would yield potT = 270.92 K, not the dumped
  270.34 K) and carries no dump-based pin.
- Per-category flux inputs for the export aggregation (fsensn, flatn, evapn,
  swabs, flwoutn, fhocnn, freshn, fsaltn, fswthrun) are not dumped — only the
  merged `cice_flux.csv` values are. The tests use the uniform per-category
  value that makes the aicen-weighted sum reproduce the dumped `exit` value
  (`flwoutn = exit flwout + (1−emissivity)·flw` for the outgoing-longwave
  emissivity correction). The albedo and state pins use real per-category data:
  `cice_albedo.csv` outputs × the coupling-time state of `cice_thermo_state.csv`
  entry(S+1) — in the nuopc driver the export at the end of step S runs after
  step S's dynamics (step_radiation at CICE_RunMod.F90:353, coupling_prep at
  :363), so exports(S) = entry(S+1) to round-off.
- `strairxT` / `strairyT` = dumped `Faii_taux` / `Faii_tauy` × aice / areacor
  (inversion of the scale_fluxes + export chain; the pre-scale T-point stress
  is not dumped).
- Similarity-test inputs: `Sa_u`/`Sa_v`/`Sa_ptem`/`Sa_shum`/`Sa_dens`/`Sa_z`
  are the dumped imports of `cice_run_imports.csv`; `uvelT`/`vvelT` are the
  uarea-weighted averages of the four surrounding U-point uvel/vvel of
  `cice_dyn_state.csv` (phase `d_entr`, the uvelT_icep reconstruction of
  ice_step_mod.F90:114-128 with the gx1v7 HTN/HTE cell-edge lengths);
  `Tsfcn`/`aicen_init` are the entry(S) per-category states of
  `cice_thermo_state.csv` (the surface temperature and the aicen saved before
  the thermodynamics, which merge_fluxes weights Trefn/Qrefn/Urefn by) and
  `aicen` is entry(S+1) (the aice of the scale_fluxes division).

### Transcription fixes found while writing the tests

1. **scale_fluxes `/aice` division missing.** `scale_fluxes`
   (cicecore/cicedyn/general/ice_flux.F90:1565-1592) divides the aggregated
   albedos, stresses, and fluxes by the ice area aice, so the coupled exports
   are per unit ice area; the skeleton's aggregation equations omitted the
   division. Dump-proven: Si_avsdr = 0.97303 = sum(albedo·aicen)/ailohi =
   0.967465/0.994260. Fixed in `cice_coupling.esm`: the 13 aggregation
   equations and Faii_taux/tauy now divide by `aice`, with the ice_flux.F90
   line citations.
2. **Faii_lwup ifelse arity.** The outgoing-longwave fix
   (ice_import_export.F90:1149-1164) was transcribed with a malformed
   3-argument `>` (only caught at test-compile time:
   `invalid_operator_arity: operator '>' takes exactly 2 argument(s)`); rebuilt
   as `ifelse(flwout > -1e-11, -stefan_boltzmann·(Tf+Tffresh)^4/ailohi, flwout)
   · mod2med_areacor`.
3. **Tf law.** The `sea_freezing_temperature` template transcribed the icepack
   mushy liquidus, but the coupled import path (CESMCOUPLED) computes
   `Tf = shr_frz_freezetemp(sss) = max(sss/(-18.48+0.01848·sss), -2.0)`
   (share/shr_frz_mod.F90.in:94-105, ice_import_export.F90:818-821; bugs.md
   #6). Updated the template; all 384 dumped rows match the coupled law to
   round-off (the −2 °C clamp is not exercised: the run's sss < 35.64 psu).
4. **Nuopc-cap constant overrides (bugs.md #11).** The CESM nuopc cap
   (`ice_comp_nuopc.F90:415-450`) overwrites ~20 icepack parameters with
   SHR_CONST values at init: cp_ice = 2117.27 (SHR_CONST_CPICE), cp_ocn =
   3996 (SHR_CONST_CPSW), cp_air = 1004.64 (SHR_CONST_CPDAIR), zvir =
   0.6077930728241565, Lsub = 2,834,700, Lfresh = 333,700, snowpatch =
   0.005, Tocnfrz = −1.836 (rhoi/rhow/rhos/gravit/vonkar/stefan/cp_wv/
   Tffresh unchanged). The first zTin transcription used the icepack
   defaults (cp_ice = 2106, Lfresh = 334000, …) and was off by 0.3–0.4%
   against the new dump columns; with the run's constants all 8 mushy layers
   reproduce the dumped zTin to round-off. The 16 affected templates, the
   coupling model's inline zvir/cp_air (eqs 68, 69, 75), and the `Control`
   Lheat binding (`cice.esm` eq 6, Lsub) now carry the run values, with the
   cap override cited in the template descriptions.
5. **Dynamics probe pins — findings.** (a) The `cice_evp_stress.csv` corner
   quantities are area-weighted (strain rate × corner area, m²/s); the dump
   header's `zetax2 [kg/s]` is the Fortran variable's nominal unit — the
   corner value is strength/max(Delta, DminTarea) = kg/(m²·s), and the area
   factors cancel in the stress update so `zetax2·divu − rep_prs` closes as
   N/m (the probe declares the honest kg/m²/s dimension, annotated). The
   D-class constants were verified by fitting, not assumed: e_factor =
   e_yieldcurve²/e_plasticpot⁴ = 0.25 reproduces the dumped Delta from
   divu/tension/shear on all 184,320 rows, epp2i = 1/e_plasticpot² = 0.25
   reproduces etax2/zetax2, and capping = 1/Ktens = 0 reproduces
   zetax2/rep_prs — no constant mismatch found. (b) The
   `cice_evp_momentum.csv` taux/tauy columns carry the **stepu partial
   stresses** `taux = vrel·waterx`, `tauy = vrel·watery` (ice_dyn_shared.F90:
   933-935, "NOTE this is not the entire ocn stress term"), not the
   dyn_finish full form `vrel·((uocn−uvel)·cosw − …)` of F90:1345-1348 that
   `cice_dynamics.esm` eqs 26-27 transcribe (the difference is the
   −vrel·uvel term); the probe transposes the stepu form, with the
   dyn_prep2 waterx/watery rotation (F90:817-818) and `sign(1,fm)` carried
   as the `ocn_hemis` parameter (bound from the dumped fm; unobservable in
   this run since sinw = 0). (c) `TbU = 0` in every traced row (no landfast
   grounding in the traced region), so the Cb pins check the wiring with a
   zero numerator; a nonzero-Cb pin needs a landfast cell from a later run.
6. **emissivity 0.985 → 0.95.** The `surface_heat_flux` template body (both
   the `ε·flw` absorption and the `ε·σ·TsfK⁴` emission terms) and the model's
   `flwoutn` equation transcribed the icepack default emissivity = 0.985
   (icepack_parameters.F90:156), but the D-class forcing_nml overrides it to
   0.95 — the value the coupling aggregation already assumed (the cap's
   flwout sum uses `(1−ε)·flw` = 0.05·flw, icepack_flux.F90:209-211).
   Dump-proven: with ε = 0.985 the step-1 (3,5) flwout aggregate misses the
   dump by rel 3.3e-3; with ε = 0.95 it reproduces it to 6e-11. Fixed in the
   template body and description, the model's flwoutn, and the model
   metadata (the template library's constants note already recorded the
   D-class 0.95 override).
7. **freshn/fsaltn corrected rate forms.** The skeleton's `freshn` carried
   the ktherm = 1 form `2·evapn − ρi·(congel+meltb−meltt) − fsnow − ρs·melts`
   (per-step, with /dt), but the F90:510-513 law cancels the sublimation and
   new-snowfall terms exactly: `freshn = evapn − (ρi·dhi + ρs·dhs)/dt` with
   `dhi = hin − worki` (which includes −evapin/ρi) and `dhs = hsn − works −
   hsn_new` (the hsn_new snowfall cancels), leaving the rate form
   `freshn = −(ρi·(congeln + meltbn + melttn) + ρs·meltsn)`, and
   `fsaltn = −ρi·dhi·ice_ref_salinity·0.001` with the same dhi (F90:513, 521;
   no separate /dt — the per-step dt divides out). Dump-proven: at step-1
   (3,5), `Σ aicen·freshn + frain·aice` reproduces the dumped fresh to
   3e-17 (the rain term is the "let rain drain through to the ocean" line of
   icepack_therm_itd.F90:2043, with the model's m/s frain × rhofresh), and
   `Σ aicen·fsaltn` reproduces the dumped fsalt exactly (ratio 1.0;
   evapin = 0 at step 1, all sublimation from snow). The per-category
   `dsnown = fsnow/ρs − meltsn + evapsn/ρs` rate law is added alongside.
8. **prep_radiation scale factor.** The step-S (S ≥ 2) thermo consumes the
   STORED step-S−1 dEdd outputs rescaled by
   `scale_factor = netsw(S)/netsw(S−1)` (icepack_shortwave.F90:3893-3922),
   with `netsw = Σ sw·(1 − al_band_ai)` computed from the grid-box-mean
   albedos `al_band_ai = Σ aicen·al_bandn` **per unit grid area** (the quirky
   aggregation of F90:3902-3922 — NOT divided by aice) and with netsw(S−1)
   the step's coupling_prep "Save net shortwave" value (the cice_flux exit
   scale_factor of CICE_RunMod.F90:618-624). Only the fswsfcn/fswintn/
   fswthrun (+ fswthrun_*) outputs are rescaled — the incoming bands and the
   albedos are not. The step-2 radiation pins reproduce the merged
   fswabs/fswthru to 1.5e-16 by feeding the model the step-1 bands × ratio
   with the unchanged step-1 albedos (the dEdd solve is linear in the input
   bands) and aggregating with the entry aicen(2) = albedo_dump(1) aicen.

### Coverage limits (why the other files have no tests)

- **Thermodynamics re-transcribed to the D-class configuration.** The first
   stage-2 skeleton transcribed the icepack namelist defaults (BL99 7/1,
   bubbly, ccsm3 albedos, cesm pond), which the run does not exercise: the
   D-class run (stage1/cime/cases/D_cice_inst/run/ice_in) used ktherm = 2
   (mushy) with nilyr = 8/nslyr = 3, conduct = MU71, shortwave = dEdd,
   frzpnd = hlid, kstrength = 1 (Rothrock), and atmbndy = similarity. The
   thermodynamics model and templates now transcribe the D-class laws, and the
   dEdd surface laws are dump-pinned. Pinned in this file: the dEdd surface
   laws (snow fraction, effective pond area, broadband albedos, radiation
   conservation), the mushy enthalpy↔temperature inversions `zTin =
   icepack_mushy_temperature_mush(qicen, zSin)` (eq 1) and `zTsn =
   snow_temperature_from_enthalpy(qsnon)` (eq 2) against the dumped
   `cice_thermo_state.csv` zTin_1..8/zTsn_1..3 columns, the step-1 surface
   energy-block aggregates flwout/fsens (test
   `thermo_surface_fluxes_step1_cell3_5`), and the step-2 radiation
   aggregates fswabs/fswthru through the prep_radiation rescale (tests
   `thermo_radiation_scale_step2_cell3_5`, `..._cell4_5`). Remaining
   unpinnable from the stage-1 dumps: the delta-Eddington 3-band solver
   itself (compute_dEdd_3bd, deferred — its outputs
   alvdrn/alidrn/alvdfn/alidfn, fswsfcn, fswthru_* and Iswabs are carried as
   parameters), the vertical conduction ODEs and the growth/melt rates
   (depend on the Tsfc implicit solve, see the divergences), and the merged
   ocean flux aggregates fhocn/fresh/fsalt/evap and fsurf/fcondtop (the
   step-1 fsurf residual is ~6e-8 and the step-2 surface fluxes carry the
   ~1e-6..1e-3 growth/melt residuals of the unpinned closures; the fresh
   rain term of icepack_therm_itd.F90:2043 is verified numerically to 3e-17
   but the freshn/fsaltn laws are pinned only as documented rate forms —
   their dump agreement at step 1 is exact by construction of the
   cancellation, see transcription-fix 7).
- **Tsfc tracer (dump bug, now fixed).** The first instrumented run read
   the `nt_Tsfc`/`nt_qice`/`nt_qsno` tracers with the j and tracer indices
   swapped (`trcrn(i,nt,j,n,blk)` vs the declared `(i,j,nt,n,blk)`), so
   Tsfc/Tsfcn/qice/qsno read 0 in every `cice_thermo_state.csv` row and
   Tsfcn/alvl/apnd/hpnd/ipnd read 0 in `cice_albedo.csv`. Fixed in the
   instrumented fork (`ced43eb`, superproject `c8a3bef`; rerun job
   11864874) — the columns now carry physical values (Tsfc ≈ −14.7 °C,
   qice_1 ≈ −3.03e8 J/m3) and the Si_t export is pinned against them
   (test `icecap_export_sit_step1_cell2_4`).
- **Spatial-D lowering blocks dynamics/transport/assembly tests.** The
   internal-stress divergence (strintx = D(sigma11,x) + D(sigma12,y)) and the
   transport continuity D(aicen·u) are spatial rewrite targets (esm-spec 4.2);
   an unlowered spatial D makes the whole-model compile fail
   (`unlowered_operator: rewrite-target operator D reached evaluation without
   being lowered`), even for `{0,0}` tests asserting only non-D variables, so
   `cice_dynamics.esm`, `cice_transport.esm`, and `cice.esm` carry `tests: []`.
   The dynamics **pointwise** laws are pinned in the coupling cap instead:
   the EVP stress maps and viscosities (`cice_evp_stress.csv`) and the stepu
   ocean-ice stress (`cice_evp_momentum.csv`), plus the Rothrock strength
   (`cice_ice_strength.csv`) — all inputs are dumped, so the pins bind them
   via `parameter_overrides` and assert the outputs to round-off. What stays
   deferred to the stage-3 discretization lowering: the spatial-D-dependent
   ODEs and stencils — the strain-rate corner stencils that feed
   divu/tension/shear from the velocity field, the stress divergence
   strintx/strinty, the momentum balance ODE, and the stress relaxation ODEs
   (the stepu implicit-solve machinery cca/ccb/ab2/cc1/cc2 of
   `cice_dynamics.esm` is documented there rather than transcribed).
   `cice_dyn_state.csv` continues to feed the similarity tests through the
   uvelT/vvelT reconstruction.
- **fswthru merge timing.** `merge_fluxes` weights fswthrun by the mid-thermo
   aicen (icepack_therm_vertical.F90:3169), which is not a dumped state; a
   per-category fswthru pin misses the dump by 0.1-18%, so Fioi_swpen uses the
   uniform construction. The per-category fswthrun of `cice_albedo.csv` stays
   available for the stage-3 dEdd tests.
- **Si_t, the Sa_pbot branch, and the lwup fix branch** are unpinned (Tsfc
   tracer not dumped; pbot branch not exercised by the run; flwout < −puny on
   the ice-covered cell so the outgoing-longwave substitution never fires).
- **Per-band Fioi_swpen_vdr/vdf/idr/idf** are dumped, but the stage-2 skeleton
   carries them as deferred zero parameters; `cice_budget.csv` has a single
   diagnostics row (hemispheric budgets are a stage-3 assembly concern).

`esm validate` passes on every file; `esm test stage2/` runs the 115 MOSART
assertions plus the 91 CICE assertions (69 coupling + 22 thermodynamics) with
0 failures and 0 errors.

## Scheme gaps (continuous ODE vs Fortran discrete)

- The Fortran's backward-Euler tridiagonal, the stepu implicit solve, and the
  EVP subcycle iteration become ODEs; the discrete-equivalent trajectories are
  matched at stage 3 through the integrator and the discretization templates.
- The incremental-remapping advection (departure points, transport integrals,
  the category-0 open-water mass) is a discretization, not a pointwise law;
  the skeleton carries its continuous continuity system.
- The ridge redistribution (ridge_shift, icepack_mechred.F90:1068+), the
  delta-Eddington shortwave solver, the salinity-profile module (zSin is a
  prescribed parameter with the 'constant' profile), the melt-pond
  tracers hpnd/apnd/alvl (prescribed parameters), and the wave/aerosol/isotope
  tracers are deferred. (The Rothrock strength law and the similarity atmo
  boundary layer are transcribed and pinned in the coupling cap; the assembly
  `Control` atmo scheme remains the constant branch.)
- The `Si_imask` ocean mask is a grid field (deferred on the single-cell
  skeleton).

## Next steps (stage 2.5 / stage 3)

1. Instrument the CICE run (stage 1) and extract per-component tuples
   (thermo column, EVP subcycle, remap transport, cap aggregation).
2. Fill the `tests` blocks; reconcile any transcription mismatches.
3. Submit the discretization templates (spatial-D lowering, incremental remap)
   to EarthSciDiscretizations; wire the Tsfc implicit solve; mount the full
   B-grid stencil.
4. Couple CICE with the other components in the top-level `eqesm.esm` model
   and integrate with `github.com/earthsciml/earthsciintegrator`.
