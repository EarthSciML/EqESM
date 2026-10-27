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

Default configuration transcribed (namelist defaults): `shortwave = 'dEdd'`
(albedo/shortwave laws of the ccsm3 branch — the delta-Eddington solver itself
is deferred), `kdyn = 1` (EVP), `grid_ice = 'B'`, `ndtd = 1`, `ndte = 120`,
`revised_evp = .false.`, `elasticDamp = 0.36`, `Ktens = 0`,
`e_yieldcurve = e_plasticpot = 2`, `capping = 1`, `u0 = 5e-5`, `cosw = 1`,
`sinw = 0`, `kstrength = 1` (Rothrock ridging deferred; the Hibler strength is
transcribed), `atmbndy = 'similarity'` (deferred; the constant scheme is
transcribed), `ktherm = 1` (BL99), `conduct = 'bubbly'`,
`tfrz_option = 'mushy'`, `saltflux_option = 'constant'`, `frzpnd = 'cesm'`,
`ncat = 5`, `nilyr = 7`, `nslyr = 1`, `l_brine = .true.`, `dt = 3600 s`.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `cice_templates.esm` | Expression-template library (41 templates): snow cover / albedo / absorbed-penetrating shortwave | `icepack_shortwave.F90`: compute_albedos 466-597, absorbed_solar 701-896 |
| | saturation humidity, surface heat flux, liquidus (linear + mushy), enthalpy↔temperature (BL99 brine quadratic, fresh, snow) | `icepack_therm_shared.F90`: surface_heat_flux 147-165, liquidus 383-418; `icepack_mushy_physics.F90`: 227-261; calculate_Tin_from_qin 68-101; snow temperature 470-484 |
| | conductivity (bubbly Pringle 2007 / MU71) and interface conductances, congelation enthalpy, ocean heat flux, melt-pond volume change | `icepack_therm_bl99.F90`: conductivity 874-961; `icepack_therm_vertical.F90`: qbot 1432-1443, fbot 662-678; `icepack_meltpond_lvl.F90`: 154-174 |
| | constant atmo boundary layer (tau, shcoef, lhcoef) | `icepack_atmo.F90`: atmo_boundary_const 398-496 (460-464, 493-494) |
| | Hibler strength, Delta invariant, viscosities/replacement pressure, classic-EVP stress maps, ocean-ice stress, seabed stress | `icepack_mechred.F90`: 1720; `ice_dyn_shared.F90`: strain_rates 2081-2161, visc_replpress 2444-2473, set_evp_parameters 451-484; stepu 923-964 |
| `cice_thermodynamics.esm` | Model `Thermodynamics` — prep_radiation shortwave/albedo, step_therm1 vertical heat conduction + surface energy balance, step_therm2 growth/melt/ponds and ocean heat/fresh/salt fluxes; 67 equations over ncat × nilyr/nslyr | `icepack_step_mod.F90`/`icepack_therm_bl99.F90` (tridiagonal 1038-1210 as semi-discrete flux-divergence ODEs), `icepack_therm_vertical.F90` (thickness changes, freshn/fsaltn 510-522, fhocnn 1625-1630, sublimation 1355-1403, 1502-1559) |
| `cice_dynamics.esm` | Model `Dynamics` — B-grid momentum balance (continuous limit of stepu) and classic-EVP stress relaxation ODEs, Hibler strength, viscosities; 34 equations over one cell | `ice_dyn_shared.F90`: stepu 923-964, dyn_prep2 815-833, dyn_finish 1345-1348, Tmass 540; `ice_dyn_evp.F90`: stress 1584-1609 |
| `cice_transport.esm` | Model `Transport` — per-category conservative advection d(a)/dt + div(a·u) = 0 (and volume forms) as spatial-D rewrite targets; enthalpy-tracer statement documented, deferred to stage 3 | `ice_transport_driver.F90`: 251-700; `ice_transport_remap.F90`: 317-831 |
| `cice_coupling.esm` | Model `IceCap` — med↔ice cap renames, areacor scaling, Sa_pbot potT/rhoa branch, export aggregation (ailohi, Tsrf, merge_fluxes sums), Faii_lwup fix; 62 equations | `ice_import_export.F90`: ice_import 414-700, ice_export 1030-1317; `icepack_flux.F90`: merge_fluxes 180-260 |
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
| Constants | `puny = 1e-11`, `Tffresh = 273.15`, `rhoi = 917`, `rhos = 330`, `rhow = 1026`, `cp_ice = 2106`, `cp_ocn = 4218`, `cp_air = 1005`, `Lfresh = 3.34e5`, `Lsub = 2.835e6`, `Lvap = 2.501e6`, `depressT = 0.054`, `emissivity = 0.985`, `stefan_boltzmann = 5.67e-8`, `gravit = 9.80616` | (icepack_parameters.F90:36, 82, 62, 61, 65, 74-75, 281, 83-84, 77, 156, 159, 56) |

**Unit notes**

- `"C"` parses as Coulomb in the esm unit checker, so all Celsius fields use
  `degC` (identical to K in the checker's algebra).
- `lhcoef` is W/m² (per specific-humidity difference, not per K): the Fortran
  `lhcoef = latncoef·Lheat·rhoa·wind` (icepack_atmo.F90:494) carries the latent
  heat inside the coefficient.
- `evapn` requires the ρ in the numerator: the Fortran
  `evapn += dhs·rhos` (icepack_therm_vertical.F90:1502-1521) is transcribed as
  the continuous rate `evapn = esub·rhos/((qsn − rhos·Lvap)·dt)`; `fsaltn`
  divides by dt (`fsaltn = −rhoi·dhi·ice_ref_salinity·0.001/dt`, F90:521).
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
| Ridging (kstrength = 1 Rothrock) | TBD | — deferred (the Hibler kstrength ≠ 1 strength is transcribed so the dynamics closes) |
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
| `ice_temperature_from_enthalpy` | brine quadratic `cp_ice·T² + ((cp_ocn−cp_ice)·Tmltk − qin/ρi − L)·T + L·Tmltk = 0` solved by the negative root (F90:68-101) | `ifelse(l_brine > 0.5, min((−bb1−√(bb1²−4·cp_ice·cc1))/(2·cp_ice), Tmltk), (L+qin/ρi)/cp_ice)`; the Fortran abort on a negative discriminant is a noted guard |
| `ice_conductivity_bubbly` | `k = max((2.11 − 0.011·T + 0.09·S/min(−puny,T))·ρi/917, kimin)` (icepack_therm_bl99.F90:918-924) | same with `min(−1e-11, T)` |
| `congelation_enthalpy` | `qbot = −ρi·(cp_ice·(Tmlts−Tbot) + L·(1−Tmlts/Tbot) − cp_ocn·Tmlts)`, capped at `−0.5·ρi·L` (icepack_therm_vertical.F90:1432-1443) | same, `ifelse` brine/fresh, `min` cap |
| `ocean_heat_flux` | `fbot = max(−cp_ocn·rhow·0.006·max(sst−Tf,0)·ustar, frzmlt)` (F90:662-678) | same; ustar floored at 0.005 in the model |
| `meltpond_cesm_dvn` | `dvn = rfrac/rhofresh·(meltt·ρi + melts·ρs + frain·dt)·aicen − volpn·(1−exp(rexp·dTs/Tp))` (icepack_meltpond_lvl.F90:154-174) | same with `Tp = −2`, `rexp = 0.01`, `dTs = max(−2−Tsfcn, 0)` |
| `atmo_constant_tau/shcoef/lhcoef` | `tau = ρa·0.0012·wind`; `shcoef = 0.0012·1005·ρa·wind`; `lhcoef = 0.0015·Lheat·ρa·wind` (icepack_atmo.F90:460, 493-494) | identical, `Lheat = Lsub = 2.835e6` for ice |
| `ice_strength_hibler` | `strength = Pstar·vice·exp(−Cstar·(1−aice))` (icepack_mechred.F90:1720) | identical (`Pstar = 2.75e4`, `Cstar = 20`) |
| `viscosity_scale` | `tmpcalc = capping·P/max(Δ,Dmin) + (1−capping)·P/(Δ+Dmin)` (ice_dyn_shared.F90:2467-2468) | identical; `Dmin = deltaminEVP = 1e-11` |
| `evp_stressp_update` | `stressp_new = (stressp·(1−arlx1i·revp) + arlx1i·(zetax2·divu − rep_prs))·denom1` (ice_dyn_evp.F90:1584-1591) | transcribed as the relaxation ODE `D(stressp)/t = (stressp_target − stressp)/Td` with `Td = 2·elasticDamp·dt` (see divergences) |

## Documented divergences (all annotated in file metadata)

1. **Heat conduction as semi-discrete ODEs.** The Fortran assembles and solves
   an implicit backward-Euler tridiagonal for Tsfc and the layer temperatures
   (icepack_therm_bl99.F90:1038-1210); the skeleton writes the equivalent
   energy-conserving method-of-lines form
   `d(q_k)/dt = (kh_up·(T_up−T_k) + kh_dn·(T_dn−T_k) + Iswabs_k)/hilyr`
   (volumetric-enthalpy balance, no ρ·cp factor), with the interface
   conductances kh of the conductivity subroutine (F90:874-961) and the
   per-layer absorbed shortwave of absorbed_solar (icepack_shortwave.F90:
   862-880). The stage-3 integrator reproduces the backward-Euler trajectory.
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
5. **Hibler strength in place of Rothrock.** The default `kstrength = 1`
   Rothrock potential-energy form (icepack_mechred.F90:1666-1712) depends on
   the ridge distribution state; the `kstrength ≠ 1` Hibler form (F90:1720) is
   transcribed so the dynamics closes.
6. **Constant atmo scheme in place of similarity.** The default
   `atmbndy = 'similarity'` Monin-Obukhov iteration (icepack_atmo.F90:200-444)
   is deferred; the constant scheme (F90:398-496) is transcribed. Consequently
   the reference-height exports Si_u10/Si_tref/Si_qref are zero (the constant
   scheme does not compute Tref/Qref/Uref).
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

## Testing

The stage-2 `tests` blocks pin the stage-1 instrumentation dumps
(`stage1/dumps/cice/`) with instantaneous laws (`time_span {0,0}`, MOSART
pattern). Run `./esm test stage2/` for the full suite (115 MOSART + CICE):

| file | tests | assertions | what is pinned |
|---|---|---|---|
| `cice_coupling.esm` | 29 | 32 | the import laws (areacor scaling, frain/rhofresh, potT/rhoa branch wiring), the freezing-temperature law, the export state fields, the per-category albedo aggregation, the scale_fluxes `/aice` flux aggregation, and the air-ice stress exports |
| `cice_thermodynamics.esm` | 0 | 0 | none — the run's configuration (ktherm = 2 mushy, MU71, dEdd, nilyr = 8/nslyr = 3, hlid, similarity) exercises none of the transcribed laws (see coverage limits) |
| `cice_dynamics.esm` | 0 | 0 | none — unlowered spatial-D stencils block compilation of any test (see coverage limits) |
| `cice_transport.esm` | 0 | 0 | none — same spatial-D compile block |
| `cice.esm` | 0 | 0 | none — the assembly inherits the dynamics/transport compile block |

All 32 coupling assertions reproduce the dumped values to round-off
(`rel 1e-12`); every expected value traces to a specific dump row (step, cell,
ncat) cited in the test description.

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

### Coverage limits (why the other files have no tests)

- **Thermodynamics configuration mismatch.** The run
   (stage1/cime/cases/D_cice_inst/run/ice_in) used ktherm = 2 (mushy) with
   nilyr = 8/nslyr = 3, conduct = MU71, shortwave = dEdd, frzpnd = hlid,
   kstrength = 1 (Rothrock), and atmbndy = similarity (the default); the
   skeleton deliberately transcribes the icepack defaults (BL99 7/1, bubbly,
   ccsm3 albedos, cesm pond, constant atmo, Hibler strength). The dEdd
   albedos of `cice_albedo.csv` differ from the ccsm3 analytic laws by up to
   ~10% (ncat 2, step 2: alvdfn = 0.9527 vs 0.8598), the Rothrock strength
   P = 3.21e5 differs from the Hibler value 5.47e5, and the similarity scheme
   gives nonzero Si_u10/Si_tref/Si_qref (1.104, 266.36, 1.95e-3) where the
   transcribed constant scheme gives zero — so no dump value pins the
   transcribed thermodynamics equations.
- **Tsfc tracer not dumped.** The nt_Tsfc tracer reads 0 in every
   `cice_thermo_state.csv` row (aggregate and per-category), so the Si_t
   export and any surface-temperature law are unpinnable; the real value is
   visible only in the exports (Si_t = 258.49 K at step 0, Tsfc ≈ −14.7 °C).
- **Spatial-D lowering blocks dynamics/transport/assembly tests.** The
   internal-stress divergence (strintx = D(sigma11,x) + D(sigma12,y)) and the
   transport continuity D(aicen·u) are spatial rewrite targets (esm-spec 4.2);
   an unlowered spatial D makes the whole-model compile fail
   (`unlowered_operator: rewrite-target operator D reached evaluation without
   being lowered`), even for `{0,0}` tests asserting only non-D variables, so
   `cice_dynamics.esm`, `cice_transport.esm`, and `cice.esm` carry `tests: []`.
   The EVP momentum/stress and ice-strength dumps (`cice_evp_momentum.csv`,
   `cice_evp_stress.csv`, `cice_ice_strength.csv`, `cice_dyn_state.csv`,
   `cice_transport_state.csv`) will be pinned once the EarthSciDiscretizations
   lowering lands (stage 3).
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
assertions plus the 32 CICE coupling assertions with 0 failures and 0 errors.

## Scheme gaps (continuous ODE vs Fortran discrete)

- The Fortran's backward-Euler tridiagonal, the stepu implicit solve, and the
  EVP subcycle iteration become ODEs; the discrete-equivalent trajectories are
  matched at stage 3 through the integrator and the discretization templates.
- The incremental-remapping advection (departure points, transport integrals,
  the category-0 open-water mass) is a discretization, not a pointwise law;
  the skeleton carries its continuous continuity system.
- The ridging (kstrength = 1), the delta-Eddington shortwave solver, the
  similarity atmo boundary layer, the salinity-profile module (zSin is a
  prescribed parameter with the 'constant' profile), the melt-pond
  tracers hpnd/apnd/alvl (prescribed parameters), and the wave/aerosol/isotope
  tracers are deferred.
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
