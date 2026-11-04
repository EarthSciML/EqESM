# CAM6 stage 2 — `.esm` implementation

Stage-2 implementation of the CAM6 atmospheric physics slice for the
EqESM project: the saturation specific humidity `qsat`, the energy
fixer, the Rayleigh friction, and the surface exchange (cam_out
exports + cam_in imports), as esm 1.3.0 files pinned by the
stage-1 instrumentation dumps of the SCAM reference run
(`stage1/dumps/cam/cam_qsat.csv`, `stage1/dumps/cam/cam_phys_tend.csv`,
`stage1/dumps/cam/cam_surface_exchange.csv`,
`stage1/dumps/cam/cam_surface_cflx.csv`).
Sources: `components/cam/src/atmos_phys/to_be_ccppized/wv_saturation.F90`
and `wv_sat_methods.F90`, `components/cam/src/atmos_phys/schemes/
conservation_adjust/check_energy/check_energy_fix.F90` +
`components/cam/src/physics/cam/check_energy.F90`,
`components/cam/src/atmos_phys/schemes/rayleigh_friction/
rayleigh_friction.F90`, and `components/cam/src/control/camsrfexch.F90`
(cam_export) + `components/cam/src/utils/cam_diagnostic_utils.F90`
(cpslec) of the instrumented cam6_4_212 fork
(`stage1/EqESM-cesm-cam`, branch `instrument/cam-dumps`); site
numbering per `stage1/runs/CAM_LANE.md` (energy_fixer = tphysbc site 1,
dump call physpkg.F90:2399; rayleigh = tphysac site 16, dump call
physpkg.F90:1757) and the boundary table in
`stage1/runs/cam_instrument.json`.

## Files

| file | content |
|---|---|
| `cam_templates.esm` | expression-template library of the shared pointwise saturation laws: `ice_weight`, `goff_gratch_water`, `goff_gratch_ice`, `svp_trans`, `svp_to_qsat`, `hltalt_no_ip`, `hltalt_mixed`, `tterm_transition`; plus the surface-exchange laws `density_dry_air`, `exner_surface`, `identity` |
| `cam_wv_saturation.esm` | the `estbl` function table (251 entries, Goff–Gratch mixed phase) + the `WvSaturation` model (the run's qsat law) + the `SvpTransLaw` model (pins the table-generation formula) |
| `cam_rayleigh_friction.esm` | the `RayleighFriction` model (otau profile + Euler-backward tendencies) |
| `cam_energy_fixer.esm` | the `EnergyFixer` model (uniform heat_glob tendency + implied surface flux) |
| `cam_surface_exchange.esm` | the `CamOutExport` model (cam_export's cam_out export laws: bottom-level state copies, dry-air density, Exner potential temperature, cpslec sea-level pressure, the precipitation split, the radiative and trace-gas pass-throughs) + the `CamInImports` model (the 20 cam_in coupler imports with their consumption sites, and the flx_net accumulation) |

## The qsat law (wv_saturation.F90:739-784, as this run evaluates it)

The run's `qsat_line` with ip = 0 evaluates the **tabulated** form:

```
es0 = estblf(T)                       ! linear interp of the estbl table,
                                      ! t_tmp = clamp(T - 127.16 K, 0, 248 K),
                                      ! flat clamps at both ends (F90:378-393)
qs  = svp_to_qsat(es0, p)             ! epsilo·es/(p - omeps·es), or 1 when
                                      ! p - es <= 0 (wv_sat_methods.F90:185-197)
es  = min(es0, p)
hltalt = calc_hltalt(T)               ! mixed-phase latent heat (F90:521-573)
tterm  = transition d(es)/dT term     ! pcf polynomial / ttrice (F90:548-573)
desdt = hltalt·es/(rh2o·T²) + tterm   ! d(es)/dT (F90:671-673)
dqsdt = qs·p·desdt/(es·(p - omeps·es)), 0 when qs == 1   (F90:676-679)
gam   = dqsdt·hltalt/cpair            ! the dumped gam (dqsdT), (F90:683-686)
```

**estbl**: `wv_sat_init` (F90:191-254) builds `plenest =
ceiling(tmax - tmin) + 2 = 251` entries at T = 127.16..377.16 K from
`svp_trans` (wv_sat_methods.F90:479-513): the Goff–Gratch water
(F90:569-580) and ice (F90:611-620) formulas blended by the ice weight
over the 20 K band below 273.15 K (pure ice below 253.15 K, pure water
at and above 273.15 K). The 251st knot (377.16 K) is never reachable
through `estblf` (T clamps at tmax = 375.16 K, so the last used entry
is 376.16 K); it is kept for table fidelity. In the esm the table is a
`function_tables` entry and `estblf` is `table_lookup` with
`interpolation: linear`, `out_of_bounds: clamp` on `min(T, 375.16 K)`
— the esm blend `tᵢ + w·(tᵢ₊₁ − tᵢ)` differs from the Fortran's
`(1−w)·tᵢ + w·tᵢ₊₁` by at most 1 ulp, far inside the 1e-12 test
tolerance. Constants (shr_const_mod.F90:24-55): rh2o = AVOGAD·BOLTZ/MWWV
= 461.5046398201599 J/(kg·K), epsilo = MWWV/MWDAIR =
0.6219705862045155, omeps = 1 − epsilo = 0.3780294137954845, cpair =
1.00464e3 J/(kg·K), latvap = 2.501e6, latice = 3.337e5 J/kg, ttrice =
20 K, pcf = [5.04469588506e-01, −5.47288442819e+00,
−3.67471858735e-01, −8.95963532403e-03, −7.78053686625e-05] (Pa,
Pa/K, …, Pa/K⁴).

The `SvpTransLaw` model pins the analytic generation formula itself
(`goff_gratch_water`, `goff_gratch_ice`, `svp_trans` of
`cam_templates.esm`) at five temperatures spanning the pure-ice,
transition, and water branches, including the estbl first and last
reachable knots (T = 127.16 and 375.16 K) as table-consistency pins.
The Goff–Gratch curve-fit coefficients are dimensionless as written
(the fit is dimensionally incoherent: a log10 of an hPa-valued number,
exponentiated and ×100); the ×100 factor carries the Pa dimension so
the analytic laws read as Pa.

## The energy fixer (check_energy_fix.F90:37-44)

```
ptend_s = heat_glob                     ! uniform dry-energy tendency, J/kg/s
eshflx  = heat_glob·(ps − ptop)/gravit  ! = heat_glob·(pint(pver+1) − pint(1))/gravit
```

`heat_glob` is the global-mean heating from `check_energy_cam_fix`
(check_energy.F90:637-702); in the SCAM run it is 0 because the IOP
forcing (ARM97_4scam.nc) is not CAM-generated (`use_camiop =
.false.`), so the dumped `ptend_s` of the chkenerg rows are all 0. `ps`
(state%ps) is dumped; `ptop` = pint(1) is **not** dumped and is derived
from the dump rows as `ps − Σ pdel = 225.5239523947239 Pa`, constant
across the window (documented in the model). The nonzero-branch pin
(heat_glob = 0.25 J/kg/s → eshflx = 2361.3585380814457 W/m²) is a
compiled-Fortran reference (see below).

## The Rayleigh friction (rayleigh_friction.F90)

```
krange = raykrange, or (rayk0 − 1)/2 when raykrange = 0     (F90:89-92)
otau(k) = otau0·(1 + tanh((rayk0 − k)/krange))/2,  otau0 = 1/(86400·raytau0)
c2 = 1/(1 + otau·ztodt);  c1 = −otau·c2;  c3 = 0.5·(1 − c2²)/ztodt
dudt = c1·u;  dvdt = c1·v;  dsdt = c3·(u² + v²)              (F90:149-156)
```

When `raytau0 = 0` the run returns zeros before the profile (F90:113-
115); the model reproduces that branch exactly (otau = 0 → c2 = 1, c1 =
c3 = 0). The dump window runs raytau0 = 0, so the dump pins exercise
the zero branch; the enabled-branch pins are compiled-Fortran
references. Note the glibc tanh saturation: with the default profile
(krange = 0.5) the argument (2 − k)/0.5 reaches −20 at k = 12 and
tanh(−20) = −1 exactly, so otau and the tendencies vanish exactly at
k ≥ 12 even when enabled — pinned as `enabled_krange0_k12_saturation`;
a raykrange = 10 pin set exercises nonzero damping of the O(1) winds.
The esm c3 writes `/ztodt` where the Fortran multiplies by rztodt =
1/ztodt (≤ 1 ulp).

## The surface exchange (camsrfexch.F90 cam_export, cam_in)

`cam_surface_exchange.esm` covers the coupler-facing boundary of the
physics. Model `CamOutExport` transcribes the cam_out export laws of
`cam_export` (camsrfexch.F90:405-577):

```
tbot   = t(pver)                          ! state copies (F90:495, 497, 499-507)
zbot   = zm(pver);  ubot/vbot/qbot/pbot = u/v/q/pmid(pver)
rho    = pbot/(rair·tbot)                 ! dry-air density (F90:503),
                                          ! rair = AVOGAD·BOLTZ/MWDAIR (shr_const)
thbot  = tbot·(ps/pmid)^(rair/cpair)      ! exner (F90:496, dp_coupling.F90:731-733)
psl    = cpslec(t, pmid, ps, phis)        ! ECMWF sea-level pressure
                                          ! (cam_diagnostics.F90:1080-1081,
                                          !  cam_diagnostic_utils.F90:20-81)
precc  = max(0, prec_dp + prec_sh)        ! pbuf PREC/Snow fields (F90:537-575)
precl  = max(0, prec_sed + prec_pcw)
precsc = min(precc, max(0, snow_dp + snow_sh))
precsl = min(precl, max(0, snow_sed + snow_pcw))
co2diag = co2vmr·1e6 (ppmv)               ! chem_surfvals CO2VMR (F90:506-511)
co2prog = 0                               ! co2_transport() false here (F90:512-516)
ozone = SRFOZONE pbuf;  netsw/flwds/soll/sols/solld/solsd = radiation outputs
```

Laws vs pass-through imports:

| cam_out field | status in the esm |
|---|---|
| tbot, zbot, ubot, vbot, qbot, pbot | **laws**: direct copies of the bottom-level state (F90:495-507) |
| rho, thbot | **laws**: `density_dry_air` and `exner_surface` templates (cam_templates.esm), constants SHR_CONST_RDAIR/CPDAIR |
| psl | **law**: the cpslec chain (see below) |
| precc/precl/precsc/precsl | **law**: the pbuf-field split + clamps (F90:537-575), unpinned (see gaps) |
| co2diag | **law** co2diag = CO2VMR·1e6 with the input **imported** (chem_surfvals); provenance open question (see gaps) |
| co2prog | **law**: identically 0 in this configuration (co2_transport off) |
| ozone, netsw, flwds, soll, sols, solld, solsd | **pass-through imports** from the radiation component (rrtmg/radiation.F90:1272-1316,1429) and the SRFOZONE pbuf — modeled as `identity` passes and pinned as bookkeeping |

The psl law runs in `diag_phys_writeout` **before** this step's
radiation update (verified bit-exact against the dumped psl using the
pre-radiation state), so the model takes two bottom-level temperatures:
`t_bot` (at cam_export) and `t_psl` (at psl time). The dump pins
reconstruct `t_bot` as T + ptend_s·dt/cpairv with cpairv = cpair =
1004.64 J/(kg·K) (this run has `dry_air_species = ''`; verified bit-
exact against the dumped tbot at steps 1/2/24/48).

Model `CamInImports` documents the 20 cam_in fields as **pure coupler
imports** (camsrfexch.F90:86-125) with their physics consumption sites:
albedos asdir/asdif/aldir/aldif → SW radiation (physpkg.F90:3086-3087);
lwup → LW surface skin temperature (rrtmg_state.F90:115); shf/lhf/cflx
→ qneg4 (physpkg.F90:1572-1574) and the vdiff/CLUBB energy checks
(1768-1769, 2775-2779); wsx/wsy → CLUBB surface momentum flux
(clubb_intr.F90:3902, 3916-3917); ts/sst → CLUBB (2742-2744, 2637-2698);
landfrac → convect_deep (2493); tref/qref → diagnostics only
(cam_diagnostics.F90:1801-1811); u10 → CARMA dust/sea-salt emissions
(off in this run); snowhland/snowhice/icefrac/ocnfrac → CLUBB/CARMA.
The per-constituent cflx array (34 constituents, all coupler-provided)
is dumped in cam_surface_cflx.csv and documented on the `cflx_Q`
parameter. CAM has **no law** relating cam_in%lhf to cam_in%cflx: the
coupler computes both consistently (dumped lhf/cflx_Q = LATVAP exactly
on every row), so the esm keeps them as independent inputs — the
consistency is a coupler property, not implemented. The one CAM-side
combination of coupler fluxes is the flx_net accumulation of tphysac
(physpkg.F90:1533):

```
flx_net = net_flx_rad + shf
         + (precc + precl)·latvap·rhoh2o + (precsc + precsl)·latice·rhoh2o
```

with net_flx_rad the radiation-side net surface flux set by
radiation_tend (physpkg.F90:3090-3092, radheat.F90:348-364); neither
term is dumped (gap).

## Testing

`./esm validate stage2/cam_*.esm` and `./esm test stage2/` (the full
stage-2 suite). All tests are instantaneous (`time_span {0,0}`), rel
1e-12, citing the dump file and row.

| model | tests | pins |
|---|---|---|
| `WvSaturation` | 10 | cam_qsat.csv rows spanning both phases and steps 1/24/48, T = 209.0..301.9 K, p = 364..92152 Pa: 4 pure-ice (T < 253.15), 3 transition (253.15-273.15), 3 water (≥ 273.15); es/qs/gam expected = dumped values, hltalt/tterm = Fortran-formula values (validator replica, reproduces the dump to 1.8e-12 Pa on es) |
| `SvpTransLaw` | 5 | analytic formula at 127.16/220.16/260.16/300.16/375.16 K (ice, transition, water branches; knots 1 and 250 of estbl) |
| `RayleighFriction` | 9 | 3 dump pins (site rayleigh, steps 1/48: all-zero branch), 3 compiled-Fortran pins (enabled, default profile: k = 1, 11, tanh saturation), 3 compiled-Fortran pins (enabled, raykrange = 10: k = 12, 16, 32) |
| `EnergyFixer` | 3 | 2 dump pins (site energy_fixer, steps 1/48: heat_glob = 0), 1 compiled-Fortran pin (heat_glob = 0.25 J/kg/s) |
| `CamOutExport` | 9 | 4 dump pins of the state laws + psl (steps 1/2/24/48, bc_out, t_bot reconstructed from the radiation-site rows), 4 dump pins of the pass-throughs (incl. the stale step-2 values, the night step, and the albedo-1.0 step), 1 cpslec-chain pin (branch 2, with the five intermediates) |
| `CamInImports` | 0 | all fields are coupler imports; the one law (flx_net) needs net_flx_rad, which is not dumped → gap (no data invented) |

Total for the CAM slice: 65 + 27 + 6 + 68 = 166 assertions. The compiled-
Fortran references compile the actual subroutines verbatim from the
fork with gfortran 11.5 in a driver at the dump's step-1 u/v state
(`scratch/cam_laws_driver/` in the worktree, the mosart-precedent
style; the driver source and all reference output files are kept
there). The stage-1 validator (`stage1/analysis/validate_cam_dumps.py`)
independently reproduces the qsat dump to round-off, and the dump
window was verified solution-neutral by the stage-1 lane.

`./esm units stage2/cam_*.esm` reports all equations consistent except
the single `es0 = table_lookup` equation of `WvSaturation`, which the
checker cannot dimension (a `table_lookup` node has no determinable
dimension — inherent to the checker, documented); no
`units_dimensionless_literal_mismatch` anywhere.

## Coverage gaps

- **Saturated branches** (p − es ≤ 0 → qs = 1, es = min(es,p), dqsdt =
  0) and the **estblf clamps** (T < 127.16 K, T > 375.16 K): the dump
  window never saturates (es < p on every row) and spans only T =
  208.7..302.5 K; the branches are implemented per the law but not
  pinned by data.
- **pint(1) of the energy fixer** is not dumped; supplied as the
  dump-derived ptop = ps − Σpdel (constant 225.5239523947239 Pa across
  the window). No data was invented.
- **rayleigh enabled branch**: the run's namelist has raytau0 = 0; the
  enabled-branch pins are compiled-Fortran references at the dump's
  wind state, not dump pins.
- **gam/dqsdT**: only `gam` is dumped; `dqsdt` and `desdt` are
  unpinned intermediates of the pinned chain.
- The esm c3 blend (`/ztodt` vs Fortran `·1/ztodt`) and the interp
  blend order differ by ≤ 1 ulp (documented above).
- **zbot**: zm(pver) at cam_export time is not dumped. `physics_update`
  recomputes zm after the temperature update (physics_types.F90:423-449),
  so the pre-radiation zm of the dump differs from what cam_export sees;
  reconstructing it needs the geopotential law (geopotential_temp.F90)
  plus the condensate mixing ratios at the radiation site (not dumped).
  The esm pins zbot only as a pass-through bookkeeping pin (input = the
  dumped zbot) — the geopotential law itself is a separate component.
- **precipitation split**: the eight pbuf fields (PREC_DP, SNOW_DP,
  PREC_SH, SNOW_SH, PREC_SED, SNOW_SED, PREC_PCW, SNOW_PCW) are not
  dumped, so the precc/precl/precsc/precsl laws are unpinned. The dump
  shows nonzero precc/precl on 15 of the 48 window steps, so the split
  is exercised by the run but its inputs are not observable; the laws
  are implemented per camsrfexch.F90:537-575.
- **cpslec branches**: all 48 window steps take branch 2 (Tstar >
  290.5 K and TT0 > 290.5 K → alph = 0, Tstar ← 0.5·(290.5+Tstar)). The
  flat branch (|phis/gravit| < 1e-4 → psl = ps), branch 1 (Tstar ≤
  290.5, TT0 > 290.5) and branch 3 (incl. the Tstar < 255 K sub-branch)
  are implemented per the Fortran but unpinned — the 887 m SCAM column
  never exercises them.
- **flx_net / net_flx_rad**: neither tend%flx_net nor the radiation
  side's net_flx is dumped, so the CamInImports accumulation law is
  unpinned. No data was invented (inputs would be fabricated).
- **co2diag provenance (open question, not a CESM defect)**: the dumped
  co2diag ≈ 1.127278 ppm (drifting slowly across the window) could not
  be traced to any prescribed input: the LBC climo file
  (LBC_2000climo_CMIP6_0p5degLat_c180227.nc) implies a global-mean
  CO2VMR of ~369 ppm, the IOP file has no co2vmr variable, and no
  namelist sets co2vmr. chem_surfvals_get('CO2VMR') is treated as a
  pass-through import with the discrepancy documented; resolving it
  would need a stage-1 dump of the chem_surfvals value at the cam_export
  site.
- **lhf = cflx·LATVAP** holds in the dump on every row but is a coupler
  (aoflux) property, not a CAM law — deliberately not implemented.

## Bugs found in the CAM code

None. (The 251-entry table with a never-used last knot is faithful
Fortran behavior, not a bug; the tanh saturation is glibc behavior.)
