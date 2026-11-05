# CAM6 stage 2 — `.esm` implementation

Stage-2 implementation of the CAM6 atmospheric physics slice for the
EqESM project: the saturation specific humidity `qsat`, the energy
fixer, and the Rayleigh friction, as esm 1.3.0 files pinned by the
stage-1 instrumentation dumps of the SCAM reference run
(`stage1/dumps/cam/cam_qsat.csv`, `stage1/dumps/cam/cam_phys_tend.csv`).
Sources: `components/cam/src/atmos_phys/to_be_ccppized/wv_saturation.F90`
and `wv_sat_methods.F90`, `components/cam/src/atmos_phys/schemes/
conservation_adjust/check_energy/check_energy_fix.F90` +
`components/cam/src/physics/cam/check_energy.F90`, and
`components/cam/src/atmos_phys/schemes/rayleigh_friction/
rayleigh_friction.F90` of the instrumented cam6_4_212 fork
(`stage1/EqESM-cesm-cam`, branch `instrument/cam-dumps`); site
numbering per `stage1/runs/CAM_LANE.md` (energy_fixer = tphysbc site 1,
dump call physpkg.F90:2399; rayleigh = tphysac site 16, dump call
physpkg.F90:1757) and the boundary table in
`stage1/runs/cam_instrument.json`.

## Files

| file | content |
|---|---|
| `cam_templates.esm` | expression-template library of the shared pointwise saturation laws: `ice_weight` (clamp form `min(1, max(0, (tmelt − t)/ttrice))`, bit-identical to the Fortran's three-branch if), `goff_gratch_water`, `goff_gratch_ice`, `svp_trans`, `svp_to_qsat`, `hltalt_no_ip`, `hltalt_mixed`, `pcf_transition_poly` (the shared pcf(5) Horner polynomial in tc), `tterm_transition`, and the zero-parameter constant fragments `tmelt`/`tboil`/`ttrice` (applied with empty bindings; TKTRIP = 273.16 K stays inline in `goff_gratch_ice`) |
| `cam_wv_saturation.esm` | the `estbl` function table (251 entries, Goff–Gratch mixed phase) + the `WvSaturation` model (the run's qsat law) + the `SvpTransLaw` model (pins the table-generation formula) |
| `cam_rayleigh_friction.esm` | the `RayleighFriction` model (otau profile + Euler-backward tendencies) |
| `cam_energy_fixer.esm` | the `EnergyFixer` model (uniform heat_glob tendency + implied surface flux) |

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

Total for the CAM slice: 65 + 27 + 6 = 98 assertions. The compiled-
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

## Bugs found in the CAM code

None. (The 251-entry table with a never-used last knot is faithful
Fortran behavior, not a bug; the tanh saturation is glibc behavior.)
