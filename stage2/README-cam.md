# CAM6 stage 2 — `.esm` implementation

Stage-2 implementation of the CAM6 atmospheric physics slice for the
EqESM project: the saturation specific humidity `qsat`, the energy
fixer, the Rayleigh friction, the surface exchange (cam_out exports +
cam_in imports), and the vertical diffusion (the implicit surface-stress
solve of `diffusion_solver.F90`), as esm 1.3.0 files pinned by the
stage-1 instrumentation dumps of the SCAM reference run
(`stage1/dumps/cam/cam_qsat.csv`, `stage1/dumps/cam/cam_phys_tend.csv`,
`stage1/dumps/cam/cam_vdiff_prof.csv`, `stage1/dumps/cam/cam_vdiff_q.csv`,
`stage1/dumps/cam/cam_surface_exchange.csv`,
`stage1/dumps/cam/cam_surface_cflx.csv`).
Sources: `components/cam/src/atmos_phys/to_be_ccppized/wv_saturation.F90`
and `wv_sat_methods.F90`, `components/cam/src/atmos_phys/schemes/
conservation_adjust/check_energy/check_energy_fix.F90` +
`components/cam/src/physics/cam/check_energy.F90`,
`components/cam/src/atmos_phys/schemes/rayleigh_friction/
rayleigh_friction.F90`, `components/cam/src/atmos_phys/schemes/
vertical_diffusion/vertical_diffusion_solver.F90` +
`diffusion_solver.F90`, and `components/cam/src/control/camsrfexch.F90`
(cam_export) + `components/cam/src/utils/cam_diagnostic_utils.F90`
(cpslec) of the instrumented cam6_4_212 fork
(`stage1/EqESM-cesm-cam`, branch `instrument/cam-dumps`); site
numbering per `stage1/runs/CAM_LANE.md` (energy_fixer = tphysbc site 1,
dump call physpkg.F90:2399; rayleigh = tphysac site 16, dump call
physpkg.F90:1757; vdiff = site vdiff) and the boundary table in
`stage1/runs/cam_instrument.json`.

## Files

| file | content |
|---|---|
| `cam_templates.esm` | expression-template library of the shared pointwise laws: the saturation laws `ice_weight` (clamp form `min(1, max(0, (tmelt − t)/ttrice))`, bit-identical to the Fortran's three-branch if), `goff_gratch_water`, `goff_gratch_ice`, `svp_trans`, `svp_to_qsat`, `hltalt_no_ip`, `hltalt_mixed`, `pcf_transition_poly` (the shared pcf(5) Horner polynomial in tc), `tterm_transition`, and the zero-parameter constant fragments `tmelt`/`tboil`/`ttrice` (applied with empty bindings; TKTRIP = 273.16 K stays inline in `goff_gratch_ice`); the two surface-exchange laws `density_dry_air` (ρ = p/(rair·t), camsrfexch.F90:503) and `exner_function` ((ps/p)^κ, dp_coupling.F90:731-733), both generic in their physical constants (rair, cappa taken as template parameters, bound by the caller to the named parameters of `CamOutExport`); and the vertical-diffusion laws `interface_diffusivity` (K = kv·dpidz_sq on the interfaces), `surface_layer_source` (the pver injection g·flux/dz[1]), `surface_stress_coefficient` (ksrf = ksrftms + max(tau/ws, ksrfmin)), `wind_damping_rate` (tau_damp = −g·ksrf/dz[1] at pver, 0 elsewhere) |
| `cam_wv_saturation.esm` | the `estbl` function table (251 entries, Goff–Gratch mixed phase) + the `WvSaturation` model (the run's qsat law) + the `SvpTransLaw` model (pins the table-generation formula) |
| `cam_rayleigh_friction.esm` | the `RayleighFriction` model (otau profile + Euler-backward tendencies) |
| `cam_energy_fixer.esm` | the `EnergyFixer` model (uniform heat_glob tendency + implied surface flux) |
| `cam_surface_exchange.esm` | the `CamOutExport` model (cam_export's cam_out export **laws**: dry-air density, Exner potential temperature, topo, the cpslec sea-level-pressure chain incl. the flat branch, the precipitation split, co2diag/co2prog; the pass-through exports tbot/zbot/ubot/vbot/qbot/pbot/netsw/flwds/soll/sols/solld/solsd/ozone are documented **imports** with no equations) + the `CamInImports` model (the cam_in coupler imports with their consumption sites, and the flx_net accumulation with the latent-heat constants as named parameters) |
| `cam_vertical_diffusion.esm` | the `VerticalDiffusion` model (the implicit Euler-backward vdiff solve as four PDEs, section below) |

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
rho    = pbot/(rair·tbot)                 ! dry-air density (F90:503),
                                          ! rair = AVOGAD·BOLTZ/MWDAIR (shr_const)
thbot  = tbot·exner_function(pbot, ps, κ) ! (F90:496, dp_coupling.F90:731-733),
                                          ! κ = cappa = rair/cpair (this run)
topo   = phis/gravit                      ! surface topographic height (F90:498)
psl    = cpslec(t, pmid, ps, phis)        ! ECMWF sea-level pressure
                                          ! (cam_diagnostics.F90:1080-1081,
                                          !  cam_diagnostic_utils.F90:20-81)
precc  = max(0, prec_dp + prec_sh)        ! pbuf PREC/Snow fields (F90:537-575)
precl  = max(0, prec_sed + prec_pcw)
precsc = min(precc, max(0, snow_dp + snow_sh))
precsl = min(precl, max(0, snow_sed + snow_pcw))
co2diag = co2vmr·1e6 [ppmv]               ! chem_surfvals CO2VMR (F90:506-511)
co2prog = 0                               ! co2_transport() false here (F90:512-516)
```

The remaining cam_out exports — tbot/zbot/ubot/vbot/qbot/pbot (direct state
copies, F90:495-507) and ozone/netsw/flwds/soll/sols/solld/solsd (the pbuf
SRFOZONE field and the radiation_tend outputs,
rrtmg/radiation.F90:1272-1316,1429) — are **documented imports**: cam_export
copies them without transformation, so the esm declares them as parameters
with no equations (an equation `tbot = t_bot` cannot fail and would only
restate the import). The physical constants are named, unit-carrying
parameters: rair/gravit/cpair (physconst.F90:89-102 =
SHR_CONST_RDAIR/G/CPDAIR), xlapse and the branch temperatures tstar0 =
290.5 K / tstar_cold = 255 K (cam_diagnostic_utils.F90:44, 64-72), and
cappa = rair/cpair (this run's cappav, air_composition.F90:354).

| cam_out field | status in the esm |
|---|---|
| tbot, zbot, ubot, vbot, qbot, pbot | **imports**: direct copies of the bottom-level state (F90:495-507), no equations |
| rho, thbot | **laws**: `density_dry_air` and `exner_function` templates (cam_templates.esm), bound to the named parameters rair and cappa |
| topo | **law**: topo = phis/gravit (F90:498), pinned as a derived value from the dumped phis (topo itself is not a dump column, see gaps) |
| psl | **law**: the cpslec chain incl. the flat branch (see below) |
| precc/precl/precsc/precsl | **law**: the pbuf-field split + clamps (F90:537-575), unpinned (see gaps) |
| co2diag | **law** co2diag = CO2VMR·1e6 [ppmv] with the input **imported** (chem_surfvals); provenance open question (see gaps) |
| co2prog | **law**: identically 0 in this configuration (co2_transport off) |
| ozone, netsw, flwds, soll, sols, solld, solsd | **imports**: pass-throughs from the radiation component (rrtmg/radiation.F90:1272-1316,1429) and the SRFOZONE pbuf, no equations |

The psl law runs in `diag_phys_writeout` **before** this step's
radiation update (verified bit-exact against the dumped psl using the
pre-radiation state), so the model takes two bottom-level temperatures:
`tbot` (at cam_export) and `t_psl` (at psl time). The dump pins
reconstruct `tbot` as T + ptend_s·dt/cpairv with cpairv = cpair =
1004.64 J/(kg·K) (this run has `dry_air_species = ''`; verified bit-
exact against the dumped tbot at steps 1/2/24/48).

The cpslec chain is transcribed with the branch structure of
cam_diagnostic_utils.F90:44-77: the shared branch predicates (b1:
Tstar ≤ 290.5 ∧ TT0 > 290.5; b2: Tstar > 290.5 ∧ TT0 > 290.5) are named
boolean intermediates used by both `psl_alph` and `psl_tstar_final`, and
the flat branch (|phis/gravit| < 1e-4 → psl = ps, F90:57-59) is written
as `ifelse(|topo| < 1e-4 m, ps, ps·exp(...))` reusing the topo law. The
pressure-reduction series reads psl = ps·exp(β·(1 − αβ/2 + (αβ)²/3))
with the shared product αβ = `psl_alphbeta` (F90:77). Note: the first
draft of this file had the two series signs flipped (1 + αβ/2 − (αβ)²/3);
all 48 window steps take branch 2 where α = 0 makes the correction
vanish, so the dump pins could not detect it — the esm now transcribes
the Fortran signs exactly, and the flip would only matter in branches 1
and 3, which remain unpinned (see gaps).

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

## Vertical diffusion (vertical_diffusion_solver.F90, diffusion_solver.F90; site vdiff)

`cam_vertical_diffusion.esm` transcribes the CAM6 vertical-diffusion solve
(the implicit Euler-backward finite-volume solve of `diffusion_solver.F90`
called by `vertical_diffusion_solver.F90`) as four PDEs on the
EarthSciDiscretizations `column_nonuniform_1d` grid:

    du/dt = d/dz(K_mom du/dz) + tau_damp u + src_u        (v mirrored)
    ds/dt = d/dz(K_th   ds/dz) + dtk/dt + src_s
    dq/dt = d/dz(K_q    dq/dz) + src_q

with `K_* = interface_diffusivity(kv_*, dpidz_sq) = kv_*·dpidz_sq` on
the interfaces (one template for the momentum/heat/tracer fields) and
the surface sources `src_* = surface_layer_source(flux, g, dz, zero)`
— the pver injection g·flux/dz[1]: `src_u = g·tauresx_enter·ramda/dz[1]`
(the residual-stress injection, F90:491-499), `src_s = g·shflux/dz[1]`
(F90:779), `src_q = g·cflx·dryfac/dz[1]` (the explicit q surface flux,
F90:974-975, with the pdel/pdeldry conversion of dry constituents,
F90:1056-1061) — with `ramda = dt/timeres` = 0.25 (the am_correction =
.false. branch, F90:491-496; timeres = 7200 s, ztodt = 1800 s in the
run). The wind damping is `tau_damp = wind_damping_rate(ksrf, g, dz,
zero)` — −g·ksrf/dz[1] at the surface layer and 0 elsewhere (F90:274-
285; nonzero only at CAM lev 32 in the window) — with the drag
coefficient `ksrf = surface_stress_coefficient(taux, tauy, u0, v0,
ksrftms) = ksrftms + max(sqrt(taux²+tauy²)/max(sqrt(u0²+v0²), wsmin),
ksrfmin)` (F90:232-239; wsmin/ksrfmin are the module parameters
F90:42-44). The `zero` argument of the templates is a units-free
declaration bound into the ifelse branches (the checker requires the
branches to share dimensions and the four sources carry different
units) — the ysu.esm idiom of the kdudz_bot/kdudz_top names below.
The diffusion term is lowered by the ESD rule
`varcoeff_face_laplacian_lev_flux_bc`, whose stencil (interface-valued
coefficient, no face averaging, layer-centre gradient spacing, prescribed
domain-face fluxes) is exactly the CAM `linear_1d_operators.F90`
finite-volume operator (spr/sub/diag, deriv_diag) — verified term by term.
The residual-stress memory is
`tauresx = tauresx_enter − tauimpx` with
`tauimpx = ∫ (u−u0)/(g·dt) dz` over the column (F90:536-589; taux,
tautmsx, taubljx, dragblj ≡ 0 in the window), lowered by the ESD rule
`integral_lev_whole`; tauimpy mirrors v. `dtk = 0.5(u0²+v0²) −
0.5(u²+v²)` is the do_beljaars kinetic-energy dissipation (F90:598-646),
bit-exact against the dump.

### Representation (esd order ↔ dump order)

All parameters and states are stored in esd order (lev 1 = surface),
reversed from the dump rows of `stage1/dumps/cam/cam_vdiff_prof.csv`
(one row per (step, lev), CAM lev 1 = model top) and
`cam_vdiff_q.csv` (one row per (step, lev, constituent)):

| esd | CAM | quantity |
| --- | --- | -------- |
| lev i        | k = 33−i     | u, v, s, q, tau_damp, u0/v0/s0, dz[i] = pdel(k) |
| interface n  | interface 34−n | ze[n] = −pint(34−n), dpidz_sq, K_* (kv = 0 in the window) |
| ze           | −pint reversed | Pa, increasing upward; zc[i] = −mid(33−i), z_sfc = ze[1], z_top = ze[NLEV+1] = −pint(1) |

The double sign flip (coordinate reversal plus the flux-direction flip of
the conservative stencil) makes the esd operator row i equal to the CAM
matrix row 33−i exactly: `L u[i] = (A·u)[33−i]`, re-derived and confirmed
against spr/sub/diag and deriv_diag by hand and by the compiled reference.

### Verified laws (dump window, all 48 steps, site vdiff)

All rel errors from the independent Python re-verification
(`stage1/analysis/vdiff_verify.py`, committed) of the instrumented run,
scanning every step of the window (the numbers quoted in the test
tolerances below):

- `tau_damp(32) = −g·ksrf/pdel32` (the wind_damping_rate +
  surface_stress_coefficient templates, F90:232-239, 274-285) to rel
  5.2e-15 across the window (1.96e-15 / 3.91e-15 / 1.82e-15 at steps
  1/2/48); τ = sqrt(taux²+tauy²), ws = max(sqrt(u0²+v0²), wsmin), with
  u0 = u1 − ptend_u·dt recovered from the dumped (u1−u0)/dt.
- `ptend_u/v(32) = tau_damp(32)·u1/v1(32) + g·taures_enter·ramda/pdel32`
  with ramda = dt/timeres = 0.25 (F90:491-499): rel 2.008e-12 (u, step
  1), 9.79e-14 (step 2), 3.89e-14 (step 48); v: 5.95e-13 (step 1),
  1.45e-14 (step 2), 1.039e-11 (step 48). Window maxima: 4.4e-12 (u,
  step 30), 9.8e-11 (v, step 21). The dump's ptend is the solver's own
  `(u1−u0)·rdt`, so the residual is dump-side rounding.
- `dtk = 0.5(u0²+v0²) − 0.5(u1²+v1²)` bit-equal (rel 0); the dumped
  `ptend_s = (s1−s0)·rdt` is `dtk/dt` to abs ≤ 6.4e-15 at the pinned
  steps (rel ≤ 6.1e-9 at the step-48 near-zero row; the window-wide
  worst row is rel 8.8e-6 at step 26 lev 32, abs 1.2e-14 — cancellation
  noise of the recovered u0 on a near-zero tendency).
- taures law `taures = taures_enter − tauimpx` (F90:536-589) to rel
  2.74e-12 (x) / 7.93e-13 (y) at step 2 (step 1: rel 0 / 2.3e-16;
  step 48: 1.92e-12 / 3.20e-14); window maxima 9.3e-12 (x, step 5) and
  8.7e-13 (y, step 5). Step-1 taures_enter is the observed pbuf initial
  (1.4478837360691117e-4 Pa) inverted from the law — the pbuf is
  zero-initialized at nstep == 0, and the run's first vdiff step already
  carries this value, so it is treated as data, not a bug.
- q: `ptend_q(32) = g·cflx·dryfac/pdel32 = g·cflx/pdeldry32` (F90:974-975,
  1056-1061) to rel ≤ 3.8e-15 for DMS, bc_a4, SOAE at steps 1/2/48 (Q is
  moist with cflx = 0); the all-row window maximum is 9.2e-11 on the
  H2SO4 rounding row (value 4.3e-25, abs 4e-34), and every dumped-zero
  row has |qsrc| ≤ 3.0e-38. dryfac = pdel/pdeldry = 1.0113874330578143 /
  1.0115819043995633 / 1.0182904760308067 for the dry constituents at
  steps 1/2/48.

### Tolerances (dump-side rounding is real and quantified)

- ptend_u/v at lev 1: rel 1e-10 (dump noise up to 9.8e-11 across the
  window, 1.04e-11 at the pinned steps, from cancellation in the dumped
  `(u1−u0)/dt`).
- tauresx/y: rel 1e-10 (law noise ≤ 9.3e-12).
- tau_damp: rel 1e-12 (noise ≤ 5.2e-15).
- dtk: rel 1e-12 (bit-exact).
- ptend_s: {rel 1e-12, abs 1e-14} — pinned-step abs noise ≤ 6.4e-15,
  principled as ≤ ulp(s1)/dt.
- ptend_q: rel 1e-12 (pinned noise ≤ 3.8e-15).
- nonzero-K operator pin: rel 1e-12 (compiled reference, no dump noise).

### The nonzero-K operator pin (compiled Fortran, not dumps)

The dump window runs all kvh/kvm/kvq ≡ 0, so no dump row exercises the
diffusion operator. The operator is pinned by a compiled-Fortran
reference driver (`stage1/analysis/cam_vdiff_driver/`, committed):
`coords_1d.F90` and `linear_1d_operators.F90` verbatim from the fork
with minimal shr_kind/shr_log/shr_sys stubs, built with the CAM gfortran
configuration (`-cpp -DUSE_CONTIGUOUS=` empty — CAM turns the F2008
contiguous attribute off for gfortran, cam/doc/ChangeLog:84560; with
`-DUSE_CONTIGUOUS=contiguous` the source `real(r8), USE_CONTIGUOUS
intent(in)` lacks the comma the attribute requires and does not compile).
The step-1 pint/dpidz_sq/u0 (u0 = u1 − ptend_u·ztodt) are transcribed
into DATA statements by `extract_step1_data.py` (repr round-trip, the
dump's 17-digit precision) with explicit `_r8` kind suffixes — a bare
literal in a DATA statement is default (single) kind and would be
rounded before conversion to real(r8). The driver builds `Coords1D`
from the dumped step-1 pint (CAM order), applies `diffusion_operator`
with K = dpidz_sq (kv ≡ 1 m²/s) and the default BoundaryFixedFlux
conditions at both ends (zero prescribed flux — exactly the esm's
kdudz_bot = kdudz_top = 0) to u0, and prints A·u0 (`build.sh` →
`reference_output.txt`). The esm pin test reverses the rows:
ptend_u(1) = 2.0633949733127771e-05, ptend_u(16) =
1.7093444110831203e-07, ptend_u(32) = 8.6780995093086661e-67 (the
model-top value with the near-vanishing top winds and the large step-1
top dpidz_sq); all 32 driver rows are bit-identical to an independent
Python replication of the stencil, and the surface and interior rows
were additionally hand-verified against the spr/sub/diag stencil
arithmetic.

### ESD/loader findings recorded for future CAM components

- **The integral rule's `where` matches only a bare layer field**: the
  CLI's shape inference does not reduce compound expressions
  (`integral((u−u0)/(g·dt), …)` fails the `f: {shape: [lev]}` check and
  the node survives to the unlowered_operator gate), so the integrand is
  bound to a named observed (`tauimpx_f`, `tauimpy_f`) first — the
  fixture precedent.
- **The rule's `kdudz_bot`/`kdudz_top` free names cannot be renamed at
  the import edge** (`rename` only touches edge exports; the names are
  free names inside the template body) and, being shared across every
  firing, one units-typed declaration can satisfy only one field (the
  momentum flux K du/dz is kg/s⁴, the heat flux kg·m/s⁵). The ysu.esm
  idiom is the fix: declare the two names **without units** (a
  units-free variable is dimension-agnostic to the checker) and share one
  rule instance across all unit classes — valid here because every
  prescribed surface flux in the window is zero and enters through the
  explicit sources.

`./esm units` on `cam_vertical_diffusion.esm` reports 28 consistent
equations, 0 mismatched, and one not checked: the grid geometry equation
`zc = column_cell_center(ze)`, whose lowered `faq` node has no
determinable dimension — the same checker artifact the WvSaturation
`table_lookup` note describes; no
`units_dimensionless_literal_mismatch` anywhere.

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
| `CamOutExport` | 6 | 4 dump pins of the state laws (steps 1/2/24/48, bc_out; rho/thbot/psl = dumped, topo = derived from the dumped phis), 1 cpslec-chain pin (step 1, branch 2, with the nine intermediates incl. the shared predicates psl_b1/psl_b2, psl_alpha and psl_alphbeta), 1 flat-branch pin (branch-exercising phis = 5e-4 m²/s², psl = ps; fails if the branch is missing, see the test description) |
| `CamInImports` | 0 | all fields are coupler imports; the one law (flx_net) needs net_flx_rad, which is not dumped → gap (no data invented) |
| `VerticalDiffusion` | 9 | 23 dump pins at steps 1/2/48 (site vdiff: tau_damp(1) at steps 1/2/48 plus tau_damp(2) = 0 and ramda = dt/timeres = 0.25 at step 1, tauresx/y, dtk, ptend_u/v/s at lev 1), 5 q dump pins (DMS/bc_a4/SOAE at step 1, DMS at steps 2/48), 3 compiled-Fortran pins (nonzero-K operator at esd lev 1/16/32) |

Total for the CAM slice: 65 + 27 + 6 + 23 + 31 = 152 assertions. The
compiled-Fortran references compile the actual subroutines verbatim from
the fork with gfortran 11.5 in a driver at the dump's step-1 u/v state
(the rayleigh/energy-fixer driver is kept in `scratch/cam_laws_driver/`
in the worktree, the mosart-precedent style; the vdiff operator driver
is committed under `stage1/analysis/cam_vdiff_driver/`). The stage-1
validator (`stage1/analysis/validate_cam_dumps.py`) independently
reproduces the qsat dump to round-off, and the dump window was verified
solution-neutral by the stage-1 lane.

`./esm units stage2/cam_*.esm` reports all equations consistent except
the single `es0 = table_lookup` equation of `WvSaturation`, which the
checker cannot dimension (a `table_lookup` node has no determinable
dimension — inherent to the checker, documented), and the `zc =
column_cell_center(ze)` grid-geometry equation of
`cam_vertical_diffusion.esm` (the same `faq`-node artifact); no
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
- **s-equation top boundary** of the vdiff solve: CAM's dse solve uses
  BoundaryFixedLayer(0)/BoundaryData(dse_top) at the top, whose flux
  `2·K(1)·(dse_top − s(1))/del(1)` differs from the rule's zero-flux
  kdudz_top whenever K(1) > 0; the window has kvh(1) = 0 and dse_top = 0,
  so the two coincide here and the gap is documented, not bridged.
- **Countergradient terms** of the vdiff solve (cgh/cgs and the qmincg
  clipping, dropmixnuc, ubc_mmr fixed upper boundary): all zero or inert
  in the window (cgh = cgs ≡ 0, cnst_fixed_ubc false for the pinned
  constituents, dropmixnuc flux zero); not carried.
- **tauresx/y at steps ≥ 2** use the previous step's updated dump value
  as the entering residual (the pbuf memory), per the verified law.
- **nonzero-K stencil** is pinned by the compiled reference, not by
  dumps (all kv ≡ 0 in the window).
- **zbot**: zm(pver) at cam_export time is not dumped. `physics_update`
  recomputes zm after the temperature update (physics_types.F90:423-449),
  so the pre-radiation zm of the dump differs from what cam_export sees;
  reconstructing it needs the geopotential law (geopotential_temp.F90)
  plus the condensate mixing ratios at the radiation site (not dumped).
  zbot is an import (no equation) — the geopotential law itself is a
  separate component.
- **topo**: cam_out%topo (camsrfexch.F90:498) is not a dump column; the
  law topo = phis/gravit is implemented and pinned as a derived value
  from the dumped phis (constant 8704.816911504207 m²/s² across the
  window). A stage-1 dump of cam_out%topo would make it a direct pin.
- **precipitation split**: the eight pbuf fields (PREC_DP, SNOW_DP,
  PREC_SH, SNOW_SH, PREC_SED, SNOW_SED, PREC_PCW, SNOW_PCW) are not
  dumped, so the precc/precl/precsc/precsl laws are unpinned. The dump
  shows nonzero precc/precl on 15 of the 48 window steps, so the split
  is exercised by the run but its inputs are not observable; the laws
  are implemented per camsrfexch.F90:537-575.
- **cpslec branches**: all 48 window steps take branch 2 (Tstar >
  290.5 K and TT0 > 290.5 K → alph = 0, Tstar ← 0.5·(290.5+Tstar)). The
  flat branch (|phis/gravit| < 1e-4 → psl = ps) is implemented and
  pinned by a branch-exercising pin (phis = 5e-4 m²/s², a near-sea-level
  column; the test description quantifies how the assertion fails if the
  branch is missing). Branch 1 (Tstar ≤ 290.5, TT0 > 290.5) and branch 3
  (incl. the Tstar < 255 K sub-branch) are implemented per the Fortran
  but unpinned — the 887 m SCAM column never exercises them. Because
  branch 2 sets alph = 0, the dump pins also cannot distinguish the
  sign of the series correction (1 − αβ/2 + (αβ)²/3); the esm transcribes
  the Fortran signs (the first draft had them flipped, see the surface-
  exchange section).
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
  site. The esm gives co2diag the unit ppmv (the Fortran's mol/mol
  CO2VMR × 1e6); the dump header labels the co2diag column mol/mol while
  the dumped values are ppmv-scale (~1.127) — an instrumentation
  labeling slip in our fork, not a data corruption (the README of the
  dump is the source of truth for the values).
- **lhf = cflx·LATVAP** holds in the dump on every row but is a coupler
  (aoflux) property, not a CAM law — deliberately not implemented.

## Bugs found in the CAM code

None. (The 251-entry table with a never-used last knot is faithful
Fortran behavior, not a bug; the tanh saturation is glibc behavior.)
The step-1 taures pbuf initial of the vdiff site (nonzero at the first
vdiff call despite the pbuf being zero-initialized at nstep == 0) is
provenance-unexplained data, not a confirmed bug.
