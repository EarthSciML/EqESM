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

## Vertical diffusion (vertical_diffusion.F90, diffusion_solver.F90; site vdiff)

`cam_vertical_diffusion.esm` transcribes the CAM6 vertical-diffusion solve
(the implicit Euler-backward finite-volume solve of `diffusion_solver.F90`
called by `vertical_diffusion_solver.F90`) as four PDEs on the
EarthSciDiscretizations `column_nonuniform_1d` grid:

    du/dt = d/dz(K_mom du/dz) + tau_damp u + src_u        (v mirrored)
    ds/dt = d/dz(K_th   ds/dz) + dtk/dt + src_s
    dq/dt = d/dz(K_q    dq/dz) + src_q

with `K_* = kv_* * dpidz_sq` on the interfaces and the surface sources
`src_u = g*tauresx_enter*ramda/dz[1]` (the pver injection of
diffusion_solver.F90:491-499), `src_s = g*shflux/dz[1]`,
`src_q = g*cflx*dryfac/dz[1]` (the explicit q surface flux, F90:974-975,
with the pdel/pdeldry conversion of dry constituents, F90:1056-1061).
The diffusion term is lowered by the ESD rule
`varcoeff_face_laplacian_lev_flux_bc`, whose stencil (interface-valued
coefficient, no face averaging, layer-centre gradient spacing, prescribed
domain-face fluxes) is exactly the CAM `linear_1d_operators.F90`
finite-volume operator (spr/sub/diag, deriv_diag) — verified term by term.
The wind damping is the solver's `tau_damp_rate` (nonzero only at CAM
lev 32 in the window) and the residual-stress memory is
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
| ze           | −pint reversed | Pa, increasing upward; zc[i] = −mid(33−i), z_sfc = ze[1], z_top = ze[33] = −pint(1) |

The double sign flip (coordinate reversal plus the flux-direction flip of
the conservative stencil) makes the esd operator row i equal to the CAM
matrix row 33−i exactly: `L u[i] = (A·u)[33−i]`, re-derived and confirmed
against spr/sub/diag and deriv_diag by hand and by the compiled reference.

### Verified laws (dump window, steps 1/2/48, site vdiff)

All rel errors from the independent Python re-verification
(`/scratch.local/ctessum/vdiff_verify.py`) of the instrumented run:

- `ptend_u/v = (u1−u0)/dt` EXACT (rel 0); the implicit identity
  `τ·u1 + g·taures_enter·ramda/pdel` matches the dumped ptend_u/v to rel
  2.0e-12 (step 1), 1.05e-13 (step 2), 1.019e-11 (step 48 v); the dump's
  ptend is the solver's own `(u1−u0)·rdt`, so the residual is dump-side
  rounding.
- `dtk = 0.5(u0²+v0²) − 0.5(u1²+v1²)` bit-equal (rel 0); `s1 = s0 + dtk`
  to ±1 ulp(s0), so the dumped `ptend_s = (s1−s0)·rdt` is `dtk/dt` ±
  ≤ 6.4e-15.
- taures law `taures = taures_enter − tauimpx` to rel 2.7e-12 (x) /
  7.9e-13 (y); step-1 taures_enter is the observed pbuf initial
  (1.4478837360691117e-4 Pa) inverted from the law — the pbuf is
  zero-initialized at nstep == 0, and the run's first vdiff step already
  carries this value, so it is treated as data, not a bug.
- q: `q1(32) − q0(32) = dt·g·cflx/pdel32` (≤ 3.8e-15) and
  `ptend_q(32) = g·cflx·dryfac/pdel32` (≤ 3.9e-15) for Q, DMS, bc_a4,
  SOAE at step 1 and DMS at steps 2/48; dryfac = pdel/pdeldry =
  1.0113874330578143 / 1.0115819043995633 / 1.0182904760308067 for the
  dry constituents (17 of them; Q is moist with cflx = 0).

### Tolerances (dump-side rounding is real and quantified)

- ptend_u/v at lev 1: rel 1e-9 (dump noise up to 1.019e-11 from
  cancellation in the dumped `(u1−u0)/dt`).
- tauresx/y: rel 1e-9 (law noise ~2.7e-12).
- dtk: rel 1e-12 (bit-exact).
- ptend_s: {rel 1e-12, abs 1e-14} — principled as ≤ ulp(s1)/dt.
- ptend_q: rel 1e-12 (noise ≤ 3.9e-15).
- nonzero-K operator pin: rel 1e-12 (compiled reference, no dump noise).

### The nonzero-K operator pin (compiled Fortran, not dumps)

The dump window runs all kvh/kvm/kvq ≡ 0, so no dump row exercises the
diffusion operator. The operator is pinned by a compiled-Fortran
reference driver (`scratch/cam_vdiff_driver/`): `coords_1d.F90` and
`linear_1d_operators.F90` verbatim from the fork with minimal
shr_kind/shr_log/shr_sys stubs, built with the CAM gfortran configuration
(`-cpp -DUSE_CONTIGUOUS=` empty — CAM turns the F2008 contiguous
attribute off for gfortran, cam/doc/ChangeLog:84560; with
`-DUSE_CONTIGUOUS=contiguous` the source `real(r8), USE_CONTIGUOUS
intent(in)` lacks the comma the attribute requires and does not compile).
The driver builds `Coords1D` from the dumped step-1 pint (CAM order),
applies `diffusion_operator` with K = dpidz_sq (kv ≡ 1 m²/s) and the
default BoundaryFixedFlux conditions at both ends (zero prescribed flux —
exactly the esm's kdudz_bot = kdudz_top = 0) to the dumped u0, and prints
A·u0. The esm pin test reverses the rows: ptend_u(1) = 2.0633949733127771e-05,
ptend_u(16) = 1.7093444110831203e-07, ptend_u(32) = 8.6780995093086661e-67
(the model-top value with the near-vanishing top winds and the large
step-1 top dpidz_sq). Surface and interior rows were additionally
hand-verified against the spr/sub/diag stencil arithmetic.

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

### Coverage gaps

- **s-equation top boundary**: CAM's dse solve uses
  BoundaryFixedLayer(0)/BoundaryData(dse_top) at the top, whose flux
  `2·K(1)·(dse_top − s(1))/del(1)` differs from the rule's zero-flux
  kdudz_top whenever K(1) > 0; the window has kvh(1) = 0 and dse_top = 0,
  so the two coincide here and the gap is documented, not bridged.
- **Countergradient terms** (cgh/cgs and the qmincg clipping,
  dropmixnuc, ubc_mmr fixed upper boundary): all zero or inert in the
  window (cgh = cgs ≡ 0, cnst_fixed_ubc false for the pinned
  constituents, dropmixnuc flux zero); not carried.
- **tauresx/y at steps ≥ 2** use the previous step's updated dump value
  as the entering residual (the pbuf memory), per the verified law.
- **nonzero-K stencil** is pinned by the compiled reference, not by
  dumps (all kv ≡ 0 in the window).

### Bugs found in the CAM code

None confirmed. The step-1 taures pbuf initial (nonzero at the first
vdiff call despite the pbuf being zero-initialized at nstep == 0) is
noted as provenance-unexplained data, not a confirmed bug.

`./esm units` on `cam_vertical_diffusion.esm` reports 25 consistent
equations, 0 mismatched, and one not checked: the grid geometry equation
`zc = column_cell_center(ze)`, whose lowered `faq` node has no
determinable dimension — the same checker artifact the WvSaturation
`table_lookup` note describes; no
`units_dimensionless_literal_mismatch` anywhere.
