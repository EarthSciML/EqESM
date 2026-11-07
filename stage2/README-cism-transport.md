# CISM stage 2 — SMB application, tracer surface/bed handling, post-transport cleanup

Stage-2 implementation of the pointwise slice of the CISM thickness/tracer
solve (`glissade_thickness_tracer_solve`) as built in CESM
(`components/cism` = cismwrap_3_00_002, source `source_cism` =
cism_main_3.00.002, cesm3.0-alphabranch checkout `stage1/cesm`) as esm 1.3.0
files: the parts of the solve that are **not** grid operations — the cap
qsmb→acab conversion and glad-level zeroing (CISM_LANE.md row 14 boundary),
the surface/basal mass-balance application with its tracer blends
(rows 6/7), and the post-transport surface/bed elevation cleanup. The
horizontal transport and the first-order vertical remap
(`libglissade/glissade_transport.F90`) are discretization operations: they
need the stage-3 assembly with the EarthSciDiscretizations stencils and are
**out of scope here**. What is pinned: the SMB application, the temperature
surface/bed handling, the thickness clip, and the cleanup, against the
stage-1 instrumentation dump `stage1/dumps/cism/cism_transport.csv`
(phase=entry at the top of the solve, before horizontal transport;
phase=exit after the SMB application and tracer cleanup) and the cap dump
`cism_cap.csv` (the qsmb import). All model logic lives in `.esm` files;
this README is the human-readable map.

Sources:
`stage1/cesm/components/cism/source_cism/libglad/glad_main.F90` (qsmb→acab),
`libglad/glad_timestep.F90` (glad-level zeroing),
`libglissade/glissade_mass_balance.F90` (the SMB application),
`libglissade/glissade.F90` (the age increment),
`libglissade/glissade_transport.F90` (the tracer surface/bed values),
`libglimmer/glimmer_utils.F90` (calc_lsrf_usrf), and the run configuration
of the T1850Gg/f09_g17_gris20 reference run (2 ice dynamics steps per cap
year, dt = 0.5 yr, upn = 11 layers, eus = 0).

Reference configuration pinned (from the dumps and the source defaults):
`which_smb = SMB_APPLY_RATE` (the reference run's mass-balance application),
`which_ho_vertical_remap = HO_VERTICAL_REMAP_FIRST_ORDER`,
`which_ho_ice_age = HO_ICE_AGE_COMPUTE`, `tstep_mbal = 8760 h` (1 yr),
`dt = 0.5 yr`, `rhoi = 917 kg/m³`, `rhoo = 1026 kg/m³` (SHR_CONST_RHOSW),
`rhow = 1000 kg/m³`, `scyr = 31536000 s`, the SIGMA_COMPUTE_GLIDE sigma
levels with n = 2 (σ₂ = 0.23140495867768607, σ₁₀ = 0.96398891966759004),
`thklim = 0` for the SMB masks, `eus = 0`.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `cism_transport.esm` | Model `GladSmb` — the cap qsmb→acab conversion acab = qsmb·tstep_mbal·hours2seconds/rhoi·windows_per_year with the named glad/glimmer constants hours2seconds = 3600 s/h, rhow = 1000 kg/m³, rhoi = 917 kg/m³, scyr = 31536000 s and the DERIVED windows_per_year = scyr/(tstep_mbal·hours2seconds) [= 1 in the reference run; the Fortran's rhow factors cancel], and the glad-level zeroing: acab = 0 where thck ≤ 0 ∧ acab < 0 (ice-free ablation) and where the cell is ice-free ocean (the `cism_ocean` predicate: topg < eus ∧ thck ≤ 0); model `MassBalance` — the SMB application: the rate law d(thck)/dt = ifelse(thck ≤ 0 ∧ net < 0, 0, net) with net = acab − bmlt, the per-step algebraic maps of the dumped entry state (melt_potential, the applied rates acab_applied/bmlt_applied, the layer-1 blend T1′ = (H₁·T1 + s·Tsurf)/(H₁ + s), the layer-10 blend T10′ = (H₁₀·T10 + b·Tbed)/(H₁₀ + b), the column-mean law, the tracer zeroing of melted columns), and the ice-age increment d(age)/dt = 1 where thck > 0; model `SurfaceUpdate` — calc_lsrf_usrf: lsrf = topg where grounded, lsrf = eus − (rhoi/rhoo)·thck where floating, usrf = lsrf + thck. The document declares a `solver` block (stiffness low, abstol 1e-10, reltol 1e-6) for the non-smooth clip — see 'The clip and the solver block' below | `glad_main.F90`: 810-824; `glad_timestep.F90`: 171-187; `glissade_mass_balance.F90`: 819-1571; `glissade.F90`: 1203-1213; `glissade_transport.F90`: 139-208; `glimmer_utils.F90`: 344-378 |
| `cism_templates.esm` | `cism_sigma` — the glide stretched-sigma coordinate sigma(x, n) = (1 − (x+1)^(−n))/(1 − 2^(−n)) of glide_calc_sigma, applied by MassBalance at σ₂ and σ₁₀; `cism_ocean` — the shared ocean predicate ocean = topg < eus ∧ thck ≤ 0 (glide_mask.inc:33-34 / glissade_masks.F90:70-96), used by GladSmb's zeroing and by the MassBalance ocean mask; `cism_tracer_mix` — the shared mixing law T′ = (H·T + s·T′)/(H + s) of the accumulation/freeze-on blends | `glide_setup.F90`: 303-308, 380-387; `glide_mask.inc`: 33-34; `glissade_masks.F90`: 70-96; `glissade_mass_balance.F90`: 1380-1386, 1458-1465 |

## The physics chain (Fortran → esm transcription)

The per-cap-step chain is: the cap accumulates qsmb (kg m⁻² s⁻¹) over the
mass-balance window `tstep_mbal`; glad converts it to a per-window thickness
and zeroes it where ice-free/ablation and in ice-free ocean
(`glad_timestep.F90:171-187`, `glad_main.F90:810-824`); each dynamics step
(0.5 yr) applies acab and bmlt to the ice column
(`glissade_apply_smb` → `mass_balance_driver` →
`add_surface_and_basal_mass_balance`, `glissade_mass_balance.F90:819-1571`),
remaps the tracers vertically, and recomputes lsrf/usrf
(`calc_lsrf_usrf`, `glimmer_utils.F90:344-378`). The stage-1 fork traces the
chain at the top of the solve (entry) and after the SMB application and
cleanup (exit).

- **GladSmb** (`glad_main.F90:810-824`): the Fortran accumulates qsmb in
  kg m⁻² s⁻¹ over `tstep_mbal` hours and converts to m of ice: acab =
  qsmb·tstep_mbal·hours2seconds/rhow·rhow/rhoi·windows_per_year with the
  named glad/glimmer constants hours2seconds = 3600 s/h
  (`glad_constants.F90:44`), rhow = SHR_CONST_RHOFW = 1000 kg/m³, rhoi =
  SHR_CONST_RHOICE = 917 kg/m³ (`glimmer_physcon.F90:41-43`) and scyr =
  31536000 s/yr (`glimmer_physcon.F90:90`) driving the DERIVED
  windows_per_year = scyr/(tstep_mbal·hours2seconds) [1/yr; = 1 in the
  reference run, tstep_mbal = 8760 h = 1 yr]. The water density cancels
  (qsmb·3600·tstep_mbal is kg/m², ÷ rhoi is m), so the esm conversion is
  acab = qsmb·tstep_mbal·hours2seconds/rhoi·windows_per_year (the dump
  column is labeled m/yr ice). The `acab = −thck_temp` cap of
  `glad_timestep.F90:233-235` ("Adjust glad acab for output") is applied
  **after** `glide_set_acab` and only modifies the glad output variable,
  never the applied field — documented, not transcribed.
- **The zeroing** (`glad_timestep.F90:171-187`): acab = 0 where thck ≤ 0
  and acab < 0 (ice-free ablation), and where GLIDE_IS_OCEAN (bed below sea
  level, no ice). Both fire independently; the esm `acab` equation nests
  them via the shared `cism_ocean` predicate template.
- **MassBalance** (`glissade_mass_balance.F90:819-1571`): the Fortran walks
  the 10 layers: accumulation adds s = acab·dt to layer 1 with the blend
  T1′ = (H₁·T1 + s·Tsurf)/(H₁ + s) at the surface temperature; ablation
  removes the top layers without changing their tracers; freeze-on adds
  b = −bmlt·dt to layer 10 with the blend at the bed temperature; basal
  melt removes the bottom layers; fully removed layers get tracer = 0; a
  column that melts out entirely accumulates the excess into
  `melt_potential` (a per-step diagnostic, zeroed at the top of every
  `mass_balance_driver` call). The esm model is the layer-free continuous
  limit with the RATE-LEVEL guard d(thck)/dt = ifelse(thck ≤ 0 ∧ net < 0,
  0, net): the derivative stops exactly when the column is gone, which
  reproduces the Fortran's layer-by-layer removal that clips the column at
  zero (F90:1398-1419) without carrying a deficit into the next step (the
  Fortran re-derives thck from the clipped layers every step; the
  state-level max(0, thck_raw) of the first review round could not express
  the per-step zeroing and was replaced by the guard). The per-step
  bookkeeping — melt_potential = max(0, −(H_pre + net·dt)), the applied
  rates acab_applied/bmlt_applied, and the tracer blends — is written as
  algebraic maps of the dumped entry state (parameter `H_pre` + entry
  tracers), evaluated at the step start, which is exact at the step
  boundaries the tests evaluate. The tracer laws keep the layer-1
  temperature, the layer-10 temperature, and the column mean through the
  shared `cism_tracer_mix` law T′ = (H·T + s·T′)/(H + s); the ablation
  mean law removes the ablated ice at the layer-1 temperature, which is
  exact while the ablation stays within layer 1. The column-mean output
  Tbar_out is derived (the Fortran has no column-mean state) and is
  UNPINNED: no dump case discriminates it (see 'Coverage gaps').
- **The ocean mask**: `glissade_get_masks` with thklim = 0
  (`glissade_masks.F90:70-96`) — ocean = topg < eus ∧ thck ≤ 0 (the same
  `cism_ocean` predicate as the glad zeroing); the SMB and BMB are not
  applied there and melt_potential stays 0 (the `effective_areafrac = 0`
  of ocean cells, F90:893-901). No traced cell is ocean (all 32 cells have
  topg ≥ 921 m): unpinned (see 'Coverage gaps').
- **The ice age** (`glissade.F90:1203-1213`): age += dt where thck > 0,
  age = 0 on ice-free cells; the tracer surface/bed values age_usrf = 0,
  age_lsrf = bottom-layer age (`glissade_transport.F90:204-208`). The esm
  ODE stops the increment on ice-free cells but cannot express the
  instantaneous reset — documented difference (unobservable in the window,
  which never refreezes). Unpinned: there is no age column in the dump
  (see 'Coverage gaps').
- **SurfaceUpdate** (`glimmer_utils.F90:344-378`): lsrf = topg where
  grounded (topg − eus ≥ −(rhoi/rhoo)·thck), lsrf = eus − (rhoi/rhoo)·thck
  where floating, usrf = lsrf + thck.

### The clip and the solver block

The rate guard is discontinuous at the clip, and the runner's default
diffsol BDF arm chokes at the crossing:

```
simulate failed: diffsol error: ODE solver error: Step size is too small at time = 0.000009337527277120702
```

(exactly the analytic crossing time thck0/|net| of the failing cell; the
same error appears at 0.07506885100532178 for the cell-29,13 pin). The
document therefore declares its own numerics in the esm-spec §2.2 `solver`
block — `stiffness: "low"` (the system is not stiff; the declaration steers
other bindings to their explicit-RK arms, which cross the kink cleanly) and
integration tolerances `abstol: 1e-10, reltol: 1e-6` (the runner's default
reltol 1e-4 / abstol 1e-6 leaves the BDF step control stuck at the kink). Under these
tolerances the BDF arm crosses the clip and freezes the state at its first
negative step, landing within ~1e-9 m of 0; the full-melt pins assert
thck = 0 with abs tolerance 1e-8 — still eight orders below the −0.35 m the
unguarded law would give, so the clip stays fully discriminating. The ERK
arm (TSIT45) lands on the clip to round-off and passes the same assertions.

## How each law is pinned

The dump rows: 32 cells × 4 steps × 2 phases, upn = 11 (temp_000 = surface,
temp_001..010 = layers 1..10, temp_011 = bed); steps 1,2 ↔ cap_step 1,
steps 3,4 ↔ cap_step 2; every row carries time = 0.5·step yr. All test
inputs come from the dump entry phases (the pre-step state) and the
cism_cap.csv import rows; all expected values come from the dump exit
phases or from the source formulas (stated in the test descriptions).
22 tests / 29 assertions.

| Law | Pin | Tolerance | Residual |
|---|---|---|---|
| qsmb→acab | cap import rows at cells 25,11 (cap steps 1, 2) and 26,11 (cap step 1) | rel 1e-12 | ≤ 1.1e-16 (exact) |
| glad zeroing (ice-free ablation) | cell 28,11 step 2 (thck = 0, raw acab = −0.7075 → 0; dump-verified); cell 28,11 step 3 (thck = 1.65e-5 > 0 → NOT zeroed; dump-verified) | exact | 0 |
| d(thck)/dt = net | accumulation: cell 26,11 steps 1/2/4 (rel 1e-6 / 1e-9 / 1e-6); ablation: cell 25,11 steps 3/4 (rel 2e-4 / 3e-4); full melt: cells 28,11 (steps 1, 3), 29,13 (step 1) → thck = 0 (abs 1e-8) | see left | advection: +9.6e-7 / +9.6e-10 / −1.7e-6 m; +1.8e-5 / +1.9e-5 m; BDF clip landing ≤ ~1e-9 m |
| full-melt clip + melt_potential | cells 28,11 (s1: 0.33379825832198695 m; s3: 0.8837069386692753 m), 29,13 (s1: 0.22629773133841058 m) — melt_potential is derived (not a dump column), from the Fortran's own conservation check | rel 1e-12 | exact |
| acab_applied | clipped branch only (cell 28,11 s1: acab_applied = −H_pre/dt = −0.039916999638080596) | rel 1e-12 | exact (source-formula: not dump columns) |
| layer-1 accumulation blend | cells 26,11 s1 (rel 1e-5), 27,11 s1 (rel 7e-5; s > H₁), 25,12 s1 (rel 2.5e-4) | see left | −1.7e-6 / −1.2e-5 / −3.6e-5 degC |
| tracer zeroing on full melt | cell 28,11 s1 (T1_out = 0; the dump exit temp_001..010 = 0 of the melted columns) | exact | 0 |
| lsrf/usrf | grounded: cell 25,11 s1; ice-free: cell 28,11 s3 (thck = 0); thin ice: cell 29,13 s4 | rel 1e-12 | exact (all 256 rows: lsrf = topg, usrf = lsrf + thck to round-off) |

The tolerance convention: 1e-12 where the esm law is exact (pure algebra
or linear-in-time integration), 5-10× the measured residual where the dump
exit value carries a physics effect the pointwise model does not contain
(horizontal advection, vertical remap) — each documented in the test
description, following the MOSART scheme-gap pattern. The full-melt abs
1e-8 covers the BDF arm's clip landing (see 'The clip and the solver
block').

## Coverage gaps (deferred or unpinnable)

- **Horizontal transport and the first-order vertical remap** are grid
  operations (`glissade_transport.F90:100-272, 1101-1441`) and belong to
  the stage-3 assembly with EarthSciDiscretizations. Their fingerprints are
  the advection residuals of the thickness and tracer pins above.
- **bmlt = 0 everywhere in the reference window** (bmlt_ground =
  bmlt_float = 0 on all 256 rows: the cold base melts nothing and the
  window has no shelf ice). The basal-melt (> 0) removal branches, the
  freeze-on (< 0) branch, and the layer-10 freeze-on blend T10′ are
  transcribed from the source but have no dump pin (the earlier
  freeze_on_blend_synthetic was dropped — synthetics pin the esm code, not
  the Fortran).
- **No ocean cells in the window** (all 32 traced cells have topg ≥ 921 m):
  the glad-level ocean zeroing branch and the glissade-level
  ocean-no-application branch (net_applied = 0, melt_potential = 0) are
  transcribed from the source but unpinned (the earlier
  zeroing_ocean_synthetic and ocean_no_application_synthetic were dropped).
- **No ice-age dump column**: the age increment is transcribed from the
  source but unpinned (the earlier age_increment_synthetic and
  age_ice_free_stops synthetics were dropped); the Fortran's ice-free reset
  (age = 0) is not representable as an ODE — documented in the equation
  description.
- **Floating lsrf/usrf**: the window has no shelf ice (the
  floating_synthetic was dropped); the floating branch
  lsrf = eus − (rhoi/rhoo)·thck is transcribed from the source and
  unpinned. All 256 dump rows satisfy lsrf = topg (grounded).
- **Tbar_out / the ablation mean law**: dropped. The column-mean output is
  derived (the Fortran has no column-mean state), and the search of all 256
  dump rows found no ablation case where the law's effect exceeds the
  advection residual (best ratios ~2.2, effects 4.7e-4 / 3.8e-5 relative —
  inside the advection-dominated tolerance), so no assertion could
  discriminate the law: the earlier ablation_mean_step3/step4 pins passed
  even with Tbar_out = Tbar and were removed. The law itself remains in the
  model (Tbar_blend → Tbar_out via the composed cism_tracer_mix law),
  documented and unpinned.
- **acab_applied / bmlt_applied / melt_potential / ocean cells** are not
  dump columns (the dump carries the input acab, not the applied values):
  pinned from the source formulas at the full-melt and applied-rate pins.
- **The layer-free tracer model** tracks the layer-1 temperature, the
  layer-10 temperature, and the column mean only; the per-layer identities
  (which layer the ablation bites into, the full 10-layer profile) are
  stage-3 vertical-discretization content (the sim1d.py reproduction of
  the exit profiles lives in `scratch/` as the verification reference).
- **The (c) 'acab = −thck_temp' output cap** of glad_timestep.F90:233-235
  is not transcribed (it adjusts the glad output variable only, after
  glide_set_acab; the applied field is unaffected) — documented.
- **melt_potential accumulation across steps**: the Fortran resets
  melt_potential to 0 at the top of every mass_balance_driver call; the
  esm form max(0, −(H_pre + net·dt)) is the per-step value of the dumped
  entry state. The two agree for any single step starting from an
  ice-covered state, which is what the tests pin.

## Bugs

No CISM bugs found in this slice; nothing added to `bugs.md`.
