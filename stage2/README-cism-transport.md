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
| `cism_transport.esm` | Model `GladSmb` — the cap qsmb→acab conversion acab = qsmb·tstep_mbal·3600/rhoi (the Fortran's rhow factors cancel) and the glad-level zeroing: acab = 0 where thck ≤ 0 ∧ acab < 0 and where the cell is ice-free ocean (topg < eus ∧ thck ≤ 0); model `MassBalance` — the SMB application: d(thck)/dt = acab − bmlt on the unclipped thickness with thck = max(0, thck_raw), the melt_potential of fully melted columns, the applied-rate bookkeeping acab_applied/bmlt_applied, the layer-1 accumulation blend T1′ = (H₁·T1 + s·Tsurf)/(H₁ + s), the layer-10 freeze-on blend T10′ = (H₁₀·T10 + b·Tbed)/(H₁₀ + b), the column-mean temperature law, the tracer zeroing of melted layers, and the ice-age increment d(age)/dt = 1 where thck > 0 with age_usrf = 0, age_lsrf = age; model `SurfaceUpdate` — calc_lsrf_usrf: lsrf = topg where grounded, lsrf = eus − (rhoi/rhoo)·thck where floating, usrf = lsrf + thck | `glad_main.F90`: 810-824; `glad_timestep.F90`: 171-187; `glissade_mass_balance.F90`: 819-1571; `glissade.F90`: 1203-1213; `glissade_transport.F90`: 139-208; `glimmer_utils.F90`: 344-378 |
| `cism_templates.esm` | `cism_sigma` appended — the glide stretched-sigma coordinate sigma(x, n) = (1 − (x+1)^(−n))/(1 − 2^(−n)) of glide_calc_sigma, applied by MassBalance at σ₂ and σ₁₀ | `glide_setup.F90`: 303-308, 380-387 |

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
  qsmb·tstep_mbal·3600/1000·rhow/rhoi. The water density cancels
  (qsmb·3600·tstep_mbal is kg/m², ÷ rhoi is m), so the esm conversion is
  acab = qsmb·tstep_mbal·3600/rhoi·windows_per_year with
  windows_per_year = 8760/tstep_mbal = 1 in the reference run (the dump
  column is labeled m/yr ice). The `acab = −thck_temp` cap of
  `glad_timestep.F90:233-235` ("Adjust glad acab for output") is applied
  **after** `glide_set_acab` and only modifies the glad output variable,
  never the applied field — documented, not transcribed.
- **The zeroing** (`glad_timestep.F90:171-187`): acab = 0 where thck ≤ 0
  and acab < 0 (ice-free ablation), and where GLIDE_IS_OCEAN (bed below sea
  level, no ice). Both fire independently; the esm `acab` equation nests
  them.
- **MassBalance** (`glissade_mass_balance.F90:819-1571`): the Fortran walks
  the 10 layers: accumulation adds s = acab·dt to layer 1 with the blend
  T1′ = (H₁·T1 + s·Tsurf)/(H₁ + s) at the surface temperature; ablation
  removes the top layers without changing their tracers; freeze-on adds
  b = −bmlt·dt to layer 10 with the blend at the bed temperature; basal
  melt removes the bottom layers; fully removed layers get tracer = 0; a
  column that melts out entirely accumulates the excess into
  `melt_potential` (a per-step diagnostic, zeroed at the top of every
  `mass_balance_driver` call). The esm model is the layer-free continuous
  limit: d(thck_raw)/dt = net_applied with thck = max(0, thck_raw) — the
  unclipped state makes the clip exact at step boundaries for the
  constant-net dump steps; melt_potential = max(0, −thck_raw) reproduces
  the Fortran's per-step value for a step starting from an ice-covered
  state (the tests pin single-step spans, so the Fortran's per-step reset
  never enters). The tracer laws keep the layer-1 temperature, the
  layer-10 temperature, and the column mean; the ablation mean law removes
  the ablated ice at the layer-1 temperature, which is exact while the
  ablation stays within layer 1 (the dump's pinned ablation cells satisfy
  this except step 4 of cell 25,11, where the bite reaches layer 2 and the
  residual is small only because T₂ ≈ T₁ there — documented in the test).
  The applied rates acab_applied/bmlt_applied use the exact per-step forms
  via the pre-step thickness H_pre = thck_raw − net_applied·dt (valid at
  step boundaries, which is where the tests evaluate them).
- **The ocean mask**: `glissade_get_masks` with thklim = 0
  (`glissade_masks.F90:70-96`) — ocean = topg < eus ∧ thck ≤ 0; the SMB
  and BMB are not applied there and melt_potential stays 0 (the
  `effective_areafrac = 0` of ocean cells, F90:893-901).
- **The ice age** (`glissade.F90:1203-1213`): age += dt where thck > 0,
  age = 0 on ice-free cells; the tracer surface/bed values age_usrf = 0,
  age_lsrf = bottom-layer age (`glissade_transport.F90:204-208`). The esm
  ODE stops the increment on ice-free cells but cannot express the
  instantaneous reset — documented difference (unobservable in the window).
- **SurfaceUpdate** (`glimmer_utils.F90:344-378`): lsrf = topg where
  grounded (topg − eus ≥ −(rhoi/rhoo)·thck), lsrf = eus − (rhoi/rhoo)·thck
  where floating, usrf = lsrf + thck.

## How each law is pinned

The dump rows: 32 cells × 4 steps × 2 phases, upn = 11 (temp_000 = surface,
temp_001..010 = layers 1..10, temp_011 = bed); steps 1,2 ↔ cap_step 1,
steps 3,4 ↔ cap_step 2; every row carries time = 0.5·step yr. All test
inputs come from the dump entry phases (the pre-step state) and the
cism_cap.csv import rows; all expected values come from the dump exit
phases or from the source formulas (stated in the test descriptions).

| Law | Pin | Tolerance | Residual |
|---|---|---|---|
| qsmb→acab | all 64 cap import rows (tests at cells 25,11 / 26,11 / 28,11 / 29,13) | rel 1e-12 | ≤ 1.1e-16 (exact) |
| glad zeroing | cell 28,11 step 2 (thck = 0, raw acab = −0.7075 → 0; dump-verified); cell 28,11 step 3 (thck = 1.65e-5 > 0 → NOT zeroed; dump-verified); synthetic ocean cell | exact | 0 |
| d(thck)/dt = net | accumulation: cell 26,11 steps 1/2/4 (rel 1e-6 / 1e-9 / 1e-6); ablation: cell 25,11 steps 3/4 (rel 2e-4 / 3e-4); full melt: cells 28,11 (steps 1,3), 29,13 (step 1) → thck = 0 exactly | see left | advection: +9.6e-7 / +9.6e-10 / −1.7e-6 m; +1.8e-5 / +1.9e-5 m; 0 |
| full-melt clip + melt_potential | cells 28,11 (s1: 0.33379825832198695 m; s3: 0.8837069386692753 m), 29,13 (s1: 0.22629773133841058 m) — melt_potential is derived (not a dump column), from the Fortran's own conservation check | rel 1e-12 | exact |
| acab_applied / bmlt_applied | full-application branch (cell 26,11 s1: acab_applied = acab); clipped branch (cell 28,11 s1: acab_applied = −H_pre/dt = −0.039916999638080596); freeze-on branch (bmlt_applied = bmlt) | rel 1e-12 | exact (source-formula: not dump columns) |
| layer-1 accumulation blend | cells 26,11 s1 (rel 1e-5), 27,11 s1 (rel 7e-5; s > H₁), 25,12 s1 (rel 2.5e-4) | see left | −1.7e-6 / −1.2e-5 / −3.6e-5 degC |
| ablation mean law | cell 25,11 s3 (rel 6e-3), s4 (rel 1.5e-4) | see left | +9.3e-4 / +1.8e-5 degC |
| freeze-on blend + mean | synthetic (bmlt = −0.5 m/yr at the cell-26,11 s1 state) | rel 1e-12 | exact |
| tracer zeroing on melt | cell 28,11 s1 (exit temp_001..010 = 0 in the dump) | exact | 0 |
| ocean mask (glissade-level) | synthetic (topg = −100 m, thck = 0: nothing applied) | exact | 0 |
| ice age | synthetic: d(age)/dt = 1 over 2 yr → age = 2; ice-free → 0 (starts at 0) | rel 1e-12 | exact |
| lsrf/usrf | grounded: cells 25,11 s1, 28,11 s3 (thck = 0), 29,13 s4 (thin); floating: synthetic | rel 1e-12 | exact (all 256 rows: lsrf = topg, usrf = lsrf + thck to round-off) |

The tolerance convention: 1e-12 where the esm law is exact (pure algebra
or linear-in-time integration), 5-10× the measured residual where the dump
exit value carries a physics effect the pointwise model does not contain
(horizontal advection, vertical remap) — each documented in the test
description, following the MOSART scheme-gap pattern.

## Coverage gaps (deferred or unpinnable)

- **Horizontal transport and the first-order vertical remap** are grid
  operations (`glissade_transport.F90:100-272, 1101-1441`) and belong to
  the stage-3 assembly with EarthSciDiscretizations. Their fingerprints are
  the advection residuals of the thickness and tracer pins above.
- **bmlt = 0 everywhere in the reference window** (bmlt_ground =
  bmlt_float = 0 on all 256 rows: the cold base melts nothing and the
  window has no shelf ice). The basal-melt (> 0) removal branches and the
  freeze-on (< 0) branch are transcribed from the source and pinned only
  synthetically (the freeze-on blend test); no dump pin exercises them.
- **No ice-age dump column**: the age increment is pinned from the source
  formula only, and the Fortran's ice-free reset (age = 0) is not
  representable as an ODE — documented in the equation description.
- **acab_applied / bmlt_applied / melt_potential / ocean cells** are not
  dump columns (the dump carries the input acab, not the applied values):
  pinned from the source formulas at the full-melt pins.
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
  esm form max(0, −thck_raw) accumulates the unclipped deficit. The two
  agree for any single step starting from an ice-covered state, which is
  what the tests pin.

## Bugs

No CISM bugs found in this slice; nothing added to `bugs.md`.
