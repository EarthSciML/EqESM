# CISM stage 2 — the DIVA velocity-solve pointwise laws (`cism_velo.esm`)

Stage-2 implementation of the **pointwise** laws of CISM's higher-order
velocity solve as built in CESM (`components/cism` = cismwrap_3_00_002,
source `source_cism` = cism_main_3.00.002, cesm3.0-alphabranch checkout
`stage1/cesm`), reference run T1850Gg/f09_g17_gris20 with `which_ho_approx =
HO_APPROX_DIVA = 4`, `which_ho_gradient = 0`, `margin = 1`,
`assemble_taud = 1`, `diva_slope_correction = .true.`, `which_ho_babc = 3`,
`max_slope = 0.5`, `thck_gradient_ramp = 50 m`, `thklim = 1 m`, `upn = 11`
(10 layers), `nquad = 4`, `rhoi = 917`, `g = 9.80616`, `dx = dy = 20000 m`.
What is implemented here is everything in `glissade_velo_higher_solve` that is
a *local law of a single column*: the vertex geometry of the driving stress
and the DIVA slope correction, the Goldberg (2011) vertical integrals of the
σ column, and the 3-D velocity recovery from the 2-D solve.  The finite-
element matrix solve itself, the vertex staggers of the per-level quantities,
and the Picard iteration are **stage-3 assembly** concerns and are left to the
full `eqesm.esm`; this file feeds their inputs from the dumps.

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/cism/source_cism/libglissade/`
(`glissade_velo_higher.F90`, `glissade_grid_operators.F90`,
`glissade_basal_traction.F90`, `glissade_calc_effecpress`), the run namelists
(`stage1/cime/cases/T1850Gg_cism_inst/run/cism.gris.config`), the stage-1
instrumentation dumps `stage1/dumps/cism/cism_velo.csv` +
`cism_velo_nbr.csv` (2026-11-05 rerun, job 12118489, fork commits
`db624ec`/`76d170d`, branch `instrument/cism-dumps`, not pushed — see
`stage1/runs/CISM_LANE.md` §4d), and `stage1/runs/cism_instrument.json`.

## File map

| File | Contents |
|---|---|
| `cism_velo.esm` | the three pointwise components `DivaGeometry`, `DivaIntegrals`, `DivaVelocity` (13 tests, 136 assertions), plus the top-level mount of `BasalTraction` and the `variable_map` edge `BasalTraction.beta → DivaVelocity.beta_raw` |
| `cism_velo_discretization.esm` | local ESD-format σ-column rules: `NLAYERS = 10`, `lev`/`lev_nodes`, the `cism_sigma_measure` / `cism_sigma_midpoint` templates, and the `integral_lev_whole` / `integral_lev_nodes_cumulative_to_top` rules imported by reference from `EarthSciDiscretizations` with `NLEV = NLAYERS` |
| `cism_templates.esm` (lines 376–742) | the DIVA expression templates `cism_hybrid_edge_gradient`, `cism_vertex_gradient`, `cism_driving_stress`, `cism_theta_slope_x`, `cism_diva_slope_factor`, `cism_beta_eff` |
| `cism_basal.esm` | the mounted `BasalTraction` (pseudo-plastic law, its own 2P/0F tests live in that file) |
| `stage1/analysis/extract_cism_velo_test_rows.py` | the committed tuple extractor: verifies every law below on **all 128 exit rows** of the dumps, then emits the numeric `parameter_overrides`/`assertions` fragments keyed by test id |

The shared DIVA templates live in `cism_templates.esm` because this file is
self-contained for review; the stage-3 assembly will import them by reference.

## DivaGeometry — the vertex geometry of the driving stress

The vertex `(i,j)` is the **NE corner** of the traced cell; the stencil is the
SW/SE/NW/NE cells around it.  Inputs are the four cells' `usrf`, `thck`,
`ice_mask`, `land_mask` plus the run constants (`thklim`, `ramp`, `max_slope`,
`dx`, `rhoi`, `g`).  Laws (all in `glissade_grid_operators.F90` /
`glissade_velo_higher.F90`):

- `z_mean_* = usrf − thck/2` per stencil cell, and the **plain vertex
  gradient** `dz_dx = (z_mean_ne + z_mean_se − z_mean_sw − z_mean_nw)/(2·dx)`
  (and `dz_dy` with the n/s pairing, `glissade_slope_angle_staggered`'s
  argument, F90:1510).
- `dusrf_dx/dy` = the **hybrid-margin surface gradient**
  (`surface_elevation_gradient`, HO_GRADIENT_MARGIN_HYBRID, F90:999-1079):
  edge gradients from the interior pairs where both stencil cells are ice,
  else the `thklim`-guarded single-cell difference, then the vertex gradient
  from the four edge means, with the `max_slope` clamp (F90:1234-1252).  The
  `ramp` smoothing (50 m) applies only in the interior of the ice sheet where
  `thck` differences are large; it is off in the pinned rows by construction
  (checked: `dusrf_dx` matches the unramped hybrid formula on every exit
  row).
- `theta_x = atan(|dz_dx|)`, `theta_y = atan(|dz_dy|)` (slope_mask = 1
  everywhere on gris20), and the DIVA slope factor
  `F = cos²θ + 4·sin²θ` (F90:2189-2193).
- `stagthck = (thck_sw + thck_se + thck_nw + thck_ne)/4` (the 4-cell mean of
  `glissade_stagger`, `ice_plus_land_mask`).
- `taudx = −rhoi·g·stagthck·dusrf_dx`, `taudy = −rhoi·g·stagthck·dusrf_dy`
  (`load_vector_gravity`, HO_ASSEMBLE_TAUD_LOCAL, F90:4353-4384).
- `effecpress_stag = rhoi·g·stagthck` — the `bwat = 0` identity of the
  effective-pressure staggering (`HO_EFFECPRESS_BWAT_BVP` with `bwat = 0`
  everywhere in this 2-year run), pinned here at the same stencil rather than
  as an input.

## DivaIntegrals — the Goldberg (2011) σ-column integrals

`compute_integrals_diva` (F90:5570-5696), all inside the
`if (active_cell(i,j))` guard (F90:5091):

- `efvs` = the **arithmetic** mean of the four quadrature-point effective
  viscosities, `efvs[l] = (efvs_q1[l] + efvs_q2[l] + efvs_q3[l] + efvs_q4[l])/4`
  (F90:5266-5272 — a mean, not a harmonic mean).
- `omega_k[k]` = the interface-to-bed cumulative integral of
  `sig_mid[l]·thck·(1/4)·Σ_p 1/efvs_qp[l]` — written as **one ESD `integral`
  operator** over the σ measure (the integrand `omega_integrand` factors out
  the quad-point mean of reciprocals by integral linearity), lowered by the
  `integral_lev_nodes_cumulative_to_top` rule (the exact suffix scan of
  F90:5631-5637).  `omega_k[11] = 0` at the surface.
- `omega_layer_avg[l] = (omega_k[l] + omega_k[l+1])/2` (the trapezoid mean,
  F90:5656) and `omega` = the whole-column `integral` of `omega_layer_avg`
  over the σ measure (`integral_lev_whole`; the exact trapezoid sum of
  F90:5655-5661 — note **no `thck` factor**: this `omega` is Goldberg's
  ω/H).
- `stag_omega = (omega + omega_se + omega_nw + omega_ne)/4` — the plain
  4-cell mean of the cell-centered omegas (`glissade_stagger`,
  `ice_plus_land_mask`, stagger_margin = 1, F90:5790-5798).
- **Inactive-cell guard**: on ice-free cells the Fortran skips the whole
  block, so `efvs` and `omega_k` are 0 there (the dumped `efvs_*` columns
  confirm 0 while the `efvs_q*` columns carry **stale** previous-step values).
  The esm equivalent is a scalar-condition `ifelse(ice_mask == 1, …, 0)`,
  which is lazy in the spec (§4.8) — the integral branch is never evaluated
  on inactive cells, so the 1/efvs_qp blow-up is avoided; `omega_layer_avg`
  and `omega` need no guard (0 propagates).  The inactive pin is the step-3
  cell (17,14) test, where `thck` collapsed 97.3 → 0.986 m between steps.

## DivaVelocity — the 2-D recovery and the 3-D profile

`compute_3d_velocity_diva` (F90:5700-5821) and the calcbeta chain:

- `beta_internal = f_ground·beta_raw`, with `beta_raw` fed by the top-level
  `variable_map` edge `BasalTraction.beta → DivaVelocity.beta_raw`
  (`param_to_var`): the mounted `BasalTraction` recomputes the pseudo-plastic
  law (φ(bed) ramp → tanφ → τ_c = tanφ·N → β = τ_c/(u₀^q·speed^(1−q))) from
  the dumped `u, v, bed, N_stag` — the same inputs the Fortran had at the
  last accepted calcbeta pass.
- `beta_eff_x = F_x·β/(F_x + β·ω)` and `beta_eff_y` analogously
  (the `cism_beta_eff` template, F90:5756-5770), with `F_x/F_y`, `omega`,
  `theta_x/theta_y`, and the per-level `stag_omega_k` fed **from the dump**:
  the vertex stagger of these per-level/vertex quantities is part of the
  stage-3 assembly, not of the pointwise column laws.
- **Basal recovery** (Goldberg 2011 eq. 34, F90:5801):
  `uvel_basal_recovered = uvel_2d·beta_eff_x/beta_internal`,
  `vvel_basal_recovered = vvel_2d·beta_eff_y/beta_internal` — pinned as
  explicit unknowns so the recovery law is tested rather than assumed.
- `btractx = β_internal·uvel_basal = beta_eff_x·uvel_2d`,
  `btracty` analogously (exit form F90:3874, solve form F90:5780); and the
  vertical profile `uvel[k] = uvel_basal + btractx·stag_omega_k[k]·cosθ_x/F_x`
  (F90:5805-5821) via the `slope_corr_x/y = cos(θ)/F` factors.
  `uvel_basal`/`vvel_basal` are **parameters** (the dumped exit
  `uvel_011`/`vvel_011`): the 3-D solve that produces them is the stage-3
  assembly.  Because `stag_omega_k_011 = 0` on every traced row, the pinned
  profiles satisfy `uvel[11] = uvel_basal` automatically — the same check the
  Fortran performs.

## Testing

`./esm test stage2/cism_velo.esm` runs 13 tests / 136 assertions
(whole suite: 2004P/0F/0E, esm 0.7.0; `./esm units --check` clean apart from
the pre-existing BasalTraction fractional-exponent notices).  Pinned rows
(from the 128 exit rows of the 2026-11-05 dumps; all tuples at full dump
precision):

| test id | row | why this row |
|---|---|---|
| `geom_step1_cell2511` | step 1, cell (25,11) | the NE-corner stencil of the standard traced column (thin ice, 1.37 m) |
| `geom_step1_cell2713` | step 1, cell (27,13) | thick ice (434.5 m) with a negative surface slope in x |
| `geom_step3_cell2414` | step 3, cell (24,14) | the fastest traced column (uvel_2d = 5.02 m/yr) |
| `geom_step3_cell2813` | step 3, cell (28,13) | the steepest pinned slope (θ_x = 0.043) |
| `int_step1_cell2511` / `cell2612` / `cell2713` | step 1 | thin / mid (83.6 m) / thick columns: the σ-column integrals across the thickness range |
| `int_step3_cell2414` | step 3 | fastest column's integrals |
| `int_step3_cell1714_inactive` | step 3, cell (17,14) | **inactive-cell guard**: `ice_mask` dropped to 0 between steps (thck 97.3 → 0.986 m) with stale nonzero `efvs_q*` inputs |
| `velo_step1_cell2511` / `cell2713` / `step3_cell2414` | — | the recovery + profile laws on the standard rows |
| `velo_step1_cell2213` | step 1, cell (22,13) | traced cell **ice-free** (thck = 0) with an active vertex — exercises the full 11-level profile through the β chain on a zero-thickness column |

Tolerances: rel 1e-12 everywhere except:

- the **β chain** (`beta_internal`, `beta_eff_x/y`, `btractx/y`, and every
  velocity derived from them) at rel 1e-5, and
- `dusrf_dx/dy` + `taudx/taudy` at rel 1e-11.

The `dusrf`/`taud` loosening is the dump round-trip of ~1800 m elevation
differences (worst deviation 6.6e-13 over the 128 rows; >15× margin).  The β
chain loosening has a real physical reason: the dumped `beta_internal` is the
last **accepted** calcbeta pass of the Picard iteration, while the exit
velocities come from the next iteration's accepted `beta_eff`; the
law-at-exit-velocity chain (β from the dumped inputs → recovery → profile)
has worst rel deviation 6.0e-7 over the 128 rows (verified by the extractor),
so 1e-5 leaves >16× margin.

## Coverage limitations (this data window)

- **The 2-D finite-element matrix solve, the vertex staggers, and the Picard
  iteration are not here** — they are the stage-3 assembly.  Consequently
  `uvel_2d`, `vvel_2d`, `uvel_basal`, `vvel_basal`, `omega`,
  `stag_omega_k_001..011`, `theta_x/y`, `F_x/y` are dump-fed inputs in this
  file (independent-model tests), and `effecpress_stag` is pinned via the
  bwat = 0 identity instead of being computed.
- **`stag_omega_k` stagger gap**: `omega_k` is computed per cell in
  `DivaIntegrals`, but its vertex stagger (`stag_omega_k`) is a dumped input
  in `DivaVelocity`; only the scalar `stag_omega` (4-cell mean) is
  implemented.  The per-level stagger is stage 3.
- **`efvs_q*` provenance**: on inactive cells the dumped `efvs_q*` columns
  are stale (the last values computed while the cell was active) — the tests
  feed them as-is and the guard must zero the outputs, exactly as the
  Fortran's `active_cell` branch does.
- **`f_ground = 1` everywhere** on gris20 in this window, so the
  `beta_internal = f_ground·beta_raw` weighting is pinned only in its
  saturated branch; the subgrid floatation weighting itself
  (F90:595-630) is not exercised by this data window.
- The reference run's velocities are in **solver-internal m/yr** (the
  `scale_output` pass runs after the exit trace); all pinned numbers are in
  those units.

## Next steps (stage 3)

Couple `DivaGeometry` + `DivaIntegrals` + `DivaVelocity` (+ `BasalTraction`)
into the top-level `eqesm.esm`: implement the 2-D DIVA matrix solve
(`glissade_velo_higher_solve`'s Picard loop), the vertex staggers
(`glissade_stagger` of the per-level `omega_k`, `theta`, `F`), and the
`effecpress` staggering, then test against full-model simulations with
`earthsciintegrator` time stepping.
