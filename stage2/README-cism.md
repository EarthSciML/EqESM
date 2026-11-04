# CISM stage 2 — basal traction, effective pressure, flow factor

Stage-2 implementation of the CISM ice-sheet component as built in CESM
(`components/cism` = cismwrap_3_00_002, source `source_cism` =
cism_main_3.00.002, cesm3.0-alphabranch checkout `stage1/cesm`) as esm 1.3.0
files: the three pointwise laws of the reference run's basal/thermodynamic
physics that CISM_LANE.md §2 marks **same law** as EqIce — the pseudo-plastic
basal drag (row 8), the Bueler–van Pelt effective pressure (row 9), and the
Paterson–Budd flow factor (row 4). The shared math is **imported by reference
from the EqIce stage-2 libraries** (`../../eqice/stage2/`, read-only, PISM
v2.3.2-derived); the CISM-specific pieces — the bed-elevation φ ramp, the
bwat=0 → N=overburden branch, the min(·, P) cap and floating zeroing, the
T* = clamp(T−T_pmp, −50, 0) + 273.15 construction, and the CISM constants
(arrmlh/arrmll/actenh/actenl/gascon, glide_types defaults) — live here. The
`tests` blocks pin the stage-1 instrumentation dumps
(`stage1/dumps/cism/`) — see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/cism/source_cism/libglissade/`
(`glissade_basal_traction.F90`, `glissade_therm.F90`,
`glissade_grid_operators.F90`), `source_cism/libglimmer/glimmer_physcon.F90`,
`source_cism/libglide/glide_types.F90`, the run namelists
(`stage1/cime/cases/T1850Gg_cism_inst/run/cism.gris.config`), and the CISM→EqIce
mapping of `stage1/runs/CISM_LANE.md` §2.

Reference configuration pinned (from `cism.gris.config` of the T1850Gg run and
the source defaults): `which_ho_babc = 3` (HO_BABC_PSEUDO_PLASTIC_OLD),
`which_ho_effecpress = 4` (HO_EFFECPRESS_BWAT_BVP), `flow_law = 2`
(FLWA_PATERSON_BUDD), `which_ho_flow_enhancement_factor = 0` (E = 1),
`pseudo_plastic_q = 0.5`, `pseudo_plastic_u0 = 100 m/yr`, the source defaults
`pseudo_plastic_phimin/phimax = 5/40 deg`, `pseudo_plastic_bedmin/bedmax =
−300/700 m` (glide_types.F90:2380-2383, not exposed in the namelist),
`bwat_till_max = 2 m`, `N_0 = 1000 Pa`, `e_0 = 0.69`, `C_c = 0.12`,
`effecpress_delta = 0.02` (glide_types.F90:2297-2330), `rhoi = 917 kg/m³`,
`grav = 9.80616 m/s²`, `pmlt = 9.7456e-8 K/Pa`, and the flow-law constants
`arrmlh = 1.733e3`, `arrmll = 3.613e-13` (Pa⁻³ s⁻¹), `actenh = 139e3`,
`actenl = 60e3` J/mol, `gascon = 8.314` J/mol/K (glimmer_physcon.F90:89-104).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `cism_templates.esm` | Expression-template library (2 match-less templates): `cism_phi_ramp` — the PSEUDO_PLASTIC_OLD bed-elevation friction-angle ramp (φ = phimin below bedmin, phimax above bedmax, linear ramp between); `cism_tempcor` — the Paterson–Budd pressure-corrected temperature clamp(temp − pmptemp, −50, 0) degC | `glissade_basal_traction.F90`: 256-277; `glissade_therm.F90`: 2596-2598 |
| `cism_basal.esm` | Models `Effecpress` (cell-centered N: overburden P = ρᵢgH, the Bueler–van Pelt N0(δP/N0)^s·10^((e0/Cc)(1−s)) law with s = clamp(bwat/bwat_till_max, 0, 1) imported from EqIce `basal_templates.esm` `mohr_coulomb_N_interior`, the bwat=0 → N = P branch, the min(·, P) cap and the floating-cell zero) and `BasalTraction` (φ ramp → tanφ, speed = √(u²+v²+0.01²), τ_c = tanφ·N_stag, β = τ_c/(u0^q·speed^(1−q)) — the drag law itself imported from EqIce `basal_templates.esm` `pseudo_plastic_drag`) | `glissade_basal_traction.F90`: calc_effecpress 668-1114 (BWAT_BVP 944-985, init 784, stagger 1096-1098, cap 1082-1086), calcbeta 69-664 (PSEUDO_PLASTIC_OLD 250-292, smallnum 136/189) |
| `cism_flowlaw.esm` | Model `FlowFactor` — flwa = E·a·exp(−Q/(RT*)) with T* = tempcor + 273.15, the (a, Q) pair switching at T* = 263.15 K; the softness law imported from EqIce `ice_rheology_templates.esm` `pb_softness`, the CISM constants (arrmlh/arrmll/actenh/actenl/gascon) bound at the call site, the T* construction CISM-specific | `glissade_therm.F90`: 2408-2640 (FLWA_PATERSON_BUDD 2585-2606, pmptemp 2170); `glimmer_physcon.F90`: 89-102 |

The cell→vertex stagger of the effective pressure (`effecpress_stag`, the
plain 4-cell average with `stagger_margin = 0`, glissade_grid_operators.F90:125-133
via glissade_basal_traction.F90:1096-1098) is a grid discretization step: it
belongs to the stage-3 assembly with the velocity solve, not to these pointwise
law components. In this lane the staggered N enters `BasalTraction` as the
input `N_stag` (the traced N column of `cism_basal_traction.csv`, identical to
the `effecpress_stag` column of `cism_velo.csv` at the same vertex).

## EqIce reuse map (CISM_LANE.md §2 rows 4, 8, 9)

| CISM law | EqIce import | CISM-specific additions here |
|---|---|---|
| Effective pressure (row 9): N = N0·(δP/N0)^s·10^((e0/Cc)(1−s)) | `basal_templates.esm` `mohr_coulomb_N_interior` (the identical functional form, Bueler & van Pelt 2015 Eq. 23 / Tulaczyk 2000; free names N0, W_till_max, e0, Cc resolved to the glide_types defaults) | s = clamp(bwat/bwat_till_max, 0, 1) (F90:960; applied to the water thickness before the law, algebraically identical since W_till_max > 0); **bwat = 0 → N = overburden** (F90:784 initialization; EqIce's law has no such branch — PISM always evaluates the formula); min(·, P) cap (F90:974-975); floating-cell zero (F90:982-985); P = ρᵢgH via EqIce `overburden_pressure` |
| Basal drag (row 8): β = τ_c/(u0^q·speed^(1−q)) | `basal_templates.esm` `pseudo_plastic_drag` (β = τ_c·(ε²+u²+v²)^((q−1)/2)·u0^(−q) with ε = 0.01 m/yr — exactly CISM's formula; the speed regularization is EqIce's `basal_speed`) | τ_c = tanφ·N with φ from the **bed-elevation ramp** (EqIce derives φ from data/inversion), using the cell-center topg (the stagtopg quirk, bugs.md #23); tanφ = tan(φ·π/180); the ramp parameters are glide_types defaults |
| Flow factor (row 4): flwa = E·a·exp(−Q/(RT*)) | `ice_rheology_templates.esm` `pb_softness` (same law, same constants family) | T* = clamp(T−T_pmp, −50, 0) + 273.15 (CISM's T* construction; PISM takes the pressure-adjusted temperature directly); the branch predicate is T* ≥ 263.15 K — CISM's `tempcor >= −10` (F90:2602) and EqIce's `T_pa < T_critical` agree at the boundary; **CISM's arrmll = 3.613e-13 differs from PISM's A_cold = 3.61e-13 in the third digit**, so the constants are bound at the CISM call site (A_cold = arrmll, Q_cold = actenl, A_warm = arrmlh, Q_warm = actenh, T_critical = 263.15 K, R_ideal_gas = gascon); E multiplies outside the law (F90:2579, 2603, 2605) |

Not in this lane (later lanes / different laws): the DIVA velocity solve
(row 2), the temperature step (row 5, incl. the pmptemp column law
pmptemp = −ρᵢg·pmlt·thck·stagsigma, F90:2170, whose value is an input here),
transport (row 6), the local-till hydrology (row 10), calving (row 12), and the
cap coupling (row 15). The φ ramp is the only basal piece without an EqIce
counterpart; everything else is a binding.

## Pointwise laws — Fortran vs esm transcription notes

- **PSEUDO_PLASTIC_OLD** (glissade_basal_traction.F90:250-292): bed =
  topg − eus with the **cell-center topg** (the code comment at F90:261-263
  admits the stagtopg bug; kept for backward compatibility, recorded in
  bugs.md #23); φ from the ramp; tanφ = tan(φ·π/180); τ_c = tanφ·N_stag;
  speed = √(u²+v²+smallnum²) with smallnum = 1.0d-2 m/yr (F90:136,189);
  β = τ_c/(u0^q·speed^(1−q)). The speed cap of F90:196-198
  (`beta_powerlaw_umax > 0`) and the β floor/f_ground weighting of
  F90:597-633 (`beta_grounded_min = 100 Pa yr/m`, `which_ho_beta_limit = 0`)
  are applied **after** the traced law: not transcribed here because they are
  solve-assembly steps (and no-ops on the traced rows: f_ground = 1 and
  β > 100 at every traced vertex; beta_powerlaw_umax = 0 in the config).
- **BWAT_BVP** (F90:944-985): relative_bwat = max(0, min(bwat/bwat_till_max,
  1)); N = N0·(δ·P/N0)^relative_bwat·10^((e0/Cc)(1−relative_bwat)); N =
  min(N, P); floating cells zeroed. The bwat > 0 branch only fires where
  bwat > 0 — **the reference run never fires it** (see How N is pinned).
  The final 0..P capping of F90:1082-1086 is a no-op on this construction and
  is documented, not transcribed.
- **FLWA_PATERSON_BUDD** (F90:2585-2606): tempcor = min(0, temp − pmptemp),
  tempcor = max(−50, tempcor); flwa = E·a·exp(−Q/(R·(tempcor + 273.15))) with
  the (a, Q) pair switching at tempcor ≥ −10; E = flow_enhancement_factor
  (1.0 in the reference run). The waterfrac enhancement of F90:2608-2610 and
  the damage branch are off in the reference config (waterfrac = 0,
  no damage).

## How N is pinned (the effecpress trace is empty)

`cism_effecpress.csv` has **no rows in this run**: the gland20 initial state
has bwat = 0 everywhere and the 2-year cold-base run produces no basal melt
(cism_thermo.csv: bwat = bmlt_ground = 0 on all 256 rows), so the bwat > 0
branch of glissade_calc_effecpress never fires. The effective pressure is
therefore pinned three ways:

1. **The bwat = 0 branch from the dumps**: with bwat = 0, N (cell) =
   overburden P = ρᵢg·thck (the F90:784 initialization). The Effecpress tests
   `bwat0_overburden_*` assert N = ρᵢg·thck at the velo-entry thck of three
   dump cells (steps 1 and 4).
2. **The staggered value from the dumps**: at every traced vertex,
   effecpress_stag = ρᵢg·stagthck to 2e-16 relative (verified across all 4
   velo-entry rows checked: cells 25,11 / 26,12 / 28,11 step 2 / 19,13 step 4)
   — the full N chain (N = ρᵢgH → stagger_margin=0 4-cell average) is
   consistent with the traced N columns; the stagger itself is the stage-3
   grid operation.
3. **The bwat > 0 branches from the source formula**: the tests
   `bwat_partial_plateau`, `bwat_saturated`, `bwat_clamped` and
   `floating_zero` evaluate the F90:960-985 law at the glide_types defaults
   with synthetic bwat (0.5, 2.0, 3.0 m) and dump-derived overburden; expected
   values are hand-computed from the source formula (documented in each test
   description), in the same spirit as the MOSART Fortran-driver pins for
   branches the dump window never reaches. The trace machinery fires whenever
   bwat > 0 (a wetter stage-1 window would replace these with dump pins).

## Testing

All `tests` blocks are `time_span {0,0}` instantaneous pins with rel 1e-12
(the dumps are written `es26.17`; the recomputations land at ≤ 4.5e-16
relative, so 1e-12 is never stretched). Dump rows are cited by
(file, cism_step, cell i,j, velo_iter/layer) in every test description.

| File | Tests | Assertions | What is pinned |
|---|---|---|---|
| `cism_templates.esm` | — | — | template library only (validated, no models) |
| `cism_basal.esm` | 15 | 41 | **Effecpress** (7 tests): the bwat=0 overburden branch at three dump cells/steps, the Bueler–van Pelt mid-range plateau (min(·,P) cap at work), full saturation, the s clamp, and the floating zero — the last four from the source formula; **BasalTraction** (8 tests): φ ramp + tanφ + speed + τ_c + β on dump rows covering both ramp-branch (φ = 34.33°, bed = 538 m) and phimax-branch (φ = 40°) cells, a 3-iteration Picard chain of one vertex (N/φ/τ_c invariant, β tracking the speed), the largest-N vertex (N = 1.79e6 Pa), a step-2 row, and the last iteration of the last step |
| `cism_flowlaw.esm` | 10 | 19 | warm-branch flwa at 6 (cell, layer, step) dump rows (incl. the largest |pmptemp| = 0.072 K traced), the cold branch, the branch boundary at T* = 263.15 K exactly, and both tempcor clamp bounds — the last four from the source formula |

Run the suite with

    esm test stage2/

(`esm validate` and `esm units` pass for all three files; the EqIce imports
are read-only and resolve from any EqESM worktree via
`../../eqice/stage2/`). Whole-suite counts (all lanes): see the commit
messages.

### Coverage limitations (this data window)

- **Cold flow-law branch**: all 1280 traced (cell, layer) rows have tempcor in
  [−5.16, −1.13] degC — the −50 floor, the 0 cap, and the tempcor < −10 cold
  branch never fire. Pinned synthetically from the source formula
  (`branch_cold`, `branch_pmp_exceeded`, `branch_floor`, `branch_boundary`).
- **bwat > 0 effective pressure**: empty trace (see How N is pinned).
- **φ = phimin branch**: no traced bed is below bedmin (traced beds span
  538..1786 m), so the ramp's lower plateau is covered only through the
  template structure (the ramp formula is pinned at both ends by the
  ramp/phimax branch tests).
- **Floating cells**: none in the window (all 256 velo rows floating_mask =
  0); the floating zeroing and the f_ground weighting are source-formula
  pins / documented no-ops.
- **Speed cap and β floor** (F90:196-198, 597-633): no-ops on the traced rows
  (beta_powerlaw_umax = 0; β > 100 Pa yr/m and f_ground = 1 everywhere);
  documented, not transcribed.
- **pmptemp column law** (F90:2170): the flowlaw rows carry pmptemp as input
  (the dump has no stagsigma column); the −ρᵢg·pmlt·thck·stagsigma law is
  part of the temperature-solve lane (the stage-1 validator's bpmp check
  pinned its bed value −ρᵢg·pmlt·thck in `cism_instrument.json`).

## Next steps (stage 2.5 / stage 3)

1. Mount `Effecpress` + `BasalTraction` (+ `FlowFactor`) into the CISM
   subassembly with the cell→vertex stagger and the velocity solve; integrate
   with earthsciintegrator against the full-CESM trajectory
   (velo-exit uvel/vvel vs the solve, `btractx/btracty` = β·u·f_ground).
2. A wetter stage-1 window (nonzero bwat / basal melt) to dump-pin the
   bwat > 0 effecpress branch and the cold flow-law branch.
3. The temperature-step lane (row 5) is implemented (see
   [Temperature step](#temperature-step-cism_lanemd-2-row-5) below): it supplies
   pmptemp/stagsigma internally; the hydrology lane (row 10) supplies bwat.

## Temperature step (CISM_LANE.md §2 row 5)

The prognostic temperature solve (`whichtemp = TEMP_PROGNOSTIC`) as
`IceTemperatureColumn` (`cism_thermo.esm`): the vertical heat equation
ρᵢcᵢ ∂T/∂t = ∂/∂z(kᵢ ∂T/∂z) + ρᵢcᵢΦ on CISM's σ coordinate, the pointwise laws of
the step, and the local EarthSciDiscretizations-format rule that lowers the
diffusion term (`cism_thermo_discretization.esm`, a candidate ESD PR). Pinned
against the extended `cism_thermo.csv` dump (256 rows × 78 columns). The dump
extension — `dttem`, `dissip_001..010`, `sigma_001..011`, `stagsigma_001..010`
(fork branch `instrument/cism-dumps`, committed at all three submodule levels,
not pushed) — was required: a dissip = 0 reproduction fails the operator
identity at 3–5% of the step-1 temperature change, and a pin needs the exit
rows and the grid geometry, not just the entry rows. Source:
`glissade_therm.F90` (glissade_therm_driver, glissade_temperature_matrix_elements
F90:1324-1507, glissade_basal_melting_ground F90:1930-2012,
glissade_basal_melting_internal F90:2032-2159, glissade_pressure_melting_point
F90:2163-2192).

### Representation: 12-node σ grid, node-centered finite volumes

CISM's temperature lives on UPN+1 = 12 nodes per column: the surface and bed
faces σ(1) = 0, σ(11) = 1 and the 10 layer midpoints stagsigma(1..10)
(glide_calc_sigma, SIGMA_COMPUTE_GLIDE). The model carries them on `lev_nodes`
with sigma_node = [0, stagsigma, 1]; the control volume of node n is the CISM
layer [σ(n−1), σ(n)] of width thck·dsigma_cv[n] (dsigma_cv[1] = dsigma_cv[12] =
0 — the boundary nodes carry no diffusion row). The diffusion term lowers, by
the rule `cism_sigma_diffusion` (priority 10, match-scoped to lev_nodes-shaped
u), to the exact dups stencil of glissade_temperature_matrix_elements:
L·T[n] = (G[n+1] − G[n])/(thck·dsigma_cv[n]) with the interface flux
G[n] = coni·(T[n] − T[n−1])/(thck·(sigma_node[n] − sigma_node[n−1])) (the
node-to-node gradient; stagsigma(0) ≡ σ(1), stagsigma(upn) ≡ σ(upn), F90:151-163).
Boundary rows: node 1 (surface) lowers to 0 — the matrix row holds the Dirichlet
value T(0) = min(0, artm) (F90:843, 1405-1408) — and node 12 (bed) lowers to 0 —
CISM's bed row is not a diffusion row but the algebraic basal BC (frozen: flux
balance; thawed: T = bpmp; F90:1473-1476, 1497-1504), carried as the `temp_bed`
law. dTdt = L·T/(ρᵢcᵢ) + dissip on the interior nodes; the surface and bed dTdt
rows are 0.

### Pin strategy: the instantaneous operator identity at the exit state

The CISM step is ONE backward-Euler solve of the tridiagonal system
(diag − dt·L)·T_exit = T_entry + dt·Φ. The harness integrators cannot reproduce
one exact stiff backward-Euler step, so the model does not integrate the
entry→exit map; instead the operator is pinned INSTANTANEOUSLY at the exit
state: for the interior rows the identity (T_exit − T_entry)/dttem = L·T_exit + Φ
holds exactly (the dump round-trips at 1e-13..6e-9 relative), and
`IceTemperatureColumn`'s dTdt is asserted against (T_exit − T_entry)/dttem with
the ICs set to the exit profile and the dumped dissip. Five such pins span the
window (step 1 and step 4 of cell (27,14), step 1 of cells (27,13) and (25,11) —
the thinnest solved column at thck = 1.3 m — and step 4 of cell (19,13));
per-node residuals ≤ 6.0e-9 relative → rel 1e-8 tolerances. The algebraic laws
of the step are pinned at rel 1e-12: temp_sfc = min(0, artm) (all 128 exit
rows), bpmp = −ρᵢg·pmlt·thck, the frozen basal row temp_bed = temp[11] +
(bfricflx + bhydroflx − bheatflx)·dsigbot·thck/coni (closes to 6.2e-17 W/m² in
flux terms), lcondflx_ground = −(coni/thck)·(temp_bed − temp[11])/dsigbot (one
thin-column assertion at rel 1e-8 — there the flux is the ratio of two cancelling
differences, temp_bed − temp[11] = 2.5e-8 K from ~1.87 degC values), and the melt
chain bmlt_raw / bmlt_ground / btemp_ground (all-zero in the window: the closed
basal balance |bfricflx + lcondflx_ground| ≤ 1.5e-15 W/m² is zeroed by the eps11
guard). The thawed-bed row, the freeze-on branch, and the internal-melt cap are
pinned synthetically from the source formulas (see gaps); the thin-cell surface
law is pinned on cell (28,11) (thck = 0.02 m ≤ thklim_temp = 1 m: the column
bypasses the solve — temp = min(0, artm) on all 12 nodes, dissip = 0).

### Files

| File | Contents | Fortran source |
|---|---|---|
| `cism_templates.esm` | +5 templates: `cism_pmp_temp` (Clapeyron pmp −ρᵢg·pmlt·depth), `cism_internal_melt_rate` (layer rate rhoi·shci·(T−Tpmp)·layer_thck/(rhoi·lhci·dt) with the max(0, ·) branch guard), `cism_basal_melt_rate` (flux balance with the eps11 zeroing applied to the RATE), `cism_basal_freeze_on` (predicate bmlt < 0 ∧ bwat = 0 ∧ btemp ≥ bpmp − eps08), `cism_basal_melt_cap` (btemp_ground = min(btemp, bpmp); freeze-on drop to min(bpmp − 2·pmp_threshold, bpmp − eps08)) | `glissade_therm.F90`: 2163-2174, 2127-2155, 1930-1967, 2007-2012, 1979 |
| `cism_thermo_discretization.esm` | Local ESD-format rule file: the `cism_sigma_diffusion` lowering rule + the `cism_interface_flux` stencil template, with the UPN/lev_nodes grid contract | `glissade_therm.F90`: 1324-1507, 151-163 |
| `cism_thermo.esm` | Model `IceTemperatureColumn` (the ODE d(temp)/dt = dTdt, the 12-node geometry parameters from the dump, the algebraic laws above) + 10 tests (5 operator pins, thin-cell law, thawed bed, internal melt, freeze-on, eps11) | `glissade_therm.F90` as above |

Constants (CCSMCOUPLED): rhoi = 917 kg/m³, shci = 2117.27 J/kg/K, lhci = 333700
J/kg, coni = 2.1 W/m/K, grav = 9.80616 m/s², pmlt = 9.7456e-8 K/Pa,
scyr = 31536000 s, pmp_threshold = 1e-3 degC, eps08 = 1e-8, eps11 = 1e-11,
dttem = 1.5768e7 s (the 0.5-yr glissade step).

### Gaps and deferrals

- **Entry→exit backward-Euler map**: not integrated (see pin strategy). A
  stage-3 coupling that must reproduce the CISM trajectory exactly would
  transcribe the tridiagonal solve itself (an algebraic map, not an ODE) or
  accept first-order-in-dt trajectory differences.
- **Surface re-pin jump**: the surface value is re-pinned inside the solve, so
  the ENTRY rows of steps 1 and 3 carry the previous value (the initial 0 and
  the previous year's artm); the surface law is pinned on the exit rows (all
  128).
- **waterfrac ≡ 0**: TEMP_PROGNOSTIC carries no internal water (the dump's
  waterfrac columns are all 0); the enthalpy scheme (waterfrac > 0) is deferred.
- **Thin-ice bypass**: cells with thck ≤ thklim_temp = 1 m skip the solve
  entirely (temp = min(0, artm) on all nodes; bpmp/btemp_ground = 0; the basal
  melt laws do not apply). 21 of 128 traced rows are such cells; cell (28,11)
  is pinned.
- **Thawed-bed / freeze-on / internal-melt branches**: no traced row exercises
  them (the window is frozen: every entry temp_011 < bpmp − pmp_threshold, worst
  margin −1.425 degC; bmlt = bwat = bheatflx = bhydroflx = 0); pinned
  synthetically from the source formulas.
- **Internal-melt temperature reset**: the melt cap's paired state update
  temp → pmptemp (F90:2148) is a post-solve mutation, documented, not
  transcribed as an equation.
- **pmptemp profile**: the dump carries only the bed value bpmp; the per-node
  profile law −ρᵢg·pmlt·thck·sigma_node is a model computation, pinned at the
  bed against the dumped bpmp and synthetically for the interior (the
  internal_melt_rate test).
- **Dump extension**: dttem, dissip_001..010, sigma_001..011, stagsigma_001..010
  were added to `cism_thermo.csv`; sigma/stagsigma are derivable from the grid
  construction, but dissip and dttem are not, and the operator pins need them.

### Next steps for this lane

Mount `IceTemperatureColumn` (+ `Effecpress`/`BasalTraction`/`FlowFactor` and
the cap lanes) into the CISM subassembly: wire artm/thck/bwat/bfricflx/
bhydroflx/bheatflx from the surface/velocity components, the dissip heating from
the velocity solve, and the temperature state into the flow factor; the
hydrology lane (row 10) supplies bwat. The basal-melt outputs (bmlt_ground,
btemp_ground) feed the thickness/geometry evolution (not in this window).
