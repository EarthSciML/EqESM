# CISM stage 2 — basal traction, effective pressure, flow factor, cap coupling boundary

Stage-2 implementation of the CISM ice-sheet component as built in CESM
(`components/cism` = cismwrap_3_00_002, source `source_cism` =
cism_main_3.00.002, cesm3.0-alphabranch checkout `stage1/cesm`) as esm 1.3.0
files: the three pointwise laws of the reference run's basal/thermodynamic
physics that CISM_LANE.md §2 marks **same law** as EqIce — the pseudo-plastic
basal drag (row 8), the Bueler–van Pelt effective pressure (row 9), and the
Paterson–Budd flow factor (row 4) — plus the CESM-only **cap coupling
boundary** (row 15, missing in EqIce). The shared math is **imported by
reference from the EqIce stage-2 libraries** (`../../eqice/stage2/`,
read-only, PISM v2.3.2-derived); the CISM-specific pieces — the bed-elevation
φ ramp, the bwat=0 → N=overburden branch, the min(·, P) cap and floating
zeroing, the T* = clamp(T−T_pmp, −50, 0) + 273.15 construction, the cap
tsfc −tkfrz conversion with the −250 °C sentinel clamp, and the CISM
constants (arrmlh/arrmll/actenh/actenl/gascon, tkfrz, glide_types defaults) —
live here. The `tests` blocks pin the stage-1 instrumentation dumps
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
| `cism_templates.esm` | Expression-template library (3 match-less templates): `cism_phi_ramp` — the PSEUDO_PLASTIC_OLD bed-elevation friction-angle ramp (φ = phimin below bedmin, phimax above bedmax, linear ramp between); `cism_tempcor` — the Paterson–Budd pressure-corrected temperature clamp(temp − pmptemp, −50, 0) degC; `cism_tsfc_import` — the cap import tsfc conversion tsfc = ifelse(tsfc_raw − tkfrz < −250, 0, tsfc_raw − tkfrz) | `glissade_basal_traction.F90`: 256-277; `glissade_therm.F90`: 2596-2598; `glc_import_export.F90`: 326-333 |
| `cism_cap.esm` | Models `GlcCapImport` (the cap import transformations: tsfc −tkfrz conversion + −250 °C sentinel clamp via `cism_tsfc_import`, qsmb pass-through) and `GlcCapExport` (the export laws: topo = usrf on the active grid, ice_covered = (thck > 0) ∧ (usrf > 0), and the coupling-interval averages rofi/rofl/hflx written as storage ODEs — the continuous limit of the Fortran's per-dynamics-step accumulation, see the cap section below) | `glc_import_export.F90`: 326-333; `glad_output_states.F90`: 48-158; `glad_output_fluxes.F90`: 77-172; `glad_main.F90`: 971-1027 |
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
transport (row 6), the local-till hydrology (row 10), and calving (row 12).
The φ ramp is the only basal piece without an EqIce counterpart; the cap
coupling (row 15) is CESM-only and is covered by `cism_cap.esm` (see the cap
section below).

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
3. The temperature-step lane (row 5) supplies pmptemp/stagsigma internally;
   the hydrology lane (row 10) supplies bwat.

## Cap coupling boundary (row 15) — `cism_cap.esm`

The mediator↔CISM coupling layer around `glad_gcm` (CISM_LANE.md §2 row 15,
"missing in EqIce — CESM-specific coupling conversions"). The stage-1 dump
`stage1/dumps/cism/cism_cap.csv` traces the boundary in `source_glc/glc_RunMod.F90`
(instrumented fork): **import** rows with tsfc [degC] and qsmb [kg/m²/s]
*exactly as consumed by CISM* (after the cap `import_fields` processing), and
**export** rows with topo [m], hflx [W/m²], rofi/rofl [kg/m²/s = mm H₂O/s] and
ice_covered [−] *exactly as produced* by `glad_set_output_fields`. Two cap
steps (one per coupling year, GLC_NCPL=1; tstep_mbal = 8760 h; each year = 2
CISM dynamics steps of dt = 0.5 yr), 32 cells × 2 phases × 2 steps = 128 rows.

### Import law (`GlcCapImport`)

`glc_import_export.F90:326-333` (`import_fields`):

- **tsfc** = ifelse(tsfc_raw − tkfrz < −250, 0, tsfc_raw − tkfrz), with
  tsfc_raw the mediator `Sl_tsrf` [K] and tkfrz = SHR_CONST_TKFRZ = 273.15 K
  (`glc_constants.F90:22-23`, `shr_const_mod.F90:37`). The −250 °C clamp is
  the "Jer hack" sentinel (F90:328-332): CLM-ocean land points send raw 0 K,
  which the conversion would drop to −273.15 °C; the clamp resets them to
  0 °C. Written as the `cism_tsfc_import` template (the subtraction appears
  twice because the template body has no let-binding; the Fortran converts
  and clamps in place).
- **qsmb** passes through unchanged (`state_getimport`, F90:679-714, applies
  no conversion).
- The mediator→model halo conversion (`vector_to_spatial` /
  `parallel_convert_nonhaloed_to_haloed`) is a grid operation, deferred to
  the stage-3 assembly.
- The model-side continuation — the mbal averaging of qsmb/tsfc over the
  coupling interval and the qsmb → acab conversion
  (acab = qsmb·tstep_mbal·3600/1000·rhow/rhoi m/yr ice,
  `glad_main.F90:800-824`, `glad_timestep.F90:200`, `glide_set_acab`) — is the
  surface-forcing/transport lane's boundary (already pinned by the stage-1
  `smb_chain` validator check and the transport dumps), not the cap's.

### Export law (`GlcCapExport`)

`glad_set_output_fields` (`glad_main.F90:971-1027`) → `set_output_states`
(`glad_output_states.F90:48-158`) + `glad_output_fluxes`
(`glad_output_fluxes.F90:77-172`):

- **topo** = usrf on the active grid (usrf > 0, `is_in_active_grid`
  F90:110-133), 0 elsewhere (F90:91, 97-99).
- **ice_covered** = 1 where thck > min_thck = 0 (`glad_constants.F90:63`,
  F90:152) **and** usrf > 0, else 0 (F90:82-86).
- **rofi** = coupling-interval average of ρᵢ·(calving_thck + removal_thck)/dt
  (F90:104-105) — calving + the other solid-ice removals (icebergs/isthmuses,
  forced retreat, tiny-thickness cleanup, F90:91-96).
- **rofl** = coupling-interval average of ρᵢ·bmlt + ρᵢ·melt_thck/dt (F90:128-129)
  — basal melt + lateral/submarine melt.
- **hflx** = coupling-interval average of −ucondflx (F90:138-139; ucondflx
  positive down, `glissade_therm.F90:669`). The sign-flip comment says
  "positive up", which contradicts the "positive down" labels of
  `glad_main.F90:625` and the CMEPS field dictionary — see `bugs.md` (TBD
  entry); unpinned by the dump because hflx = 0 throughout this window.

**Time average as a storage ODE, not `integral`.** The Fortran forms each flux
as a sum over the n_icetstep = 2 dynamics steps of the coupling interval
(reset at the interval start, `glad_timestep.F90:149`; one accumulation per
`glissade_tstep`, F90:290), then divides by the step count
(`calculate_average_output_fluxes`, F90:151-159). With piecewise-constant
integrands over the steps, Σₖ fₖ/count = (1/T_coupling)·∫₀^T_coupling f dt
with T_coupling = count·dt = 8760 h = 3.1536e7 s — so the export is the
coupling-interval average in the integral sense. It is written here as the
storage ODEs D(S_rofi)/dt = ρᵢ·(calving_rate + removal_rate) (and
analogously S_rofl, S_hflx) with the export S/T_coupling: the `integral` op of
esm 1.3.0 is a *spatial* PIDE term that the CLI cannot lower or evaluate
(probed: `unlowered_operator`, esm-spec §4.2/§9.6.8 — a file using it
validates but cannot simulate), so the differential form — the same
continuous-limit pattern the MOSART lane uses for `D(wh)/dt` — is used
instead. The ODE is exact for the Fortran's piecewise-constant integrands.
Model time is in seconds (the dump's thour column × 3600; cap_step 1 at 0 h,
cap_step 2 at 8760 h). The Fortran exports the average only at the interval
end; the .esm evaluates S/T_coupling at any t, and the reset of S at each
interval boundary is a discrete event for the stage-3 integration.

**Inputs the dump does not carry.** The integrands — calving_thck,
removal_thck, bmlt, lateral_melt_thck, ucondflx — are per-dynamics-step CISM
state; `cism_cap.csv` carries only the final interval averages (all zero in
this window: the 2-year cold-base run has no calving, no basal/lateral melt,
and ucondflx = 0). They enter as parameters. usrf/thck are likewise not in
the cap dump and are taken from the transport dumps in the topo tests
(cross-file identity topo == transport-exit usrf verified on all 64 export
rows). `ice_sheet_grid_mask` (the same active-grid gate) is not dumped and is
not modeled.

### Cap tests (14 tests, 43 assertions)

All instantaneous pins are rel 1e-12 `time_span {0,0}`; the average-law pins
integrate one coupling interval [0, 3.1536e7 s] with constant integrands
(linear ODEs — exact to round-off).

| Test group | Tests | What is pinned |
|---|---|---|
| Import conversion (dump) | 5 | tsfc = tsfc_raw − 273.15 and the qsmb pass-through at cap_step 1/2 rows (cells 25,11 / 28,11 / 19,13 / 27,14 / 28,12, incl. the negative-SMB cell 28,11, the coldest −9.678 °C and warmest −0.004 °C traced tsfc); the raw input is reconstructed as tsfc_dump + 273.15 (no dump row hits the clamp) |
| Import clamp (source formula) | 4 | the sentinel firing (raw 0 K — the Jer hack case), the strict-< boundary at exactly −250 °C, just below it, and the warm linear branch |
| Export states (dump) | 6 | topo = usrf and ice_covered = 1 at both cap steps (cells 25,11 / 26,12 / 19,13 / 28,11), incl. the thinnest traced ice (thck = 1.65e-5 m, the min_thck = 0 boundary exercised from above) |
| Export states (source formula) | 2 | the ice-free branch (thck = 0 → ice_covered = 0, topo unchanged) and the inactive-grid branch (usrf = 0 → topo = ice_covered = 0) |
| Flux averages | 4 | the dumped zeros (all 64 export rows have hflx = rofi = rofl = 0), and the averaging law with constant integrands: solid (rofi = ρᵢ·1.75 m/yr, incl. an intermediate S_rofi pin at t = T/2), liquid, heat, and all three combined |

### Coverage limitations (this data window)

- **The −250 °C clamp never fires in the dumps** (traced tsfc ∈ [−9.68,
  −0.004] °C — all CLM land/ice points); pinned synthetically from the source
  formula.
- **rofi/rofl/hflx are identically zero** in the window (no calving, no
  melt, ucondflx = 0); the dump pins the zero case, the averaging law is
  pinned with synthetic constant integrands. A wetter/dynamic stage-1 window
  would dump-pin nonzero averages (and the hflx sign).
- **ice_covered = 1 on every export row** and usrf > 0 everywhere; the
  ice-free and inactive-grid branches are source-formula pins.
- **hflx sign convention** is unpinned (all zeros) — bugs.md TBD entry.
- The qsmb → acab SMB chain is the transport lane's pin
  (`smb_chain` validator check; the cap dump's qsmb is its input).
