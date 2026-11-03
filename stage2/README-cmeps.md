# CMEPS stage 2 — `.esm` bulk flux laws

Stage-2 implementation of the CMEPS atmosphere–ocean bulk flux computation
(cesm3.0-alphabranch checkout `stage1/cesm/components/cmeps/cesm/flux_atmocn/`,
`flux_atmocn_Large.F90`) as esm 1.3.0 files: the expression-template library
and the `AoFlux` component transcribing the pointwise laws of the Large & Pond /
Large & Yeager scheme exactly, with `tests` blocks pinning the stage-1
instrumentation dump (`stage1/dumps/cmeps/cmeps_aoflux.csv`, 288 rows = 12
PET-0 ocean cells × 24 coupling steps of the `C_JRA_t62` DATM+MOM6 run). See
[Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/cmeps/cesm/flux_atmocn/flux_atmocn_Large.F90`
(the traced scheme), `stage1/cesm/share/src/shr_wv_sat_mod.F90`
(`shr_wv_sat_qsat_liquid` / `shr_wv_sat_svp_to_qsat` / `GoffGratch_svp_liquid`),
`stage1/cesm/components/cmeps/mediator/med_phases_aofluxes_mod.F90` (the
`add_gusts` / `minwind` / `seq_flux_atmocn_minwind` configuration), and the run
record in `stage1/runs/cmeps_instrument.json` / `stage1/runs/CMEPS_LANE.md`.

Configuration transcribed (the `C_JRA_t62` run's `nuopc.runconfig`):
`ocn_surface_flux_scheme = 0` (flux_atmOcn_large), `add_gusts = .false.`,
`coldair_outbreak_mod = .false.`, `aofluxes_use_shr_wv_sat = .true.` (the
shr_wv_sat Goff & Gratch path; `wv_sat_use_tables = .false.` in the traced run,
so the analytic form is used, not the lookup table), `flux_convergence = 0.01`,
`flux_max_iteration = 5`, `minwind = 0.5` m/s.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `cmeps_templates.esm` | Expression-template library (13 templates): the Goff & Gratch (1946) liquid-water saturation vapor pressure (`shr_wv_sat_mod.F90:942-956`) and the epsilo·es/(p−omeps·es) saturation-humidity conversion with the p ≤ es guard (`F90:351-363`), the Large & Yeager (2009) neutral 10-m drag coefficient cdn(U) = 0.0027/min(33,U) + 0.000142 + 0.0000764·min(33,U) − 3.14807e-13·min(33,U)⁶ (`flux_atmocn_Large.F90:153-154`), the Paulson (1970) unstable psimhu/psixhu (`F90:163-164`), the combined stable/unstable psi_m/psi_h (`F90:258-259`), the Monin-Obukhov stability parameter (unclamped, with karman/g/zvir as named parameters, `F90:252-253`) and its sign-preserving |hol| ≤ 10 clamp `stability_parameter_clamp` (`F90:254`), the exchange-coefficient height shift c/(1 + c/karman·(alz − psi)) (`F90:282-284`) and its clamped variant with max(..., −0.5) used for the wind-speed shift (`F90:272`) — both with karman as a named parameter, the neutral heat coefficient (1−stable)·0.0327 + stable·0.018 (`F90:238`), and the frozen-state selection (input-state pinning for the composed iteration) | `flux_atmocn_Large.F90`, `shr_wv_sat_mod.F90` |
| `cmeps_aoflux_step.esm` | Model `AoFluxStep` — ONE pass of the stability do-while (`F90:242-290`), authored once and mounted 6× by reference (models.step1..step5 + exit of `cmeps_aoflux.esm`): 22 equations = the 17 update states (hol_raw, the unclamped Monin-Obukhov parameter of `F90:252-253`; hol = stability_parameter_clamp(hol_raw), the |hol| ≤ 10 clamp of `F90:254`; stable, xsq, xqq, psimh, psixh, the clamped wind-shift rdw, u10n, rdn, rhn, rd, rh, re, ustar, tstar, qstar), the per-update convergence flag conv = |(ustar − ustar_prev)/ustar| ≤ 0.01, and 4 frozen-state fields pinning the INPUT state of the first-converged update (ustar_prev/tstar_prev/qstar_prev/rdn_prev, via the `frozen_state_select` template); 20 parameters = the 7 fixed pointwise fields (vmag, delt, delq, alz, zbot, thbot, qbot) + ren, the CESM constants karman/g/zvir (defaults, not coupled — the ren precedent), the previous state (ustar_prev, tstar_prev, qstar_prev, rdn_prev), conv_prev, and the 4 frozen threading inputs | `flux_atmocn_Large.F90` |
| `cmeps_aoflux.esm` | Model `AoFlux` — the full pointwise chain, composed: 55 equations = the 12 raw inputs; the derived inputs vmag (minwind-clamped relative wind), wind0 = vmag (the resolved wind of `F90:204`; with add_gusts = .false. the Fortran computes it identically to vmag, so the alias is exact), esat_raw/esat (the min(es, pslv) cap), qsat (from the RAW es, as shr_wv_sat_qsat_liquid calls svp_to_qsat before the cap), ssq = 0.98·qsat, delt, delq, alz, cp (with the CESM constants cpdair/cpvir/latvap/stebol/karman declared once as named parameters); the neutral first estimates stable0/rdn0/rhn0/ren0/ustar0/tstar0/qstar0 (`F90:236-244`); the convergence representation conv_1..4 = step1..4.conv and iter = first converged k, else 5; the 16 exit-state fields = the 15 exit.* fields + ren (the iteration of `F90:242-290` is reproduced by the 6 top-level mounts of `./cmeps_aoflux_step.esm` chained by 86 variable_map couplings: ustar0/tstar0/qstar0/rdn0 seed step1, each step's state/conv/frozen-input fields feed the next, the 7 fixed fields feed every step; the 6th mount exit is fed step5's 4 pinned input-state fields + the 7 fixed fields and recomputes the 15 exit-state fields — the exit state of the first-converged update k is a deterministic function of the pre-converged state (update k−1 outputs) plus the fixed fields); the 7 fluxes (`F90:299-311`); the 10 diagnostics (`F90:316-333`) | `flux_atmocn_Large.F90` |

The `ren = 0.0346` water coefficient is constant across the iteration and is
carried once (`ren0`, the final `ren`, and each AoFluxStep mount's `ren`
parameter bind the same 0.0346 value; the parameter is not coupled). The same
"declare once, not coupled" pattern carries the other CESM constants:
`AoFlux` declares karman = 0.4, cpdair = 1004.64 J/kg/K, cpvir = 1810/1004.64 − 1,
latvap = 2.501e6 J/kg and stebol = 5.67e-8 W/m²/K⁴ as named parameters used by
the cp/lat/lwup/evap/fac_t/fac_q equations, and each AoFluxStep mount declares
karman/g/zvir (the stability-parameter constants) as parameters with the same
defaults; none are coupled (shr_const_mod.F90:22, 32, 46, 48, 52, 21, 31).

## Units

All fields are SI; temperatures are Kelvin; humidity kg/kg (the unit checker's
algebra treats kg/kg as dimensionless, which is how qsat/ssq/delt/… close);
stresses N/m²; fluxes W/m²; water flux kg/m²/s.

| Quantity | Symbols | Units |
|---|---|---|
| Heights / pressures | `zbot`, `ztref`/`zref`, `pslv`, `esat`, `esat_raw` | m; Pa |
| Atmosphere state | `ubot`, `vbot`, `us`, `vs`, `vmag`, `wind0`, `ugust`, `u10n`, `u10res`, `duu10n` | m/s; m²/s² |
| Temperatures | `thbot`, `tbot`, `ts`, `delt`, `tstar`, `tref` | K |
| Humidity | `qbot`, `qsat`, `ssq`, `delq`, `qstar`, `qref` | kg/kg (dimensionless in the checker) |
| Coefficients | `rdn`, `rhn`, `ren`, `rd`, `rh`, `re`, `rdw`, `stable`, `stable0` | 1 |
| Stability | `hol`, `hol2`, `xsq`, `xqq`, `xsq2`, `xqq2`, `psimh`, `psixh`, `psix2`, `fac_t`, `fac_q` | 1 |
| Fluxes | `tau`, `taux`, `tauy`; `sen`, `lat`, `lwup`; `evap` | N/m²; W/m²; kg/m²/s |
| Constants | `minwind = 0.5` m/s, `zref = 10` m, `ztref = 2` m, `karman = 0.4`, `grav = 9.80616` m/s², `cpdair = 1004.64` J/kg/K, `cpvir = 1810/1004.64 − 1`, `latvap = 2.501e6` J/kg, `stebol = 5.67e-8` W/m²/K⁴, `zvir = 28.966/18.016 − 1`, `epsilo = 18.016/28.966`, `omeps = 1 − epsilo`, `tboil = 373.16` K (Goff & Gratch) | (shr_const_mod.F90; the nuopc-cap values the run uses) |

**Unit notes**

- The `tref` potential-to-actual correction `0.01·ztref` (F90:323) is a lapse
  rate: the 0.01 is declared `K/m`, so `0.01 [K/m]·2 [m] = 0.02 K` closes.
- The Goff & Gratch formula evaluates in hPa and the trailing `×100` converts
  to Pa; in the unit algebra the 100 is declared `Pa` (the dimensionless 10^()
  factor multiplies it), per the template description.
- The empirical cdn coefficients carry the SI dimensions that make each term
  dimensionless (0.0027 [m/s], 0.0000764 [s/m], 3.14807e-13 [s⁶/m⁶]) — the
  template closes as a dimensionless coefficient.
- `sign(0.5, x)` (the Fortran two-argument sign) is emulated as the 0/1
  selection `ifelse(x >= 0, 1, 0)` for both `stable0` (F90:237, on delt) and
  the per-iteration `stable` (F90:255, on hol); the sign-preserving |hol| ≤ 10
  clamp of F90:254 is transcribed as `stability_parameter_clamp` =
  `sign(hol_raw)·min(|hol_raw|,10)` (the unary esm sign returns ±1, and 0 at
  hol_raw = 0 where 0·min = 0 anyway), which coincides bit-identically with the
  conditional `ifelse(hol_raw >= 0, min(|hol_raw|,10), −min(|hol_raw|,10))` on
  every value.

## Pointwise laws — Fortran vs esm transcription notes

| Law | Fortran | esm |
|---|---|---|
| Wind magnitude | `vmag = max(minwind, sqrt((ubot-us)^2 + (vbot-vs)^2))` (F90:193-196; no gusts) | `max(sqrt((ubot−us)² + (vbot−vs)²), 0.5 m/s)` |
| Saturation vapor pressure | `es = 10^(−7.90298·(tboil/t−1) + 5.02808·log10(tboil/t) − 1.3816e-7·(10^(11.344·(1−t/tboil))−1) + 8.1328e-3·(10^(−3.49149·(tboil/t−1))−1) + log10(1013.246))·100` (shr_wv_sat_mod.F90:942-956) | `goff_gratch_liquid_svp(ts)` template; `esat = min(esat_raw, pslv)` at the call site |
| Saturation humidity | `qs = epsilo·es/(p − omeps·es)`, `qs = 1` when `p ≤ es` (F90:351-363) | `svp_to_qsat(esat_raw, pslv)` — bound to the RAW es (the Fortran caps es AFTER computing qsat, F90:439-441) |
| Stability parameter | `hol = karman·g·zbot·(tstar/thbot + qstar/(1/zvir + qbot))/ustar²`, then `hol = sign(min(|hol|,10), hol)` (F90:252-254) | `stability_parameter(ustar, tstar, qstar, thbot, qbot, zbot, karman, g, zvir)` (unclamped) + `stability_parameter_clamp(hol_raw)` templates |
| Stability functions | `psimh = −5·hol·stable + (1−stable)·psimhu(xqq)`, xqq = sqrt(sqrt(|1−16·hol|)) (F90:256-259) | `psi_momentum` / `psi_scalar` templates |
| Wind-speed shift | `rd = rdn/(1 + max(rdn/karman·(alz−psimh), −0.5))`; `u10n = vmag·rd/rdn` (F90:272-273) | `exchange_coefficient_shift_clamped` (momentum only, wind-shift step only) |
| Coefficient shift | `rd/rh/re = c/(1 + c/karman·(alz−psi))` (F90:282-284) | `exchange_coefficient_shift` |
| Fluxes | `tau = rbot·ustar²`; `taux/tauy = tau·(ubot−us)/vmag` etc.; `sen = cp·tau·tstar/ustar`; `lat = latvap·tau·qstar/ustar`; `lwup = −stebol·ts⁴`; `evap = lat/latvap` (F90:299-311) | identical |
| 2-m diagnostics | `hol2 = hol·ztref/zbot`; `psix2` at (hol2, stable, xqq2); `tref = thbot − delt·(rh/karman)·(alz+al2−psixh+psix2) − 0.01·ztref`; `qref = qbot − delq·(re/karman)·(…)`; `duu10n = u10n²`; `u10res = u10n·wind0/vmag` (F90:316-333) | identical; the diagnostics consume the exit-state rh/re/psimh/psixh/stable, exactly as the Fortran leaves them after the do-while |

## Documented divergences (all annotated in file metadata)

1. **Iteration unrolled and exit-selected, composed by reference.** The
   Fortran do-while of F90:242-290 runs until |(ustar − ustar_prev)/ustar| ≤
   flux_convergence (checked after each update, at most flux_max_iteration =
   5). esm 1.3.0 has no loop construct, so the 16-equation update law is
   authored ONCE in the one-step model `AoFluxStep` (`cmeps_aoflux_step.esm`)
   and mounted six times by reference in `AoFlux` as top-level models
   step1..step5 + exit (esm-spec §4.7), chained by 86 variable_map couplings.
   Each mount threads the previous state (ustar/tstar/qstar/rdn), the previous
   convergence flag, and 4 frozen-state fields that pin the INPUT state of the
   FIRST-converged update (ustar_prev/tstar_prev/qstar_prev/rdn_prev, via
   `frozen_state_select` in `cmeps_templates.esm`): the Fortran exits at the
   first passing update k, leaving the update-k exit-state fields, which are
   deterministic functions of the PRE-converged state (the update k−1 outputs)
   plus the fixed fields — and the iteration states do not repeat after
   convergence (ustar_2 differs from ustar_4 at the 1.8e-4 level on the pinned
   rows), so the input state must be carried forward through the remaining
   mounts. The 6th mount `exit` is fed step5's 4 pinned fields + the 7 fixed
   fields (conv_prev stays at its 0 default, so the frozen selection inside it
   is the identity) and its 16 fresh unknowns ARE the exit state; the parent
   aliases the 15 exit-state fields from `exit.*` (all except the intermediates
   hol_raw/rdw and the unused conv), which coincide bit-identically with the
   values the Fortran leaves in the iteration variables. The parent exposes
   the same external names as the pre-refactor unrolled model, so the pinned
   tests are unchanged: `conv_1..4 = step1..4.conv`, `iter` = the first
   converged k, else 5, and the 16 exit-state fields = the 15 `exit.*` fields
   + ren. The all-5-iterations case keeps conv_prev = 0 throughout, so the
   frozen inputs are the fresh 4th-update values, which the exit mount passes
   through the update law to produce the 5th-update values — exactly what the
   do-while leaves in the iteration variables after its final pass. The
   stage-1 replay (`stage1/analysis/validate_cmeps_dump.py`) mirrors this
   exact structure and reproduces the dumped `iter` and exit state on all 288
   rows, so the representation is pinned, not just the laws.
   - **Why top-level mounts and not subsystems.** The pinned CLI (esm 0.7.0)
     silently drops variable_map couplings whose targets live under
     `subsystems` in the single-top-level-model build path
     (`build_array_compiled` builds via `from_file` when the document has one
     top-level model and applies couplings only when flattening documents
     with several top-level models). Probes in `scratch/` confirmed the §4.7
     top-level form works (`probe_toplevel.esm`, 8/8) while the subsystem form
     leaves values at their defaults.
   - **Why not a coupled recurrence.** A direct coupled-recurrence formulation
     (each update referencing the previous update through the mount boundary)
     trips the §4.9.6 cycle rule; the CLI reports: `Observed dependency cycle:
     x -> y -> x. Each of these observed variables is defined in terms of the
     next, so no evaluation order satisfies every definition (esm-spec
     §4.9.6). Break the cycle by splitting one observed into a pre-value and a
     post-value.` The step model therefore splits every chained state into a
     pre-value parameter (x_prev) and a post-value unknown (x) — precisely the
     split the rule prescribes. `scratch/recur_probe2.esm` is the minimal
     failing probe (a single-variable recurrence, `scratch/recur_probe.esm`,
     passes 4/4, but the coupled multi-variable form does not).
2. **rainc unused.** With `add_gusts = .false.` the physics never reads rainc;
   the stage-1 dump sanitizes it to 0 after the mediator's 1-element dummy
   allocation was exposed as an out-of-bounds read (bugs.md #28). The
   gustiness law `ugust(min(rainc, 6.94444e-4))` of F90:197-202 is not
   transcribed (no dump coverage); `ugust = 0` is a plain equation.
3. **Coverage limits.** All 288 dump rows are unstable Antarctic-coast cells
   with low-to-moderate wind, so the stable branch (stable0 = 1, the
   `rhn = 0.018` ratio, the −5·hol stable psi of F90:258-259), the cdn 33 m/s
   wind cap, the gustiness path and the cold-air-outbreak branch carry no dump
   pins — the stable-branch and cdn laws are transcribed in the templates
   (they are exercised by construction only in the psi_momentum/psi_scalar and
   neutral_heat_coefficient branches). The minwind = 0.5 m/s floor and the
   |hol| ≤ 10 clamp ARE pinned (rows with vmag = 0.5 and hol = −10).
4. **esat vs qsat ordering.** The dump's `esat` column is the capped
   min(es_raw, pslv) while `qsat` is computed from the uncapped es_raw — the
   model carries both (`esat_raw` internal, `esat` capped output); the two
   coincide in the traced rows (es ≈ 700 Pa ≪ pslv ≈ 99.8 kPa), so the cap is
   not exercised.

## Testing

The stage-2 `tests` block pins the stage-1 instrumentation dump
(`stage1/dumps/cmeps/cmeps_aoflux.csv`) with instantaneous laws
(`time_span {0,0}`, MOSART pattern): 6 tests × 52 assertions = 312 assertions,
all `rel 1e-12`, byte-identical to the pre-refactor unrolled model (the
composition keeps the same external variable names). The full stage-2 suite —
`./esm test` on `cmeps_aoflux.esm` plus the CICE and MOSART files — reports
561 assertions (312 aoflux + 249 CICE/MOSART); this component needs no
`ESD_ROOT`.

| test | dump row (step, gid) | what is pinned |
|---|---|---|
| `iter2_typical_unstable` | (1, 41583) | the full chain at the typical traced state: convergence at iteration 2 (145/288 rows), vmag 2.72 m/s, hol −1.84; all 47 dumped columns + conv_1..4 + iter |
| `iter3_cold_cell` | (1, 42662) | the coldest traced cell: convergence at iteration 3, vmag 1.33 m/s, hol −4.12 |
| `iter4_cold_cell` | (3, 42662) | convergence at iteration 4, vmag 0.80 m/s, hol −7.20 (approaching the clamp) |
| `iter5_hol_clamp` | (8, 42662) | all 5 iterations run (no convergence) with hol pinned at the |hol| ≤ 10 clamp (−10.0 exactly), vmag 0.52 m/s near the minwind floor |
| `iter5_minwind_clamp` | (22, 42662) | all 5 iterations run with vmag pinned at the minwind = 0.5 m/s floor while hol = −9.76 stays just off the clamp |
| `iter2_late_step` | (24, 41583) | the same cell as the first test 23 coupling steps later (convergence at iteration 2 again, vmag 2.22 m/s, hol −2.58) |

Every test overrides the 12 raw inputs at full dump precision and asserts the
11 derived inputs, the 4 neutral first estimates, the 16-field exit state, the
7 fluxes, the 10 diagnostics, and the convergence representation
(`conv_1..4` from the replay, `iter` from both). The replay
(`stage1/analysis/validate_cmeps_dump.py`) reproduces all 288 dump rows to
≤ 4.4e-15 relative on every computed column, so `rel 1e-12` pins the dump to
round-off. Row identity (step, gid, local n, lat, lon) and the exercised laws
are cited in each test description.

### Coverage summary

| branch | traced rows | pinned |
|---|---|---|
| unstable psi (Paulson psimhu/psixhu, stable = 0) | 288 | ✓ |
| |hol| ≤ 10 clamp | 21 rows at hol = −10 | ✓ |
| minwind = 0.5 m/s floor | 17 rows at vmag = 0.5 | ✓ |
| iteration convergence 2 / 3 / 4 / 5 | 145 / 85 / 37 / 21 | ✓ |
| nonzero ocean currents (us/vs) | all rows (|u| ≤ 0.13 m/s) | ✓ |
| stable psi (−5·hol branch), rhn = 0.018 | 0 | ✗ gap |
| cdn 33 m/s cap | 0 | ✗ gap |
| gustiness (rainc → ugust), cold-air outbreak | 0 | ✗ gap |

### Transcription fixes found while writing the tests

1. **Exit state vs last unrolled update.** The first draft of the replay/selector
   returned the state after the final unrolled update; the dump carries the
   state at the CONVERGED iteration (the Fortran leaves the last-update values
   in the variables when the do-while exits), which for the iter-2 rows differs
   at the 1e-2 level (hol −3.683 vs −3.626 at the traced (11, 42124) row). The
   component selects the exit state by `iter`, and the validator replay mirrors
   the selection; with the fix all 288 rows reproduce to ≤ 4.4e-15.
2. **svp_to_qsat guard units.** The template's `(p − es) <= 0` guard declared
   the zero dimensionless, which only fails dimensional analysis once the
   template is bound to typed variables (p = pslv Pa, es = esat_raw Pa); the
   zero is declared Pa.
