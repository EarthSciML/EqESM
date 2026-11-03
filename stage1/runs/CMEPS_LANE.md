# CMEPS lane (mediator atmosphere–ocean bulk fluxes, stage 1 + stage 2)

Stage-1 instrumentation and stage-2 transcription of the CMEPS atmosphere–ocean
bulk flux computation `flux_atmOcn_large` (`cesm/flux_atmocn/flux_atmocn_Large.F90`,
the Large & Pond / Large & Yeager scheme of the CESM mediator, cesm3.0-alphabranch
checkout), following the EqIce/EqESM dumps convention. Machine-readable run record:
`cmeps_instrument.json` in this directory.

## Scheme and run configuration

The traced run is the `C_JRA_t62` DATM+MOM6 compset (one month, 24 hourly coupling
steps) with the flux settings of `CaseDocs/nuopc.runconfig`:

| Setting | Value | Consequence for the traced laws |
|---|---|---|
| `ocn_surface_flux_scheme = 0` | flux_atmOcn_large | the Large & Pond / Large & Yeager pointwise chain |
| `add_gusts = .false.` | no gustiness | `ugust = 0`, `vmag = wind0`, `rainc` unused |
| `coldair_outbreak_mod = .false.` | no Mahrt & Sun wind boost | `vscl` branch not exercised |
| `aofluxes_use_shr_wv_sat = .true.` | shr_wv_sat path | analytic Goff & Gratch (1946) es + epsilo·es/(p−omeps·es) qsat (`wv_sat_use_tables = .false.`) |
| `flux_convergence = 0.01`, `flux_max_iteration = 5` | do-while of F90:242-290 | 5-update fixed-point iteration, exit test |(ustar_k − ustar_{k−1})/ustar_k| ≤ 0.01 |
| `minwind = 0.5` m/s | `seq_flux_atmocn_minwind` | vmag/wind0 floor |

The instrumented case is `stage1/cime/cases/C_JRA_t62_inst`, cloned from the
reference `C_JRA_t62` case. Only deltas: `user_nl_cpl` (dump gating) and
`user_nl_mom` (`AUTO_MASKTABLE=False`, which matches the reference case's
`MOM_input` and fixes a MOM6 init FATAL in the instrumented case; `buildnml`
regenerates `MOM_override` at `case.build`).

## Fork commits (local only, never pushed)

Fork worktree `stage1/EqESM-cesm-cmeps`, branch `instrument/cmeps-dumps` (cmeps
submodule on the same branch):

| Commit (cmeps submodule) | Change |
|---|---|
| `10741c34` | dump header: `iter` column position fixed to right after `gid` (data write order); subnormal flush to 0.0 |
| `4750efc5` | `es26.17e3` float format (gfortran `es26.17` omits the exponent letter, bugs.md #21); `rainc` fallback init `aoflux_in%rainc(1) = 0.0` (bugs.md #20) |
| `bd88a843` | dump passes a sanitized `rainc` (= 0 when `add_gusts = .false.`) instead of the out-of-bounds `rainc(n)` read of the mediator's 1-element dummy array (bugs.md #20) |

Pointer commits in `EqESM-cesm-cmeps`: `5d367e0`, `0403cf3`, `735eb8b`.

## The dump

`stage1/dumps/cmeps/cmeps_aoflux.csv` — gated by the mediator attribute
`dump_cmeps` (default `.false.`), PET 0 only, 12 masked ocean cells × 24 coupling
steps = 288 rows, 67 columns: the 12 raw inputs, the 11 derived inputs, the 4
neutral first estimates, the 16-field exit state of the stability iteration
(`iter` included), the 7 fluxes and the 10 diagnostics, plus point metadata
(`n`, `gid`, `lat`, `lon`). `# description` / `# columns` headers; floats written
`es26.17e3` with subnormals flushed (both header and format are gfortran-workaround
decisions documented in the header). `hol`/`xsq`/`xqq` are the zbot-level
pre-2m-rescale values (`hol0`/`xsq0`/`xqq0` of the instrumentation), `hol2`/`xsq2`/
`xqq2` the 2-m diagnostics.

## Validation

- **Neutrality.** The instrumented run's restarts are md5-identical to the
  reference run's: `cpl.r 0f47e2f0f4f207240148d8ef333d5044`, `datm.r
  e14609a64cbba942d56a9eaa961c4bda`, `dice.r 8e1656e609b0c189b68284b940b2a19b`,
  `drof.r 002ca0b2a42c135246a5134e2517e00e`, `mom6.r
  2f19de039ced94bcd30930385cc29b76`; `ocean.stats.nc` byte-identical.
- **Full replay.** `stage1/analysis/validate_cmeps_dump.py` replays
  `flux_atmOcn_large` from the dumped inputs and reproduces all 46 computed
  columns on all 288 rows to ≤ 4.4e-15 relative (machine precision) and the
  dumped `iter` on every row. The replay mirrors the stage-2 esm structure:
  all 5 updates unrolled, `conv_k = |(ustar_k − ustar_{k−1})/ustar_k| ≤ 0.01`
  (k = 1..4), `iter` = first converged k else 5, exit state selected by `iter`.
- **Structure.** 2 header lines, 288 rows, 67 columns, zero wrong-arity rows,
  zero NaN, zero spval-like values (`rainc = 0` on every row after the #20 fix).
- **Coverage (and its limits).** All 288 rows are unstable Antarctic-coast cells
  with low-to-moderate wind: `hol ∈ [−10.0, −1.16]` (21 rows pinned at the
  |hol| ≤ 10 clamp), `vmag ∈ [0.5, 4.8]` m/s (17 rows pinned at the minwind
  floor), ocean currents small but nonzero (`us`/`vs` ∈ ±0.13 m/s, so the
  stress direction cosines bite), iteration counts 2 (145 rows), 3 (85), 4 (37),
  5 (21). NOT exercised: the stable branch (stable0 = 1, the `rhn = 0.018`
  ratio, the −5·hol stable psi), the cdn 33 m/s wind cap, gustiness, the
  cold-air-outbreak branch. The stable-branch laws are transcribed in the
  templates but carry no dump pins — a future run over a stable window (e.g. a
  daytime tropical ocean slice) would close the gap.

## Stage-2 handoff

- `stage2/cmeps_templates.esm` — 11 expression templates (Goff & Gratch svp,
  svp→qsat, Large & Yeager cdn, Paulson psimhu/psixhu, psi_momentum/psi_scalar,
  stability_parameter with the |hol| ≤ 10 clamp, the exchange-coefficient
  shift and its clamped wind-shift variant, neutral_heat_coefficient).
- `stage2/cmeps_aoflux.esm` — model `AoFlux`, 136 equations: derived inputs,
  neutral first estimates, the 5-update unrolled iteration (16 states each),
  `conv_1..4` + `iter`, the 16 exit-state selections, fluxes, diagnostics.
- `stage2/README-cmeps.md` — the human-readable map (transcription notes,
  divergences, testing table).
- 6 tests × 52 assertions = 312 assertions, all green under `./esm test stage2/`
  (561 total with the MOSART/CICE suites) at `rel 1e-12`, `time_span {0,0}`,
  pinning dump rows at iterations 2/3/4/5/5/2 including the minwind-floor row
  (vmag = 0.5) and the hol-clamp row (hol = −10).

## Bugs found in the CESM code

- bugs.md #20: `med_phases_aofluxes_mod.F90` `set_aoflux_in_pointers` allocates
  `aoflux_in%rainc(1)` uninitialized when `add_gusts = .false.` (benign in
  pristine code — rainc is never read on the no-gusts path — but the dump
  exposed it as an out-of-bounds read).
- bugs.md #21: gfortran `es26.17` omits the exponent letter for subnormals and
  3-digit exponents (unparseable CSV); `es26.17e3` fixed.

## Status

Stage 1 complete (dump validated, solution-neutral, documented); stage 2 complete
(model + templates + tests green, 312 new assertions). Next: import `cmeps_aoflux`
by reference into the top-level `eqesm.esm` coupling (stage 3) and extend the
dump window to a stable/strong-wind slice for the untested branches.
