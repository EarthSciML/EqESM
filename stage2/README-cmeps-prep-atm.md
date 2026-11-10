# CMEPS atmosphere-export merge component (stage 2)

This lane closes the last documented gap of `README-cmeps-ocnalb.md` ("the
Sx_*/Faxx_* atmosphere merges of prep_atm need a CAM case"): a second
instrumented CMEPS case, `F2000climo_f19_inst` (CAM6 + CLM + CICE + DOCN +
MOSART, f19 grid, 6 h run, 12 coupling steps at 1800 s), traced with a new
prep_atm phase dump.

| File | Model | Laws | Source |
|---|---|---|---|
| `cmeps_prep_atm.esm` | `AtmMerge` | the fraction-weighted lnd/ice/ocn -> atm merges | `mediator/med_phases_prep_atm_mod.F90:156-166` (with `med_merge_mod.F90:405-412`, `esmFldsExchange_cesm_mod.F90`) |

The model imports `fraction_weighted_sum` from `./cmeps_templates.esm` by
reference. Every atm export is a three-source merge accumulated in the Fortran
term order `dp = 0 + Sl*lfrac, dp = dp + Si*ifrac, dp = dp + So*ofrac`
(`med_merge_mod.F90:405-412`), so the .esm composes `fraction_weighted_sum`
twice: the land+ice pair weighted by `lfrac`/`ifrac` first (the inner
application), then the open-ocean source weighted by `ofrac` added on top (the
outer application with the identity `ice_frac = 1.0`, which multiplies
exactly). The 14 merged exports are `Sx_avsdr`/`Sx_avsdf`/`Sx_anidr`/
`Sx_anidf` (L875-887, 902-914, 929-941, 956-968), `Sx_tref`/`Sx_u10`/
`Sx_qref` (L988-1002, 1015-1029, 1046-1062), `Faxx_taux`/`Faxx_tauy`/
`Faxx_lat`/`Faxx_sen`/`Faxx_evap`/`Faxx_lwup` (L1085-1099, 1112-1126,
1139-1153, 1166-1180, 1197-1211, 1236-1250) and `Sx_t` (L1267-1279); the
15th export `So_t_export` is the unweighted `So_t` copy of L1282-1285
(`mrg_type = 'copy'`). The weights are the `FBFrac(compatm)` fractions:
`lfrac` for the land source (`mrg_fracname_lnd2atm_state/flux = 'lfrac'`,
`med_internalstate_mod.F90:289-290`), `ifrac` for the ice source, `ofrac` for
the ocean source (the traced case has `samegrid_atmlnd = .true.` and
`component_computes_enthalpy_flux = none`).

## The instrumentation fix

`med_dump_prep_atm_active` was originally gated on the *ogrid* aoflux init
flag (`is_init`) and the ogrid traced-cell counter (`step`). On this case
PET 0's f19 ogrid block is all-land (the Antarctica rows come first in the
local element order; `So_omask = 0` everywhere on PET 0), so the ogrid init
bails ("WARNING: no masked cells traced" in ocn.log), `is_init` stays
`.false.` and the prep_atm dump was permanently gated off — a silent
instrumentation failure that produced no rows. The fork fix (cmeps submodule
commit `17eef08f`, pointer `ea48423`):

- the atm-export dump gets its own counters `step3`/`ymd3`/`tod3`;
  `med_dump_prep_atm_tick` advances them from the mediator clock and is
  called by `med_phases_prep_atm_dump`, and `med_dump_prep_atm_active` gates
  on `is_init3`/`step3` instead of the ogrid flags;
- `med_dump_init_agrid` (the atm-grid dump init, first exercised by this
  dump) called `ESMF_MeshGet` on unallocated `elementIDs`/`elementCoords`
  (the ogrid path allocates first) — a latent segfault on the first
  execution; it now queries `elementCount` first and allocates.

The ogrid aoflux dump (`cmeps_aoflux_f2000.csv`) still produces nothing for
this case by design: PET 0 has no masked ogrid cells, and a cross-PE redesign
of the traced-cell selection is out of scope (the C_JRA_t62 lane already
pins the aoflux laws).

## The dump

`stage1/dumps/cmeps/cmeps_prep_atm.csv` (741708 B, 432 rows = 36 traced atm
cells x 12 coupling steps), written by `med_phases_prep_atm_dump` after the
merges, one row per traced atm-grid cell per step: `step,ymd,tod,n,gid,lat,lon`
(metadata; n is the PET-0 local atm index, gid the global element sequence
id), the fractions `lfrac,ifrac,ofrac`, then per field the three mapped
sources and the merged export: `Sl_/Si_/So_/Sx_avsdr`, `avsdf`, `anidr`,
`anidf`, `tref`, `u10`, `qref`, `Fall_/Faii_/Faox_/Faxx_taux`, `tauy`, `lat`,
`sen`, `evap`, `lwup`, `Sl_/Si_/So_/Sx_t`, and the `So_t` copy export
(duplicated column; the extractor/validator rename the second `So_t` to
`So_t_export`). 67 columns; floats `es26.17e3`, subnormals flushed.

Traced cells: the `dump_first_cells = 12` base window of PET 0 (n = 1..12,
present in every step) unioned with the per-step partial-ice (0 < ifrac < 1)
and partial-land (0 < lfrac < 1) extras, capped at 12 each; on PET 0's f19
slice the extras are also the same cells every step, so the traced set is the
fixed 36 cells n = 1..12, 18, 20, 22, 34, 36, 38, 40, 50, 52, 54, 66, 68, 70,
72, 88, 95, 102, 109, 111, 118, 120, 123, 136, 139 (all coastal Antarctica
cells of PET 0's slice; the stable (n, gid) pairing holds across all 12
steps). The pinned fraction regimes: fully-land
(lfrac = 1, n = 1), open-ocean (ofrac = 1, n = 5), partial-ice (ifrac =
0.99875, n = 2), the mixed three-fraction cell (n = 6: lfrac 0.0753, ifrac
0.6206, ofrac 0.3041), and partial-land (n = 95, 102). Sources on
zero-fraction cells carry the CMEPS spval 1e30 (the 13 land-source columns
`Sl_*`/`Fall_*` on the 20 lfrac = 0 cells: 240 spval-like values per column)
or 0.0 (the `Si_*` columns of the ice-free Antarctica rows) — both inert at
weight 0 (0.0 x 1e30 = 0.0 exactly, no NaN), and the tests pin the resulting
values implicitly.

## Validation (stage1/analysis/validate_cmeps_prep_atm_dump.py)

The validator replays the merge law per group (all 14 merged exports plus the
`So_t` copy) in the Fortran term order and checks every row:

- all 14 merge groups reproduce to **worst relative error <= 2.9e-16** over
  all 432 rows (bit-identical modulo the last ulp), with the spval-like
  sources checked as spval, not NaN (0 NaN anywhere);
- `lfrac + ifrac + ofrac == 1` (tol 1e-12) on every row;
- `So_t_export == So_t` everywhere (the copy law);
- fraction-regime counts: 252 partial-ice / 180 open-ocean / 0 fully-iced /
  168 partial-land rows; distinct ifrac values 80 (0.0000..0.9988);
- every step traces the base window (n = 1..12 present in every step, 36
  cells/step), and the stable (n, gid) pairing holds across all 12 steps.

`RESULT: prep_atm PASS`. Neutrality: the instrumented run's `cpl.r`, `cice.r`
and `docn.r` restarts are md5-identical to the reference `F2000climo_f19`
case (b9a9f0e2... / dc5e4704... / 0e5c3203...); `cam.r`/`clm2.r`/`mosart.r`
differ only in the filename variables (`cpath`, `hrestpath`, `nhfil`,
`op_field1/2`, `locfnh`, `locfnhr`) and the `caseid`/`case_title`/`history`
attributes, as expected.

## Transcription notes

- The 15 exports are `fraction_weighted_sum` composed twice (no dedicated
  three-source template): `((Sl*lfrac + Si*ifrac) + So*ofrac)` reproduces the
  `med_merge` accumulation order exactly, including the spval-at-zero-weight
  rows.
- `Sx_qref`/`Faxx_evap` are the bulk blank-suffix branch of the
  `water_bulk_or_tracers` split (the traced run is bulk).
- The `Faox_*` ocean sources come from the aoflux phase output
  (`FBMed_aoflux_a`; on this case's `aoflux_grid = ogrid` they are the ogrid
  values mapped to the atm grid by `med_aofluxes_map_ogrid2agrid_output`,
  called at the top of `med_phases_prep_atm`, F90:143-145); the `So_*` albedo
  sources come from the ocnalb phase output (`FBMed_ocnalb_a`, `addmap_ocnalb`).
- The `So_t` copy export is named `So_t_export` in the model because the
  dump carries both the `So_t` source and the copy export columns.

## Bugs found in the CESM code (this lane)

- bugs.md #38: the prep_atm dump gating bug (fork instrumentation; see "The
  instrumentation fix" above) and the latent `ESMF_MeshGet` crash in
  `med_dump_init_agrid`.

## Tests

5 tests x 15 assertions = 75 assertions, all green under `./esm test stage2/`
(whole suite **2737P/0F/0E**; the worktree's starting commit 5726a2d ran
2662P/0F/0E, so this lane adds exactly its 75) at `rel 1e-12`, `time_span
{0,0}`, plus `./esm validate` clean. The pinned rows exercise the fraction
regimes: `step1_fully_land` (merge reduces to the land source), `step1_open_ocean`
(reduces to the ocean source), `step1_partial_ice` (two-source ice + ocean),
`step1_mixed` (full three-way merge, every source valid) and `step12_mixed`
(the same cell at the last step, where the 05:30 local-time albedo sources
are the dark-side 1.0 from all three sources and the merge still sums to
1.0). Every merged export is asserted at full dump precision, so a wrong
weight, a wrong source, a wrong association or a wrong fraction all fail.
Full 432-row coverage is the validator's job (bit-identical replay).

## Gaps

- The ogrid aoflux dump for this case produces no rows (PET 0 all-land; a
  cross-PE traced-cell selection would be needed — out of scope; the aoflux
  laws are pinned by the C_JRA_t62 lane).
- No fully-iced rows in this window (ifrac max 0.9988); the lfrac = 1 rows
  pin the land-only reduction instead. The dark side is pinned only through
  the albedo = 1.0 sources of `step12_mixed`.
- The `lnd2atm` merges on grids with `samegrid_atmlnd = .false.` (the
  `lfrin` fallback of the `Sl_lfrac` custom copy, med_phases_prep_atm_mod.F90)
  are not exercised (the traced case is same-grid).
