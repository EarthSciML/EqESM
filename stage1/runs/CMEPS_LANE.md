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
| `d81b5cc1` | dynamic traced-cell selection (the "wide" dump): new MED attributes `dump_file`, `dump_stable_cells`, `dump_wind_cells`; when their sum is > 0 the selection switches from fixed first-N cells to a per-coupling-step rank of ALL masked PET-0 cells — the top `dump_wind_cells` by relative wind sqrt((ubot−us)²+(vbot−vs)²) and the top `dump_stable_cells` by delt = thbot − ts with delt > 0 — union, deduped, sorted by local index; ranking reads only the raw aoflux inputs (solution-neutral, verified below); fixed first-N mode kept for backward compatibility (`namelist_definition_drv.xml` entries + `med_phases_aofluxes_mod.F90` call site) |
| `10741c34` | dump header: `iter` column position fixed to right after `gid` (data write order); subnormal flush to 0.0 |
| `4750efc5` | `es26.17e3` float format (gfortran `es26.17` omits the exponent letter, bugs.md #29); `rainc` fallback init `aoflux_in%rainc(1) = 0.0` (bugs.md #28) |
| `bd88a843` | dump passes a sanitized `rainc` (= 0 when `add_gusts = .false.`) instead of the out-of-bounds `rainc(n)` read of the mediator's 1-element dummy array (bugs.md #28) |

Pointer commits in `EqESM-cesm-cmeps`: `bb957c8` (wide run), `5d367e0`, `0403cf3`, `735eb8b`.

## The dump

Two dump files in `stage1/dumps/cmeps/`, both gated by the mediator attribute
`dump_cmeps` (default `.false.`), PET 0 only, 67 columns each: the 12 raw inputs,
the 11 derived inputs, the 4 neutral first estimates, the 16-field exit state of
the stability iteration (`iter` included), the 7 fluxes and the 10 diagnostics,
plus point metadata (`n`, `gid`, `lat`, `lon`). `# description` / `# columns`
headers; floats written `es26.17e3` with subnormals flushed (both header and
format are gfortran-workaround decisions documented in the header).
`hol`/`xsq`/`xqq` are the zbot-level pre-2m-rescale values (`hol0`/`xsq0`/`xqq0`
of the instrumentation), `hol2`/`xsq2`/`xqq2` the 2-m diagnostics. Note: the
column name `lat` is used for both the cell latitude (first occurrence) and the
latent-heat flux (second); the replay/validators resolve the collision by name
(the latitude is metadata, the latent-heat flux is asserted via evap = lat/LATVAP).

- `cmeps_aoflux.csv` — fixed selection (first `dump_first_cells` masked cells):
  12 masked ocean cells × 24 coupling steps = 288 rows, local ogrid indices
  n 6933..7117, gid 41583..42667, all ~69 S on the Antarctic coast (PET 0's slice
  of the tx0.66v1 ocean grid), lon 74.0–77.3 E (ESMF [−360, 0] convention).
- `cmeps_aoflux_wide.csv` — dynamic selection (`dump_stable_cells=8`,
  `dump_wind_cells=8`, `dump_first_cells=0`): per coupling step the strongest-
  wind 8 of ALL masked PET-0 cells plus the most-stable 8 (delt > 0), union,
  deduped, sorted by n. 192 rows = 8 rows/step × 24 steps, all wind-selected
  (the stable selection found no cell with delt > 0 on PET 0 in this window).
  9 distinct cells, a subset of the fixed window's 12 (the 3 weakest-wind cells
  n 7023/7113/7114 drop out of every step's top-8).

## Validation

- **Neutrality.** Both instrumented runs' restarts are md5-identical to the
  reference run's: `cpl.r 0f47e2f0f4f207240148d8ef333d5044`, `datm.r
  e14609a64cbba942d56a9eaa961c4bda`, `dice.r 8e1656e609b0c189b68284b940b2a19b`,
  `drof.r 002ca0b2a42c135246a5134e2517e00e`, `mom6.r
  2f19de039ced94bcd30930385cc29b76`; `ocean.stats.nc` byte-identical. The
  dynamic per-step ranking of the wide run reads the aoflux inputs read-only,
  so it is provably solution-neutral too (the run job 12055973, resubmitted
  after the SLURM_MEM_* env conflict of 12047777; the fork code is d81b5cc1).
- **Full replay.** `stage1/analysis/validate_cmeps_dump.py` replays
  `flux_atmOcn_large` from the dumped inputs and reproduces all 46 computed
  columns on all rows to machine precision and the dumped `iter` on every row:
  ≤ 4.4e-15 relative on the 288 fixed-window rows, ≤ 3.4e-15 on the 192
  wide-window rows (worst columns fac_t/fac_q). The replay mirrors the stage-2
  esm structure: all 5 updates unrolled, `conv_k = |(ustar_k −
  ustar_{k−1})/ustar_k| ≤ 0.01` (k = 1..4), `iter` = first converged k else 5,
  exit state selected by `iter`.
- **Structure.** Both dumps: 2 header lines, 67 columns, zero wrong-arity rows,
  zero NaN, zero spval-like values (`rainc = 0` on every row after the #28 fix);
  the wide dump passes the per-step row-count check (8 rows/step, cell set
  varying between steps as expected of the dynamic selection).
- **Coverage — fixed window.** All 288 rows are unstable Antarctic-coast cells
  with low-to-moderate wind: `hol ∈ [−10.0, −1.16]` (21 rows pinned at the
  |hol| ≤ 10 clamp), `vmag ∈ [0.5, 4.8]` m/s (17 rows pinned at the minwind
  floor), ocean currents small but nonzero (`us`/`vs` ∈ ±0.13 m/s, so the
  stress direction cosines bite), iteration counts 2 (145 rows), 3 (85), 4 (37),
  5 (21).
- **Coverage — wide window (and its limits).** The dynamic selection ranks ALL
  masked PET-0 cells every step, so its coverage is exhaustive over the domain:
  the strongest-wind rows the run offers — `vmag ∈ [0.63, 4.81]` m/s (domain
  max 4.81, still far below the cdn 33 m/s cap) — at iteration counts 2 (126
  rows), 3 (57), 4 (9); `hol ∈ [−9.90, −1.16]` (no clamp rows). NOT exercised
  even by the exhaustive ranking: the stable branch (delt > 0 nowhere on PET
  0's 12 masked cells in this January window, so the stable selection is empty
  every step) and the cdn 33 m/s wind cap (no cell of PET 0 exceeds 4.8 m/s).
  The stable-branch and cdn-cap laws are transcribed in the templates but carry
  no dump pins; closing those gaps needs a different window (a PET slice or
  season with warm-air-over-cold-water cells and/or storm-force winds), not
  more cells of this domain. Gustiness and the cold-air-outbreak branch are
  namelist-off in this case and were not turned on (documented, not run).

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
- 9 tests × 52 assertions = 468 assertions (6 from the fixed dump, 3 from the
  wide dump), all green under `./esm test stage2/` (1327P/0F/0E whole suite) at
  `rel 1e-12`, `time_span {0,0}`, pinning dump rows at iterations 2/3/4/5/5/2
  (fixed window, including the minwind-floor row vmag = 0.5 and the hol-clamp
  row hol = −10) and 2/3/4 (wide window: `wide_maxwind_cell` (4, 42667) the
  domain-max vmag 4.81 m/s at iter 2, `wide_iter3_coldest_cell` (13, 42126) at
  iter 3, `wide_iter4_holclamp_adj` (19, 41583) at hol −9.901, a tenth off the
  clamp). The wide tests pin the wind-ranked selection's rows and iteration
  variants; the stable branch and the cdn cap remain gaps (see above) and are
  documented, not pinned.

## Bugs found in the CESM code

- bugs.md #28: `med_phases_aofluxes_mod.F90` `set_aoflux_in_pointers` allocates
  `aoflux_in%rainc(1)` uninitialized when `add_gusts = .false.` (benign in
  pristine code — rainc is never read on the no-gusts path — but the dump
  exposed it as an out-of-bounds read).
- bugs.md #29: gfortran `es26.17` omits the exponent letter for subnormals and
  3-digit exponents (unparseable CSV); `es26.17e3` fixed.

## Status

Stage 1 complete: fixed dump (288 rows) plus wide dump (192 rows, dynamic
per-step stable/wind selection), both validated to machine precision and
solution-neutral (restarts md5-identical; wide run 12055973). The wide window
provably exhausts PET 0's domain — the strongest-wind rows available (max
4.81 m/s) are pinned; the stable branch and the cdn 33 m/s cap are not reachable
in this window and stay documented gaps (closing them needs a different
window/domain, and gustiness/cold-air-outbreak would need namelist changes).
Stage 2 complete (model + templates + tests green, 468 cmeps assertions, whole
suite 1327P/0F/0E). Next: import `cmeps_aoflux` by reference into the top-level
`eqesm.esm` coupling (stage 3).

## 2026-11-05: ocean albedo and ocean-export merge dumps (ocnalb/merge lane)

Second instrumentation round on the same fork (branch `instrument/cmeps-dumps`,
commits `cb613113` in the cmeps submodule, pointer `43ccf07` in
EqESM-cesm-cmeps): two more phase traces of the same C_JRA_t62_inst run, gated
by the new mediator attributes `dump_ocnalb`/`dump_merge` (default .false.).

- `cmeps_ocnalb.csv` (288 rows = 12 cells × 24 steps): the
  `med_phases_ocnalb_run` albedo phase — orbital inputs/outputs
  (nextsw_cday, eccen/obliqr/lambm0/mvelpp, delta/eccf), the albedo settings
  (albdif/albdir/min_albedo/use_min_albedo/use_briegleb/flux_albav), rlat/rlon,
  cosz, the FBFrac fractions (ifrac/ofrac/ifrad/ofrad) and the four albedos
  (anidr/avsdr/anidf/avsdf), 31 columns.
- `cmeps_merge.csv` (288 rows): the `med_phases_prep_ocn_accum` ocn-export
  merges plus the `med_phases_prep_ocn_custom` shortwave computation — merge
  inputs (atm/aoflux/ice/runoff), the 4 ocean albedos (`FBMed_ocnalb_o` So_*,
  the link to the albedo phase), the derived columns (fswabsv/fswabsi, the
  scaled fractions), the 13 merged outputs and the 4 shortwave bands, 60
  columns.

Trace selection: the same fixed `dump_first_cells = 12` window (n 6933..7117)
unioned with the per-step partial-ice cells; all 12 window cells are
partial-ice (ifrac 0.0683..0.6246), so the traced set is the fixed 12 every
step.

Runs (all solution-neutral: restarts md5-identical to the reference — cpl.r
0f47e2f0…, mom6.r 2f19de03… — and ocean.stats.nc byte-identical):

| Job | Build | Result |
|---|---|---|
| 12119693 | first ocnalb/merge attempt | hung, killed (scancel), superseded |
| 12126699 | writers v1 | merge dump complete; ocnalb dump headers-only (shared selection cursor bug) |
| 12127523 | cursor fix (`inext2o`) | both dumps 288 rows, trailing-comma artifact |
| 12134714 | + trailing-comma fix | **definitive** dumps; also regenerated cmeps_aoflux_wide.csv (md5 166a90fc…, stats identical to the recorded 12055973 stats) |

Validation (`stage1/analysis/validate_cmeps_ocnalb_merge_dumps.py`, replay in
Fortran term order): **all 7 ocnalb columns and all 23 merge columns reproduce
bit-identically** (0.000e+00 relative error over all 288 rows of each dump),
including the Berger 1978 declination/cosz/Briegleb chain and the swnet chain
with the scaled fractions and the band split. Two findings pinned by the
dumps: (1) the custom shortwave band coefficients `c1..c4 = 0.285/0.215`
(med_phases_prep_ocn_mod.F90:830-834) are single-precision literals promoted
to r8 — every row has vdr/swnet = 0.28499999642372131, not 0.285; (2) the
radiation fractions lag the current fractions by one coupling step (the
ocnalb phase copies ifrac/ofrac into ifrad/ofrad after its albedo loop:
ifrad(k) = ifrac(k−1), verified on all 288 rows).

Bugs: bugs.md #34 (ocnalb init else branch sets the dead local `use_min_ocnalb`
instead of the module `use_min_albedo`, benign in the traced run); two
instrumentation bugs fixed during the lane (shared `inext2` cursor between the
new writers starved the ocnalb dump; trailing empty CSV field from the
`55(es26.17e3,",")` format — fixed in 12134714).

Stage-2 handoff (this round): 12 new expression templates in
`stage2/cmeps_templates.esm` (earth_true_longitude, solar_declination,
earth_sun_distance_factor, solar_zenith_cosine, briegleb_direct_albedo,
taylor_direct_albedo, ocean_albedo_scheme, ocean_albedo_visible,
ocean_albedo_nearinfrared, scaled_fraction, fraction_weighted_sum,
shortwave_band); `stage2/cmeps_ocnalb.esm` (model `OcnAlb`, 7 equations, 9
tests × 7 assertions) and `stage2/cmeps_merge.esm` (model `OcnMerge`, 23
equations, 9 tests × 23 assertions), both importing the shared templates by
reference; `stage2/README-cmeps-ocnalb.md`. Whole stage2 suite
**1685P/0F/0E** under `./esm test stage2/` at `rel 1e-12`, `time_span {0,0}`.
Gaps (documented, unpinned): Taylor scheme, flux_albav branch, use_min_albedo
clamp and dark side (all 288 ocnalb rows are Briegleb day-side with
flux_albav = use_min_albedo = 0); the import_swpen_by_bands per-band branch
(SICE exports no Fioi_swpen_* bands); the enthalpy merges (not dumped); the
prep_atm atmosphere merges (needs a CAM case).

## 2026-11-09: prep_atm atmosphere-export merge dump (F2000climo_f19_inst case)

Third instrumentation round on the same fork (branch `instrument/cmeps-dumps`):
the `med_phases_prep_atm` lnd/ice/ocn -> atm merges, gated by the new mediator
attribute `dump_prep_atm` (default .false.). This closes the last
prep_atm gap of the ocnalb/merge round ("needs a CAM case").

- **Case:** `stage1/cime/cases/F2000climo_f19_inst` (cloned from the reference
  `F2000climo_f19`; deltas: `user_nl_cpl` dump gating only), CAM6+CLM+CICE+
  DOCN+MOSART, f19 grid, 6 h run, 12 coupling steps at 1800 s, 36 PEs. Built
  and run in the agent allocation with `./case.build` + `python3 .case.run`
  (run 22:28-22:35 local, logs `*.log.12284979.261109-222841.gz`; the
  `aoflux_grid = ogrid` setting is the case default).
- **Dump:** `cmeps_prep_atm.csv` (741708 B, 432 rows = 36 traced atm cells x
  12 steps, fixed cell set every step: n = 1..12, 18, 20, 22, 34, 36, 38, 40,
  50, 52, 54, 66, 68, 70, 72, 88, 95, 102, 109, 111, 118, 120, 123, 136, 139
  — the `dump_first_cells = 12` base window plus the per-step partial-ice /
  partial-land extras, capped at 12 each). Columns: step,ymd,tod,n,gid,lat,lon,
  lfrac,ifrac,ofrac, then per field the three mapped sources and the merged
  export (Sl_/Si_/So_/Sx_ for the 4 albedos, tref/u10/qref, the
  Fall_/Faii_/Faox_/Faxx_ 6 fluxes) and Sl_/Si_/So_/Sx_t, plus the duplicated
  So_t copy-export column; 67 columns, `es26.17e3`, subnormals flushed.
- **Root cause found (fork instrumentation bug, bugs.md #38):**
  `med_dump_prep_atm_active` was gated on the *ogrid* aoflux init flag
  (`is_init`) and counter (`step`). On this case PET 0's f19 ogrid block is
  all-land (Antarctica rows first in local order; `So_omask` = 0 everywhere on
  PET 0), the ogrid init bails ("WARNING: no masked cells traced" in ocn.log),
  `is_init` stays .false., and the prep_atm dump never fired — silently (no
  rows, no error). Fixed (cmeps submodule `17eef08f`, fork pointer `ea48423`):
  the atm-export dump gets its own counters `step3`/`ymd3`/`tod3`, advanced by
  `med_dump_prep_atm_tick` from the mediator clock inside
  `med_phases_prep_atm_dump`; `med_dump_prep_atm_active` gates on
  `is_init3`/`step3`. A second latent crash fixed in the same pass:
  `med_dump_init_agrid` (first exercised by this dump) called `ESMF_MeshGet`
  on unallocated `elementIDs`/`elementCoords` — it now queries `elementCount`
  first, then allocates.
- **Validation** (`stage1/analysis/validate_cmeps_prep_atm_dump.py`):
  **all 14 merge groups reproduce to worst relative error <= 2.9e-16 over
  all 432 rows** (0 NaN anywhere; the spval 1e30 land sources on the 20
  lfrac = 0 cells — 240 spval-like values per Sl_*/Fall_* column — are
  checked as spval, not NaN, and are inert at weight 0); `lfrac+ifrac+ofrac
  == 1` (tol 1e-12) everywhere; `So_t_export == So_t` everywhere (the copy
  law of esmFldsExchange L1282-1285); regime counts 252 partial-ice / 180
  open-ocean / 168 partial-land rows; the base window n = 1..12 present in
  every step and the (n, gid) pairing stable across steps. `RESULT: prep_atm
  PASS`.
- **Neutrality:** `cpl.r` / `cice.r` / `docn.r` md5-identical to the
  reference F2000climo_f19 run (b9a9f0e229d2fd94b97244f2f9879690 /
  dc5e470490799c48bda968fb62b6d850 / 0e5c32030dc324addbe661e3273df90d);
  `cam.r`/`clm2.r`/`mosart.r` differ only in the filename variables (cpath,
  hrestpath, nhfil, op_field1/2, locfnh, locfnhr) and the caseid/case_title/
  history attributes — the dump code changes nothing when enabled.
- **Existing dump md5s preserved:** cmeps_aoflux.csv 6c2428e5…, cmeps_
  aoflux_wide.csv 166a90fc…, cmeps_merge.csv 7ca39ed3…, cmeps_ocnalb.csv
  cef4d00d… (untouched by this case).
- **Ogrid aoflux dump skipped by design:** `cmeps_aoflux_f2000.csv` produces
  nothing for this case (PET 0 has no masked ogrid cells); a cross-PE
  redesign of the traced-cell selection is out of scope (the aoflux laws are
  pinned by the C_JRA_t62 lane).
- **Stage-2 handoff (this round):** `stage2/cmeps_prep_atm.esm` (model
  `AtmMerge`, 15 equations — the 14 merged exports as `fraction_weighted_sum`
  composed twice in the Fortran accumulation order plus the `So_t` copy
  export — importing the shared template by reference; 5 tests x 15
  assertions pinning the fraction regimes fully-land / open-ocean /
  partial-ice / mixed three-fraction, plus the last-step mixed cell with the
  dark-side albedo 1.0 sources) and `stage2/README-cmeps-prep-atm.md`.
  Whole stage2 suite **2737P/0F/0E** under `./esm test stage2/` (ESD_ROOT set;
  the worktree's starting commit 5726a2d ran 2662P/0F/0E — this lane adds
  exactly its 75 assertions) at `rel 1e-12`, `time_span {0,0}`.
