# CMEPS ocean albedo and ocean-export merge components (stage 2)

This lane extends the CMEPS stage-2 work (`README-cmeps.md`) with two more
pointwise components transcribed from the instrumented C_JRA_t62 run:

| File | Model | Laws | Source |
|---|---|---|---|
| `cmeps_ocnalb.esm` | `OcnAlb` | solar geometry + ocean albedo | `mediator/med_phases_ocnalb_mod.F90` (with `share/src/shr_orb_mod.F90`) |
| `cmeps_merge.esm` | `OcnMerge` | fraction-weighted ocn-export merges + custom shortwave | `mediator/med_phases_prep_ocn_mod.F90` (with `med_merge_mod.F90`, `esmFldsExchange_cesm_mod.F90`) |

Both import their shared math from `./cmeps_templates.esm` by reference (the
library now carries, in addition to the aoflux templates, the Berger 1978
solar-geometry templates `earth_mean_anomaly` / `earth_true_longitude_series` /
`earth_true_longitude` / `solar_declination` / `earth_sun_distance_factor` /
`solar_zenith_cosine`, the albedo templates `briegleb_direct_albedo` /
`taylor_direct_albedo` / `ocean_albedo_scheme` / `ocean_albedo_visible` /
`ocean_albedo_nearinfrared`, and the merge templates `scaled_fraction` /
`fraction_weighted_sum` / `open_ocean_weighted_sum` /
`ocean_absorbed_shortwave` / `shortwave_band`). The mean anomaly
`lambm0 + (calday−80.5)·2π/365 − mvelpp` is factored into `earth_mean_anomaly`
and bound once by `earth_true_longitude` (which passes it to the
equation-of-center series `earth_true_longitude_series`); `fswabsv` and
`fswabsi` share the single `ocean_absorbed_shortwave` law
`sw_dr·(1−α_dr) + sw_df·(1−α_df)` (the same law as the CLM `absorbed_shortwave`
of clm_templates.esm, named distinctly because both template libraries are
imported together and the esm-spec deep-equality rule rejects same-name
non-deep-equal definitions); the two-source merges whose sources share the
open-ocean weight (Faxa_rain/Faxa_snow/Foxx_lwnet) use
`open_ocean_weighted_sum` = `(a + b)·ofrac` rather than `fraction_weighted_sum`
(which is kept for the different-weight Foxx_taux/Foxx_tauy merges).

## The dumps

`stage1/dumps/cmeps/cmeps_ocnalb.csv` (288 rows = 12 traced cells × 24 steps)
pins the albedo phase; `stage1/dumps/cmeps/cmeps_merge.csv` (same 12 cells × 24
steps) pins the merges and the custom shortwave computation. Both come from the
definitive instrumented run **12134714** (fork commits
`cb613113`/`43ccf07` on `instrument/cmeps-dumps`), which is solution-neutral:
its restarts are md5-identical to the reference C_JRA_t62 run and to the two
earlier instrumented runs (12126699, 12127523 — three independent runs, same
restart md5s), and `ocean.stats.nc` is byte-identical. The traced cells are the
fixed `dump_first_cells = 12` window of PET 0 (local ogrid indices n
6933..7117, all ~69 S, January); all 12 are partial-ice cells
(ifrac 0.0683..0.6246), so the traced set is the same 12 cells every step.

Row contents:

- **ocnalb** (31 columns): nextsw_cday (1.0417..2.0 day), the orbital inputs
  eccen/obliqr/lambm0/mvelpp and outputs delta/eccf, the albedo settings
  albdif/albdir/min_albedo/use_min_albedo/use_briegleb/flux_albav, rlat/rlon,
  cosz, the FBFrac fractions ifrac/ofrac/ifrad/ofrad, and the four albedos
  anidr/avsdr/anidf/avsdf.
- **merge** (60 columns): the 4 fractions, the 9 atm inputs (rain/snow
  convective+large-scale, lwdn, 4 shortwave bands), the 5 aoflux inputs
  (lwup/sen/taux/tauy/evap), the 6 ice imports (taux/tauy/melth/meltw/salt/
  swpen), the 2 runoff inputs, the 4 ocean albedos (So_*, the `FBMed_ocnalb_o`
  link to the ocnalb phase), the 6 derived columns (fswabsv/fswabsi, the four
  scaled fractions), the 13 merged outputs and the 4 shortwave bands.

## Validation (stage1/analysis/validate_cmeps_ocnalb_merge_dumps.py)

The validator replays both phases from the dumped inputs in the Fortran term
order and compares every dumped output column:

- **ocnalb**: all 7 columns (delta, eccf, cosz, anidr, avsdr, anidf, avsdf)
  reproduce **bit-identically** (0.000e+00 relative error over all 288 rows).
  The Berger 1978 declination series (ve = 80.5, dayspy = 365), the cosz chain
  and the Briegleb et al. (1986) law with its cubic correction term
  `0.026/(cosz^1.7 + 0.065) + 0.15·(cosz−0.1)(cosz−0.5)(cosz−1.0)` are exact.
- **merge**: all 23 columns reproduce **bit-identically** (0.000e+00), including
  the swnet chain with the scaled fractions and the band split.

Two findings surfaced during validation, both pinned by the dumps:

1. **The band coefficients are single-precision literals.** The Fortran
   `c1 = 0.285` (med_phases_prep_ocn_mod.F90:830-834) carries no `_r8` suffix,
   so it is a float32 literal promoted to r8: every dumped row has
   `Foxx_swnet_vdr/Foxx_swnet = 0.28499999642372131` (= float32 0.285), not
   0.285. The `shortwave_band` template documents this and the components bind
   the exact promoted values (0.2849999964237213 / 0.2150000035762787).
2. **The radiation fractions lag by one coupling step.** The ocnalb phase
   copies ifrac/ofrac into ifrad/ofrad *after* its albedo loop, so at the
   merge phase's dump point (which runs after ocnalb in the next step's prep)
   the ifrad/ofrad used by the swnet computation are the *previous* step's
   fractions: ifrad(k) = ifrac(k−1), verified for all 288 rows. The merge
   component treats ifrad/ofrad as inputs (that is what the Fortran does); the
   lag is documented in both components.

## Transcription notes

- The ocnalb phase skips the whole computation when `flux_albav`; the traced
  run has `flux_albav = .false.`, so that branch is transcribed but unpinned,
  as are the Taylor scheme (`ocean_albedo_scheme = 0` in the run), the
  `use_min_albedo` clamp (`min_albedo = 0` unset) and the dark-side branch
  (all 288 rows are day-side, cosz 0.036..0.695). The `nextsw_cday >= -0.5`
  guard is always true here.
- `Foxx_swnet` weights the *open-ocean absorbed* shortwave by the lagged
  radiation fraction `ofracr_scaled` and the ice-penetrating shortwave by the
  current `ifrac_scaled` — the asymmetry is the Fortran's (F90:814), not a
  transcription artifact.
- The `Fioi_*_ice` dump labels map to `FBImp(compice,compocn)` fields whose
  names drop the `_ice` suffix; the components use the dump labels.
- The `import_swpen_by_bands` per-band branch (F90:824-827) is not exercised
  (SICE exports no `Fioi_swpen_vdr/vdf/idr/idf`) and is documented, not
  transcribed as a separate template.
- The enthalpy merges (Foxx_hlat etc.) are not dumped (documented gap); the
  Sx_*/Faxx_* atmosphere merges of prep_atm need a CAM case (no prep_atm phase
  in this run sequence).

## Bugs found in the CESM code (this lane)

- bugs.md #34: `med_phases_ocnalb_mod.F90` init: the `ocean_albedo_limit` else
  branch sets the dead local `use_min_ocnalb` instead of the module
  `use_min_albedo` (benign — the namelist default leaves the limit unset in
  the traced run, and the albedo loop reads the module variable).
- Instrumentation bugs fixed during this lane (not CESM bugs): the ocnalb
  writer initially shared the merge writer's selection cursor and wrote
  headers only (fixed with a dedicated `inext2o`), and both new writers
  emitted a trailing empty CSV field (fixed in run 12134714).

## Tests

9 tests × 4 assertions (ocnalb) and 9 tests × 15 assertions (merge), 18 tests
total, all green under `./esm test stage2/` (whole suite 2297P/0F/0E) at
`rel 1e-12`, `time_span {0,0}`. The rows are the four (step, cell) corners
plus step 2 cell 6933 (fraction-lag context), mid-run step 12 cell 7024, and
the min/max/median cosz (ocnalb) and ifrac (merge) rows; the informative
outputs are asserted per test. Pins that cannot fail are dropped: the ocnalb
anidf/avsdf (copies of the overridden albdif) and avsdr (duplicate of anidr),
the merge scaled fractions (equal to their inputs in every row), the
Foxx_rofl/Foxx_rofi copies of Forr_rofl/Forr_rofi, and the identically-zero
Faxa_rain and Fioi_salt (see Gaps). Full 288-row coverage is the validator's
job (bit-identical replay).

## Gaps

- Taylor scheme, flux_albav branch, use_min_albedo clamp, dark side (ocnalb);
  import_swpen_by_bands branch, enthalpy merges, prep_atm atmosphere merges
  (merge). All documented in the component descriptions.
- Faxa_rain and Fioi_salt (merge): both are identically zero in all 288 traced
  rows — the traced DATM forcing has no rain (Faxa_rainc/Faxa_rainl = 0) and
  the SICE data-ice exports no salt (Fioi_salt_ice = 0) — so their merge laws
  `(Faxa_rainc + Faxa_rainl)·ofrac` and `Fioi_salt_ice·ifrac` are transcribed
  but unpinned; a rain-bearing or salt-exporting case would pin them.
