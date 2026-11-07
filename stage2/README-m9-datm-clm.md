# M9: datm -> CLM atmospheric-forcing subassembly (review round)

`stage2/m9_datm_clm.esm` mounts `cdeps_datm_gswp3.esm:Clmncep`,
`clm_atm_forcing.esm:ImportFields` and `clm_atm_forcing.esm:Downscale` and
couples them exactly as the CESM mediator hands the datm exports to the land
in the I2000Clm60Sp f19 reference run. The 12 tests per file pin the coupled
chain on the overlapping (cell, step) rows of the stage-1 dumps at rel 1e-12
(dumps are es26.17): datm stream/export rows (step s-1, cell n=1) in, CLM
forcing rows (step s, column c=1) of `clm_bareground.csv`,
`clm_canopyflux.csv` (stage=ent) and `clm_surfrad.csv` out. The one-step
offset is the CESM data-atmosphere sequence itself (the datm runs after the
land within each coupling interval, so the land is forced by the previous
interval's atmosphere; every forc_* dump value equals its one-step-earlier
datm export bit-exact — verified by `stage1/analysis/extract_m9_test_rows.py`,
which also computes the two values that are not dump columns,
forc_wind = sqrt(u^2+v^2) and forc_solar = solad_vis+solad_nir+solai_vis+solai_nir).

## Pinned (per test, step 1..12, cell 1)

| variable | law | dump column(s) |
|---|---|---|
| Downscale.forc_t | import copy + ifelse (pass-through at needs_downscaling=0); kept once (step1) as coupling-wiring check | clm_bareground.csv forc_t |
| Downscale.forc_q | datm Sa_shum qsat-at-tbotmax clamp + import copy (m9); import copy (clm_atm_forcing) | clm_bareground.csv forc_q |
| Downscale.forc_pbot | pass-through; kept once (step1) as coupling-wiring check | clm_bareground.csv forc_pbot |
| Downscale.forc_rho | density law (pbot - 0.378 vp)/(rair t) | clm_bareground.csv forc_rho |
| ImportFields.forc_u/forc_v | datm wind split Sa_u = Sa_v = wind/sqrt(2) (m9); import copy (clm_atm_forcing) | clm_bareground.csv forc_u/forc_v |
| ImportFields.forc_wind | sqrt(u^2+v^2); kept once (step1), derived (not a dump column) | = datm strm_wind |
| Downscale.forc_rain | rain channel, identically 0 (all-snow); kept once (step1) as coupling-wiring check | clm_canopyflux.csv forc_rain |
| Downscale.forc_snow | datm all-snow precip split + partition_precip (frac_rain = 0) | clm_canopyflux.csv forc_snow |
| Downscale.forc_solad_vis/nir, ImportFields.forc_solai_vis/nir | datm direct/diffuse ratios + copies | clm_surfrad.csv forc_solad_vis/nir/solai_vis/nir |
| Downscale.forc_solar | sum of the four solar components; kept once (step1), derived (not a dump column) | = datm strm_swdn |

## Coverage gaps (implemented but not pinned here)

- **ImportFields copies**: forc_hgt (= Sa_z, incl. the datamode 30 m fallback
  when strm_z = spval), forc_topo (= Sa_topo), forc_th (= Sa_ptem),
  forc_lwrad (= Faxa_lwdn), forc_hgt_u/t/q, forc_rainc/rainl/snowc/snowl
  (only the totals are pinned; the traced Faxa_rainc/rainl are 0).
- **ImportFields derived**: forc_vp (vapor pressure), forc_po2
  (o2_molar_const factor), forc_rh (qsat_qs + percent) — the traced columns
  are all at T ~ 249 K, so forc_rh would exercise only the ice branch of the
  QSat law, and even that is unpinned.
- **Glacier topo downscaling** (needs_downscaling = 1): the whole tbot_c /
  Hbot / pbot_c / thbot_c / qbot_c / rhos_c branch of downscale_forcings
  (atm2lndMod.F90:202-267) — the traced columns are non-glacier soil columns
  at the South Pole with topo_col == forc_topo_grc (2824.921 m), so
  hsurf_diff = 0 and the branch is inactive in every test.
- **Glacier precip ramp** (landunit_glc = 1): the glc all_snow_t / frac_rain_slope
  ramp endpoints (ctsm60_params.c260518.nc) — no glacier landunit is traced.
- **Ramp interior**: frac_rain is 0 on all 12 traced steps (T ~ 249 K below
  the nonglc all-snow temperature 273.15 K); 0 < frac_rain < 1 and the
  rain/snow conversion bookkeeping (rain_to_snow_conversion,
  snow_to_rain_conversion, eflx_sh_precip_conversion = mm_to_m*denh2o*hfus)
  never trigger, so only the all-snow branch is pinned.
- **Longwave downscaling** (glcmec_downscale_longwave = 1): lwrad_down (the
  Van Tricht lapse_rate_longwave = 0.032 W/m2/m with the
  longwave_downscaling_limit clamp) — the flag is 0 in the traced run.
- **DATM branches inactive in the traced window**: tbotmax < 50 K tkfrz
  addition, pbot = 0 -> 101325, the rain branch of the precip split
  (frac > 0), the swnet albedo branch, Faxa_lwdn end-to-end (asserted
  nowhere), and the bias-correction/anomaly-forcing block (inactive for
  GSWP3v1).
- **Cells/steps**: only the South Pole row (c = 1..5, value-identical) and
  the first 12 half-hourly steps of 2000-01-01 (one geographic point, one
  day, no warm-season or glacier columns). The datm-side laws themselves are
  pinned more broadly by the cdeps_datm_gswp3.esm component tests.

Closing these gaps needs stage-1 dump extension: a traced glacier column
(with needs_downscaling = 1 and landunit_glc = 1) in a run that also visits
the ramp interior (T near 273.15 K) and glcmec_downscale_longwave = .true.,
plus dump columns for forc_vp/forc_po2/forc_rh/forc_lwrad (and the conversion
fluxes) on those rows.

## Review-round notes (2026-11-06)

- The 13 ImportFields -> Downscale couplings are restated verbatim in
  `m9_datm_clm.esm` because the CLI cannot mount a file's internal couplings
  by reference; the exact validation error is recorded in the m9 file
  description (probed under `scratch/`).
- The test tuples (including the two derived values) are re-verified against
  the dumps by `stage1/analysis/extract_m9_test_rows.py` at rel 1e-12.
