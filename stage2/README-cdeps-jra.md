# CDEPS DATM JRA55 stage 2 — `.esm` implementation

Stage-2 implementation of the CDEPS DATM data-atmosphere component for the
JRA55 v1.5 daily TL319 forcing with the JRA datamode (cdeps1.0.112, fork
`stage1/EqESM-cesm-cdeps` branch `instrument/cdeps-dumps`) as esm 1.3.0
files: the stream time-interpolation laws (the same `shr_strdata_advance`
layer as the GSWP3 slice, factored into the shared `stream_interp`
template), the JRA55 forcing-derivation laws, the eight JRA55 data sources
read through the esm netcdf reader, and the tests pinning the stage-1
instrumentation dumps (`stage1/dumps/cdeps/datm_jra_streams.csv`,
`datm_jra_exports.csv`, `cdeps_strdata_streams_jra.csv`) — see
[Testing](#testing). The shared pointwise laws live in
`cdeps_templates.esm` (this lane appends `stream_interp`, `density_jra`,
`rain_frac_jra`, `jra_albedo`; the GSWP3-lane laws are untouched).

Sources: `components/cdeps/datm/datm_datamode_jra_mod.F90` (the JRA55
datamode), `components/cdeps/streams/dshr_strdata_mod.F90` +
`dshr_tinterp_mod.F90` (the stream layer, identical to the GSWP3 slice),
`components/cdeps/share/shr_const_mod.F90` (the constants), and the
stream configuration of the reference C_JRA_t62 run
(`stage1/cime/cases/C_JRA_t62/Buildconf/datmconf/datm.streams.xml`);
boundary/derivation tables in `stage1/runs/CDEPS_LANE.md` (JRA55 slice
section) + `stage1/runs/cdeps_instrument.json` (`jra_run`).

## Files

| file | content |
|---|---|
| `cdeps_templates.esm` | expression-template library of the shared pointwise laws; the JRA lane appends `stream_interp` (value = flb·lb + fub·ub, 0 at night), `density_jra` (pbot/(rdair·tbot·(1 + qcoef·shum)) with the r4-rounded qcoef, see bugs.md), `rain_frac_jra` (the tbot < TKFRZ step), `jra_albedo` (the LY08 latitude law with the r4-rounded coefficients, see bugs.md) |
| `cdeps_datm_jra.esm` | the component: `data_sources` (the 8 JRA55 v1.5 daily TL319 files), the `Streams` interpolation model, the `Jra` datamode transformation model, and the `JRAFiles` netcdf-reader model; tests from the dump rows |

## Model structure

```
JRA55 v1.5 daily TL319 netCDF (2920 records of 1958; 8 vars: prec, lwdn,
   swdn, q_10, slp, t_10, u_10, v_10; netcdf3 classic float32)
   │  esm netcdf reader (metadata.esio_format = "netcdf", CONST sources,
   │  time axis at "all", lat/lon fixed at cell (0,0), extent -> N_REC)
   ▼
JRAFiles.*  (raw record values, shaped over index_sets.record)
   │  [the bilinear spatial map to the T62 model mesh is out of scope for
   │   this slice — the trace pins the post-map values (see Risks)]
   ▼
Streams.*   value = flb·lb + fub·ub   (coszen: lb·flb by day, 0 at night)
   │  per stream: prec/lwdn/q_10/slp/t_10/u_10/v_10 linear, swdn coszen
   ▼
Jra.*       the datamode laws (below)
   ▼
Sa_*, Faxa_*  export state consumed by CAM/MOM6 via the mediator
```

`Streams` carries the same time-interpolation law as the GSWP3 slice
(`shr_strdata_advance`, dshr_strdata_mod.F90:1292-1391) — factored here
into the shared `stream_interp` template instead of per-model equations.
`Jra` carries the JRA55 datamode derivations of
`datm_datamode_jra_advance` (datm/datm_datamode_jra_mod.F90:271-322):

- copies: Sa_pslv = Sa_pbot = strm_pslv; Sa_tbot = Sa_ptem = strm_tbot;
  Sa_u/Sa_v = strm_u/strm_v; Sa_shum = strm_shum (**no qsat cap** — the
  JRA datamode has no humidity law); Faxa_swdn/lwdn = strm_swdn/lwdn;
  Sa_z = 10 m; Sa_u10m/v10m = Sa_u/Sa_v.
- Sa_dens = pbot/(rdair·tbot·(1 + qcoef·shum)), rdair =
  SHR_CONST_RDAIR = 287.0423113650487 J/kg/K, qcoef = **r4(0.608)** =
  0.6079999804496765 (the Fortran literal is unsuffixed → single
  precision; see bugs.md).
- precip split: when the streams separate rain/snow (prrn/prsn present),
  Faxa_rainl/snowl = prrn/prsn passthrough; otherwise a step at
  tbot < TKFRZ: rainl = prec, snowl = 0 above freezing, the reverse
  below. Faxa_rainc = Faxa_snowc = 0. The C_JRA_t62 case has no separate
  streams (separates_rain_snow = 0, prrn = prsn = spval).
- shortwave: Faxa_swvdr/ndr/vdf/ndf = swdn·0.28/0.31/0.24/0.17 (exact
  `_R8` literals); avg_alb = **r4(0.069)** − **r4(0.011)**·cos(2·lat·degtorad)
  (LY08, Large & Yeager 2008; the coefficients are unsuffixed → single
  precision, see bugs.md); Faxa_swnet = swdn·(1 − avg_alb). lat is the
  model-mesh element latitude yc — an ESMF 8.9.1 probe confirms
  `ownedElemCoords` equals the mesh file's `centerCoords` (the dumped
  `lat` column), so no mesh-coordinate correction is needed.

`JRAFiles` pins the netcdf reader path into the eight JRA55 v1.5 files:
2920 daily records of 1958 (first record at 1.5 h after 1958-01-01 00Z),
lat ascending from −89.57°, netcdf3 classic float32; the C_JRA_t62 case
aligns model year 1 to data year 1958 (taxmode = cycle, year_first =
1958, year_last = 2023, year_align = 1) — the reader tests use record 1
(= 1958-01-01) and record 2920 (= 1958-12-31).

## The data_sources recipe (esm 0.7.0 CLI, esio enabled)

Identical recipe to the GSWP3 lane: **CONST** sources (no `temporal`
block — DISCRETE sources fail the static checks in this CLI),
`metadata.esio_format = "netcdf"`, a `file://` `url_template` (the
`JRA.v1.5.<var>.TL319.1958.210504.nc` files in
`stage1/inputdata/ocn/jra55/v1.5_noleap/`), `select.axes` = `["all",
{fixed: 0}, {fixed: 0}]` (0-based indices; the time axis MUST stay "all"),
and `extent.metaparameter` binding the delivered record count (2920) to
the `N_REC` metaparameter that sizes the `index_sets.record` interval.
Each source needs its OWN extent metaparameter when counts differ; all
eight JRA sources agree on 2920, so they share `N_REC`. The netcdf
reader delivers the file's native float32 values widened to float64; the
expected values in the tests are the exact widenings. The reader tests
(records 1 and 2920 at cell lat0/lon0) take noticeably longer than the
GSWP3 lane's because 2920-record sources are read per test — if the
suite runtime ever becomes a problem, the JRAFiles tests can be dropped
without losing law coverage (the interpolation laws are pinned by the
Streams tests).

## Testing

`./esm validate stage2/cdeps_templates.esm stage2/cdeps_datm_jra.esm`
(clean) and `./esm test stage2/cdeps_datm_jra.esm` (18 tests, 122
assertions, all passing; `./esm units` reports SI base units for the
equations). `${ESD_ROOT}`-dependent files are not used by this lane.

- **Streams tests** (6): the four synthetic branches (linear, coszen day,
  night zero, nearest) plus two traced rows of
  `cdeps_strdata_streams_jra.csv` — the step 0 coszen swdn row (flb IS
  cosz/tavCosz = 0.996590623955609711, fub = 0) and the step 11 linear
  t_10 row (flb = 1/3, fub = 2/3). Expected values are the dumped `value`
  column to 1e-12 relative.
- **Jra tests** (10): seven synthetic branch tuples (all-snow step,
  all-rain step, separated-streams passthrough, density, polar albedo +
  sw split, equatorial albedo, copies/constants) plus three traced
  (step, point) rows exercising the full stream→export chain: step 0
  point 1 (all-snow branch, 23 export assertions), step 11 point 1
  (zero-precip row — both branches of the tbot step are exercised by the
  traced rows), step 0 point 16 (per-point variation, density/swnet at
  the 8th significant digit). Inputs are the dumped `strm_*` values of
  `datm_jra_streams.csv`, expected values the dumped exports of
  `datm_jra_exports.csv`.
- **JRAFiles tests** (2): raw record 1 and record 2920 of all eight
  variables through the netcdf reader at cell (0,0).

The dumps were verified solution-neutral (bit-identical restart files vs
the reference C_JRA_t62 run, all five components) and law-exact on every
row (`stage1/analysis/check_datm_jra_laws.py`: 4608 streams rows + 192
datamode rows, 0 failures).

## Deferred / notes

- The **spatial map** (bilinear ESMF regrid from the TL319 source mesh to
  the T62 model mesh) is applied before the traced layer; the trace pins
  post-map values, so the stage-2 model consumes the same values. The
  regrid itself is an EarthSciDiscretizations concern (PR later).
- The **temporal cadence**: the dumps pin the interpolated values and
  factors per step; the stream clock semantics (why the step-11 row sits
  at flb = 1/3) follow from `shr_strdata_advance`'s interval bookkeeping
  and are not re-expressed here — a stage-3 integrator concern
  (EarthSciIntegrator).
- The JRA datamode's **dead code** (`cosfactor`/`phs_c0`/`dLWarc` in
  `datm_datamode_jra_advance`) is documented in bugs.md and omitted from
  the Jra model.
- The **r4-rounded literals** (0.608, 0.069, 0.011) are a genuine CDEPS
  finding — the laws differ from a naive transcription at the 2e-8
  (albedo) / 2e-11 (density) level. The templates and tests carry the
  exact r4-rounded constants with provenance notes (bugs.md), and the
  check script reproduces them.
- The traced 16 points are all on the southernmost T62 row (polar day in
  January): the coszen **night branch** is covered only synthetically
  (`night_zero`), as in the GSWP3 slice.
- The aerosol/deposition streams (ns 9/10 of the streams dump) are pure
  pass-throughs with unit scales, documented but not modeled (inactive
  in the C_JRA_t62 compset's datm exports consumed by the models).

## Units hygiene

Every literal in the JRA-lane templates and in the `Streams`/`Jra`
equations of `cdeps_datm_jra.esm` carries an explicit `units` field, so
`./esm units` reports every equation consistent (SI). Dimensional
readings: rdair carries J/(kg·K); the density qcoef and the albedo
coefficients are dimensionless; degtorad carries rad/deg so the albedo's
argument reads rad; the stream fields' own units vary per stream (the
Streams law is unit-agnostic, as in the GSWP3 slice).
