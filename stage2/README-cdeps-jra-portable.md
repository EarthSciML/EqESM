# CDEPS JRA55 — portable data sources and shared `stream_interp` (M8 follow-up)

Slice of the CDEPS DATM lane (branch `wt/cdeps-jra-portable`): makes the
JRA55 data sources portable (public NCAR URLs instead of cluster-local
`file://` paths), drops the slow JRAFiles reader tests, and moves the GSWP3
`Streams` model onto the shared `stream_interp` template. The lane README
(`README-cdeps-jra.md`) is managed by the lane owner; this file records the
slice's changes, gaps, and the CLI behavior found.

## Changes

| file | change |
|---|---|
| `stage2/cdeps_datm_jra.esm` | the eight `JRA_*` `data_sources` now point at the public NCAR CESM inputdata server with `[year]` url templates, `temporal` blocks, and `reference` blocks; the two `JRAFiles` netcdf-reader tests are removed (`tests: []`); metadata description updated |
| `stage2/cdeps_datm_gswp3.esm` | the `Streams` model's inline time-interpolation law is replaced by one `apply_expression_template` call to the `stream_interp` template of `cdeps_templates.esm` (imported by reference); tests unchanged |
| `stage2/cdeps_templates.esm` | untouched (append-only; the `stream_interp` template already exists) |

### JRA data sources (public NCAR archive)

The case's `datm.streams.xml` (`C_JRA_t62_datm_inst`) reads the
`JRA.v1.5.<var>.TL319.<year>.210504.nc` files of `ocn/jra55/v1.5_noleap/`.
Every source now resolves against

`https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/ocn/jra55/v1.5_noleap/`

with the template `JRA.v1.5.<var>.TL319.[year].210504.nc`, a `temporal`
block (start 1958-01-01, end 2020-12-31 — the archive's holdings with the
`.210504` naming, verified by directory listing for all eight variables;
2021+ files carry different suffixes and are **not** on the server —
`file_period P1Y`, `frequency PT3H`, 2920 records/file, time variable
`time`, the file's own `hours since 1958-01-01` axis), and a `reference`
block (Kobayashi et al. 2015, doi:10.2151/jmsj.2015-001, archive URL).
HTTP 200 verified via `curl -sI` for all eight 1958 files and the listing.

**URL-template syntax note (CLI finding):** the esm 0.7.0 CLI's linked
reader (earthsciio 0.1.4) expands only the time-crate `[year]`-style
placeholders in `url_template`, **not** the esm-spec `{date:...}` token
form: with `{date:%Y}` the transport fetches the literal template and the
server answers HTTP 400 (observed), ERRORing every assertion in the
document at build (extent discovery). The EarthSciModels port's
`{date:%Y-%m}` templates therefore cannot be ingested by this CLI either
(their "network access at test time" claim is unverified for the https
form). This document uses the `[year]` form — the one the repo's CLI
expands — and says so in the metadata description. `[year]` remains
expandable by the newer EarthSciIO (its `[`-form branch is retained).

### JRA frequency correction (2026-11 review round)

The eight `JRA_*` `data_sources` declared `frequency: "P1D"` next to
`records_per_file: 2920` — inconsistent (2920 records/yr is 3-hourly in a
365-day noleap year).  Verified against the local copies under
`$M/stage1/inputdata/ocn/jra55/v1.5_noleap/`: the `time` variable of
`JRA.v1.5.lwdn.TL319.1958.210504.nc` runs `hours since 1958-01-01` in
steps of 3.0 h (1.5, 4.5, 7.5, …; `time_bnds[0] = [0, 3]`), and the file
history names the source `prra_input4MIPs_atmosphericState_OMIP_
MRI-JRA55-do-1-5-0` — the 3-hourly JRA-55-do means (Tsujino et al. 2018),
not daily values.  Corrected `frequency` to `PT3H` in all eight sources,
and the "daily TL319 files"/"daily mean" wording of the metadata
descriptions, citations, tags, `N_REC` description and the
`JRAFiles`-reader descriptions to "3-hourly".  The CDEPS `datm.streams.xml`
of the C_JRA_t62 case carries no frequency (its `dtlimit = 1.5` guards the
time-interpolation search only), so the record mapping is driven by the
file's own time axis and is unchanged; the JRA `Streams`/`Jra` assertions
still pass at the same values (the reader picks the record by the file
time axis, which `frequency` never fed).

### JRAFiles tests removed — coverage gap

The two `JRAFiles` reader tests (`first_record_cell_00`,
`last_record_cell_00`) read all 2920 records of the eight files (16
assertions, ~27 min of suite time for 2 tests) and were **removed**; the
`JRAFiles` model (variables, data-fed `update` wiring, extent binding to
`N_REC`) is unchanged and `tests` is now empty. Coverage gap: the raw
netcdf record values at the fixed cell are no longer asserted by the
suite. The values were verified during development (record 1 and 2920 of
every variable, exact float64 widenings of the file float32s — the same
numbers the removed tests asserted; the probe in `scratch/` reproduced
record 1 of `lwdn` = 142.08123779296875 from the remote URL). The
interpolation and datamode laws remain pinned by the `Streams` (6) and
`Jra` (10) tests. The reader path stays exercised downstream by consumers
with data access (as in the EarthSciModels port).

### GSWP3 `Streams` on the shared template

The GSWP3 `Streams` model now invokes `stream_interp` (one law, imported
by reference) instead of inlining the `ifelse` law — identical to the JRA
`Streams` model. All 58 assertions pass unchanged (the template body is
byte-identical to the inlined law). The `stream_interp` template
description in `cdeps_templates.esm` still says the GSWP3 model "carries
the identical inline equation" — now stale; left untouched (append-only
rule), the lane owner can reword it.

GSWP3 data sources are kept as they were (cluster-local `file://` + the
live `GSWP3Files.netcdf_records` test are EqESM-only), per the M8 task.
Its live test is not slow (the whole file's tests, including the 248-record
reader test, run in ~2m08s) — reported, left as-is.

## `esm test` with remote data URLs (task item 3 — observed behavior)

- `esm validate` and `esm units` are offline and pass.
- `esm test` builds every data-fed parameter up front (extent discovery
  samples each consumed source, even when no test reads it — the
  `JRAFiles.*` parameters), so running the JRA document **requires network
  access to the NCAR archive at build time**. No network (or a
  non-expandable template) ⇒ every assertion in the document errors at
  compile with the transport error. This is a property of the CLI build,
  not of the document (the old `file://` form had the same build path).
- First run downloads all eight 1958 files in full (~2.39 GB each, ~19 GB)
  into EarthSciIO's content-addressed cache. **The default cache root is
  `$TMPDIR/earthsci-esm-cache` — RAM-backed on this cluster; set
  `ESM_CACHE_DIR=/scratch.local/...` to avoid OOM.** Timing measured:
  ~21.6 s for one file end-to-end; the full JRA document (8 files +
  106 assertions) 2m25s cold; the whole `stage2/` suite 2m26s warm.
- Re-runs reuse the cache (conditional GET), no re-download.

## Suite numbers

Baseline on the starting commit (6884c43, recorded at da0a151):
2662P/0F/0E. Final: **2646P/0F/0E** — exactly the 16 assertions of the two
removed JRAFiles tests — validate-failures=0, units-issues=0.
