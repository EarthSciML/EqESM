# CDEPS lane — data atmosphere model (DATM, milestone M8, data_sources part)

Parallel to the MOSART/CICE/CLM lanes. Scope: **CDEPS data models**, first
slice = **DATM with the GSWP3v1 forcing and the CLMNCEP datamode** — the
atmospheric forcing for the I-class land reference runs (and, in the coupled
model, the boundary condition CLM consumes). The deliverable is a stage-2
`data_sources`-based DATM component (`stage2/cdeps_datm_gswp3.esm` +
`stage2/cdeps_templates.esm`) whose equations are pinned to instrumented
CESM DATM traces, plus this lane record and `stage1/runs/cdeps_instrument.json`.

DATM is the natural first CDEPS slice: the GSWP3v1→CLMNCEP path is pure
**forcing production** (time interpolation of observed data + a handful of
diagnostic derivations), it has no prognostic state, and it is already
exercised end-to-end by the `I2000Clm60Sp` reference run
(`stage1/cime/cases/I2000Clm60Sp_f19`, job 11826959).

## Ownership / collision guards

- **Ours:** `stage2/cdeps_*.esm`, `stage2/README-cdeps.md`, fork
  `stage1/EqESM-cesm-cdeps` (branch `instrument/cdeps-dumps`),
  `stage1/dumps/cdeps/`, `stage1/runs/cdeps_instrument.json`, this file.
- **Peer lanes:** MOSART (`stage2/mosart_*.esm`), CICE
  (`stage2/cice_*.esm`), CLM (`stage2/clm_*.esm`). The DATM dumps feed the
  CLM lane's forcing boundary, but the two lanes are independent.
- **Coordinate:** `PROGRESS.md` (contended — add only a short pointer),
  `lib/cesm_constants.esm`, ESD templates, the shared cluster queue and
  `stage1/build`.

## DATM data flow (what we pin)

```
GSWP3v1 netCDF files (0.5°, 3-hourly)
   │  shr_stream_findBounds / readLBUB  (streams/dshr_stream_mod.F90,
   │                                     streams/dshr_strdata_mod.F90)
   ▼
raw stream values on the source grid, LB + UB in time
   │  spatial map (mapalgo=bilinear, ESMF regrid) → model grid
   ▼
fldbun_data(LB), fldbun_data(UB)                     ── traced (cdeps_strdata_streams.csv)
   │  time interpolation, per stream tintalgo:
   │    coszen  : v = LB·cosz/tavCosz  (0 at night)  (streams/dshr_strdata_mod.F90:1292)
   │    linear  : v = flb·LB + fub·UB                 (streams/dshr_strdata_mod.F90:1384)
   │    nearest : v = LB or UB                        (dshr_tinterp_mod.F90:36)
   │    lower   : v = LB
   ▼
fldbun_model (model grid, model time)                ── traced (datm_gswp3_streams.csv)
   │  datm_datamode_clmncep_advance derivations (datm/datm_datamode_clmncep_mod.F90)
   ▼
export state (Sa_*, Faxa_*)                          ── traced (datm_gswp3_exports.csv)
   │  mediator → CLM
```

Entry points: `atm_comp_nuopc.F90 datm_comp_run` (namelist, advance loop,
restart) → `shr_strdata_advance` (all streams advanced together) →
`datm_datamode_clmncep_advance` (per-point derivations).

## Stream configuration (I2000Clm60Sp, from Buildconf/datmconf/datm.streams.xml)

| # | stream | tintalgo | readmode | mapalgo | fields (file var → model var) |
|---|---|---|---|---|---|
| 1 | CLMGSWP3v1.Solar | **coszen** | single | bilinear | FSDS → Faxa_swdn |
| 2 | CLMGSWP3v1.Precip | **nearest** | single | bilinear | PRECTmms → Faxa_precn |
| 3 | CLMGSWP3v1.TPQW | **linear** | single | bilinear | TBOT → Sa_tbot, WIND → Sa_wind, QBOT → Sa_shum, PSRF → Sa_pbot, FLDS → Faxa_lwdn |
| 4 | presaero.clim_2000 | linear | single | bilinear | 14 BCDEPWET/…/DSTX04DD → Faxa_bcphiwet…Faxa_dstdry4 |
| 5 | presndep.clim_2000_cmip6 | linear | single | bilinear | NDEP_NHx_month → Faxa_ndep_nhx (×1e-3), NDEP_NOy_month → Faxa_ndep_noy (×1e-3) |
| 6 | preso3.clim_2000 | linear | single | bilinear | O3 → Sa_o3 |
| 7 | topo.observed | **lower** | single | bilinear | TOPO → Sa_topo |

Streams 1–3 are the GSWP3 forcing (the stage-2 scope); 4–6 are the
prescribed aerosol/ndep/ozone climatologies (stage-2 **noted as separate**
data_sources; their fields pass through unchanged — `datm_pres_*_mod.F90`
is a pure pass-through with unit-scale application — and are not traced
beyond the streams layer); 7 supplies Sa_topo with a `lower` tintalgo
(constant within a year).

All GSWP3 files: `clmforc.GSWP3.c2011.0.5x0.5.*` 1991–2010,
`taxmode=cycle`, `year_align=1991`, `dtlimit=1.5`, 3-hourly
(`time` units "days since <month>-01-01" per file; 248 records/month),
on the shared TPQWL SCRIP/ESMF mesh (720×360). The datm model mesh is
`fv1.9x2.5_141008_ESMFmesh.nc` (144×96), so every stream value the
datamode sees is the bilinear-map of the source field onto the model
grid — the stage-2 `data_sources` must reproduce this mapping
(`EarthSciData` regrid/ESMF routehandle or an equivalent discretized
map, PRs to `EarthSciML/EarthSciDiscretizations` if needed).

## Field derivations (CLMNCEP datamode) — source of truth

All in `datm/datm_datamode_clmncep_mod.F90` (pristine cdeps1.0.112;
line numbers from that file). Stream values are already model-grid,
model-time interpolants. `tbotmax` = global max of `strm_Sa_tbot` on the
**first** call (units detector), `anidrmax` = global max of the `Sx_anidr`
import (existence detector).

| export | law | F90 ref |
|---|---|---|
| Sa_tbot | = strm_Sa_tbot; **if tbotmax < 50 K add TKFRZ (273.15)**; clamp `max(180,·)` | :328,:384,:386 |
| Sa_ptem | = Sa_tbot | :387 |
| Sa_topo | = strm_Sa_topo (1e30 spval, SHR_CONST_SPVAL, if stream absent) | :330 |
| Sa_z | = strm_Sa_z if present else **30.0 m** | :377 |
| Sa_pbot | = strm_Sa_pbot; **0 → pstd (101325 Pa)**; pstd if stream absent | :390 |
| Sa_pslv | = Sa_pbot | :399 |
| Sa_u, Sa_v | = wind/√2 each | :402 |
| Sa_shum | from Sa_shum stream: qsat = 0.622·e/(pbot − 0.378·e), e = esat(T), **Sa_shum = min(qsat, stream)**; RH stream: e = RH·0.01·esat(T); tdew stream: esat(Tdew) (Tdew + TKFRZ if tdewmax<50) | :408,:416,:420 |
| Sa_dens | vp = q·p/(0.622+0.378·q); ρ = (p − 0.378·vp)/(T·Rd), Rd = 287.042 J/kg/K | :428 |
| Faxa_lwdn | = strm_Faxa_lwdn (formula fallback ea·σ·T⁴ only if stream absent) | :432 |
| Faxa_swdn | = strm_Faxa_swdn | :329 |
| Faxa_swndr/vdr | 0.5·swdn·ratio_rvrf; ratio_rvrf = clamp(0.29548 + 0.00504·s − 1.4957e-5·s² + 1.4881e-8·s³, 0.01, 0.99), s = 0.5·swdn | :450 |
| Faxa_swndf/vdf | 0.5·swdn·(1 − ratio_rvrf) (vis coefficients: 0.17639 + 0.00380·s − 9.0039e-6·s² + 8.1351e-9·s³) | :454 |
| Faxa_swnet | Σ(1−alb)·sw per band (Sx_anidr/Sx_avsdr/Sx_anidf/Sx_avsdf imports); 0 if anidrmax < 1e-8 or ≥ 1e30·0.9 (spval 1e30 = SHR_CONST_SPVAL); plain sum if imports absent | :466 |
| Faxa_rainc/l | = 0.1/0.9 · strm_Faxa_precn | :479 |
| Faxa_snowc/l | = rainc/l · (1 − frac); frac = clamp((T−TKFRZ)·0.5, 0, 1) (`shr_precip_mod.F90:17`); then rainc/l ×= frac (all `max(0,·)`) | :483 |
| esat(T) | water poly a0…a6 / ice poly b0…b6, T clamped [−50, 50] °C, ×100 Pa; water if Tbot ≥ TKFRZ else ice | :545 |

Bias-correction / anomaly-forcing blocks (:498–:539) are inactive for the
GSWP3 case (no `bias_correct`/`anomaly_forcing` streams), so the stage-2
component covers the unmodified path; the blocks are documented but not
traced.

### Time interpolation laws (streams layer)

- **coszen** (`dshr_strdata_mod.F90:1292`): v(n) = LB(n)·cosz(n)/tavCosz(n)
  if cosz(n) > 0.001 else 0. cosz from `shr_tInterp_getCosz`
  (`dshr_tinterp_mod.F90:276`, clamped to ≥ 0.001; declination/eccentricity
  from the orbital params set by `shr_strdata_setOrbs`); tavCosz = mean of
  cosz over the LB→UB interval (`shr_tInterp_getAvgCosz`, :162).
- **linear** (:1384): v = flb·LB + fub·UB, factors from
  `shr_tInterp_getFactors` (`dshr_tinterp_mod.F90:36`).
- **nearest** (same call, `algo='nearest'`): f = (1,0) or (0,1) by which
  bound is closer in time (ties → LB).
- **lower/upper**: (1,0) / (0,1).
- Calendar handling: model = `noleap` (CLM), stream = `noleap` — no
  Feb-29 adjustment for this case (`ymdmod = ymd`).

## Instrumentation (fork `stage1/EqESM-cesm-cdeps`, branch `instrument/cdeps-dumps`)

Fork from pristine `stage1/cesm` (never modified); cdeps submodule
`@f80e1f2` (1 commit on `cdeps1.0.112`), superproject `@a94e749`. All
dumps are **additive** (pure reads + CSV writes on PET 0, no state
change) and **namelist-gated** (`dump_datm`, default `.false.`), so a
default build/run is bit-identical to un-instrumented DATM.

| file | change |
|---|---|
| `datm/datm_dump_mod.F90` (new) | datamode-layer dumps: `datm_gswp3_streams.csv` (interpolated stream values as consumed), `datm_gswp3_exports.csv` (export state + Sx_* imports, spval 1e30 = SHR_CONST_SPVAL when absent). PET 0, first `dump_first_cells` mesh points × first `dump_first_steps` calls; `datm_dump_init` from `init_pointers` (captures lon/lat), `datm_dump_finalize` from `ModelFinalize`. |
| `datm/datm_datamode_clmncep_mod.F90` | `datm_dump_init` at end of `init_pointers`; `advance` signature gains `target_ymd, target_tod`; dump calls at end of `advance` (after the bias/anomaly block), gated on `mainproc`. |
| `streams/dshr_strdata_mod.F90` | `shr_strdata_dump_set` / `shr_strdata_dump_finalize` (public); per-stream dump `shr_strdata_dump_stream` called at the end of the time-interp loop in `shr_strdata_advance` (gated on `istr=='datm'` + mainproc): `cdeps_strdata_streams.csv` with LB/UB values, flb/fub factors and the interpolated value, per (stream, field, point). coszen factor recomputed with the same `shr_tInterp_getCosz` call. |
| `datm/atm_comp_nuopc.F90` | `dump_datm/dump_dir/dump_first_cells/dump_first_steps` in the existing `datm_nml` group (read+broadcast: bcasttmp 10→12 + `dump_dir` CX/CL broadcast); forwarded to both dump layers; `datm_dump_finalize` + `shr_strdata_dump_finalize` in `ModelFinalize`. |
| `datm/cime_config/namelist_definition_datm.xml` | 4 `stage1_instrumentation` entries, group `datm_nml`. |
| `datm/CMakeLists.txt` | add `datm_dump_mod.F90`. |

Dump CSVs land in `stage1/dumps/cdeps/` (gitignored, regenerable; set by
`dump_dir` in `user_nl_datm`).

## Run (regeneration recipe)

```sh
source stage1/runs/env.sh
python3 stage1/cesm/cime/scripts/create_newcase \
  --case stage1/cime/cases/I2000Clm60Sp_cdeps_inst \
  --compset I2000Clm60Sp --res f19_f19_mg17 --machine illinois \
  --run-unsupported --srcroot $PWD/stage1/EqESM-cesm-cdeps
cd stage1/cime/cases/I2000Clm60Sp_cdeps_inst
printf '\ndump_datm = .true.\ndump_dir = "<repo>/stage1/dumps/cdeps/"\ndump_first_cells = 16\ndump_first_steps = 12\n' >> user_nl_datm
./xmlchange STOP_OPTION=ndays,STOP_N=1,JOB_WALLCLOCK_TIME=02:00:00,TOTALPES=36
./case.setup && ./case.build
env -u SLURM_MEM_PER_NODE -u SLURM_MEM_PER_CPU -u SLURM_MEM_PER_GPU ./case.submit
```

Notes: `--srcroot` must point at the fork worktree (CIME scripts come from
the pristine checkout via `stage1/cesm/cime`); the fork's submodules
(`git submodule update --init`, plus cdeps's `fox`/`share/genf90` and
clm's `src/fates` nested submodules) must be checked out. Same compset /
grid / namelists as the reference `I2000Clm60Sp_f19` run, so restart files
are comparable for the solution-neutrality check.

## Waves

- **W0 — recon** ✅: DATM code study (all derivations above), reference
  case `I2000Clm60Sp_f19` COMPLETED (job 11826959), GSWP3 inputdata on
  disk, stream config extracted from the reference case's
  `Buildconf/datmconf/datm.streams.xml`.
- **W1 — instrumentation + run** ✅. W1a: fork + dump module + run — two
  runs of `I2000Clm60Sp_cdeps_inst` (jobs 11979948, 11982594), both
  COMPLETED with CSVs in `stage1/dumps/cdeps/` (4800 + 192 + 192 data
  rows; no blank lines; spval text corrected to 1e30 = SHR_CONST_SPVAL in
  cdeps commit 9d6dc95 — the run-2 CSVs are the deliverable). W1b:
  `stage1/runs/cdeps_instrument.json` + this file + the check scripts
  `stage1/analysis/check_datm_neutral.py` / `check_datm_laws.py` ✅.
  Verification: restart files bit-identical to the reference run
  (all 5 components; the instrumented case ran on a different pes
  decomposition, 80 tasks/2 nodes vs 36/1 node, and is still
  bit-identical); all 4800 streams-layer rows and 192 datamode rows
  verified to round-off against the Fortran laws (0 failures).
- **W2 — dumps → stage-2 component** ✅. `stage2/cdeps_templates.esm`
  (shared laws: esat polynomials, qsat, vapor_pressure, density,
  ratio_rvrf NIR/VIS, rain/snow ramp, wind split) +
  `stage2/cdeps_datm_gswp3.esm` (`data_sources` for the 3 GSWP3 streams +
  the topo stream, the `Streams` interpolation model, the `Clmncep`
  transformation model, the `GSWP3Files` netcdf-reader model) +
  `stage2/README-cdeps.md`. Tests use the dumped LB/UB+factor rows so the
  interpolation and the derivations are both pinned (7 + 7 + 1 tests =
  58 assertions, all passing); the netCDF-reading tests exercise the
  `data_sources` path against the actual GSWP3 files (CONST sources,
  `metadata.esio_format = netcdf`, `select` with the time axis at "all"
  and fixed lat/lon, `extent.metaparameter` + `index_sets` + `coords`
  assertions 1-based; each source needs its own extent metaparameter —
  topo delivers 1 record, the 3-hourly files 248). `./esm validate` clean
  for both files, `./esm test stage2/` green apart from the three
  pre-existing `${ESD_ROOT}` cice template-import errors (CICE lane's
  pending EarthSciDiscretizations PR).
- **W3 — PRs** into EarthSciModels (`cdeps_datm_gswp3.esm` as a new
  component; `data_sources`/regrid PRs as needed).

## Risks

- `data_sources` in esm 0.7.0 is young: the netcdf reader rejects
  time-axis `select` (provider owns cadence), DISCRETE sources fail
  static checks, and only CONST (build-time) sources with `extent` +
  `index_sets` + `coords` assertions are proven; the temporal
  interpolation of the stage-2 model is therefore expressed in the
  model equations (LB/UB from the file's records, factors from the
  dumped laws) rather than in `data_sources` `temporal` blocks.
- The dump traces the **model-grid** values (bilinear map already
  applied). A stage-2 test of the spatial map needs the regrid weights —
  out of scope for this slice (noted; ESD regrid PRs later).
- Shared cluster queue: serialize submissions with the peer lanes.
- The datamode `tbotmax`/`tdewmax` units detectors are global-max logic
  — a test must reproduce the max over the traced points, which is a
  property of the whole field, not per-point; keep `tbotmax`/`anidrmax`
  columns in the dumps (done).

---

# JRA55 datamode slice (2026-11-05) — C_JRA_t62_datm_inst

Same stream layer, different datamode: DATM%JRA-1p5-2023 with the JRA55
v1.5 daily TL319 forcing (8 streams: prec, lwdn, swdn, q_10, slp, t_10,
u_10, v_10 → strm_prec/lwdn/swdn/shum/pslv/tbot/u/v) on the T62 datm mesh
(grid a%T62_l%null_oi%tx2_3v2, mesh T62_040121_ESMFmesh.nc), compset
2000_DATM%JRA-1p5-2023_SLND_DICE%SSMI_MOM6_DROF%JRA-1p5-2023_SGLC_SWAV_SESP.
All streams taxmode=cycle, year_first=1958, year_last=2023, year_align=1
(model year 1 = data year 1958; the JRA datamode has no climatology
option), readmode=single, mapalgo=bilinear.

## JRA datamode laws (datm_datamode_jra_mod.F90:271-322) — source of truth

- Copies: Sa_pslv = Sa_pbot = strm_pslv; Sa_tbot = Sa_ptem = strm_tbot;
  Sa_u/Sa_v = strm_u/strm_v; Sa_shum = strm_shum (NO qsat cap — the JRA
  datamode has no humidity law); Faxa_swdn/lwdn = strm_swdn/lwdn.
- Sa_z = 10 m; Sa_u10m/v10m = Sa_u/Sa_v.
- Sa_dens = pbot/(rdair·tbot·(1 + qcoef·shum)), rdair = SHR_CONST_RDAIR =
  287.0423113650487, qcoef = **r4(0.608)** = 0.6079999804496765 (see
  bugs.md: the 0.608 literal is unsuffixed → single precision).
- Precip split: when the streams separate rain/snow (prrn/prsn present),
  Faxa_rainl/snowl = strm_prrn/prsn passthrough; otherwise a step:
  Faxa_rainl = prec, Faxa_snowl = 0 when tbot ≥ TKFRZ, else the reverse.
  Faxa_rainc = Faxa_snowc = 0 (no convective precip in the datm).
- Shortwave: Faxa_swvdr/ndr/vdf/ndf = swdn·0.28/0.31/0.24/0.17 (exact
  _R8 literals); avg_alb = **r4(0.069)** − **r4(0.011)**·cos(2·lat·degtorad)
  (LY08; the coefficients are unsuffixed → single precision, see bugs.md);
  Faxa_swnet = swdn·(1 − avg_alb). lat = the model-mesh element latitude
  yc (ESMF ownedElemCoords of the T62 mesh — verified by a direct ESMF
  8.9.1 probe to equal the file centerCoords, i.e. the dumped lat column).
- Dead code: cosfactor/phs_c0/dLWarc computed but never used (bugs.md).

## Instrumentation delta (fork branch instrument/cdeps-dumps)

The datm_dump_mod tag (GSWP3 slice) extended with a `_jra` suffix variant
and the JRA55-specific columns: datm_jra_streams.csv (strm_* values +
separates_rain_snow flag) and datm_jra_exports.csv (avg_alb + all
exports) from datm_datamode_jra_mod.F90 (mainproc arg + allocatable
avg_alb added), plus the cdeps_strdata_streams_jra.csv tag suffix in
dshr_strdata_mod.F90. Commits 2de56ba (signature fix), 2c40df2 (JRA
dumps), adf3f96 (superproject pointer).

## Run (2026-11-05)

Case `C_JRA_t62_datm_inst` (job 12132709, COMPLETED 7:32; earlier
attempts failed on the SLURM_MEM_* env conflict — fixed with
`MPI_RUN_COMMAND='env -u SLURM_MEM_PER_CPU -u SLURM_MEM_PER_GPU -u
SLURM_MEM_PER_NODE srun --mpi=pmix -n 72'`), STOP_N=1 ndays, ATM_NCPL=24,
TOTALPES=72. user_nl_datm: dump_datm=.true., dump_first_cells=16,
dump_first_steps=12. user_nl_mom: AUTO_MASKTABLE = False (required on
this grid — MOM6 init FATALs otherwise; the reference case has it too).
Dumps: datm_jra_streams.csv (192 rows), datm_jra_exports.csv (192 rows),
cdeps_strdata_streams_jra.csv (4608 rows) in stage1/dumps/cdeps/; GSWP3
dumps untouched. **Neutrality: all 5 restart components bit-identical to
the reference C_JRA_t62 run** (check_datm_neutral.py pattern).

## Stage 2 (worktree branch wt/cdeps-datm-jra)

`stage2/cdeps_datm_jra.esm` + the `stream_interp`/`density_jra`/
`rain_frac_jra`/`jra_albedo` templates appended to
`stage2/cdeps_templates.esm`; `stage1/analysis/check_datm_jra_laws.py`
verifies every dump row to round-off (0 failures). Reader path: the 8
JRA55 v1.5 TL319 files (2920 daily records, lat ascending from −89.57,
netcdf3 classic float32) via the esm netcdf reader (CONST sources,
record 1 = first day of 1958). The r4-literal finding (0.608/0.069/0.011)
is the reason the albedo/density laws differ from a naive transcription
at the 2e-8/2e-11 level — the templates and tests carry the r4-rounded
constants with provenance notes.

## Risks

- The traced 16 points are all on the southernmost T62 row (polar day in
  January): the coszen night branch is covered only synthetically.
- Same data_sources caveats as the GSWP3 slice (CONST sources only,
  cadence owned by the Provider, spatial map out of scope).
- The JRA files' 2920-record cadence (daily, first record at 1.5 h)
  makes the reader tests heavier than GSWP3's; if the full `esm test`
  run stays slow, the JRAFiles tests can be dropped without losing the
  law coverage (the interpolation laws are pinned by the Streams tests).
