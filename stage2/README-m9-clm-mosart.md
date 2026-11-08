# M9: CLM runoff -> MOSART river routing — gate result: not feasible on current dumps

`stage2/m9_clm_mosart.esm` **was not written**. The Step-1 gate of the task
fails: the CLM-side runoff exports and the MOSART imports do **not** overlap
on any (cell, step) in the current stage-1 dumps, so there is nothing yet to
pin a coupled land+river test against. This README records the gate analysis
and the exact dump fields/points a stage-1 job would need to unblock Step 2.

## Gate question and answer

**Q:** Do the CLM-side runoff exports (qflx_surf/qflx_drain/qflx_runoff etc.
as passed to `lnd_export`) and the MOSART imports (`mosart_run_imports.csv`)
overlap on any (cell, step), and are the two runs bit-identical where they
overlap (same compset, same inputs)?

**A:** **No overlap.** The CLM lane (`I2000Clm60Sp_clm_inst`) dumped **no
runoff fields at all**: none of the 11 files in `stage1/dumps/clm/`
(`clm_qsat`, `clm_urbalb`, `clm_urban_lw`, `clm_urban_rad`, `clm_surfalbedo`,
`clm_surfrad`, `clm_bareground`, `clm_soilfluxes`, `clm_canopyflux`,
`clm_run_imports`, `clm_soiltemp`) carries `qflx_rofliq_*`, `qflx_surf`,
`qflx_drain`, `qflx_qrgwl`, `qirrig`, or any `Flrl_*`/`lnd_export` field
(verified by grepping all dump CSVs; the lane traced the urban / radiation /
surface-energy / soil-temperature path only, per
`stage1/runs/clm_instrument.json`). The MOSART imports exist on the **rof
grid** (r05 reaches `nt,nr,gindex`) at 8 coupling steps. Nothing in the CLM
dumps is on the rof grid, and nothing in the MOSART dumps records the
land-side values that produced them. The bit-identical question is therefore
vacuous — there is no shared dump point to compare.

Both runs are nonetheless the same configuration and were each verified
solution-neutral against the same reference case, so the underlying
trajectories agree where they overlap in time (see next section).

## The two lane runs (same configuration, different traces)

Both lanes ran compset `2000_DATM%GSWP3v1_CLM60%SP_SICE_SOCN_MOSART_DGLC%NOEVOLVE_SWAV_SESP`
at `f19_f19_mg17`, 2000-01-01, 1 day, 36 PEs, and the same inputdata
(verified in `stage1/cime/cases/{I2000Clm60Sp_clm_inst,I2000Clm60Sp_mosart_inst,I2000Clm60Sp_f19}/env_case.xml`:
identical COMPSET; both inst cases and the reference carry the `_SESP`
suffix). Coupling cadence is identical in both: `ATM_NCPL=48` (30 min;
`LND_NCPL=$ATM_NCPL`), `ROF_NCPL=8` (3 h = 10800 s coupling period)
(`env_run.xml` of both cases).

| | CLM lane (`I2000Clm60Sp_clm_inst`) | MOSART lane (`I2000Clm60Sp_mosart_inst`) |
|---|---|---|
| traces | 11 CSV files, f19 subgrid points (column/patch/landunit), PET 0, 12 steps of 1800 s | 9 CSV files, r05 reach points (`nt,nr,gindex`), PET 0, 8 coupling steps |
| runoff fields | **none** (urban/radiation/surface-energy/soil-temp path only) | `mosart_run_imports.csv`: qsur/qsub/qgwl/qirrig/qglc_liq/qglc_ice [m3/s] + volr/wt/wr/wh state at `mosart_run` entry (`mosart_driver.F90:500-507`) |
| neutrality | restart vs `I2000Clm60Sp_f19`: 290/292 vars bit-identical, 2 differ only in embedded case-name path strings (`clm_instrument.json` checks) | restart/history vs `I2000Clm60Sp_f19`: bit-identical modulo timestamps/case attributes (`mosart_instrument.json` checks) |

So: the two runs are effectively the same trajectory (each bit-identical to
the same reference on all shared state), and the mosart imports were
produced from CLM runoff values that are deterministic and reproducible in
the same config — but **no dump records those CLM values**, and the dump
point sets do not intersect (different grids, different points, different
cadence).

## The lnd->rof chain as implemented (what the future .esm would model)

Cited against the stage-1 forks (`stage1/EqESM-cesm-ctsm`,
`stage1/EqESM-cesm-mosart`, `stage1/EqESM-cesm-cmeps`; commit pins in the
lane `*_instrument.json` files):

1. **Production (column level, mm/s):** `HydrologyDrainageMod.F90:112-119`
   computes `qflx_surf_col` (surface runoff), `qflx_drain_col` (subsurface
   drainage), `qflx_drain_perched_col`, `qflx_qrgwl_col` (glacier/wetland/
   lake runoff); `IrrigationMod` computes `qflx_sfc_irrig_col`; snow capping
   produces `qflx_ice_runoff_col`.
2. **Gridcell accumulation** (`main/lnd2atmMod.F90:370-470`, `subroutine
   lnd2atm`): `qflx_rofliq_qsur_grc = c2g(qflx_surf_col, urbanf, unity)`,
   `qflx_rofliq_qsub_grc = c2g(qflx_drain_col)`, `qflx_rofliq_drain_perched_grc`,
   `qflx_rofliq_qgwl_grc = c2g(qflx_qrgwl_col) + qflx_liq_dynbal_grc`,
   `qirrig_grc = c2g(qflx_sfc_irrig_col)`, `qflx_rofice_grc = c2g(qflx_ice_runoff_col)
   + qflx_ice_dynbal_grc`. The `waterlnd2atmbulk` fields are time-averaged
   over the 30-min export interval by the standard accumulator machinery
   (`clm_driver.F90:1480-1500` `UpdateAccVars`).
3. **CLM export** (`cpl/nuopc/lnd_import_export.F90:895-935`, the "output to
   river" block of `export_fields`): `Flrl_rofsur = qflx_rofliq_qsur_grc`,
   `Flrl_rofgwl = qflx_rofliq_qgwl_grc`, `Flrl_rofi = qflx_rofice_grc`,
   `Flrl_irrig = -qirrig_grc` (negative demand), `Flrl_rofsub =
   qflx_rofliq_qsub_grc + qflx_rofliq_drain_perched_grc` (+
   `qflx_rofliq_stream_grc` when `use_hillslope_routing`). Units kg/m^2/s of
   water = mm/s.
4. **Mediator** (`components/cmeps/mediator/med_phases_prep_rof_mod.F90`):
   accumulates the 6 lnd exports per 3-h rof window into `FBlndAccum2rof_l`
   (L261), averages (L332-340), maps lnd->rof with `maptype_lnd2rof =
   mapconsf` (conservative flux mapping, L351, L514; `mapfcopy` fallback),
   with a normalized-map special case for irrigation using `lfrac`
   (L625-665), producing `x2r_Forr_rof`, `x2r_Forr_rof_subsurface`,
   `x2r_Forr_rof_gwl`, `x2r_Forr_rof_irrig` on the r05 grid (+
   `Fgrg_rofl`/`Fgrg_rofi` from glc when `rof_from_glc`).
5. **ROF import conversion** (`rof_import_export.F90:589-595`
   `state_getimport`): `q [m3/s] = field [mm/s] * med2mod_areacor * area [m2]
   * 0.001` with `med2mod_areacor = mesh_areas/model_areas` (F90:210). This
   conversion is already implemented and tested standalone in
   `stage2/mosart_coupling.esm` (`RofCap`, pins the step-1 cell-100312 rows
   at rel 1e-12).
6. **MOSART run entry** (`mosart_driver.F90:500-507`): the dump call that
   produced `mosart_run_imports.csv` (the only existing dump of the chain).

## What mosart_run_imports.csv contains (the reference side that exists)

256 data rows = 8 coupling steps x 16 cells x 2 tracer rows (nt=1 LIQ, nt=2
ICE; scalar columns repeat on both). Traced cells: nr 1..7 land (mask=1:
gindex 100312, 100313, 100314, 101032, 101033, 101034, 101035), nr 8..16
ocean (mask=2: gindex 101515..101523). Runoff is nonzero **only at coupling
step 1** (ymd 20000101, tod 10800, i.e. the 00:00-03:00 window covering CLM
steps 1..6), LIQ rows only:

| (step, nt, nr, gindex) | qsur [m3/s] |
|---|---|
| (1, 1, 1, 100312) | 7.51365591037279353e-05 |
| (1, 1, 2, 100313) | 3.75682795518607151e-05 (= exactly half of nr 1) |
| (1, 1, 4, 101032) | 7.53755999886219255e-05 |
| (1, 1, 5, 101033) | 3.76877999943059754e-05 (= exactly half of nr 4) |

`qsub = qgwl = qirrig = qglc_liq = qglc_ice = 0` on all 256 rows (so the
subsurface-runoff, irrigation-demand and glc->rof paths are unexercised in
this run). The qsur value of cell 100312 is the one already pinned by the
`mosart_coupling.esm` test (which derives `Flrl_rofsur =
qsur/(area*mm_to_m) = 2.59072497089951852e-11 mm/s` over `area =
2.900213644740529e9 m2`).

## Exactly which dump fields/points a stage-1 job would need

To unblock Step 2 (a `stage2/m9_clm_mosart.esm` mounting the CLM export side
+ `mosart_coupling.esm:RofCap` with `variable_map` couplings implementing the
mediator lnd->rof mapping, tested against `mosart_run_imports.csv`), a
stage-1 job in the same I2000Clm60Sp f19 configuration must produce **two
new dumps**:

**A. CLM export side (clm fork, `instrument/clm-dumps` — owned by the CLM
lane; new subroutine in `main/clm_dump_mod.F90` following the existing
`dump_*` pattern, e.g. `clm_rofliq_exports.csv`):**
- fields (exactly the six exported in `lnd_import_export.F90:895-935`,
  [mm/s]): `qflx_rofliq_qsur_grc`, `qflx_rofliq_qsub_grc`,
  `qflx_rofliq_drain_perched_grc`, `qflx_rofliq_qgwl_grc`, `qirrig_grc`,
  `qflx_rofice_grc`; optionally `qflx_rofliq_grc` (total) for the
  qsur+qsub+qgwl balance;
- points: PET 0, **all** land gridcells (`begg:endg`, ~384 of 13824 with 36
  PEs), at **each lnd export step s = 1..48** (all six 30-min exports per 3-h
  rof window are needed to reproduce the mediator average), carrying the
  gridcell global index, lat/lon, area [m2] and land fraction. The traced
  r05 cells (gindex 100312/100313/101032/101033) are a mid-latitude basin —
  **not** the South-Pole first-16 cells the CLM lane currently traces — so a
  first-16-cells filter would miss the contributing f19 cells; the
  contributing stencil can be read from the lnd2rof `mapconsf` mapping file;
- optionally (for a future CLM hydrology component test, and identical to
  the CTSM lane's planned Infiltration/TotalSurfaceRunoff/Drainage dump,
  `stage1/runs/CTSM_LANE.md`): the column-level production fluxes
  `qflx_surf_col`, `qflx_drain_col`, `qflx_drain_perched_col`,
  `qflx_qrgwl_col`, `qflx_sfc_irrig_col`, `qflx_ice_runoff_col` on the same
  filter as `clm_soilfluxes` (first 16 non-lake columns of PET 0), to pin
  the `c2g` urbanf/unity aggregation.

**B. The regridded values MOSART actually received (one of two points):**
- (a) **rof side (mosart fork — owned by the MOSART lane):** the raw
  `Flrl_*` fields as received, i.e. the `fldptr` values in
  `state_getimport` (`rof_import_export.F90:589-595`) after the
  `med2mod_areacor` correction and before `* area * 0.001`, for the same 16
  cells x 8 steps x 2 tracers filter as `dump_run_imports` — this pins the
  mediator mapping and the `RofCap` conversion in one place; or
- (b) **cmeps side:** `x2r_Forr_rof`, `x2r_Forr_rof_subsurface`,
  `x2r_Forr_rof_gwl`, `x2r_Forr_rof_irrig` on the rof grid right after the
  map in `med_phases_prep_rof_mod.F90` (L351-359), same points.

With (A) + (B) the coupled test pins the full chain: f19 gridcell values at
lnd export steps -> (mediator average of 6 + `mapconsf` mapping, written as
a `variable_map` with the mapping-file weights; irrigation via the
normalized lfrac mapping) -> `RofCap` conversion -> qsur/qsub/qgwl/qirrig
against `mosart_run_imports.csv` at rel 1e-12, on the rows of the table
above.

## Coverage caveats for the future test

- Only 4 (cell, step) rows are nonzero; everything else in
  `mosart_run_imports.csv` is identically zero. The nonzero qsur rows pin
  the chain; the zero rows can only serve as wiring checks (they cannot
  fail on value), and the qsub/qgwl/qirrig/qglc_* paths stay unpinned
  unless a longer/different-season run (or a `rof_from_glc` config)
  produces nonzero values.
- The mediator irrigation branch (`med_map_field_normalized` with `lfrac`,
  `med_phases_prep_rof_mod.F90:625-665`) is unexercised (`qirrig = 0`).
- No CESM defect was found during this gate analysis; the missing CLM
  runoff dump is a dump-coverage gap, so nothing is appended to `bugs.md`.

## Step 2 status

Not started (gated). The pieces that already exist for it:
`stage2/mosart_coupling.esm` (RofCap import/export conversions, merged) and
its test pins the same cell-100312 step-1 rows referenced above; the CLM
export law (items 1-3 of the chain) has no stage-2 component yet and is the
CLM lane's planned hydrology slice; the mediator mapping is M8/M9 work.
