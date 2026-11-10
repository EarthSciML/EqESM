# CLM land→river runoff exports — stage 2 data & validator map

Stage-1 instrumentation lane for the CLM land→river (MOSART) coupling
(milestone M9 part A; ctsm5.4.063, cesm3.0-alphabranch; instrumented fork
commits `ada2cb03b` + fixes `dccc7bbb9` on branch `instrument/clm-dumps`):
the gridcell runoff exports at the `lnd_import_export` export call, with the
column-level production fluxes and c2g weights pinned alongside so the
`urbanf`/`unity` aggregation is reproducible exactly. This lane is at the
stage-1 → stage-2 boundary: the instrumented dumps exist, the validator
(`stage1/analysis/verify_clm_rofliq_dumps.py`) pins the Fortran laws on them,
and the `.esm` coupling components will be built against the same dumps (no
`.esm` files yet — see [Status](#status)).

The dumps come from the instrumented rerun of case `I2000Clm60Sp_clm_inst`
(in-node `case.submit --no-batch` in slurm allocation 12285825), model LID
**261110-045347** (run log `/scratch.local/ctessum/clmrofliq/run_rofliq2.log`;
the first attempt LID 261110-040834 crashed in the final step's export — see
[Bring-up bugs](#bring-up-bugs)). The run is solution-neutral: the
instrumented restart is bit-identical to the uninstrumented baseline
`I2000Clm60Sp_f19` restart (292/292 common numeric variables, 0 differ,
mask-aware compare), and all 15 pre-existing `clm_*.csv` are bit-identical
(md5) to the pre-lane backup (`/scratch.local/ctessum/clmrofliq/dumps_backup/`).

New namelist: `dump_rofliq` (default `.false.` = off; a default build/run is
bit-identical to un-instrumented CTSM) + `dump_rofliq_steps` (0 = all steps),
`mpi_bcast`'d in `control_spmd` together with the other `dump_*` variables
(the W2d deadlock lesson). Both dump routines are read-only traces gated by
`dump_rofliq`, called from serial call sites only, and write via PET 0 after
an all-PE `MPI_Gatherv` (so **all** 5666 f19 land gridcells / 28,494 columns
are dumped, not just PET 0's slice — the segmented decomp spreads PET-0
cells across the globe).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `stage1/dumps/clm/clm_rofliq_exports.csv` | Gridcell land→river exports at the export call, one row per gridcell per step (48 steps × 5666 cells): `Flrl_rofsur` = qflx_rofliq_qsur_grc, `Flrl_rofsub` = qflx_rofliq_qsub_grc + qflx_rofliq_drain_perched_grc, `Flrl_rofgwl` = qflx_rofliq_qgwl_grc, `Flrl_irrig` = −qirrig_grc (the mediator negates the demand), `Flrl_rofi` = qflx_rofice_grc, all [mm/s] = [kg/m²/s] positive downward; plus the total qflx_rofliq_grc and the dynamic-land-cover-change corrections qflx_liq/ice_dynbal_grc (0 in this run); and gindex/lat/lon/area [m²] = ldomain%area [km²]·1e6/lfrac | `lnd_import_export.F90`: export_fields "output to river" block :890-935 (dump hook after the Flrl_rofsub export) |
| `stage1/dumps/clm/clm_rofliq_cols.csv` | Column-level production fluxes and c2g weights at the end of the runoff section of lnd2atm, one row per column per step (48 × 28,494): qflx_surf_col, qflx_drain_col, qflx_drain_perched_col, qflx_qrgwl_col (already including qflx_liq_from_ice_col), qflx_sfc_irrig_col, qflx_ice_runoff_col + wtgcell, ctype, canyon_hwr, active + gindex_grc/gindex_col | `lnd2atmMod.F90`: runoff section :381-470 (dump hook :470-487); `subgridAveMod.F90`: c2g_1d :761-818, set_c2l_scale :104-165 |
| `stage1/analysis/verify_clm_rofliq_dumps.py` | Validator: pins the c2g aggregation laws on both CSVs at round-off, reproduces the first MOSART coupling-step runoff imports from the exports + mesh geometry + mediator accumulation, verifies the 15 pre-existing dumps bit-identical (md5) against a backup, and the restart neutrality (mask-aware compare) | — |

## The traced laws (validated)

All c2g identities hold to round-off (max rel 4e-16, es26.17 dump
precision). For every step and gridcell the validator recomputes

```
sumwt = Σ_c wtgcell_c                      (active columns; == 1.0 ± 2e-15 here)
exported_grc = Σ_c (scale_c · wtgcell_c · flux_c) / sumwt + dynbal_grc
```

from the dumped columns, with the `urbanf` law (`set_c2l_scale 'urbanf'`,
urbpoi derived from the column type — only urban landunits carry column
types 71..75):

```
scale = 3·canyon_hwr   sunwall (72) / shadewall (73)
scale = 3              road_imperv (74) / road_perv (75)
scale = 1              roof (71) and all non-urban types (1, 5, 215, 401..410)
```

applied to the six dumped column fluxes (qsur→qflx_rofliq_qsur_grc,
qflx_drain→qflx_rofliq_qsub_grc, qflx_drain_perched→
qflx_rofliq_drain_perched_grc, qflx_qrgwl→qflx_rofliq_qgwl_grc,
qflx_sfc_irrig→qirrig_grc, qflx_ice_runoff→qflx_rofice_grc), plus the total:

```
qflx_rofliq_grc = c2g(qflx_runoff_col) + qflx_liq_dynbal_grc
qflx_runoff_col = qsur_col + qdrain_col + qperch_col + qrgwl_col
```

(checked as the c2g of the sum; the lnd2atm `qflx_liq_from_ice_col`
adjustment at :456-473 is already inside the dumped qflx_qrgwl_col). The
dynbal terms are added to the qgwl/total/rofice gridcell values
(lnd2atmMod.F90 :445-452, :465-469) and are exactly 0 in this run; the
validator still adds them to pin the law shape.

Sign/unit conventions: mm/s = kg/m²/s positive downward; qsur/qperch/qirrig/
rofice and the column production fields are ≥ 0, but qflx_drain_col (bottom
drainage/capillary rise, min −1.1e-4 mm/s here) and qflx_qrgwl_col (glacier
columns accumulate; min −2.5 mm/s) — and hence qsub, qgwl and the total —
can be negative (`mosart_in` sets `qgwl_runoff_option = "negative"`).
`Flrl_irrig` is exported negated (`minus = .true.` at
lnd_import_export.F90:916-922); qirrig_grc = 0 everywhere on this cold
start's first day, so the irrig identity pins the law shape but not a
nonzero value.

## The lnd→rof (MOSART) coupling, reproduced

The mediator accumulates the lnd exports between rof coupling steps
(CMEPS `med_phases_prep_rof_accum` from `med_phases_post_lnd`), time-averages
them (`fldbun_average` by `lndAccum2rof_cnt` = 6 exports over the 3 h window
[0, 10800] in this run), and maps them with ESMF `mapconsf`, built **online
at run time** from `mesh_lnd = fv1.9x2.5_141008_ESMFmesh.nc` and `mesh_rof =
r05_nomask_c110308_ESMFmesh.nc` (`lnd2rof_map` unset in `nuopc.runconfig`).
The four nonzero MOSART imports (cells 100312/100313/101032/101033 = a 2×2
block of 0.5° r05 cells at lon 295.5..296.5°E, lat −20.5..−19.5, fully inside
land gridcell 1981) reproduce as

```
import [m³/s] = A_dest · avg6(qsur_src) · W · 1e-3
```

with `A_dest` the spherical cell area (R = 6371.229 km, exact lat/lon quad),
`avg6 = (qsur_step1 + qsur_step2)/6 = 2.590714e-11` mm/s (gridcell 1981's
qsur is nonzero only at steps 1-2 of this cold start — the initial runoff
burst), and `W = A_overlap/A_dest = 1.0` (100312, 101032) / `0.5` (100313,
101033). Worst relative residual **1.35e-6**, which is the ESMF-internal
spherical-polygon intersection / sphere-radius convention — the exact
weights are computed online at run time and never persisted, so they are not
reproducible from the mesh files alone (documented limitation; the validator
checks the reproduction at rel 1e-5).

## Coverage gaps

- **Exact r05 map weights** are ESMF-internal (online `mapconsf`, great-circle
  spherical-polygon intersection, area/radius convention), never persisted;
  the analytic reproduction is structural at ~1e-6. If bit-exact weights are
  wanted later, the run would need the `lnd2rof_map` file (or an ESMF
  RegridWeightGen invocation with the same options) persisted from the case.
- **MOSART dump window** (`mosart_run_imports.csv`, the mosart lane's dump)
  is partial: first 16 local cells of PET 0 × 8 coupling steps, so only the
  4 step-1 nonzero cells are pinnable; and its imports are **exactly zero at
  coupling steps 2..8** although the lnd exports are nonzero all day (~1150
  cells with qsur > 0 at every step) — an open question for the mosart lane
  (dump semantics or coupling state; the step-1 values reproduce exactly).
- **qirrig = 0** all day on this cold start (no irrigation demand on the
  first day), so the −qirrig export negation pins the law shape but no
  nonzero value; the mosart-side irrigation handling is unpinned.
- **dynbal terms = 0** this run (no dynamic land cover change on day 1);
  they are dumped and added in the identities, but unpinned in value.

## Bring-up bugs

Two bugs in this lane's own dump instrumentation (not CLM physics), fixed in
fork commit `dccc7bbb9` (superproject `1664928`):

1. **Final-step export crash** ("Fortran runtime error: Unit number is
   negative"): `dump_finalize` (called at `nlend` in clm_driver before the
   last step's export call) closed the NEWUNIT handles; gfortran NEWUNITs are
   negative and a closed negative unit cannot be reopened in place. Fix:
   `dump_ensure_open` reopens with a fresh NEWUNIT in append mode
   (`status='old', position='append'`) after finalize, so the final step's
   rows land in the same CSVs.
2. **Garbage gridcell index in the cols dump**: the cols dump was passed the
   gridcell-indexed `gindex_grc` where a column-indexed array was needed
   ("\*\*\*\*\*\*\*\*" in the first data row). Fix: lnd2atm builds
   `col_gidx_grc(begc:endc) = gindex_grc(col%gridcell(begc:endc))` as a temp.

Plus one benign CLM quirk logged to `bugs.md`: `lun%canyon_hwr` is never
initialized for non-urban landunits (NaN in the dump; never read by
`set_c2l_scale` for non-urban, so harmless).

## Status

Stage-1 instrumentation complete and verified: `verify_clm_rofliq_dumps.py`
— ALL CHECKS PASSED (structure, c2g identities at round-off, mosart
reproduction, `--existing` 15 CSVs bit-identical, `--restart` neutrality
292/292). The next step (stage 2, a separate PR) is the `.esm` coupling
components: the gridcell export identity feeds the `mosart_coupling.esm`
`RofCap` pin (Flrl_rofsur = 2.59072497089951852e-11 mm/s at r05 cell 100312
over 2.900213644740529e9 m² at coupling step 1) and the land→river flux
export component; the column-level identities pin the `urbanf`/`unity` c2g
aggregation used by every gridcell flux.
