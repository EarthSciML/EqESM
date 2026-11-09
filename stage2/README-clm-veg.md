# CLM vegetated-canopy physics — stage 2 data & validator map

Stage-2 implementation lane for the CTSM/CLM vegetated-canopy physics
(ctsm5.4.063, cesm3.0-alphabranch; instrumented fork commit `afbae49dd` on top
of `65f12f0d8`, branch `instrument/clm-dumps`): the CanopyFluxes vegetated
branch, the TwoStream canopy radiation law, canopy interception/throughfall,
and the soil-water solve. This lane is at the stage-1 → stage-2 boundary: the
instrumented dumps exist, the validator (`stage1/analysis/verify_clm_veg_dumps.py`)
pins the Fortran laws on them, and the `.esm` files will be built against the
same dumps (no `.esm` files yet — see [Status](#status)).

The dumps come from the instrumented rerun of case `I2000Clm60Sp_clm_inst`
(in-node `case.submit --no-batch`), model LID **261108-072738** (run log
`/scratch.local/ctessum/clmveg/run_veg_final.log`; diagnostic LIDs
261108-064127 / 261108-071237 with marker edits, all reverted). The run is
solution-neutral: the instrumented restart is bit-identical to the
uninstrumented baseline `I2000Clm60Sp_f19` restart (292/292 common variables,
0 differ, mask-aware compare), and all 11 pre-existing `clm_*.csv` are
bit-identical (md5) to the pre-veg dumps.

## The vegetated-patch selection

New namelist `dump_veg_patches` (default 0 = off). When > 0, every PET calls
`dump_select_veg` at the top of every step (before the OpenMP clump loops;
serial, reads canopy state only, so it never feeds back into the physics):

1. **Writer PET, chosen once at the first call (step 1)** by two integer
   allreduces over the CLM communicator: each PET counts its qualifying
   patches (`frac_veg_nosno = 1`, `elai+esai > 0`, `liqcan+snocan > 0`),
   reduces the max, and the PET(s) at the max reduce the lowest PET rank.
   On this cold start (no `finidat`) every PET has 0 qualifying patches at
   step 1, so all PETs tie at the max and the tie-break picks **PET 0**.
2. **Per traced step**, the writer PET ranks its qualifying patches by
   `elai+esai` (descending), then `liqcan+snocan` (descending), then patch
   index (ascending), keeps the top `dump_veg_patches` (16 in this run) and
   stores them sorted by patch index together with their distinct columns.

Traced steps are **2..12** (11 steps): step 1 of the cold start has
`liqcan + snocan == 0` on every patch, so nothing qualifies. All four files
share the same per-step selection. Row counts: `clm_canopyflux_veg.csv` 176
(16 patches × 11 steps), `clm_twostream.csv` 92 (the `coszen > 0` daylit
subset), `clm_interception.csv` 176 entry + 176 hyd rows, `clm_soilwater.csv`
318 (entry/exit pairs for 15 distinct selected columns on steps 2..5 and 12,
14 on steps 6..11 — the number of distinct columns of the top-16 patches
varies with the selection).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `stage1/dumps/clm/clm_canopyflux_veg.csv` | CanopyFluxes vegetated-branch exit boundary: one row per selected patch per step with the leaf-temperature-iteration outputs (`num_iter`, `det`, `dele`), the canopy microclimate (`t_veg`, `dt_veg`, `ram1`, `rah1`, `rb1`, `sabv`, `thm`, `z0mv`, `displa`, `cgrnds`, `cgrndl`, `taux`, `tauy`), the fluxes (`eflx_sh_veg`, `qflx_evap_veg`, `qflx_tran_veg`) and the canopy storage **after** the end-of-routine update (liqcan/snocan) | `CanopyFluxesMod.F90`: CanopyFluxes :191-1765 (leaf-temp iteration :1305-1370, 1440-1458; storage update :1618-1635) |
| `stage1/dumps/clm/clm_twostream.csv` | TwoStream canopy radiation law: one row per selected daylit (`coszen > 0`) patch per radiation step with inputs (`coszen`, `elai`, `esai`, `t_veg`, `fwet`, `fcansno`, `rho`/`tau` vis/nir, `albgrd`/`albgri` vis/nir) and outputs (`albd`/`albi` vis/nir, `fabd`/`fabi` vis/nir, `fabd_sun`/`fabd_sha`/`fabi_sun`/`fabi_sha` **vis only**, `ftdd`/`ftid`/`ftii` vis/nir) | `SurfaceAlbedoMod.F90`: TwoStream :1311-1823 (Dickinson 1983 / Sellers 1985 / Bonan 2011) |
| `stage1/dumps/clm/clm_interception.csv` | Canopy interception/throughfall: one row per selected patch per step at stage `ent` (entry: storage before interception) and `hyd` (after `UpdateState_AddInterceptionToCanopy`/`BulkFlux_CanopyExcess`: storage + throughfall + intercepted fluxes + `fwet`/`fcansno`/`fdry`) | `CanopyHydrologyMod.F90`: BulkFlux_CanopyInterceptionAndThroughfall :479-566, UpdateState_AddInterceptionToCanopy, BulkFlux_CanopyExcess, BulkDiag_FracWet :1220-1273 |
| `stage1/dumps/clm/clm_soilwater.csv` | Soil-water solve boundary: one row per selected **column** per step at stage `entry` (state before the solve) and `exit` (after), with the 10-layer `h2osoi_liq` and the aquifer `wa`. The hk-law columns (`vwc_liq`..`dsmpdw`) are NaN in this run — see [Coverage gaps](#coverage-gaps) | `SoilWaterMovementMod.F90`: SoilWater :1347-1506 (zengdecker2009 hk law :1412-1458, tridiagonal solve :1460-1506) |
| `stage1/analysis/verify_clm_veg_dumps.py` | Validator: pins the Fortran laws on the four CSVs, checks the selection consistency, verifies the 11 pre-existing dumps bit-identical (md5) against a backup, and the restart neutrality (mask-aware compare of the instrumented vs baseline `clm2.r` restart) | — |

## The traced laws (validated)

All checks pass to round-off / tiny epsilon (es26.17 dump precision). The
validator recomputes the laws from the dumped inputs and compares against the
dumped outputs; the laws below are the ones it encodes.

**Interception** (stage `hyd` rows; `use_clm5_fpi = .true.` in this run):

```
fpiliq  = interception_fraction * tanh(elai + esai)     ! interception_fraction = 1.0
fpisnow = 1 - exp(-0.5*(elai + esai))
qflx_through_liq  = (forc_rain + irrig)*(1 - fpiliq), etc.
liqcan(hyd) = min(liqcanmx, max(0, liqcan(ent) + dtime*qflx_intercepted_liq))
              with liqcanmx = liq_canopy_storage_scalar*(elai+esai)  ! scalar 0.1
snocan(hyd) in [0, max(0, snocan(ent) + dtime*qflx_intercepted_snow)]  (snow
              unloading removes water but forc_t/forc_wind are not dumped)
fwet    = min((h2ocan/((elai+esai)*liq_canopy_storage_scalar))^(2/3),
              maximum_leaf_wetted_fraction)            ! 0.05, h2ocan = liqcan+snocan > 0
fcansno = min((snocan/((elai+esai)*snow_canopy_storage_scalar))^0.15, 1)  ! scalar 6.0
fdry    = (1 - fwet)*elai/(elai + esai)
```

**CanopyFluxes exit** — the leaf-temperature iteration converges to
`det < dtmin = 0.01 K`, `dele < dlemin = 0.1 W/m²` within `num_iter ≤
itmax_canopy_fluxes + 1 = 41` (the convergence check runs at `itlef = 3..41`;
`num_iter = 41` rows may be at the iteration cap, so `det`/`dele` are not
checked there). `t_veg ∈ (150, 400) K`, `qflx_tran_veg ≥ 0` (`qflx_evap_veg`
may be below `qflx_tran_veg` when canopy condensation occurs). The step order
within a model step is **CanopyHydrology → CanopyFluxes → SurfaceAlbedo**
(empirically verified), so the `canopyflux_veg` rows carry the **post-evap**
storage from the end-of-CanopyFluxes update (:1618-1635):

```
if t_veg > tfrz (273.15):
  if (evap - tran)*dtime > liqcan(hyd):     ! all liquid evaporated
     snocan(cf) = max(0, snocan(hyd) + liqcan(hyd) + (tran - evap)*dtime)
  liqcan(cf) = max(0, liqcan(hyd) + (tran - evap)*dtime)
else: (mirror with the roles of liqcan/snocan swapped)
```

verified against the interception `hyd` rows of the same step/patch (e.g.
p = 135, step 2: `0.48338 = 0.47969 + (0 − (−2.05e-6))·1800`), with
`truncate_small_values` possibly zeroing sub-epsilon values afterwards.

**TwoStream** — the full Dickinson/Sellers/Bonan two-stream computation is
recomputed per row from the dumped inputs plus `xl(pft)` read from
`ctsm60_params.c260518.nc` (`chil = clamp(xl, −0.4, 0.6)`, |chil| ≤ 0.01 →
0.01; `phi1`/`phi2`/`gdir`/`twostext`/`avmu`, `asu`/`betadl`/`betail`, the
`snowveg_affects_radiation = .true.` fcansno blend of `omega = (1−fcansno)·
omegal + fcansno·0.8` (vis) / `0.4` (nir) with `betads = betais = 0.5`, the
`h1..h10`/`p1..p4`/`u1..u3`/`tmp2..tmp9` chain, and the output formulas for
`albd`/`albi`/`fabd`/`fabi`/`fabd_sun`/`fabd_sha`/`fabi_sun`/`fabi_sha`/
`ftdd`/`ftid`/`ftii`). All match to round-off. The `cosz = max(0.001, coszen)`
floor applies inside TwoStream; the sun/shade columns exist for vis only in
the dump.

**Soil water** — with `soilwater_movement_method = 1` (moisture_form) the
zengdecker2009 hk loop never runs, so only structure is pinned: entry/exit
row pairing per (step, column), exit `h2osoi_liq ≥ 0`, `wa` unchanged by the
solve, and the hk-law columns NaN.

## Run configuration

`user_nl_clm` deltas on top of the W1-W3 case config: `dump_veg_patches = 16`,
`itmax_canopy_fluxes = 40`. Kept from the existing case: `use_clm5_fpi =
.true.` (fpiliq = `interception_fraction·tanh`), `snowveg_affects_radiation =
.true.` (TwoStream fcansno omega blend), `paramfile = ctsm60_params.c260518.nc`,
`dump_first_steps = 12`, `dtime = 1800 s`, 36 MPI tasks. Params read from the
paramfile by the validator: `interception_fraction = 1.0`,
`liq_canopy_storage_scalar = 0.1`, `snow_canopy_storage_scalar = 6.0`,
`maximum_leaf_wetted_fraction = 0.05`, `xl` indexed by pft 0..78.

## Coverage gaps

- **`soilwater_movement_method = 1` (moisture_form)** — the zengdecker2009
  hk law (`SoilWaterMovementMod.F90:1412-1458`) and its tridiagonal solve
  never run in this case, so `clm_soilwater.csv` has entry/exit rows only and
  the hk-law columns are NaN. A moisture_form solve trace could be added
  later without touching this lane's files. (The `layer, h2osoi_liq: 182 2
  -20.878` debug print at :1422 and the step-1 balance warnings are
  pre-existing baseline behavior, deterministic and identical in both runs —
  see `bugs.md` for the underlying negative-h2osoi_liq candidate entry.)
- **Interception** — snow unloading is not pinned exactly (needs
  `forc_t`/`forc_wind`, not dumped); `fpisnow` and the `fcansno` law are
  pinned. The `irrig` term never binds (no irrigation in this run).
- **TwoStream** — sun/shade absorption columns exist for vis only; the
  `lSFonly = .true.` branch and `use_SSRE` are not exercised.
- **CanopyFluxes** — `num_iter = 41` rows are at the iteration cap
  (`itmax_canopy_fluxes = 40`); the vegetation below/at `tfrz` storage-update
  mirror branch is exercised only if a traced patch is frozen (check the
  data; `t_veg` is above tfrz on most rows).
- **Selection** — step 1 of a cold start traces nothing (no qualifying
  patches); the traced steps are 2..12. The writer-PET choice is
  deterministic (PET 0 here) but data-dependent, so the `.esm` tests must
  read the (step, patch) keys from the CSVs rather than assuming which PET
  or which patches are traced.

## Validator usage

```
python3 stage1/analysis/verify_clm_veg_dumps.py \
    --existing <backup-dir-with-11-md5s> \
    --restart <instrumented clm2.r> <baseline clm2.r>
```

- `--existing`: md5-compare the 11 pre-existing `clm_*.csv` in
  `stage1/dumps/clm/` against a backup (bit-identical neutrality of the
  re-dumps).
- `--restart`: mask-aware bit-identity of the instrumented restart vs the
  uninstrumented baseline (the naive compare's few "diffs" are masked/NaN
  artifacts).
- Without flags: pins the laws on the four veg CSVs (selection consistency,
  interception, canopyflux incl. the storage-update identity, full TwoStream
  recomputation, soilwater structure).

Result of the final run: **ALL CHECKS PASSED**.

## Status

Stage-1 complete for this lane (instrumentation, runs, neutrality proof,
validator). Stage-2 `.esm` implementation of the traced laws (CanopyFluxes
vegetated branch, TwoStream, interception, soil water) is planned against
these CSVs; the validator's recomputations are the reference for the
`.esm` `tests` blocks. No `.esm` files exist for this lane yet.
