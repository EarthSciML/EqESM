# CLM stage 2 — vegetated canopy hydrology & storage update

Stage-2 slice for the CTSM/CLM vegetated-patch canopy physics
(ctsm5.4.063, cesm3.0-alphabranch; dump merge `eee5c6e` of the
`instrument/clm-dumps` fork): the canopy interception/throughfall hydrology
of `CanopyHydrologyMod.F90` and the end-of-CanopyFluxes canopy-storage
update of `CanopyFluxesMod.F90:1615-1643`, pinned on the vegetated-patch
dumps `stage1/dumps/clm/clm_interception.csv` (176 ent + 176 hyd rows) and
`stage1/dumps/clm/clm_canopyflux_veg.csv` (176 rows; the pre-flux interface
and the rest of `CanopyHydrologyMod` were pinned on the W2 dumps in
`stage2/clm_canopy.esm` — see `README-clm.md`/`README-clm-veg.md` for the
stage-1 lane map). The leaf-temperature iteration of CanopyFluxes
(`num_iter`/`det`/`dele`) is a Newton loop and is **out of scope** — see
[Coverage gaps](#coverage-gaps).

Sources: `stage1/cesm/components/clm/src/biogeophys/CanopyHydrologyMod.F90`
(BulkFlux_CanopyInterceptionAndThroughfall :479-566,
UpdateState_AddInterceptionToCanopy :770-787, BulkFlux_CanopyExcess
:793-838, BulkDiag_FracWet :1220-1273),
`CanopyFluxesMod.F90:1614-1643` (dew-accumulation update and the
`truncate_small_values` call), `utils/NumericsMod.F90:50-99`
(truncate_small_values), the parameter file `ctsm60_params.c260518.nc`
(`interception_fraction = 1.0`, `liq_canopy_storage_scalar = 0.1`,
`snow_canopy_storage_scalar = 6.0`, `maximum_leaf_wetted_fraction = 0.05`),
and the traced case's `user_nl_clm` (`use_clm5_fpi = .true.`,
`itmax_canopy_fluxes = 40`, dtime = 1800 s).

## Files

| File | Contents |
|---|---|
| `stage2/clm_canopy.esm` | Model `CanopyHydrology` (existing): 5 new tests pinning the vegetated interception rows of `clm_interception.csv` (see below); a new `interception_fraction` parameter (ctsm60_params.c260518.nc) feeding the corrected `canopy_intercept_liq` template. Model `CanopyFluxInterface` (existing): untouched. |
| `stage2/clm_canopy_storage.esm` | **New.** Model `CanopyStorageUpdate`: the end-of-CanopyFluxes canopy-storage update (`CanopyFluxesMod.F90:1615-1643`) as a pure function of the post-hydrology storage and the CanopyFluxes outputs — above/below freezing branches, the all-liquid/all-snow-evaporated transfer, the `max(0,.)` clamps, and the `truncate_small_values` pass (`abs(snocan) < 1e-10*abs(snocan_hyd)` → 0, with the intermediate `snocan_raw` to keep the truncation a one-way equation). 6 tests pin the canopyflux_veg rows against the same-step interception hyd rows. |
| `stage2/clm_templates.esm` | **Corrected in place:** `canopy_intercept_liq` was the `use_clm5_fpi = .false.` branch `0.25*(1 - exp(-0.5*(elai+esai)))`; the vegetated dumps show this run's law is the `use_clm5_fpi = .true.` branch `interception_fraction*tanh(elai+esai)` (CanopyHydrologyMod.F90:538) — `qflx_intercepted_liq = forc_rain*tanh(elai+esai)` on all 101 intercepting rows, ruling out the old form (which is 4.4× smaller at the traced L ≈ 4.5-6.2). The template now takes `interception_fraction` and uses the tanh body; the description records the correction. All other canopy templates (`canopy_intercept_snow`, `canopy_storage_max`, `canopy_excess_flux`) were verified against the Fortran and the dumps and are unchanged. |
| `stage1/analysis/clm_canopy_veg_test_tuples.py` | **New.** Tuple extractor/verifier for the two files' `tests` blocks: for a fixed list of (step, patch) keys it recomputes the interception and storage-update laws from the dumped inputs (mirroring `verify_clm_veg_dumps.py`) and prints the es26.17 tuples that the hand-authored tests pin. Run: `python3 stage1/analysis/clm_canopy_veg_test_tuples.py`. |

## The pinned laws

All pins are at `rel = 1e-12` (es26.17 dump precision) and every assertion
can fail. The dumps carry no `fpiliq`/`fpisnow` columns (the interception
CSV writes the fluxes), so the interception coefficients are pinned through
`qflx_intercepted_liq = (forc_rain+irrig)*fpiliq` and
`qflx_through_liq = (forc_rain+irrig)*(1-fpiliq)`.

**Interception (`clm_interception.csv` hyd rows; `CanopyHydrologyMod.F90`):**

| Test | Row | Pins |
|---|---|---|
| `clm_canopy_interception_veg_step2_p135` | step 2, p 135 | the tanh interception law (fpiliq derived + qflx_intercepted_liq/qflx_through_liq), the fwet = 0.05 clamp, fdry, h2ocan/vegt, liqcanmx/snocanmx (derived) — this row's hyd liqcan saturates at liqcanmx (the discrete min-clamp of AddInterceptionToCanopy/CanopyExcess), so the derived liqcanmx equals the dumped hyd liqcan to round-off |
| `clm_canopy_interception_veg_step7_p714` | step 7, p 714 | interception fluxes at light rain; **unsaturated fwet** (raw power law, no clamp) and fdry; derived liqcanmx |
| `clm_canopy_interception_veg_step2_p250` | step 2, p 250 | fwet just below its clamp (4.3958e-2, the largest unsaturated value of the 49 unsaturated hyd rows), fdry, h2ocan/vegt; no precipitation this step (check_point = 0) |
| `clm_canopy_interception_veg_ode_step3_p223` | step 3, p 223 | **one-step ODE trace** (time_span {0, 1800}): from the ent-row storage, d(liqcan)/dt = qflx_intercepted_liq (the step is unsaturated so the excess-fall flux is identically 0 and the RHS is constant — any consistent integrator integrates it exactly) reproduces the hyd-row storage to round-off, matching the Fortran's discrete update `liqcan = max(0, liqcan + dtime*intercepted)`; plus the flux split and unsaturated fwet/fdry at the hyd state |
| `clm_canopy_interception_veg_step12_p264` | step 12, p 264 | the fwet power law at a 1e-7-mm argument (4.5461e-5), fdry, h2ocan |

**Storage update (`clm_canopy_storage.esm`; `CanopyFluxesMod.F90:1615-1643`):**

| Test | Row | Pin |
|---|---|---|
| `clm_canopy_storage_update_step2_p135` | step 2, p 135 | the validator's canonical example: condensation (evap = −2.05e-6 mm/s, tran = 0) adds `(tran−evap)·dtime = +3.69e-3` mm to the hyd liqcan 4.7969e-1 → 4.8338e-1 as dumped; the all-evaporated branch does not bind, so snocan stays at snocan_hyd (a snocan assertion that fails if the branch logic is wrong) |
| `clm_canopy_storage_update_step2_p198` | step 2, p 198 | evaporation drawdown (net = −1.12e-1 mm), no clamp, snocan unchanged |
| `clm_canopy_storage_update_step2_p50` | step 2, p 50 | the all-liquid-evaporated boundary: `(evap−tran)·dtime` exceeds liqcan_hyd by 1 ulp, `liqcan_hyd + (tran−evap)·dtime = −3.47e-18`, and the `max(0,.)` clamp at :1623 gives 0.0 exactly as dumped — pins the clamp and the branch condition to the last bit of f64 arithmetic |
| `clm_canopy_storage_update_step7_p714` | step 7, p 714 | second independent razor-edge pin of the same clamp (margin −9.1e-18) |
| `clm_canopy_storage_update_step3_p155` | step 3, p 155 | drawdown of a saturated canopy (hyd liqcan = liqcanmx), mid-value result |
| `clm_canopy_storage_update_step12_p264` | step 12, p 264 | smallest-storage pin: liqcan 1.450e-7 → 1.0614e-7 mm (sub-relative-epsilon move reproduced to round-off) |

## Coverage gaps

- **Snow on the canopy: unpinned.** All 176 interception rows have
  `forc_snow = 0` and `snocan = 0` (the traced selection never sees
  snowfall), so `fpisnow` (1 − exp(−0.5·(elai+esai))), the snow
  throughfall/interception fluxes, the snow-unloading fluxes
  (`qflx_snotempunload`/`qflx_snowindunload`/`qflx_snow_unload`,
  :917-925), `fcansno` and the `truncate_small_values` pass are all
  identically zero on every row and cannot be pinned (the laws are in
  `clm_canopy.esm`/`clm_canopy_storage.esm` and match the Fortran, but
  nothing fails if they are wrong). `forc_t`/`forc_wind` are not in the
  veg dumps either.
- **The all-evaporated transfer into snocan is unobservable:** all rows
  have `snocan_hyd = 0`, so the `:1622` transfer `max(0, snocan_hyd +
  liqcan_hyd + net)` and the untaken-branch `snocan_hyd` both give 0; only
  the liqcan side of the boundary rows (which do take the branch) is
  pinned. The **below-freezing mirror branch** (`:1628-1635`) is not
  exercised at all — all 176 rows have `t_veg > 273.15` (min 288.68 K);
  `clm_canopy_storage.esm` transcribes it as written (including the
  unclamped `:1630` transfer — see `bugs.md` entry TBD).
- **The discrete one-step storage clamp at saturation** (`min(liqcanmx,
  ·)` in the Fortran's update chain) is not pinable through the ODE form
  (the continuous limit approaches liqcanmx asymptotically); the saturated
  rows pin `liqcanmx = 0.1*(elai+esai)` instead, via equality of the
  dumped hyd liqcan with the derived liqcanmx (48 of 176 rows saturate).
- **Ground fluxes:** the dumped `qflx_liq_grnd`/`qflx_snow_grnd` are
  column aggregates (p2c over all patches of the column), so the
  patch-level ground flux law `throughfall + canfall + unload +
  irrig_drip` (`SumFlux_FluxesOntoGround :1063-1067`) cannot be pinned
  exactly; the validator checks only a lower bound. Irrigation is zero in
  this run.
- **Leaf-temperature iteration (out of scope):** `num_iter`/`det`/`dele`,
  `t_veg` (beyond its role as an input to the storage update),
  `ram1`/`rah1`/`rb1`, `sabv`, `cgrnds`/`cgrndl`, `taux`/`tauy`,
  `eflx_sh_veg`/`qflx_evap_veg`/`qflx_tran_veg` production: the
  `CanopyFluxesMod.F90:1305-1370,1440-1458` leaf-temperature Newton loop
  (with `det < dtmin = 0.01 K`, `dele < dlemin = 0.1 W/m²` convergence,
  `num_iter ≤ itmax+1 = 41`; 13 rows exit at the iteration cap) is a
  candidate for a later stage-2 slice; the validator's structural checks
  on it are in `verify_clm_veg_dumps.py::check_canopyflux`.
- **Selection:** step 1 of the cold start traces nothing (no qualifying
  patches); traced steps are 2..12, 16 patches/step.

## Verification

```
# tuple extractor + law verification for the pinned rows
python3 stage1/analysis/clm_canopy_veg_test_tuples.py

export ESD_ROOT=.../EarthSciDiscretizations
export ESM_CACHE_DIR=/scratch.local/ctessum/earthsci-esm-cache
./esm test stage2/
```

Suite counts on this branch: 3073P/0F/0E at the starting commit →
**3107P/0F/0E** (clm_canopy.esm 103 assertions, clm_canopy_storage.esm 8,
both 0F/0E). The stage-1 lane validator (`verify_clm_veg_dumps.py`) is
unchanged and still passes ("ALL CHECKS PASSED").
