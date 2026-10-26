# MOSART stage 2 — `.esm` skeleton

Stage-2 implementation of the CESM MOSART river-routing component (mosart1.1.13,
cesm3.0-alphabranch checkout `e59f2512`) as esm 1.3.0 files: component structure,
imports/exports, solver block, the pointwise flow laws of `mosart_physics.F90`
transcribed exactly, subcomponent routing steps, and top-level assembly. The
`tests` blocks are intentionally empty pending the stage-1 instrumentation dumps
(`stage1/dumps/mosart/`) — see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/mosart/src/riverroute/mosart_physics.F90`,
`stage1/cesm/components/mosart/src/riverroute/mosart_driver.F90`,
`stage1/cesm/components/mosart/src/riverroute/mosart_tstatusflux_type.F90`,
`stage1/cesm/components/mosart/src/riverroute/mosart_tspatialunit_type.F90`,
`stage1/cesm/components/mosart/src/cpl/nuopc/rof_import_export.F90`,
and the boundary/subassembly plan in `stage1/boundaries.md` (§1.5, §2.6).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `mosart_templates.esm` | Expression-template library: `manning_velocity`, `overland_discharge`, `wetted_area`, `rectangular_depth`, `rectangular_perimeter`, `hydraulic_radius`, `channel_depth`, `channel_perimeter` | `mosart_physics.F90`: CRVRMAN_nosqrt 406-422, CREHT_nosqrt 426-437, GRMR 441-449, GRHT 453-465, GRPT 469-481, GRRR 485-497, GRHR 501-532, GRPR 536-572 |
| `mosart_hillslope.esm` | Model `Hillslope` — Euler step 1: hillslopeRouting + updateState_hillslope, `etin` lateral inflow | `mosart_physics.F90`: Euler 74-80, hillslopeRouting 225-239, updateState_hillslope 352-360 |
| `mosart_subnetwork.esm` | Model `Subnetwork` — Euler step 2: subnetworkRouting + updateState_subnetwork, `erlateral` | `mosart_physics.F90`: Euler 114-132, subnetworkRouting 243-269, updateState_subnetwork 364-381 |
| `mosart_channel.esm` | Model `Channel` — Euler step 4: Routing_KW + updateState_mainchannel, `flow` | `mosart_physics.F90`: Euler 175-209, Routing_KW 273-348, updateState_mainchannel 385-402 |
| `mosart_coupling.esm` | Model `RofCap` — med↔rof cap unit conversions only | `rof_import_export.F90`: state_getimport 589-595, state_setexport 636-643, fields 88-117, areacor 209-210, tdepth_max 401 |
| `mosart_routing.esm` | Top-level assembly: mounts RofCap/Hillslope/Subnetwork/Channel by reference + inline `Control` (mosart_run steps 2-8, subtimestep logic, aggregation), 29 `variable_map` coupling entries, enums, solver block | `mosart_driver.F90`: mosart_run 375-873; `mosart_tspatialunit_type.F90`: set_subtimesteps 624-684 |

The routing steps are written as storage ODEs (`D(wh)/dt`, `D(wt)/dt`, `D(wr)/dt`)
with the pointwise laws as observed variables via `apply_expression_template` —
the continuous limit of the Fortran's forward-Euler cascade; at stage 3 the ODEs
are integrated with `github.com/earthsciml/earthsciintegrator`. Dry-cell clamps
and the Fortran if/else branches are transcribed as `ifelse` expressions.

## Units

Per `mosart_tstatusflux_type.F90` and `mosart_tspatialunit_type.F90` geometry:

| Quantity | Symbol | Units |
|---|---|---|
| Hillslope storage / depth | `wh`, `yh` | m |
| Subnetwork storage | `wt` | m³ |
| Main-channel storage | `wr` | m³ |
| Hillslope fluxes | `qsur`, `qsub`, `ehout`, `dwh` | m/s |
| Subnetwork / channel fluxes | `etin`, `etout`, `erlateral`, `erin`, `erout`, `dwr`, `eroutUp`, `flow`, `qirrig_actual`, `flood` | m³/s |
| Velocities | `vt`, `vr` | m/s |
| Cross-section areas | `mt`, `mr` | m² |
| Depths / perimeters / radii | `yt`, `pt`, `rt`, `yr`, `pr`, `rr` | m |
| Slopes (sqrt) | `hslpsqrt`, `tslpsqrt`, `rslpsqrt` | 1 |
| Drainage density | `Gxr` | 1/m |
| Lengths / widths | `tlen`, `hlen`, `twidth`, `rlen`, `rwidth`, `rwidth0`, `rdepth` | m |
| Areas | `area`, `areaTotal2` | m² |
| Manning roughness | `nh`, `nt`, `nr` | s/m^(1/3) |
| Times | `delt`, `localdelt`, `coupling_period`, `delt_coupling`, `delt_mosart` | s |
| Cap fields (CMEPS) | `Flrl_*`, `Fgrg_*`, `Forr_*`, `Flrr_flood` | mm/s |
| Cap depth fields | `Flrr_volr`, `Flrr_volrmch`, `Sr_tdepth`, `Sr_tdepth_max` | m |

**Unit notes**

- Manning's n is declared with its SI dimension `s/m^(1/3)` (the frivinp network
  file stores bare numbers). This keeps the equations dimensionally consistent:
  `(rr²)^(1/3)·sqrt(slope)/n` has units `m^(2/3)/(s/m^(1/3)) = m/s`. The esm 1.3.0
  unit checker treats a mismatch as a hard error, so this declaration is required.
- Cap fields are imported in `mm/s` (= CMEPS kg/m²/s for water) and converted
  with `mm_to_m` [m/mm] = 0.001 exactly as `state_getimport`
  (`q [m³/s] = field [mm/s] · med2mod_areacor · area · 0.001`).
- `Flrr_volr = volr/area` and `Flrr_volrmch = wr/area` are depths [m] with **no**
  mm scaling, per `state_setexport` (F90:396-397).
- `TINYVALUE = 1.0e-14` and `SLOPE1 = 0.1` are inlined as bare literals in the
  templates (the unit checker skips comparisons against undeterminable operands).

## §2.6 mapping, as implemented

| Subassembly (boundaries.md §2.6) | Proposed | Implemented |
|---|---|---|
| Cap boundary (kg/m²s↔m³/s) | `mosart_coupling.esm` | ✓ `mosart_coupling.esm` (RofCap) — pure unit conversions |
| Direct transfer (mosart_run steps 2-4) | `mosart_coupling.esm` | ⚠ moved into inline `Control` in `mosart_routing.esm`: the qgwl threshold option and the flood rule need `Channel.wr`, so the orchestration lives with the systems (annotated deviation from the §2.6 table) |
| Irrigation/flood (steps 2-3) | `mosart_routing.esm` | ✓ `Control` (rate rules) |
| Routing hillslope+subnetwork+channel | `mosart_routing.esm` | ✓ mounted subassemblies, coupled in Euler order |
| Flow laws / geometry (D, first test targets) | `mosart_routing.esm` templates | ✓ factored into `mosart_templates.esm`, imported by reference |
| Network topology (Init, set_routehandles, set_subtimesteps) | data source + discretization | ◐ set_subtimesteps (CFL phi/numDT/localdelt) is in `Control`; set_routehandles / frivinp network data source deferred (stage 2.5) |
| Lakes/reservoirs/groundwater | TBD | — (not present in this MOSART; §2.5) |
| State/budget | harness | — (stage 3, with `mosart_budget_type.F90` as a conservation test) |

## Pointwise laws — Fortran vs esm transcription notes

All templates are match-less, invoked by name, and transcribed exactly. Fortran
formula → esm formula:

| Template | Fortran (mosart_physics.F90) | esm |
|---|---|---|
| `manning_velocity` | CRVRMAN_nosqrt: `v = ((rr²)^(1/3))·sqrtslp/n`, zeroed when `rr <= 0` | `ifelse(rr <= 0, 0, ((rr²)^(1/3))·sqrtslp/n)` |
| `overland_discharge` | CREHT_nosqrt: `eht = Gxr·yh·CRVRMAN_nosqrt(sqrthslp, nh, yh)` | same, calling `manning_velocity` by reference (the one shared-law factoring) |
| `wetted_area` | GRMR: `m = w/len` | `w/len` |
| `rectangular_depth` | GRHT: `h = m/twid`, zeroed when `m <= TINYVALUE` | `ifelse(m <= 1e-14, 0, m/twid)` |
| `rectangular_perimeter` | GRPT: `p = twid + 2·h`, zeroed when `h <= TINYVALUE` | `ifelse(h <= 1e-14, 0, twid + 2·h)` |
| `hydraulic_radius` | GRRR: `r = m/p`, zeroed when `p <= TINYVALUE` | `ifelse(p <= 1e-14, 0, m/p)` |
| `channel_depth` | GRHR three-part cross section (bankfull rectangle → trapezoid of slope SLOPE1=0.1 → flood plain): `h=0` if `m<=TINYVALUE`; `m/rwidth` if `m - rdepth·rwidth <= TINYVALUE`; trapezoid `rdepth + (-rwidth + sqrt(rwidth² + 4·deltam/SLOPE1))·SLOPE1/2`, `deltam = m - rdepth·rwidth`; flood plain `rdepth + SLOPE1·((rwidth0-rwidth)/2) + deltam/rwidth0`, `deltam = m - rdepth·rwidth - (rwidth+rwidth0)·SLOPE1·((rwidth0-rwidth)/2)/2` | nested `ifelse` with `SLOPE1 = 0.1` inlined; trapezoid threshold test `m > rdepth·rwidth + (rwidth+rwidth0)·0.1·((rwidth0-rwidth)/2)/2 + 1e-14` as in F90:517-518 |
| `channel_perimeter` | GRPR: `p=0` if `h<TINYVALUE`; bankfull `rwidth + 2·h` if `h <= rdepth+TINYVALUE`; flood plain `rwidth + 2·(rdepth + ((rwidth0-rwidth)/2)·SLOPE1·sinatan + deltah)`, `deltah = h - rdepth - ((rwidth0-rwidth)/2)·SLOPE1`; trapezoid `rwidth + 2·(rdepth + (h-rdepth)·sinatan)`; `sinatan = sinatanSLOPE1defr = 1/sin(atan(SLOPE1def))` (F90:554) | same structure; `sinatan` spelled as the identical expression `1/sin(atan(0.1))` |

Substep / routing transcriptions:

- **Hillslope** (Euler step 1): `ehout_manning = -CREHT_nosqrt(...)`;
  dry-cell clamp fires when `ehout_manning < 0` **and**
  `wh + (qsur + ehout_manning)·delt < TINYVALUE`, replacing `ehout` by
  `-(qsur + wh/delt)`; `dwh = qsur + ehout`; `yh = wh`;
  `etin = (-ehout + qsub)·area·frac`.
- **Subnetwork** (Euler step 2): no tributaries (`tlen <= hlen`) → `etout = -etin`
  (pass-through); otherwise `vt_manning = manning_velocity(tslpsqrt, nt, rt)`,
  `etout_manning = -vt_manning·mt`, dry-cell clamp `etout = -(etin + wt/localdelt)`
  when `wt + (etin + etout_manning)·localdelt < TINYVALUE`; `vt = -etout/mt` when
  the clamp fired and `mt > 0` (F90:257-259); `dwt = etin + etout`;
  `erlateral = -etout`. Geometry (`mt, yt, pt, rt`) from GRMR/GRHT/GRPT/GRRR when
  `tlen > 0` and `wt > 0`, else all zero (updateState_subnetwork).
- **Channel** (Euler step 4, RoutingMethod=1 only, F90:280-284): `erin = -eroutUp`;
  `rlen <= 0` → `vr = 0`, `erout = -(erin + erlateral)`;
  `areaTotal2/(rwidth·rlen) > 1e6` → `erout = -(erin + erlateral)` (F90:310-312);
  else `vr_manning = manning_velocity(rslpsqrt, nr, rr)`, `erout_manning = -vr_manning·mr`,
  dry-cell clamp fires when `-erout_manning > TINYVALUE` **and**
  `wr + (erlateral + erin + erout_manning)·localdelt < TINYVALUE`, then
  `erout = -(erlateral + erin + wr/localdelt)` and `vr = -erout/mr` when `mr > 0`;
  `dwr = erlateral + erin + erout + qgwl·area·frac`; `flow = -erout`. Geometry
  (`mr, yr, pr, rr`) from GRMR/GRHR/GRPR/GRRR when `rlen > 0` and `wr > 0`, else zero.
- **Control** (mosart_run + set_subtimesteps): `nsub = ceil(coupling_period/delt_mosart)`
  (F90:786-789), `delt = delt_coupling/nsub`; CFL diagnostics
  `phi_r = areaTotal2·sqrt(rslp)/(rlen·rwidth)`, `phi_t = area·sqrt(tslp)/(tlen·twidth)`;
  `numDT = max(ifelse(phi >= 10, DLevelR·log10(phi) + 1, DLevelR + 1), 1)`;
  `localdelt_t = delt/DLevelH2R/numDT_t`, `localdelt_r = delt/DLevelH2R/numDT_r`.
  Irrigation: `irrig_volume = -qirrig_m3s·coupling_period`,
  `qirrig_actual = -min(irrig_volume, wr)/coupling_period`, overage returned to
  subsurface: `qsub_i = qsub_m3s + ifelse(irrig_volume > wr, (wr - irrig_volume)/coupling_period, 0)`.
  Flood (mask==1 only): `volr = wt + wr + wh·area`, `flood = ifelse(volr > fthresh, (volr - fthresh)/delt_coupling, 0)`.
  Direct transfer: negative `qsur_f`/`qsub_i` and mask==0 cells to `direct_liq`;
  `qgwl` under `qgwl_runoff_option` (all/negative/threshold: `qgwl_m3s·coupling_period + wr < river_depth_minimum·area`),
  active only for `bypass_routing_option` direct_in_place/direct_to_outlet; frozen
  runoff `direct_ice = qsur_ice_m3s`; glc streams under `rof_from_glc`.
  Aggregation: `runoffocn_liq = flow` for mask>=2; `runofftot_liq = direct_liq + runoffocn_liq + direct_glc_liq`;
  `runofftot_ice = direct_ice + direct_glc_ice`; the ice_runoff branch of
  export_fields (rof_import_export.F90:362-382): with `ice_runoff` on,
  `rofl_m3s = runofftot_liq - direct_glc_liq`, `rofi_m3s = direct_ice`;
  with it off, `rofl_m3s = runofftot_liq + direct_ice - direct_glc_liq`, `rofi_m3s = 0`.

## Documented divergences (all annotated in file metadata)

1. **Subnetwork `vt` when no tributaries** (`tlen <= hlen`): the Fortran leaves
   `vt` stale (unused in the pass-through branch); here `vt` always takes the
   Manning value (the `ifelse` on the clamp condition handles the other branch).
   Values agree whenever `vt` is used downstream.
2. **Channel `vr` in the `areaTotal2/(rwidth·rlen) > 1e6` branch**: the Fortran
   leaves `vr` at the Manning value (unused); here the same is expressed (the
   `vr` ifelse only zeroes for `rlen <= 0` and overrides on the dry-cell clamp).
3. **`volr` uses the running `wr`** (the F90:854 sum `wt + wr + wh·area` is
   computed *after* the irrigation/flood rules in mosart_run; the esm Control
   computes `volr` from the coupled `wr` — consistent because the flood rule is
   a rate rule in the same algebraic system).
4. **Irrigation and flood as rate rules**: the Fortran removes the irrigation
   volume and flood volume from storage once per coupling period (storage
   jumps); the esm writes them as continuous rate rules
   (`qirrig_actual = -min(irrig_volume, wr)/coupling_period`, `flood = (volr - fthresh)/delt_coupling`).
   Same coupling-period totals; the intra-period trajectory differs. Stage-3
   integration may replace these with discrete events.
5. **Negative-storage warning** (`write_msg` in updateState_*, F90:361, 384)
   not transcribed (diagnostic only).
6. **`eroutUp` upstream gather and the `rh_direct` outlet gather** are plain
   parameters (`Channel.eroutUp`, `Control.direct_*`) — placeholders until the
   frivinp network data source lands (stage 2.5, boundaries.md §2.4); the
   ID/dnID ESMF SMM topology is not expressible without it.
7. **Fortran logicals** (`rof_from_glc`, `flds_r2l_stream_channel_depths`,
   `ice_runoff`) are 0/1 parameters compared with `== 1.0`.
8. **`D(wh)/dt` continuous limit**: the Fortran advances storage with forward
   Euler at `delt` (hillslope) / `localdelt` (subnetwork, channel) with
   DLevelH2R/numDT averaging; the esm ODEs integrate continuously with the
   subtimestep counts entering only through `localdelt`. Instantaneous-vs-average
   differences are expected and will be characterized against the dumps.

## Bugs found in the C++/Fortran code

Centralized list (per AGENTS.md): see `bugs.md` in the repo root; nothing new
found during this transcription. Notable observations worth a `bugs.md` entry:
`dwt`/`dwr` type comments say [m3] but the code assigns flux sums and advances
storage by `·localdelt` (rate reading adopted, noted in variable descriptions);
`qgwl` is added to `dwr` without an `ifelse` when `bypass_routing_option = none`
(F90:327 is unconditional).

## Testing

- `tests` blocks are intentionally empty (`"tests": []`) in every model until
  the stage-1 instrumentation dumps (`stage1/dumps/mosart/`) provide numeric
  tuples (TODO in each file's metadata description).
- Validation: `esm validate <file>` passes for all six files (esm 1.3.0,
  EarthSciAST CLI v0.4.0); `esm test <file>` parses cleanly with 0 assertions.
- First test targets, per boundaries.md §4: the pointwise laws (D boundaries)
  `CRVRMAN_nosqrt`, `CREHT_nosqrt`, `GRMR/GRHR/GRPR/GRRR/GRHT/GRPT` as
  instantaneous tuples, then `updateState_*`, `Routing_KW`, and the
  irrigation/flood rate rules.

## Next steps (stage 2.5 / stage 3)

1. Stage-1 instrumentation dumps → fill the `tests` blocks (per-cell traces
   with 1-based cell indices pinned to the dump order).
2. frivinp network `data_sources` (static) + ID/dnID gather → replace the
   `eroutUp`/`direct` placeholders.
3. Mount into `eqesm.esm` and integrate with earthsciintegrator; compare
   against full-CESM simulations (conservation via `mosart_budget_type.F90`).
