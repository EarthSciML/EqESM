# MOM6 stage 2 — penetrating shortwave absorption (`mom6_shortwave.esm`)

Implements the MOM6 penetrating shortwave absorption path of
`MOM_opacity.F90` (MOM6 dev/ncar_260909 @301879336, fork
`stage1/EqESM-cesm-mom6` branch `instrument/mom6-dumps`; in this MOM6
version the two routines live in `MOM_opacity.F90`, not in a separate
`MOM_shortwave_abs.F90`):

- `absorbRemainingSW` (`MOM_opacity.F90:602-869`) — heats the water column
  by the Beer-Lambert decay of the penetrating bands;
- `sumSWoverBands` (`MOM_opacity.F90:894-1041`) — accumulates the
  penetrating shortwave at the interfaces for the KPP buoyancy-flux path;
- the per-band surface split `extractFluxes1d`
  (`MOM_forcing_type.F90:619`).

Run configuration: `C_JRA_t62_mom6_inst`, `OPACITY_SCHEME=OHLMANN_03`
(chlorophyll-driven opacity lookup, `opacity_from_chl` :247-479),
`N_SW_BANDS=2`, `RHO_0=1035 kg m-3`, `C_P=3992 J kg-1 K-1`, all unit
rescale factors 1.0 (internal == MKS, so the dumped values are the
physical-unit law values).

## Files

- `stage2/mom6_shortwave.esm` — the `ShortwaveAbsorption` component
  (11 equations, 6 tests, 78 assertions at rel 1e-12);
- `stage2/mom6_sw_discretization.esm` — local discretization rules that
  lower the optical-depth integrals: `sw_optical_depth` (heating path,
  bounds `z_sfc`, measure `h`) and `sw_optical_depth_buoy` (buoyancy
  path, bounds `z_sfc_buoy`, measure `h_buoy`), sharing the flat
  `sw_optical_depth_scan` prefix faq (a candidate EarthSciDiscretizations
  PR: the ESD integral rules constrain the integrand to a bare lev-shaped
  reference, so the per-band opacity needs these band-lifted match
  rules);
- `stage1/analysis/make_mom6_sw_abs_tests.py` — extracts the test tuples
  from the dumps and self-checks every pin against the law at rel 1e-12
  (run `python3 stage1/analysis/make_mom6_sw_abs_tests.py [dump_dir]`);
- `bugs.md` — one entry: the `I_Habs` length-scale inconsistency between
  `sumSWoverBands` and `absorbRemainingSW` (see below).

## The law (as math, not as the Fortran loop)

With `penSW_bnd(n) = sw_pen_band(n)·dt/(RHO_0·C_P)` [degC m] (the surface
split, `MOM_forcing_type.F90:619`) and the cumulative optical depth at
interface `i` as the Beer-Lambert integral
`tau(n,i) = ∫₀^{z_i} opacity(n,z′) dz′` (esm-spec §4.2), written with the
ESD `integral` op and lowered by the local match rules of
`./mom6_sw_discretization.esm` — the cumulative-from-bottom integral of
`integral_lev_nodes_cumulative_from_bottom`, lifted to the `[band, layer]`
integrand because the imported EarthSciDiscretizations integral rules
constrain the integrand to a bare lev-shaped reference; the rule body is
the per-interface prefix `faq` over `j < i` with the layer thickness as
the measure.  The per-interface band flux is named once:

```
F_pen(n,i)     = penSW_bnd(n)·exp(−tau(n,i))                     [degC m]
F_pen_buoy(n,i)= penSW_bnd(n)·exp(−tau_buoy(n,i))                [degC m]
pen_sw_rem(n,k)= F_pen(n,k+1)                                    [degC m]
dT(k)          = Σ_n F_pen(n,k)·(1 − exp(−opacity(n,k)·h(k)))/h(k) [degC]
netPen(i)      = Σ_n F_pen_buoy(n,i)/dt                          [degC m s-1]
absorbed_total = Σ_k dT(k)·h(k)                                  [degC m]
```

`tau` integrates with the `h` measure = h_after (the post-flux thickness
the heating path sees, `absorbRemainingSW :748-750`); `tau_buoy` with
`h_buoy` = h_before (the pre-flux thickness of the KPP buoyancy path,
`thickness_to_dz` via `MOM_forcing_type.F90:1102`, `sumSWoverBands
:980-982`); both integrals share the surface bound `z_sfc` = `z_sfc_buoy`
= 0 m and differ only in the measure.  The netPen profile decays with
h_before (verified to 1e-16), the pen_sw_rem/dT profiles with h_after.

## Dump rows pinned

`mom6_sw_abs.csv` rows joined to `mom6_bflux.csv` rows of the same
(gi, gj, cpl_step, tod) for `h_after`/`h_before`/`dt`.  In this run the
coupling steps map to time of day as `tod = 3600 + 3600·(cpl_step−1)`
(step 1 dumped twice: tod 3600 and 7200); `dt` is 3600 s on every traced
row.

| test id | (gi, gj, cpl_step, tod) | (i, j) | lon, lat | depth | sw_pen_band (b1, b2) |
|---|---|---|---|---|---|
| sw_g3j78_step1_t3600 | (3, 78, 1, 3600) | (7, 82) | −285.33°, −69.50° | 555.0 m | 95.61, 34.92 W m-2 |
| sw_g7j80_step1_t3600 | (7, 80, 1, 3600) | (11, 84) | −282.67°, −69.03° | 501.8 m | 109.19, 39.79 W m-2 |
| sw_g4j78_step1_t3600 | (4, 78, 1, 3600) | (8, 82) | −284.67°, −69.50° | 330.8 m (shallowest) | 79.30, 28.90 W m-2 |
| sw_g3j80_step1_t3600 | (3, 80, 1, 3600) | (7, 84) | −285.33°, −69.03° | 764.1 m (deepest) | 139.92, 50.60 W m-2 |
| sw_g3j78_step5_t21600 | (3, 78, 5, 21600) | (7, 82) | −285.33°, −69.50° | 553.3 m | 173.95, 63.54 W m-2 |
| sw_g4j79_step12_t46800 | (4, 79, 12, 46800) | (8, 83) | −284.67°, −69.27° | 561.7 m | 67.14, 24.38 W m-2 |

Each test pins 13 values: `pen_sw_rem` at band 1 layers 1/3/8 and the
deepest pinnable layer, band 2 layers 1 and the deepest pinnable layer;
`dT` at layers 1/3 and the deepest layer where **both** bands are still
unclamped; `netPen` at interfaces 1/2 and the deepest interface with
|netPen| > 1e-6; and `absorbed_total` (a derived conservation pin, not a
dump column).  All dump precision is es26.17, so the pins are at rel
1e-12; every pin is self-checked by the extraction script before it is
emitted.

## Pinnable region and the min_SW_heat clamp gap

`absorbRemainingSW` rewrites the decay when the remaining SW of a band
falls below a floor: for `answer_date >= 20190101` the clamp fires at
layer `k` for band `n` when `nsw·pen_sw_rem(n,k) < 2.5e-11·dt` degC m
(`PEN_SW_FLUX_ABSORB = 2.5e-11 degC m s-1`, `PEN_SW_ABSORB_MINTHICK =
1.0 m`, `MOM_opacity.F90:711-730`), zeroing the layer's absorption or
setting it to the min-SW-heat branch.  The clamp is a **stateful
recurrence** (the clamped remaining feeds the next layer) and also
references the model's `I_Habs = 1e3·GV%H_to_m` length scale — see the
`bugs.md` entry — so it is not expressed as pure math in the component;
the pinnable region is defined exactly by the clamp condition:

- layer `k` is pinnable iff `nsw·pen_sw_rem(n,k) ≥ 2.5e-11·dt` for every
  band (band 2's much larger opacity, 0.6058 vs 0.1663 m-1, makes its
  clamp fire much shallower — the dT deep pin must require **both** bands
  unclamped);
- `netPen` IS clamped by `sumSWoverBands`' own min_SW_heat zero branch
  (`MOM_opacity.F90:989`, with the 1000×-wrong `I_Habs` of the `bugs.md`
  entry): the zero branch fires at `nsw·Pen_SW_bnd ≤ 2.5e-11·1e3·
  (h−h_min_heat)` degC m s-1 (no dt factor: `sumSWoverBands` is called
  with dt = 1.0, `MOM_forcing_type.F90:1102`), so the dump writes exact
  zeros below the first clamped interface — e.g. row (3,78,1,3600) has
  netPen = 1.67e-11 degC m s-1 at interface 21 and exactly 0.0 at
  interfaces 22+ — not exp underflow; only interfaces with |netPen| >
  1e-6 are pinned (all far above the clamp region).

Inside the pinnable region the component reproduces the dumps at
rel 3.6e-15 (pen_sw_rem, 3431 layer values), 1e-13 (dT), and 1e-16
(netPen).  `absorbed_total` pins the conservation identity across the
whole column including the clamped layers: the clamp shifts only
~1e-24 degC m of the ~0.11 degC m total in this window (the unclamped
sum matches the dumped `Σ_k dT_k·h_after_k` to ~1e-15 rel).

## Coverage gaps (this data window)

- **min_SW_heat clamp region**: `pen_sw_rem`/`dT` below the pinnable
  layers are pinned only through `absorbed_total`; the clamp branches
  (zero and min-SW-heat) are not pinned as equations (a stateful
  recurrence, see above; an exact Python simulation of the clamped
  routine reproduces the dumps to 3e-14, 402 clamp firings across the
  156 rows, so the unpinned gap is understood, not unexplained).
- **Thin-column redistribution**: the `h_heat < H_limit_fluxes`
  redistribution path of `absorbRemainingSW` (`T_chg` spreading) is not
  exercised (all traced columns are deep).
- **OHLMANN_03 surface split**: the chlorophyll-lookup split
  `A_n(chl)·SW_tot` of `set_pen_shortwave` is unpinnable here (chl and
  the downwelling SW are not dumped), so `sw_pen_band` is an input of the
  component; the split law is left for a future dump.
- **max(0,·) guards**: the `max(0.0, ...)` floor on the min-branch
  `SW_trans` and the `Pen` guard are not exercised in the pinned region.
- **netPen bottom zeros**: below the first clamped interface the dump
  writes exact zeros — the `sumSWoverBands` min_SW_heat zero branch
  (with the 1000×-wrong `I_Habs`, see the `bugs.md` entry), not exp
  underflow — unpinned (the |netPen| > 1e-6 pins stay above the clamp
  region).
- Only `nsw=2` (OHLMANN_03) bands are covered.

## bugs.md

`## TBD.` — `sumSWoverBands` hardcodes `I_Habs = 1e3*GV%H_to_m`
(`MOM_opacity.F90:950`) while `absorbRemainingSW` uses
`optics%PenSW_absorb_Invlen ≈ 1 m-1` (:711, 1211): a 1000× inconsistency
in the min_SW_heat clamp length scale.  Only affects the unpinned clamp
region in this window; documented, not corrected (the component pins the
unclamped law exactly as the model computes it).
