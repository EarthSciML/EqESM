# CLM snow/soil temperature — stage 2 `.esm` skeleton

Stage-2 implementation of the CTSM/CLM snow/soil temperature physics
(ctsm5.4.063, cesm3.0-alphabranch; instrumented fork commit `ab94ab484`,
dumps regenerated at fork clm `1c0d375dd`) as esm 1.3.0 files: the
`SoilThermProp` pointwise thermal-property laws and the implicit
`SoilTemperature` snow/soil temperature solve of
`components/clm/src/biogeophys/SoilTemperatureMod.F90`, with the shared
pointwise laws factored into an expression-template library
(`clm_templates.esm`). The `tests` blocks pin the stage-1 instrumentation
dumps (`stage1/dumps/clm/clm_soiltemp.csv`, W3 run LID 261103-111900,
16 traced columns × 34 layers × 12 steps; every pin is column 1 at step 1)
— see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `components/clm/src/biogeophys/SoilTemperatureMod.F90`
(`SoilThermProp` :604-903, `ComputeHeatDiffFluxAndFactor` :1825-1936,
`SetRHSVec_Snow` :2081-2174, `SetRHSVec_Soil` :2256-2379,
`SetMatrix_Snow` :2617-2710, `SetMatrix_Soil` :2713-2859), the constants of
`clm_varcon.F90`/`shr_const_mod.F90`, and the instrumented dump module
`stage1/EqESM-cesm-ctsm/components/clm/src/main/clm_dump_mod.F90`.

Default configuration transcribed (the D-class reference run's namelist):
`cnfac = 0.5` (Crank-Nicolson), `capr = 0.34` (top snow layer surface-flux
weight), `dtime = 1800 s`, `snow_thermal_cond_method = 'Sturm1997'`,
`nlevgrnd = 25`, `snl = -9` on the traced columns, `use_excess_ice = .true.`
(excess ice 0 on every traced row), `frac_sno_fluxes = 1` exactly,
`eflx_bot = 0` (no geothermal flux), `imelt = 0` on every traced row
(no phase change).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `clm_templates.esm` | Expression-template library (26 templates): snow bulk density, Sturm1997 snow thermal conductivity (with the float32-rounded coefficients of bugs.md #31), Jordan1991 snow thermal conductivity (transcribed, untested — see [Coverage gaps](#coverage-gaps)), Johansen (1975) soil thermal conductivity and its Kersten-number/liquid-fraction/saturated-conductivity factors, de Vries (1963) soil heat capacity, snow heat capacity with the thin-layer floor, harmonic-mean interface conductivity with the zero bottom boundary | `SoilTemperatureMod.F90`: SoilThermProp :604-903 |
| `clm_soiltemp.esm` | Model `SoilThermProp` — the per-layer pointwise laws: snow bulk density (:744), Sturm1997/Johansen/bedrock layer conductivity (:774-778, :712-727), snow/soil/bedrock heat capacity (:892-895, :852-854), harmonic interface conductivity (:825-828); 4 equations over the 34-layer CLM top-down column (j = -8..25); 68 pinned outputs | `SoilTemperatureMod.F90`: SoilThermProp :604-903 |
| | Model `SoilTemperature` — the implicit Crank-Nicolson snow/soil temperature solve as a residual-form finite-volume PDE on the EarthSciDiscretizations `column_nonuniform_1d` grid: `L_cn = D(K·D(u_cn, lev), lev)` lowered by the rule `varcoeff_face_laplacian_lev_flux_bc` (imported by reference, `${ESD_ROOT}`), with the rule's face closures `kdudz_top = hs_top_snow + dhsdT·(u_cn[NLEV] − t_entry[NLEV])` and `kdudz_bot = eflx_bot`; unknowns `u_cn/du/zc/A_capr/L_cn/fn_esd/fact_esd/kdudz_top/kdudz_bot/cv_dudt/resid`; 6 tests × 144 pinned outputs (see [The SoilTemperature model](#the-soiltemperature-model)) | `SoilTemperatureMod.F90`: ComputeHeatDiffFluxAndFactor :1825-1936, SetRHSVec_Snow :2081-2174, SetRHSVec_Soil :2256-2379, SetMatrix_Snow :2617-2710, SetMatrix_Soil :2713-2859 |

## The SoilTemperature model

CLM advances the snow/soil temperatures by solving the tridiagonal system

```
t_new(j) − t_old(j) = fact(j)·[ (fn_cn(j) − fn_cn(j−1)) + sabg(j) ]          (interior rows)
```

with `fact(j) = dtime/cv(j)` and the interface fluxes

```
fn(j) = tk(j)·(t(j+1) − t(j)) / (z(j+1) − z(j))
```

evaluated at the Crank-Nicolson temperature `u_cn = cnfac·t_entry +
(1 − cnfac)·t_exit` (the fluxes are linear in temperature, so
`cnfac·fn(entry) + (1 − cnfac)·fn(exit) = fn(u_cn)` and one evaluation at
`u_cn` reproduces the system exactly). The `.esm` model writes the CLM system
as a finite-volume diffusion PDE on the ESD `column_nonuniform_1d` grid:

- **Mapping.** ESD lev n = CLM j = 26 − n (lev 1 = deepest soil layer
  j = 25, lev 34 = top snow layer j = −8), with the interface coordinates
  `ze[n] = −zi(26 − n)` (ze[1] = −49.556 m the column bottom, ze[35] =
  +22.132 m the snow surface) and the layer centres `zc = column_cell_center(ze)`
  (exactly −z(26 − n)). The grid's `lev` axis runs bottom → top as the rule
  requires; the CLM top-down arrays are reversed at the component boundary.
- **Operator.** `L_cn = D(K·D(u_cn, lev), lev)` is lowered by the
  `varcoeff_face_laplacian_lev_flux_bc` rule to the CLM flux differences:
  `L u[n] = (G[n+1] − G[n])/dz[n]` with the interface fluxes
  `G[n] = K[n]·(u[n] − u[n−1])/(zc[n] − zc[n−1])`, so `dz[n]·L_cn[n] =
  fn_cn(26 − n) − fn_cn(27 − n)`. The rule's two face closures carry the CLM
  boundary fluxes:
  - `kdudz_top = hs_top_snow + dhsdT·(u_cn[NLEV] − t_entry[NLEV])` — the
    surface flux at the CN temperature. The CLM top snow row is *fully
    implicit* in the surface flux: SetMatrix_Snow's diagonal carries
    `−fact·dhsdT` (not `(1 − cnfac)·fact·dhsdT`), i.e. the row's surface-flux
    variation is `dhsdT·(t_exit − t_entry)`. The top-layer residual's
    `(1/(2A))` weighting doubles the `dhsdT·du` term and restores exactly that
    fully implicit form (verified against the dump to round-off).
  - `kdudz_bot = eflx_bot` — the bottom boundary flux (`fn(nlevgrnd) =
    eflx_bot` in the Fortran; 0 in the traced run).
- **Top-row factor.** The top snow layer's tridiagonal factor is the
  capr-weighted `fact(−8) = dtime/(cv·A)` with
  `A = 0.5·((ze[NLEV+1] − zc[NLEV]) + capr·(ze[NLEV+1] − zc[NLEV−1]))/dz[NLEV]`
  (ComputeHeatDiffFluxAndFactor :1927-1928), carried as the scalar unknown
  `A_capr = 0.643828751513956` in the traced column.
- **Residual form.** The model is *not* time-integrated: the explicit harness
  integrator cannot reproduce the implicit CN solve, so the equations are the
  algebraic statement that the dumped exit state solves the CLM system. With
  the dumped entry/exit temperatures, conductivities and heat capacities as
  parameters, the per-row residuals
  - interior: `cv·du/dt = 0.5·(dz·L_cn + sabg)`
  - bottom: `cv·du/dt = 0.5·(dz·L_cn + sabg) + eflx_bot` (the explicit
    `−eflx_bot` doubles the one already inside `dz·L_cn[1]` through
    `kdudz_bot`, restoring the row's `eflx_bot`)
  - top: `cv·du/dt = (1/(2A))·(dhsdT·du + dz·L_cn)` (no sabg: the top-layer
    solar absorption is inside `hs_top_snow`)
  evaluate to round-off only — max |resid| = 3.2e-10 W/m² over the 34 rows —
  pinning that the dumped `t_exit` is the exact solution of the CLM
  tridiagonal system built on the dumped properties. `fn_esd` (the dump's
  `fn` column, entry-state fluxes, `fn(nlevgrnd) = eflx_bot`) and `fact_esd`
  are pinned as well, matching the dump to rel 4.8e-14 / 1.1e-16.
- **dz parameter.** The rule's finite-volume measure `dz` is fed the dumped
  CLM dz (reversed) as a parameter rather than derived from `ze`: the top snow
  layer's (z, dz, zi) triple is internally inconsistent in the dump
  (bugs.md #32), and the CLM rows use the dumped dz — deriving dz from ze
  would shift the top-layer `L` by ~0.2 W/m³ and break the row reproduction.
  The layer centres `zc` *are* ze-derived (`column_cell_center`), which holds
  exactly.

## Coverage gaps

The traced W3 run exercises only a corner of the CLM soil-temperature code;
each gap below is a documented, deliberate non-coverage with its reason:

- **imelt = 0 everywhere** — no phase change occurs in the 12 traced steps,
  so the `Phasechange` laws of SoilTemperatureMod.F90 (the latent-heat
  tridiagonal re-solve and the freeze/melt temperature clamps) are not
  exercised and not transcribed. The SoilTemperature model covers the
  `imelt = 0` branch exactly (the dump's rows carry no latent-heat terms).
- **frac_sno_fluxes = 1 exactly** — the j = 1 snow/soil interface row of
  SetRHSVec_Soil/SetMatrix_Soil is the general
  `frac_sno_fluxes`-weighted split surface row; with `frac_sno_fluxes = 1`
  the `(1 − frac_sno_fluxes)·(hs_soil − dhsdT·t)` term vanishes, `fn(0)`
  enters unweighted and the row collapses to the standard interior form
  (plus sabg), which is what the model carries. The 0 < frac_sno_fluxes < 1
  branch is untested.
- **Jordan1991 snow conductivity untested** — `snow_thermal_cond_method =
  'Sturm1997'` in the traced run; the Jordan1991 template exists in
  `clm_templates.esm` but has no pins.
- **h2osno_no_layers = 0** — the `cv(c,1) += cpice·h2osno_no_layers`
  addition of SoilThermProp :880-882 never binds.
- **excess_ice = 0, snow h2osoi_liq = 0** — the excess-ice and snow-liquid
  terms of the conductivity/heat-capacity laws never bind.
- **eflx_bot = 0** — the bottom boundary is insulated in the traced run; the
  nonzero geothermal-flux branch of the bottom row is untested (the model
  carries it: `kdudz_bot = eflx_bot` and the bottom residual's `−eflx_bot`).
- **frac_h2osfc = 0** — the surface-water correction of SetRHSVec_Soil
  (:2373-2377) never binds (no ponded water on the traced column).
- **Snow-soil flux split** — the `j = nlevurb` urban branch and the
  `IsSimpleBuildTemp` building branches of ComputeHeatDiffFluxAndFactor are
  not transcribed (non-urban traced column).

## Dump defects

- **bugs.md #31** — the Sturm1997 conductivity literals lack `_r8` and
  execute as float32-rounded values; reproduced exactly in
  `snow_thermal_cond_sturm1997` (rel 1.6e-16 pins).
- **bugs.md #32** — the top snow layer's (z, dz) is inconsistent with its zi
  interfaces (snow hydrology leaves the center/thickness updated but not the
  interface depths); absorbed exactly by feeding the dumped `zi_bot` through
  the harmonic-mean template and the dumped dz through the PDE rule.

## Testing

`ESD_ROOT` must point at a checkout of
`github.com/earthsciml/EarthSciDiscretizations` (the SoilTemperature model
imports the column diffusion rule by reference; the whole stage-2 suite
requires it for the CICE B-grid rules as well). Then:

```
ESD_ROOT=/path/to/EarthSciDiscretizations ./esm test stage2
```

`clm_soiltemp.esm` pins 212 assertions (68 SoilThermProp + 144
SoilTemperature) from `stage1/dumps/clm/clm_soiltemp.csv` column 1 step 1:
every pin is a direct reproduction of a dumped value (rel 1e-12, abs 1e-8 for
the round-off residuals), with no fitted parameters.
