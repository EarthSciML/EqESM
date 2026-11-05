# CLM snow/soil temperature — stage 2 `.esm` skeleton

Stage-2 implementation of the CTSM/CLM snow/soil temperature physics
(ctsm5.4.063, cesm3.0-alphabranch; instrumented fork commit `1c0d375dd`,
dumps from the W3 run, LID 261103-111900) as esm 1.3.0 files: the
`SoilTemperature` heat-conduction PDE

```
cv·∂T/∂t = dz·D(K·D(T, lev), lev) + sabg
```

on the EarthSciDiscretizations `column_nonuniform_1d` grid, with the
`SoilThermProp` thermal-property laws of
`components/clm/src/biogeophys/SoilTemperatureMod.F90` as plain array
equations, the shared pointwise laws factored into an expression-template
library (`clm_templates.esm`), and the CLM top-row finite-volume measure in a
local discretization rule (`clm_soiltemp_discretization.esm`, a candidate
EarthSciDiscretizations PR). The `tests` blocks pin the stage-1 instrumentation
dump (`stage1/dumps/clm/clm_soiltemp.csv`, 16 traced columns × 34 layers × 12
steps; every pin is column 1 at step 1) — see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `components/clm/src/biogeophys/SoilTemperatureMod.F90`
(`SoilThermProp` :604-903, `ComputeGroundHeatFluxAndDeriv` :1545-1798,
`ComputeHeatDiffFluxAndFactor` :1801-1936, `SetRHSVec_Snow` :2081-2174,
`SetRHSVec_Soil` :2256-2379, `SetMatrix_Snow` :2617-2710, `SetMatrix_Soil`
:2713-2859), the constants of `clm_varcon.F90`/`shr_const_mod.F90`, and the
instrumented dump module `stage1/EqESM-cesm-ctsm/components/clm/src/main/clm_dump_mod.F90`.

Default configuration of the traced run: `cnfac = 0.5` (Crank-Nicolson),
`capr = 0.34` (top snow layer surface-flux weight), `dtime = 1800 s`,
`snow_thermal_cond_method = 'Sturm1997'`, `nlevgrnd = 25`, `snl = -9` on the
traced columns (9 snow layers j = -8..0), `use_excess_ice = .true.` (excess
ice 0 on every traced row), `frac_sno_fluxes = 1` exactly, `eflx_bot = 0` (no
geothermal flux), `imelt = 0` on every traced row (no phase change). All soil
layers are frozen (225.7 K < tfrz), so the Johansen law runs in its
frozen branch and the thermal properties are temperature-independent in the
pinned regime.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `clm_templates.esm` | Expression-template library (31 templates: 11 soil-temperature laws + 15 shared pointwise laws + 5 albedo): snow bulk density, Sturm1997 snow thermal conductivity (with the float32-rounded coefficients of bugs.md #33), Jordan1991 snow thermal conductivity (transcribed, untested — see [Coverage gaps](#coverage-gaps)), Johansen (1975) soil thermal conductivity in the named-intermediate form (`satw`/`dke`/`dksat` params, the caller computes each once), its factor laws (saturation water content, Kersten number, liquid fraction, saturated conductivity), de Vries (1963) soil heat capacity, snow heat capacity with the thin-layer floor, harmonic-mean interface conductivity with the zero bottom boundary | `SoilTemperatureMod.F90`: SoilThermProp :604-903 |
| `clm_soiltemp_discretization.esm` | ESD-format rewrite rule `varcoeff_face_laplacian_lev_flux_bc_capr` (priority 11): the conservative variable-coefficient vertical diffusion `D(K·D(u, lev), lev)` of `varcoeff_face_laplacian_lev_flux_bc` with the CLM top-row measure `M = 0.5·(ze[NLEV+1] − zc[NLEV] + capr·(ze[NLEV+1] − zc[NLEV−1]))` in place of `dz[NLEV]`; imports the two ESD stencils by reference; candidate EarthSciDiscretizations PR | `SoilTemperatureMod.F90`: ComputeHeatDiffFluxAndFactor :1895, SetMatrix_Snow :2617-2710 |
| `clm_soiltemp.esm` | Model `SoilTemperature` — the heat-conduction PDE `cv·∂T/∂t = dz·D(K·D(T,lev),lev) + sabg` with the state T, the SoilThermProp property laws (consuming `t_entry`, as in CLM), the harmonic-mean interface conductivities on the face axis, the linearized surface flux `kdudz_top = hs + dhsdT·(T[NLEV] − t_entry[NLEV])`, the bottom flux `kdudz_bot = eflx_bot`, the observeds `fn_esd`, `dTdt`, `cn_balance`; 3 tests × 205 pinned/derived assertions (see [The SoilTemperature model](#the-soiltemperature-model)) | `SoilTemperatureMod.F90`: SoilThermProp :604-903, ComputeGroundHeatFluxAndDeriv :1545-1798, the implicit solve :1825-2859 |

## The SoilTemperature model

CLM advances the snow/soil temperatures with the implicit Crank-Nicolson
solve of SetMatrix_Soil/SetMatrix_Snow on the tridiagonal row system. The
`.esm` model is a genuine PDE — the equation `cv·∂T/∂t = dz·D(K·D(T, lev),
lev) + sabg` with the state `T` — whose spatial operator, evaluated at the
pinned states, reproduces the CLM row system; it is **not** a transcription
of the tridiagonal matrix rows.

- **Mapping.** ESD lev n = CLM j = 26 − n (lev 1 = deepest soil layer
  j = 25, lev 34 = top snow layer j = −8), with the interface coordinates
  `ze[n] = −zi(26 − n)` (ze[1] = −49.556 m the column bottom, ze[35] =
  +22.132 m the snow surface) and the layer centres `zc = column_cell_center(ze)`
  (exactly −z(26 − n)). The grid's `lev` axis runs bottom → top as the rule
  requires.
- **Operator.** `D(K·D(T, lev), lev)` is lowered by the local rule
  `varcoeff_face_laplacian_lev_flux_bc_capr` to the CLM flux differences:
  interior rows `L u[n] = (G[n+1] − G[n])/dz[n]` with the interface fluxes
  `G[n] = K[n]·(u[n] − u[n−1])/(zc[n] − zc[n−1])` (the ESD stencil
  `interface_flux_K_lev`), the bottom row `(G[2] − kdudz_bot)/dz[1]`, and the
  top row `(kdudz_top − G[NLEV])/M` with the **capr-weighted measure**
  `M = 0.5·(ze[NLEV+1] − zc[NLEV] + capr·(ze[NLEV+1] − zc[NLEV−1]))` of
  ComputeHeatDiffFluxAndFactor :1895 (M = 0.01002979145516516 m in the traced
  column; A ≡ M/dz[NLEV] = 0.643828751513956). The model multiplies the
  lowered divergence by the *dumped* layer thickness `dz` (a parameter — the
  top layer's dumped dz = 0.015578352833079012 differs from the grid spacing
  ze[NLEV+1] − ze[NLEV] = 0.015517897173545236, bugs.md #34) and adds `sabg`.
- **Boundary fluxes.** The rule's free names carry the CLM boundary fluxes:
  `kdudz_top = hs + dhsdT·(T[NLEV] − t_entry[NLEV])`, the surface flux
  linearized about the entry state (ComputeGroundHeatFluxAndDeriv :1545-1798,
  SetMatrix_Snow), and `kdudz_bot = eflx_bot`, the prescribed bottom flux
  (fn(nlevgrnd) = eflx_bot in the Fortran; 0 in the traced run). The top
  snow layer's absorbed solar (sabg = 57.79 W/m² in the dump) is inside
  `hs_top_snow` (ComputeGroundHeatFluxAndDeriv :1779-1783 adds no sabg to the
  top row), so the model's `sabg` parameter carries 0 at the top layer.
- **Properties.** The SoilThermProp laws are plain array equations over `lev`
  (no faq/clamped gathers): the layer masks come from metaparameters —
  `is_snow` = makearray over N_SNOW_MAX = 9 snow layers, `is_bedrock` =
  makearray over the nbedrock = 20 bottom layers — and the layer laws are
  `ifelse(is_snow, …, ifelse(is_bedrock, …, …))` with the discarded branches
  never evaluated. `thk` composes the Sturm1997 snow law and the Johansen
  (1975) soil law with the named intermediates `satw`, `dke`, `fl`, `dksat`
  (each computed once as an observed; the ESM template format admits
  params + body only, so `soil_thermal_cond` takes the intermediates as
  params). All laws consume the entry state `t_entry`, as in CLM
  (SoilThermProp runs at the start of the step). The interface conductivities
  `K` are the harmonic mean of `clm_interface_conductivity` gathered over the
  face index set (`faq` over `lev_nodes` — the library's own adjacent-gather
  idiom, cf. `nonuniform_cell_center`): `K[n] = tk(26−n)` for n = 2..34,
  K[1] = 0 (the zero bottom boundary of :828), K[35] = 0 dummy above the top
  layer.
- **The three pins.** The state `T` is pinned by initial condition at both
  dumped stages with `time_span {0, 0}` (instantaneous evaluation), and the
  properties are frozen from `t_entry` exactly as in the CLM step:
  1. *Entry pin* (`clm_soiltemp_entry`, 137 assertions): `T = t_entry`;
     `K` (35 faces) and `cv` (34 layers) reproduce the dumped `tk`/`cv`
     through the full property chain, `fn_esd[n] = −G[n]` reproduces the
     dumped `fn` column (which is exactly the entry-state interface flux
     `−G(t_entry)` written at the exit stage, verified numerically to 1e-13),
     and `dTdt` reproduces the *derived* instantaneous PDE tendency
     `(dz·D(K·D(t_entry,lev),lev) + sabg)/cv`.
  2. *Exit pin* (`clm_soiltemp_exit`, 35 assertions): `T = t_exit` under the
     same frozen properties; `dTdt` reproduces the derived tendency at the
     exit state and `kdudz_top = hs + dhsdT·(t_exit[NLEV] − t_entry[NLEV])`
     the derived surface flux.
  3. *Derived CN identity* (`clm_soiltemp_cn_balance`, 33 assertions): the
     Crank-Nicolson energy balance
     `0.5·(F(t_entry) + F(t_exit)) = cv·(t_exit − t_entry)/dt` with
     `F(u) = dz·D(K·D(u,lev),lev) + sabg` — the flux operator is linear in
     the state (K frozen from `t_entry`), so `0.5·(F(t_entry)+F(t_exit)) =
     dz·D(K·D(t_avg,lev),lev) + sabg` with `t_avg = 0.5·(t_entry + t_exit)`,
     and the model's `cn_balance` observed evaluates the residual. With the
     dumped states it lands on the CLM solver's round-off for the interior
     and bottom rows (lev 1..33): max |resid| = 6.3e-10 W/m², pinned with
     `abs` tolerance 1e-9 W/m². This pins that the dumped `t_exit` is the
     exact solution of the CLM tridiagonal system built on the dumped
     properties.
- **Top-row caveat.** The identity holds for rows 1..33. The top row of
  SetMatrix_Snow is *fully implicit* in the surface flux — its diagonal
  carries `−fact·dhsdT` (not `(1 − cnfac)·fact·dhsdT`) — so the top-row
  balance is `cv·(M/dz[NLEV])·ΔT/dt = 0.5·(kdudz_top(t_entry) +
  kdudz_top(t_exit)) − 0.5·(G(t_entry) + G(t_exit))` with the +0.5·dhsdT·ΔT
  fully-implicit term in place of the plain CN average; equivalently the plain
  identity overshoots the top row by `0.5·dhsdT·ΔT/A` (A = M/dz[NLEV]),
  verified numerically against the dumped row. The model's top row (the
  `kdudz_top` linearization in the rule's face closure) reproduces this
  exactly; the top row of `cn_balance` is therefore not pinned.
- **eflx_bot sign.** The F90 comment at the eflx_bot declaration says
  "+upward" while the `fn(nlevgrnd) = eflx_bot` slot is +downward
  (SetRHSVec_Soil :2379); with eflx_bot = 0 in the dump the sign convention is
  unpinned and the model documents the ambiguity rather than guessing.

## Coverage gaps

The traced W3 run exercises only a corner of the CLM soil-temperature code;
each gap below is a documented, deliberate non-coverage with its reason:

- **imelt = 0 everywhere** — no phase change occurs in the 12 traced steps,
  so the `Phasechange` laws of SoilTemperatureMod.F90 (the latent-heat
  tridiagonal re-solve and the freeze/melt temperature clamps) are not
  exercised and not transcribed. The model covers the `imelt = 0` branch
  exactly.
- **frac_sno_fluxes = 1 exactly** — the j = 1 snow/soil interface row of
  SetRHSVec_Soil/SetMatrix_Soil is the general
  `frac_sno_fluxes`-weighted split surface row; with `frac_sno_fluxes = 1`
  the `(1 − frac_sno_fluxes)·(hs_soil − dhsdT·t)` term vanishes and the row
  collapses to the standard interior form. The 0 < frac_sno_fluxes < 1
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
  carries it: `kdudz_bot = eflx_bot`).
- **frac_h2osfc = 0** — the surface-water correction of SetRHSVec_Soil
  (:2373-2377) never binds (no ponded water on the traced column).
- **Snow-soil flux split** — the `j = nlevurb` urban branch and the
  `IsSimpleBuildTemp` building branches of ComputeHeatDiffFluxAndFactor are
  not transcribed (non-urban traced column).

## Dump defects

- **bugs.md #33** — the Sturm1997 conductivity literals lack `_r8` and
  execute as float32-rounded values; reproduced exactly in
  `snow_thermal_cond_sturm1997` (rel 1.6e-16 pins).
- **bugs.md #34** — the top snow layer's (z, dz) is inconsistent with its zi
  interfaces (suspected DivideSnowLayers/AddNewSnow: snow hydrology leaves the
  center/thickness updated but not the interface depths); absorbed exactly by
  feeding the dumped `zi_bot` through the harmonic-mean template, the dumped
  dz through the PDE rule's measure multiplication, and the capr measure from
  ze/zc (never the dumped dz).

## Testing

`ESD_ROOT` must point at a checkout of
`github.com/earthsciml/EarthSciDiscretizations` (the model imports the column
grid, the two stencils and the flux-BC rule by reference; the whole stage-2
suite requires it for the CICE B-grid rules as well). Then:

```
ESD_ROOT=/path/to/EarthSciDiscretizations ./esm test stage2
```

`clm_soiltemp.esm` pins 205 assertions (137 entry + 35 exit + 33 derived CN
identity) from `stage1/dumps/clm/clm_soiltemp.csv` column 1 step 1: every
direct pin reproduces a dumped value at rel 1e-12; the dTdt and kdudz_top
values are *derived* (computed in Python from the dumped quantities by
`stage1/analysis/extract_soiltemp_testdata.py`, which also verifies the
embedded tuples with `--check`); the CN identity residual is expected 0 within
abs 1e-9 W/m² (the CLM solver's round-off). No fitted parameters.
