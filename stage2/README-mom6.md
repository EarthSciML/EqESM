# MOM6 stage 2 — equation of state, surface forcing, and CVMix KPP

Stage-2 implementation of three MOM6 boundaries of the MOM6 lane
(`stage1/runs/MOM6_LANE.md`): the Wright 1997 full-range equation of state
(boundary B1, exercised by `KPP_compute_BLD`), the CVMix KPP planetary
boundary layer (boundary B2, `KPP_compute_BLD` + `KPP_calculate`), and the
surface forcing import-to-fluxes/forces mapping (boundary I), as esm 1.3.0
files with `tests` blocks pinned to the stage-1 instrumentation dumps
(`stage1/dumps/mom6/mom6_eos.csv`, `mom6_kpp_obl.csv`, `mom6_kpp.csv`,
`mom6_sfc_forcing.csv`; 12 traced PET-0 cells × 12 coupling steps of the
reference run `C_JRA_t62_mom6_inst`, MOM6 `dev/ncar_260909 @301879336`, fork
`stage1/EqESM-cesm-mom6` branch `instrument/mom6-dumps`).

Sources:
`components/mom/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90`,
`components/mom/MOM6/config_src/drivers/nuopc_cap/mom_surface_forcing_nuopc.F90`,
`components/mom/MOM6/src/core/MOM_forcing_type.F90`,
`components/mom/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90`,
`components/mom/MOM6/pkg/CVMix-src/src/shared/cvmix_kpp.F90`,
`components/mom/MOM6/pkg/CVMix-src/src/shared/cvmix_math.F90`,
`components/mom/MOM6/src/framework/mom6_dump_mod.F90`.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `mom6_templates.esm` | Expression-template library: `wright_alpha0`, `wright_p0`, `wright_lambda` (the three shared polynomials), `wright_density`, `wright_drho_dT`, `wright_drho_dS`, the 15 constant fragments `wright_coeff_a0..c5` (the single shared source of the fit coefficients), the `wright_rho` composition (the whole in-situ law in one application, used by KPPBLD), `masked_pair_average` (the AGRID h-point stress average) | `MOM_EOS_Wright_full.F90`: density_elem_Wright_full 75-91, calculate_density_derivs_elem_Wright_full 179-204, constants 23-37; `mom_surface_forcing_nuopc.F90`: 943-951 |
| `mom6_eos.esm` | Model `WrightEOS` — rho(T, S, p), drho/dT, drho/dS via the templates (imported by reference), `I_denom2`; 15 dump-pinned tests | `MOM_EOS_Wright_full.F90` as above; the EOS calls themselves are `calculate_density`/`calculate_density_derivs_scalar` in `MOM_CVMix_KPP.F90:1437` |
| `mom6_kpp_templates.esm` | Expression-template library for the KPP laws: `kpp_turbulent_scale` (the MXL w_s/w_m stability functions), `kpp_unresolved_shear` (Vt² with the Danabasoglu Cv), `kpp_bulk_richardson`, `kpp_shape_cubic` (G(σ), evaluate_cubic order), `kpp_quad_interp_c/b/a` and `kpp_quad_root` (the OBL-depth quadratic), `kpp_override_nonzero` (the KPP_IS_ADDITIVE=False replacement) | `cvmix_kpp.F90`: turbulent_scales 2416-2604, unresolved_shear 2611-2903, bulk_Richardson 2031-2132, shape_function_coeffs 3079-3112; `cvmix_math.F90`: linear_through_cubic_poly_interp 103-138, cubic_root_find 260-294, evaluate_cubic 300-340 |
| `mom6_kpp.esm` | Models `KPPBLD` (KPP_compute_BLD down to the pre-smoothing OBLdepth/kOBL: surface-layer averages, the Wright EOS as three `wright_rho` applications at the layer-top pressure, deltaBuoy/N2/N, w_s/Vt²/BulkRi, OBL_depth_low with the Monin-Obukhov limit and the MOM6 bounds) and `KPPCoefs` (cvmix_coeffs_kpp_low: the no-match shape G(σ), the 2..floor(kOBL) overwrite, the KPP_IS_ADDITIVE=False replacement via `kpp_override_nonzero` with `Kd_salt = Kd_heat`, the CVMix non-local term, all driven by the *smoothed* OBLdepth/kOBL from the dump); 8 dump-pinned tests (5 BLD rows × 34..39 + 3 coefficient rows × 21..28 = 262 assertions) | `MOM_CVMix_KPP.F90`: KPP_compute_BLD 1195-1575, KPP_calculate 740-1067, KPP_smooth_BLD 1715-1834; `cvmix_kpp.F90`: OBL_depth_low 1585-1814, kOBL 1821-1867, coeffs_kpp_low 701-1298 |
| `mom6_surface_forcing.esm` | Model `SurfaceForcing` — the 25-field IOB import → fluxes/forces mapping, 39 outputs, 6 dump-pinned tests | `mom_surface_forcing_nuopc.F90`: convert_IOB_to_fluxes 258-734, convert_IOB_to_forces 739-1050; `MOM_forcing_type.F90`: set_net_mass_forcing/get_net_mass_forcing 2748-2800 |

## Units

| Quantity | Symbol | Units |
|---|---|---|
| Potential temperature | `T` | degC |
| Salinity | `S` | ppt |
| Pressure | `p`, `p0`, `p_surf_full` | Pa |
| Density | `rho` | kg/m³ |
| Density derivatives | `drho_dT`, `drho_dS` | kg/m³/degC, kg/m³/ppt |
| Specific volume (alpha_0) | `al0` | m³/kg |
| Sound-speed-squared offset | `lam` | m²/s² |
| EOS coefficients | `a0..a2`, `b0..b5`, `c0..c5` | as in the Fortran's units comments (F90:23-37) |
| Heat fluxes | `sw`, `lw`, `latent`, `sens`, `seaice_melt_heat_o`, `heat_added`, `heat_content_*`, `latent_*_diag` | W/m² |
| Water fluxes | `evap`, `lprec_o`, `fprec_o`, `vprec`, `lrunoff_o`, `frunoff_o`, `lrunoff_glc`, `frunoff_glc`, `seaice_melt_o`, `salt_flux_o`, `salt_flux_in`, `net_mass_src` | kg/m²/s |
| Wind stresses | `u_flux`, `v_flux`, `taux`, `tauy`, `tau_mag`, `GUST_CONST` | Pa |
| Friction velocity | `ustar` | m/s |
| Wind angle | `omega_w2x` | rad |
| Wind speed squared | `u10_sqr`, `u10_sqr_o` | m²/s² |
| Masks / fractions | `mask_east`, `mask_north`, `ice_fraction`, `ice_fraction_o`, `u10_sqr_allocated` | 1 |
| Latent heats | `Lf`, `Lv` | J/kg |
| Reference density | `RHO_0` | kg/m³ |
| Layer/interface heights | `z_cell`, `z_inter`, `SLdepth`, `delH`, `OBLdepth`, `OBLdepth_s`, `root` | m |
| Layer thicknesses | `h`, `dz` | m |
| Velocities | `U_H`, `V_H`, `surfU`, `surfV`, `Uk`, `Vk`, `uStar`, `uStar_s`, `w_s`, `w_m` | m/s |
| Shear | `deltaU2`, `Vt2` | m²/s² |
| Buoyancy | `buoyFlux`, `surfBuoyFlux_s`, `deltaBuoy`, `B0` | m²/s³ (buoyFlux), m/s² (deltaBuoy) |
| Buoyancy frequency | `N2`, `N` | s⁻², s⁻¹ |
| Bulk Richardson number | `BulkRi`, `quad_a2` | 1 |
| Quadratic coefficients | `quad_b2`, `quad_c2`, `slope` | m⁻¹, m⁻², m⁻¹ |
| Fractional OBL index | `kOBL`, `kOBL_s`, `kwup`, `kfirst`, `kwfirst` | 1 |
| Diffusivities/viscosity | `Kd_heat`, `Kd_salt`, `Kv`, `Kt_in`, `Ks_in`, `Kv_in`, `Kt_out`, `Ks_out`, `Kv_out` | m²/s |
| Non-local transport | `NLT_heat`, `NLT_salt`, `nonlocal_coeff`, `sigma`, `G` | 1 |
| EOS pressure | `pRef`, `p_surf` | Pa |
| KPP constants | `GORHO` = g/ρ₀ (m⁴ s⁻² kg⁻¹), `PRESS_STEP` = g·ρ₀ (Pa/m), `MIN_OBL_DEPTH` (m), `MIN_VTSQR` (m² s⁻²), `RI_CRIT`, `VON_KARMAN`, `SURF_LAYER_EXT`, `ZETA_M` = −0.2, `ZETA_S` = −1.0, `A_M`, `C_M`, `A_S`, `C_S`, `CSTAR` (dimensionless LMD94 constants) | as noted |

`degC` and `ppt` are accepted unit strings in esm 1.3.0; ppt is treated as
dimensionless, matching the model's EOS usage (the a2/b5/c5 coefficients carry
the salinity dimension instead).

## The Wright 1997 full-range EOS

`rho(T, S, p) = (p + p0) / (lambda + al0·(p + p0))` with the three polynomials

```
al0    = a0 + (a1·T + a2·S)
p0     = b0 + (b4·S + T·(b1 + (T·(b2 + b3·T) + b5·S)))
lambda = c0 + (c4·S + T·(c1 + (T·(c2 + c3·T) + c5·S)))
```

and the analytic first derivatives (F90:199-202)

```
drho/dT = I_denom2·(lambda·(b1 + (T·(2.0·b2 + 3.0·b3·T) + b5·S))
                    − (p+p0)·((p+p0)·a1 + (c1 + (T·(c2·2.0 + c3·3.0·T) + c5·S))))
drho/dS = I_denom2·(lambda·(b4 + b5·T) − (p+p0)·((p+p0)·a2 + (c4 + c5·T)))
I_denom2 = 1/(lambda + al0·(p+p0))²
```

`WrightEOS` binds the 15 Fortran coefficients as named parameters with units
(values at F90:23-37) and computes `al0`/`p0`/`lam` as intermediate unknowns
via the shared polynomial templates, then `rho` and the derivatives via
`wright_density`/`wright_drho_dT`/`wright_drho_dS` — the same operation order
as the Fortran, so the dump pins hold to round-off (rel 1e-12).

**Derivatives are transcribed formulas, not symbolic `D(rho, T)`.** In esm
1.3.0 the `D` operator is a rewrite-target for discretization stencils
(esm-spec §4.2 / §9.6.8) and the CLI 0.7.0 evaluator rejects it for algebraic
expressions (`unlowered_operator: rewrite-target operator 'D' reached
evaluation without being lowered to a stencil by a rewrite rule`), so the
Fortran's hand-coded derivative formulas are reproduced exactly instead (same
operation order; pinned by the drho_dT/drho_dS dump columns).

**Transcription bug caught by the pins.** The first hand-authored nesting of
the shared polynomials wrote `T·(b2·T + b3·T²)` where the Fortran has
`T·(b2 + b3·T)` — an extra factor of T in `wright_p0`, `wright_lambda`, and
the b-side of `wright_drho_dT`. The CLI evaluates nested `+`/`*`/`^`
expressions exactly as written (verified with isolated probes), so the dump
pins failed on every row until the nesting was corrected. This is a
documentation-grade reminder that the .esm transcription, not the CLI, was
at fault.

## CVMix KPP boundary layer

`mom6_kpp.esm` transcribes the two MOM6 KPP call sites of the reference
configuration (`C_JRA_t62_mom6_inst`: RI_CRIT=0.3, VON_KARMAN=0.4,
SURF_LAYER_EXTENT=0.1, MINIMUM_OBL_DEPTH=5.0, MINIMUM_VT2=1e-10,
INTERP_TYPE=quadratic, INTERP_TYPE2=LMD94, MATCH_TECHNIQUE=MatchGradient,
ENHANCE_DIFFUSION=False, NLT_SHAPE=CVMix, N_SMOOTH=3, KPP_IS_ADDITIVE=False,
KPP_SHORTWAVE_METHOD=MXL_SW, l_LMD_ws=false, USE_KPP_LT_K/VT2=False,
COMPUTE_MONIN_OBUKHOV=True, COMPUTE_EKMAN=False, all unit rescale factors
1.0).

- **`KPPBLD`** reproduces `KPP_compute_BLD` (`MOM_CVMix_KPP.F90:1195-1575`)
  down to the *pre-smoothing* OBLdepth and kOBL: the surface-layer means
  (F90:1342-1405), the per-layer hydrostatic EOS pressure
  (`term(i) = PRESS_STEP·h(i−1)`, PRESS_STEP = g·ρ₀ = 10143.0 Pa/m), the
  Wright 1997 EOS applied once per state (surface-layer mean, layer value,
  layer−1 value, all at the layer-top pressure pRef) through the shared
  `wright_rho` composition of `./mom6_templates.esm` (whose 15 fit
  coefficients are the `wright_coeff_*` constant fragments — the same values
  `mom6_eos.esm` declares as parameters, so neither model redeclares them),
  deltaBuoy/N2/N (F90:1444-1457), the MXL turbulent scale w_s, the
  unresolved shear Vt² (with the Danabasoglu et al. 2006 Cv = 2.1 − 200·N)
  and the bulk Richardson numbers (F90:1509-1550), then `OBL_depth_low`
  (`cvmix_kpp.F90:1585-1814`) and the kOBL fraction (1821-1867). The
  75/76-sized sentinels (the i = 76 zero branch of N, the N_LAYER-indexed
  OBL_limit, the kfirst/kwfirst = 76 no-crossing values and the 75.75
  bottom kOBL) are spelled with the N_LAYER/N_IFACE metaparameters.
- **`KPPCoefs`** reproduces `cvmix_coeffs_kpp_low` (701-1298, not
  lStokesMOST) as called by `KPP_calculate` (740-1067): with MatchGradient +
  ENHANCE_DIFFUSION=False the shape is the no-match cubic
  G(σ) = σ·(1−σ)², the overwrite loop covers exactly interfaces
  2..floor(kOBL), KPP_IS_ADDITIVE=False keeps the incoming value where the
  KPP value is zero (F90:1023-1038) — the same `kpp_override_nonzero` law
  shared by the Kt_out/Ks_out/Kv_out triple, with Kd_salt = Kd_heat — and
  the non-local term is the CVMix law
  `nonlocal_coeff·G(σ)` (NLT_SHAPE=CVMix: MOM6's override block at
  960-1001 is inactive, so NLT_heat = NLT_salt, F90:1002-1003) with
  nonlocal_coeff = Cstar·κ·(κ·ε·c_s)^(1/3) = 6.327399015078757. Inputs are
  the *smoothed* OBLdepth/kOBL (mom6_kpp.csv), z_inter, uStar, the MXL_SW
  surfBuoyFlux = buoyFlux(1) − buoyFlux(floor(kOBL)+1) (F90:848-853) and the
  incoming Kt/Ks/Kv.

All KPP laws live in `mom6_kpp_templates.esm` (imported by reference) so the
same pointwise laws can serve other branches (e.g. the match-gradient shape
coefficients, the Ekman/Monin-Obukhov OBL limits) and the stage-3 coupled
model.

### Transcription notes

- **OBL depth via the per-layer quadratic, in closed form.** The Fortran
  builds the crossing-layer quadratic with `cvmix_math_poly_interp` (3×3
  matrix inverse) and finds its root by Newton iteration
  (`cubic_root_find`, `cvmix_math.F90:260-294`). The .esm writes the
  quadratic coefficients c/b/a (templates `kpp_quad_interp_c/b/a`) and the
  closed-form root `(−b ± √(b²−4ac))/(2c)` selected to fall in the bracket
  (`kpp_quad_root`); against all 156 traced columns the closed form
  reproduces Newton to ≤ 4.6e-13 relative (the Fortran's TOL is 1e-12, and
  the tests pin at rel 1e-12). The degenerate branches (negative
  discriminant, c = 0, no root in the bracket) return the Newton starting
  point 0.5·(x1+x2) and are never the selected crossing.
- **Surface-layer means in closed form, bit-identical.** The Fortran
  accumulates delH(k, k′) = min(max(SLdepth(k) + z_inter(k′), 0), dz(k′))
  layer by layer; the .esm writes the closed form over all 75 layers. The
  crossing-layer sum reaches SLdepth exactly (Sterbenz), so the closed form
  is bit-identical to the running hTot and delH ≡ 0 beyond the surface
  layer. The `integral` operator would be the natural spelling (the means
  are ∫₀^SLdepth T dz / SLdepth) but the CLI 0.7.0 evaluator cannot evaluate
  it for algebraic expressions, so the equivalent closed form is used.
- **Parameter arrays are only gathered at provably in-range indices.** The
  esm build checks const-node array bounds statically
  (`E_TREEWALK_CONSTARRAY_OOB`), so `h/dz/T/S/U_H/V_H` are read only at
  plain or `max(i−1, 1)`-clamped indices and `buoyFlux` at 1 and i+1 (i ≤
  75); every computed-index gather (the quadratic arrays, `root(kfirst)`,
  z_cell/z_inter at symbolic k, ...) is on an unknown array. This shapes the
  spelling: N2 is computed on the 75 layers and N on the 76 interfaces with
  the i = N_IFACE zero branch (`ifelse(i == N_IFACE, 0,
  sqrt(max(N2(i), 0)))`), which avoids an out-of-range `index(h, 76)`
  inside an untaken branch.
- **The Fortran's operation order is preserved wherever the pins can see
  it**: the bulk-Richardson numerator/denominator, the evaluate_cubic
  summation order of G(σ) = (σ + (−2)·σ²) + σ³, and uStar³ spelled
  (uStar·uStar)·uStar. The MXL stability functions and the shape coefficients
  a_m/c_m/a_s/c_s enter as exact-double parameters (CVMix init constants),
  not recomputed.
- **Smoothing is a boundary, not a law.** KPP_smooth_BLD (1715-1834) is the
  N_SMOOTH=3 horizontal Laplacian of the OBL depth — an inherently
  horizontal operation that a single-column component cannot reproduce. The
  pins therefore cover the *pre-smoothing* OBLdepth/kOBL (mom6_kpp_obl.csv)
  and the tests consume the smoothed values as inputs (mom6_kpp.csv).

### Pinned rows and coverage

| Test (BLD) | mom6_kpp_obl.csv row | kfirst | OBLdepth [m] | kOBL |
|---|---|---|---|---|
| `cell782_step3600_quad_k3` | (7, 82, 3600) | 3 | 6.613560999075406 | 3.75 |
| `cell883_step3600_quad_k3` | (8, 83, 3600) | 3 | 5.974325344745279 | 3.25 |
| `cell784_step10800_two_point_k2` | (7, 84, 10800) | 2 | 5.0 (floored) | 2.75 |
| `cell1184_step3600_kfirst1` | (11, 84, 3600) | 1 | 5.0 (floored) | 3.25 |
| `cell684_step10800_kfirst1_shallow` | (6, 84, 10800) | 1 | 5.0 (floored) | 2.75 |

Each BLD test pins OBLdepth/kOBL, BulkRi/deltaU2/N2 at layers 1..kOBL+3
(6 or 5 layers), spot z_inter/z_cell values at the same range plus the
three deep layers 73..75, and kfirst (34..39 assertions, rel 1e-12).
`kfirst` is *not a dump column*: it comes from the Python reference
reimplementation (`stage1/analysis/kpp_ref.py`) of the OBL_depth_low
branch selection and is labeled "derived (not a dump column)" in the test
descriptions. The kfirst = 1 and 2 rows exercise the `−z_cell(1)` branch
and the two-point (zero-slope) quadratic; the kfirst = 3 rows exercise the
full slope-preserving quadratic. The kOBL fraction is pinned on both sides
of the crossing: .75 (OBL below the layer-3 center, kfirst = 3 row) and
.25 (OBL above it), plus the shallow kfirst = 1 rows. The test tuples are
extracted by `stage1/analysis/extract_mom6_kpp_test_rows.py` (dump
directory as its argument; the pinned rows and the keep-ranges live in the
script's docstring).

| Test (KPPCoefs) | mom6_kpp.csv row | OBLdepth_s [m] | kOBL_s | kwup vs ktup |
|---|---|---|---|---|
| `cell782_step3600_kOBL325` | (7, 82, 3600) | 6.1872106097070425 | 3.25 | kwup = ktup + 1 |
| `cell882_step3600_kOBL375` | (8, 82, 3600) | 6.309964392417778 | 3.75 | kwup = ktup |
| `cell684_step10800_kOBL275` | (6, 84, 10800) | 5.001656526745087 | 2.75 | kwup = 2 |

Each coefficient test pins the 7 asserted columns
Kd_heat/Kd_salt/Kv/NLT_heat/Kt_out/Ks_out/Kv_out at interfaces
2..kwup+1 — the KPP-active region plus the first restoration interface —
and at the deep interface 76 (21..28 assertions, rel 1e-12), covering both
the .25 and .75 kOBL branches of the zero restoration at kwup and the
kwup = 2 minimal case. (The 0-below-kwup and Kt_out = Kt_in bulk of the
76-interface profiles is not asserted; the deep interface pins the
fallback-to-incoming law once per row.)

### Coverage limitations (this data window)

| Unpinned branch | Why |
|---|---|
| N_SMOOTH=3 horizontal smoothing | inherently horizontal (see above); pre-smoothing values pinned, smoothed values consumed as inputs |
| kfirst = N_LAYER+1 (OBL_limit) branch of OBL_depth_low | no traced column has BulkRi ≤ RI_CRIT in every layer |
| kfirst = 2 interpolated value | below the MINIMUM_OBL_DEPTH = 5 m floor in every traced column; only the floored OBLdepth can be pinned |
| unstable (surfBuoyFlux < 0) branches | all 156 traced columns are stable; the stability functions' convective branch and the nonzero non-local term (NLT ≡ 0 in all pins) are unpinned |
| Langmuir paths | l_LMD_ws = false, USE_KPP_LT_K/VT2 = False: the lStokesMOST and Langmuir-entrainment branches of coeffs_kpp_low are not transcribed |
| MINIMUM_VT2 floor | not reached in the traces (Vt² stays above 1e-10 m² s⁻²) |
| Ekman / Monin-Obukhov OBL limits | COMPUTE_EKMAN=False, so only the Monin-Obukhov limit enters OBL_limit; it is transcribed (the min(−root, OBL_limit) cap), and no traced OBLdepth reaches the cap |
| interface-indexed N profile | not asserted (only the layer-indexed N2); N is an intermediate of the BulkRi/Vt² chain, and the i = N_IFACE zero branch is the only interface-specific behavior |
| KPPCoefs below-kwup zeros and the Kt_out = Kt_in fallback bulk | the informative core (interfaces 2..kwup+1 + the deep interface 76) is pinned; the identically-zero and copy-of-input stretches are not |

## Surface forcing mapping

All fluxes/forces are pass-throughs of the 25 imported fields except:

| Output | Law | Fortran |
|---|---|---|
| `sw` | `sw_vis_dir + sw_vis_dif + sw_nir_dir + sw_nir_dif` | :597 |
| `latent` | `(latent_evap_diag + latent_fprec_diag) + latent_frunoff_diag − frunoff_glc·Lf` — assembled from the diagnostics below (bit-identical to the flat `q_flux·Lv − fprec·Lf − frunoff·Lf − frunoff_glc·Lf`; IEEE `x − y ≡ x + (−y)`) | :572-596 |
| `latent_evap_diag` | `q_flux·Lv` | :594-595 |
| `latent_fprec_diag` | `−fprec·Lf` | :575-576 |
| `latent_frunoff_diag` | `−frunoff·Lf` | :581-582 |
| `vprec` | `vprec_restore − net_fw_global_adjust` | :455-486, :694-714 |
| `salt_flux_o` | `salt_flux + salt_restore_flux` | :430-447, :686-689 |
| `taux` | `masked_pair_average(u_flux, u_flux_east, mask_east) = (u_flux + mask_east·u_flux_east)/(1 + mask_east)` — AGRID u-point average of the h-point stresses of cells (i,j) and (i+1,j) | :943-946 |
| `tauy` | `masked_pair_average(v_flux, v_flux_north, mask_north) = (v_flux + mask_north·v_flux_north)/(1 + mask_north)` — AGRID v-point average of (i,j) and (i,j+1) | :948-951 |
| `tau_mag` | `GUST_CONST + sqrt(u_flux² + v_flux²)` — built from the h-point stresses, NOT from taux/tauy | :956-957 |
| `ustar` | `sqrt(tau_mag/RHO_0)` | :960-961 |
| `omega_w2x` | `atan2(v_flux, u_flux)` | :962 |
| `net_mass_src` | the nine-term sum `lprec_o + fprec_o + vprec + lrunoff_o + frunoff_o + lrunoff_glc + frunoff_glc + evap + seaice_melt_o` | `MOM_forcing_type.F90:2748-2800` |

The dump column names for the mechanical fields are `tau_mag`, `ustar`, and
`p_surf_full` (the model records `fluxes%tau_mag = forces%tau_mag` etc. on
every row), so the `_o`-suffixed duplicate outputs of an earlier draft
(`tau_mag_o`/`ustar_o`/`p_surf_full_o`, values identical to the unsuffixed
names on every row) were removed in the factoring cleanup; the tests pin the
unsuffixed names.

**Inputs that are not IOB imports.** Three fields are ocean-state- or
global-driven and cannot be recovered from the local import tuple; they enter
as parameters, and the tests set them from the dump:

- `salt_restore_flux` — the salinity-restoring salt flux of the
  `SRESTORE_AS_SFLUX=True` path (:430-447), driven by `SSS` vs the restoring
  field; the dump records it as `salt_flux_o − salt_flux_in` (the mediator
  sends `salt_flux = 0`, so `salt_flux_o` is the restoring flux itself).
- `net_fw_global_adjust` — the area-mean net fresh-water flux subtracted by
  `ADJUST_NET_FRESH_WATER_TO_ZERO` (:694-714), a global sum over the ocean;
  the dump records it as `−vprec`, constant across the cells of each step.
- `vprec_restore` — the restoring-as-freshwater flux of the
  `SRESTORE_AS_SFLUX=False` path; identically 0 here (the sflux path skips
  it), included so the general `vprec` law reads as in the Fortran.
- `heat_added_override` — the `hflx_adj` data override (:1074-1077); absent
  in the reference run, so `heat_added ≡ 0` in the dump.

**Configuration artifacts pinned by the tests.** `u10_sqr_o ≡ 0` although
the import `u10_sqr` is nonzero: the reference run never allocates
`fluxes%u10_sqr` (:644-645), modelled as `u10_sqr_allocated = 0`.
`lrunoff_glc`/`frunoff_glc ≡ 0`: the mediator does not send glacier runoff in
this configuration.

**tauy uses the north neighbor, not the east.** The first draft averaged
`tauy` with the *east* neighbor's `v_flux`; the dump pins showed the Fortran's
stencil is the v-point average with the *north* neighbor (i, j+1) (:948-951).
The j=84 traced rows have ocean north neighbors that are not traced (their
`tauy ≠ v_flux` with no available `v_flux_north`), so no j=84 row can pin
`tauy`; all six test rows have traced north neighbors.

## Testing

Run with

    esm validate stage2/mom6_templates.esm stage2/mom6_kpp_templates.esm \
        stage2/mom6_eos.esm stage2/mom6_kpp.esm stage2/mom6_surface_forcing.esm
    esm test stage2/mom6_eos.esm stage2/mom6_kpp.esm stage2/mom6_surface_forcing.esm

All `tests` blocks are `time_span {0,0}` instantaneous pins; every expected
value is a real dump row, identified in the test descriptions by
step/cell/layer. The dump writes `es26.17` and the reference run is r8, so
the transcriptions are pinned to round-off at rel 1e-12 (the measured
f64-vs-dump residuals are ≤ 1.4e-16 for the EOS and forcing transcriptions;
the KPP closed-form quadratic reproduces the Fortran's Newton iteration to
≤ 4.6e-13 relative, still inside the 1e-12 tolerance). The run configuration
is documented in the file metadata (SRESTORE_AS_SFLUX,
ADJUST_NET_FRESH_WATER_TO_ZERO, ADJUST_NET_SRESTORE_TO_ZERO, wind_stagger =
AGRID, GUST_CONST = 4.0e-4 Pa, RHO_0 = 1035 kg m-3, Lf = 3.337e5,
Lv = 2.501e6 J/kg, and the KPP configuration listed in the KPP section, all
unit rescale factors 1.0).

| File | Tests | Assertions | What is pinned |
|---|---|---|---|
| `mom6_eos.esm` | 15 | 45 | rho/drho_dT/drho_dS at 12 surface/mid/deep states of 4 cells over steps 1, 3, 5, 6, 12: the T/S/p range extremes of the dump (T −1.99..2.04 °C, S 33.25..34.66 ppt, p 9.97e4..7.85e6 Pa), min/max rho, drho_dT (most- and least-negative), drho_dS (min/max) — every template exercised on every row |
| `mom6_kpp.esm` | 8 | 262 | See the KPP section: 5 OBL rows (BulkRi/deltaU2/N2 at layers 1..kOBL+3, z_inter/z_cell spots + 3 deep layers, OBLdepth/kOBL, kfirst derived) and 3 coefficient rows (the 7 asserted KPP output columns Kd_heat/Kd_salt/Kv/NLT_heat/Kt_out/Ks_out/Kv_out at interfaces 2..kwup+1 + the deep interface 76; NLT_salt and the Kt/Ks/Kv_in inputs are not asserted), rel 1e-12 |
| `mom6_surface_forcing.esm` | 6 | 234 | All 39 output columns on 6 rows covering: ocean/land east neighbors (taux averaged vs taux = u_flux) with ocean north neighbors (tauy averaged); the runoff/enthalpy branches (rows with lrunoff/frunoff/hrofl/hrofi nonzero); the sea-ice melt branches (seaice_melt_heat/seaice_melt nonzero); evaporation and snow on every row; vprec via net_fw_global_adjust (steps 1, 3, 5, 12); salt_flux_o via salt_restore_flux; the smallest and largest tau_mag and salt_restore_flux of the traced rows |

### Coverage limitations (this data window)

- **EOS**: the dump window only spans −2..2.04 °C, 33.25..34.66 ppt,
  9.97e4..7.85e6 Pa (an Antarctic cell block in early summer); the law's
  validity range (−2..40 °C, 0..40 PSU, 0..1e8 Pa) is not exercised at the
  extremes, and no density-anomaly/spec-vol variants
  (`density_anomaly_elem_Wright_full`, `spec_vol_elem_Wright_full`) are
  dumped — only the in-situ law used by `KPP_compute_BLD`.
- **Forcing**: `lprec ≡ 0` in the whole dump (only snow falls in the
  window), so the nonzero `lprec_o` branch and the `heat_content_lprec`
  (hrain) / `heat_content_fprec` (hsnow) / `heat_content_cond` (hcond)
  columns are pinned only at 0; `frunoff_glc_import ≡ 0` so the
  `−frunoff_glc·Lf` latent term is pinned only in its zero branch;
  `u10_sqr_o ≡ 0` is a configuration artifact (see above); `taux`/`tauy`
  need the east/north neighbor imports, and the j=84 rows' north neighbors
  are untraced ocean cells, so `tauy` is unpinned there (all six test rows
  use traced neighbors); the vprec and salt_flux_o laws are pinned through
  the state-driven inputs above, not through their full Fortran sources.
- The B-grid/C-grid `wind_stagger` alternatives and the
  `ADJUST_NET_FRESH_WATER_BY_SCALING`/`ADJUST_NET_SRESTORE_BY_SCALING` paths
  are not transcribed (only the AGRID + area-mean-adjustment configuration
  of the reference run is implemented).
- **KPP**: the unpinned branches of the KPP section's coverage table — the
  unstable stability-function and non-local branches, the kfirst = 2
  interpolation below the 5 m floor, the kfirst = N_LAYER+1 OBL_limit
  branch, the Langmuir paths, the MINIMUM_VT2 floor, and the horizontal
  N_SMOOTH smoothing (see the KPP section for how the smoothed values enter
  as inputs).

## Next steps

1. A wider stage-1 window (warmer/more saline cells, unstable
   surfBuoyFlux, liquid precipitation, glacier runoff) to cover the
   untested branches above.
2. Mount into `eqesm.esm` with `earthsciintegrator` (stage 3) and compare
   against full-CESM runs.
