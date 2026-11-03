# MOM6 stage 2 — equation of state and surface forcing

Stage-2 implementation of two MOM6 boundaries of the MOM6 lane
(`stage1/runs/MOM6_LANE.md`): the Wright 1997 full-range equation of state
(boundary B1, exercised by `KPP_compute_BLD`) and the surface forcing
import-to-fluxes/forces mapping (boundary I), as esm 1.3.0 files with `tests`
blocks pinned to the stage-1 instrumentation dumps
(`stage1/dumps/mom6/mom6_eos.csv`, `mom6_sfc_forcing.csv`; 12 traced PET-0
cells × 12 coupling steps of the reference run `C_JRA_t62_mom6_inst`, MOM6
`dev/ncar_260909 @301879336`, fork `stage1/EqESM-cesm-mom6` branch
`instrument/mom6-dumps`).

Sources:
`components/mom/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90`,
`components/mom/MOM6/config_src/drivers/nuopc_cap/mom_surface_forcing_nuopc.F90`,
`components/mom/MOM6/src/core/MOM_forcing_type.F90`,
`components/mom/MOM6/src/framework/mom6_dump_mod.F90`.

## File map

| File | Contents | Fortran source |
|---|---|---|
| `mom6_templates.esm` | Expression-template library: `wright_alpha0`, `wright_p0`, `wright_lambda` (the three shared polynomials), `wright_density`, `wright_drho_dT`, `wright_drho_dS` | `MOM_EOS_Wright_full.F90`: density_elem_Wright_full 75-91, calculate_density_derivs_elem_Wright_full 179-204, constants 23-37 |
| `mom6_eos.esm` | Model `WrightEOS` — rho(T, S, p), drho/dT, drho/dS via the templates (imported by reference), `I_denom2`; 15 dump-pinned tests | `MOM_EOS_Wright_full.F90` as above; the EOS calls themselves are `calculate_density`/`calculate_density_derivs_scalar` in `MOM_CVMix_KPP.F90:1437` |
| `mom6_surface_forcing.esm` | Model `SurfaceForcing` — the 25-field IOB import → fluxes/forces mapping, 42 outputs, 6 dump-pinned tests | `mom_surface_forcing_nuopc.F90`: convert_IOB_to_fluxes 258-734, convert_IOB_to_forces 739-1050; `MOM_forcing_type.F90`: set_net_mass_forcing/get_net_mass_forcing 2748-2800 |

## Units

| Quantity | Symbol | Units |
|---|---|---|
| Potential temperature | `T` | degC |
| Salinity | `S` | ppt |
| Pressure | `p`, `p0`, `p_surf_full`, `p_surf_full_o` | Pa |
| Density | `rho` | kg/m³ |
| Density derivatives | `drho_dT`, `drho_dS` | kg/m³/degC, kg/m³/ppt |
| Specific volume (alpha_0) | `al0` | m³/kg |
| Sound-speed-squared offset | `lam` | m²/s² |
| EOS coefficients | `a0..a2`, `b0..b5`, `c0..c5` | as in the Fortran's units comments (F90:23-37) |
| Heat fluxes | `sw`, `lw`, `latent`, `sens`, `seaice_melt_heat_o`, `heat_added`, `heat_content_*`, `latent_*_diag` | W/m² |
| Water fluxes | `evap`, `lprec_o`, `fprec_o`, `vprec`, `lrunoff_o`, `frunoff_o`, `lrunoff_glc`, `frunoff_glc`, `seaice_melt_o`, `salt_flux_o`, `salt_flux_in`, `net_mass_src` | kg/m²/s |
| Wind stresses | `u_flux`, `v_flux`, `taux`, `tauy`, `tau_mag`, `tau_mag_o`, `GUST_CONST` | Pa |
| Friction velocity | `ustar`, `ustar_o` | m/s |
| Wind angle | `omega_w2x` | rad |
| Wind speed squared | `u10_sqr`, `u10_sqr_o` | m²/s² |
| Masks / fractions | `mask_east`, `mask_north`, `ice_fraction`, `ice_fraction_o`, `u10_sqr_allocated` | 1 |
| Latent heats | `Lf`, `Lv` | J/kg |
| Reference density | `RHO_0` | kg/m³ |

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

## Surface forcing mapping

All fluxes/forces are pass-throughs of the 25 imported fields except:

| Output | Law | Fortran |
|---|---|---|
| `sw` | `sw_vis_dir + sw_vis_dif + sw_nir_dir + sw_nir_dif` | :597 |
| `latent` | `q_flux·Lv − fprec·Lf − frunoff·Lf − frunoff_glc·Lf` | :572-596 |
| `latent_evap_diag` | `q_flux·Lv` | :594-595 |
| `latent_fprec_diag` | `−fprec·Lf` | :575-576 |
| `latent_frunoff_diag` | `−frunoff·Lf` | :581-582 |
| `vprec` | `vprec_restore − net_fw_global_adjust` | :455-486, :694-714 |
| `salt_flux_o` | `salt_flux + salt_restore_flux` | :430-447, :686-689 |
| `taux` | `(u_flux + mask_east·u_flux_east)/(1 + mask_east)` — AGRID u-point average of the h-point stresses of cells (i,j) and (i+1,j) | :943-946 |
| `tauy` | `(v_flux + mask_north·v_flux_north)/(1 + mask_north)` — AGRID v-point average of (i,j) and (i,j+1) | :948-951 |
| `tau_mag`, `tau_mag_o` | `GUST_CONST + sqrt(u_flux² + v_flux²)` — built from the h-point stresses, NOT from taux/tauy | :956-957 |
| `ustar`, `ustar_o` | `sqrt(tau_mag/RHO_0)` | :960-961 |
| `omega_w2x` | `atan2(v_flux, u_flux)` | :962 |
| `net_mass_src` | the nine-term sum `lprec_o + fprec_o + vprec + lrunoff_o + frunoff_o + lrunoff_glc + frunoff_glc + evap + seaice_melt_o` | `MOM_forcing_type.F90:2748-2800` |

`tau_mag_o`/`ustar_o`/`p_surf_full_o` equal `tau_mag`/`ustar`/`p_surf_full`
(the dump records `fluxes%tau_mag = forces%tau_mag` etc. on every row).

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

    esm validate stage2/mom6_templates.esm stage2/mom6_eos.esm stage2/mom6_surface_forcing.esm
    esm test stage2/mom6_eos.esm stage2/mom6_surface_forcing.esm

All `tests` blocks are `time_span {0,0}` instantaneous pins; every expected
value is a real dump row, identified in the test descriptions by
step/cell/layer. The dump writes `es26.17` and the reference run is r8, so
the transcriptions are pinned to round-off (rel 1e-12; the measured f64-vs-
dump residuals are ≤ 1.4e-16). The run configuration is documented in the
file metadata (SRESTORE_AS_SFLUX, ADJUST_NET_FRESH_WATER_TO_ZERO,
ADJUST_NET_SRESTORE_TO_ZERO, wind_stagger = AGRID, GUST_CONST = 4.0e-4 Pa,
RHO_0 = 1035 kg m-3, Lf = 3.337e5, Lv = 2.501e6 J/kg, all unit rescale
factors 1.0).

| File | Tests | Assertions | What is pinned |
|---|---|---|---|
| `mom6_eos.esm` | 15 | 45 | rho/drho_dT/drho_dS at 12 surface/mid/deep states of 4 cells over steps 1, 3, 5, 6, 12: the T/S/p range extremes of the dump (T −1.99..2.04 °C, S 33.25..34.66 ppt, p 9.97e4..7.85e6 Pa), min/max rho, drho_dT (most- and least-negative), drho_dS (min/max) — every template exercised on every row |
| `mom6_surface_forcing.esm` | 6 | 252 | All 42 output columns on 6 rows covering: ocean/land east neighbors (taux averaged vs taux = u_flux) with ocean north neighbors (tauy averaged); the runoff/enthalpy branches (rows with lrunoff/frunoff/hrofl/hrofi nonzero); the sea-ice melt branches (seaice_melt_heat/seaice_melt nonzero); evaporation and snow on every row; vprec via net_fw_global_adjust (steps 1, 3, 5, 12); salt_flux_o via salt_restore_flux; the smallest and largest tau_mag and salt_restore_flux of the traced rows |

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

## Next steps

1. `stage2/mom6_kpp.esm` — the CVMix KPP boundary-layer component
   (boundary B2, `MOM_CVMix_KPP.F90` + `MOM6_KPP`), which consumes the EOS
   component above (`calculate_density` calls in `KPP_compute_BLD`).
2. A wider stage-1 window (warmer/more saline cells, liquid precipitation,
   glacier runoff) to cover the untested branches above.
3. Mount into `eqesm.esm` with `earthsciintegrator` (stage 3) and compare
   against full-CESM runs.
