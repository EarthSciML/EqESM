# MOM6 lane — stage-1 instrumentation (W1)

Status: **run complete, dumps validated** (2026-11-02). Fork
`stage1/EqESM-cesm-mom6` (branch `instrument/mom6-dumps`), case
`C_JRA_t62_mom6_inst`, dumps in `stage1/dumps/mom6/` (gitignored).
Run record: `stage1/runs/mom6_instrument.json`.

## Run configuration (reference case `C_JRA_t62`, rerun as `C_JRA_t62_mom6_inst`)

Chosen case: **`C_JRA_t62`** (not `C_JRA_tx10`). Reasons: (1) `C_JRA_tx10`
cannot run at all — the `mom/tx10deg` testmod references six inputdata files
(`ocean_hgrid_tx10deg_20260219.nc`, `ESMF_mesh_tx10deg_20260219.nc`, …) that
exist on no inputdata server (bugs.md #2); (2) even if that testmod were
fixed, the only alternative `C_JRA_tx10`-style grid is the full `tx0.1v3`
10-km ocean (tens of thousands of PEs-hours/day), while `C_JRA_t62` at the
2° `tx2_3v2` tripole ocean completed in 7:53 wall on 2 nodes (72 PEs) —
about two orders of magnitude cheaper, and already the project's validated
C-class reference (PROGRESS.md "C-class run — SUCCESS").

| Item | Value |
|---|---|
| Compset | `C_JRA` = `2000_DATM%JRA-1p5-2023_SLND_DICE%SSMI_MOM6_DROF%JRA-1p5-2023_SGLC_SWAV_SESP` |
| Grid | `T62_t232` (atm T62 1.9x2.5; ocn/ice **tx2_3v2** tripole, NIGLOBAL=540 × NJGLOBAL=480, **NK=75** z\* layers) |
| Run | startup, 1 day (STOP_OPTION=ndays, STOP_N=1), 2000-01-01 |
| Coupling | ATM/OCN/ICE_NCPL=24 → coupler dt 3600 s; THERMO_SPANS_COUPLING=False |
| MOM6 dt | **DT=1800 s** (baroclinic dynamics), **DT_THERM=3600 s** (thermodynamics), 2 dynamics steps per coupling step |
| PEs | 72 total (NTASKS_OCN=36, NTASKS_CPL=36, 2 nodes, NTHRDS=1), `srun --mpi=pmix -n N` launcher |
| EOS | **`EQN_OF_STATE = "WRIGHT_FULL"`** (Wright 1997 full-range fit to UNESCO; the MOM6 default — not set explicitly in MOM_input), `USE_PSURF_IN_EOS=True`, TFREEZE_FORM="LINEAR" (S0=0, dTf/dS=-0.054 °C/ppt) |
| Density ref | RHO_0 = 1035 kg/m3 (Boussinesq), C_P = 3992 J/kg/°C, G_EARTH = 9.8 m/s2 |
| KPP/CVMix | `USE_KPP=True` (CVMix KPP, Large et al. 1994): MATCH_TECHNIQUE=MatchGradient, INTERP_TYPE2=LMD94, N_SMOOTH=3 (OBL Laplacian smoothing), RI_CRIT=0.3, VON_KARMAN=0.4, CS=98.96, CS2=6.32739901508, SURF_LAYER_EXTENT=0.1, MINIMUM_OBL_DEPTH=5.0 m, COMPUTE_MONIN_OBUKHOV=True, ENHANCE_DIFFUSION=False, KPP_IS_ADDITIVE=False (KPP replaces other diffusivity wherever non-zero), NLT_SHAPE=CVMix (per MATCH_TECHNIQUE), APPLY_NONLOCAL_TRANSPORT=True, KPP_SHORTWAVE_METHOD=MXL_SW. Other mixing: USE_LMD94=True (shear, NU_ZERO=0.0025, RI_ZERO=0.5, N_SMOOTH_RI=1), USE_CVMix_TIDAL=True (Jayne), USE_CVMix_DDIFF=True (convection). Langmuir/Stokes: all off (STOKES_MIXING=False, LT_K_ENHANCEMENT=False). |
| Misc | FRAZIL=True, AUTO_MASKTABLE=False (workaround, bugs.md #3), MAXTRUNC=1e6, RESTART_CONTROL=3, STOCH_EOS=True (stochastic EOS perturbations used in the pressure-gradient force; the EOS law itself is deterministic — see EOS boundary note) |
| Units | **All MOM6 rescale factors are 1.0 in this run** (no `*_RESCALE` parameters in MOM_input; defaults in `MOM_unit_scaling.F90`): internal units = SI (T °C, S ppt, p Pa, rho kg m-3, thickness m via `H_TO_M=1.0`). Dumps are written in physical units by explicitly multiplying by the conversion factors (`US%Z_to_m`, `US%L_T_to_m_s`, `GV%H_to_m`, `US%R_to_kg_m3`, `EOS%C_to_degC`, `EOS%S_to_ppt`, `EOS%RL2_T2_to_Pa`), so the dump code stays correct if a future run rescales. |

## Boundary map (3 boundaries, first bounded set)

All dumps: PET 0 only, first `dump_first_cells` mask2dT cells of PET 0's
compute domain in (j,i) loop order (the same order MOM6's (i,j) physics
loops visit), first `dump_first_steps` steps, gated by `dump_mom6=.false.`
default, OMP-critical writes, one CSV per boundary with
`# description =`/`# columns =` headers naming every column with units and
the F90 file:line it comes from.

### B1 — Equation of state (D boundary, pointwise law) → `mom6_eos.csv`

- **Law:** Wright 1997 full-range (Wright, JPO 1997, "A new salinity
  equation of state"; fit to UNESCO EOS-80 for -2<θ<40 °C, 0<S<40 psu,
  0<p<1e8 Pa):
  ρ = (p + p0) / (λ + α0 (p + p0)), with
  α0 = a0 + a1T + a2S;  p0 = b0 + b4S + T(b1 + T(b2 + b3T) + b5S);
  λ = c0 + c4S + T(c1 + T(c2 + c3T) + c5S),
  and the 15 constants a0..c5 from `MOM_EOS_Wright_full.F90:24-38`.
- **Dump site:** `MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90:1437`
  (inside `KPP_compute_BLD`, the per-column loop) — the model's in-situ
  density evaluation for the bulk Richardson number: per layer k of each
  traced column, T (potential, °C), S (ppt), p (Pa, = surface pressure +
  hydrostatic weight of layers above, i.e. pressure at the top interface of
  layer k) → ρ (kg m-3) plus drho/dT, drho/dS (kg m-3 °C-1, kg m-3 ppt-1)
  from `calculate_density_derivs_scalar` (MOM_EOS.F90:855).
  Also dumps `p_surf` and the column of interface pressures.
- **Scale factors:** `EOS%C_to_degC`, `EOS%S_to_ppt`, `EOS%RL2_T2_to_Pa`,
  `EOS%kg_m3_to_R` (all 1.0 in this run; applied explicitly).
- **Note:** `STOCH_EOS=True` only perturbs the density used in the
  pressure-gradient force (`MOM_stoch_eos.F90`); the law itself (and this
  dump) is the deterministic Wright_FULL EOS.

### B2 — KPP boundary layer (CVMix, D boundary) → `mom6_kpp.csv` + `mom6_kpp_obl.csv`

- **Entry/exit per traced column:**
  - `KPP_compute_BLD` (`MOM_CVMix_KPP.F90:1071`): computes the OBL depth.
    Dump `mom6_kpp_obl.csv` after `CS%kOBL` is set (line 1574, before the
    smoothing step): inputs per column — surface friction velocity u* (m/s),
    surface buoyancy flux B0 (m2 s-3) and its profile buoyFlux(:,k), layer
    thicknesses h (m), layer T/S (°C/ppt), layer u/v (m/s), N2 profile
    (s-2), Coriolis f, surface pressure; outputs — OBL depth (m), kOBL
    (fractional level), surface-layer averages (T, S, u, v), the Monin-
    Obukhov limit flag.
  - `KPP_calculate` (`MOM_CVMix_KPP.F90:740`): applies CVMix to produce
    diffusivities/viscosities and non-local transport. Dump
    `mom6_kpp.csv` at the end of the per-column loop (line ~1050, after the
    Kt/Ks/Kv update): inputs — same column state + pre-KPP interior Kd/Ks/Kv
    (m2/s, captured at column entry), OBL depth (m, smoothed), kOBL, u*,
    surface buoyancy flux (m2 s-3, per KPP_SHORTWAVE_METHOD), z_inter/z_cell
    (m, negative in ocean); outputs — CVMix Kdiffusivity heat/salt (m2 s-1),
    Kviscosity (m2 s-1), nonLocalTrans heat/salt (nondim) and the final
    Kt/Ks/Kv after the KPP update, plus Langmuir diagnostics (La_SL,
    LangEnhK; zero/1.0 with Langmuir off).
- **Units:** all dumped in MKS via `US%Z_to_m`, `US%s_to_T`, `GV%H_to_m`,
  `GV%HZ_T_to_m2_s`, `buoy_scale`, `US%L_T_to_m_s`.

### B3 — Surface forcing import (I boundary, per coupling step) → `mom6_sfc_forcing.csv`

- **Path:** CMEPS mediator fields → `mom_import` (`mom_cap_methods.F90:80`,
  ESMF importState → `ice_ocean_boundary_type`, incl. `med2mod_areacor`
  scaling and the cos_rot/sin_rot rotation of taux/tauy, lines 158-175) →
  `update_ocean_model` (`mom_ocean_model_nuopc.F90:570-608`:
  `convert_IOB_to_forces` + `convert_IOB_to_fluxes`) → MOM6 `forces`/`fluxes`.
- **Dump site:** `mom_ocean_model_nuopc.F90` right after
  `set_net_mass_forcing` (line ~614), i.e. once per coupling step with the
  final fluxes+forces: imported coupler fields as consumed by MOM6 —
  p (Pa, Sa_pslv), sw_flux_vis_dir/dif/nir_dir/dif (W/m2, Foxx_swnet_*),
  lw_flux (Foxx_lwnet), t_flux sens (Foxx_sen), q_flux evap (Foxx_evap),
  lprec/fprec (Faxa_rain/snow, kg m-2 s-1), lrunoff/frunoff (Foxx_rofl/rofi,
  kg m-2 s-1), u_flux/v_flux (Foxx_taux/tauy, Pa, rotated to local
  coordinates), ice_fraction, plus enthalpy fields (hrain/hsnow/hrofl/hrofi/
  hevap/hcond, W/m2) — and the resulting MOM6 forcing: fluxes%sw/sens/
  latent/lw (W/m2), evap/lprec/fprec/lrunoff/frunoff (kg m-2 s-1),
  salt_flux (kg m-2 s-1), vprec, seaice_melt, p_surf_full (Pa),
  heat_content_* (W/m2), forces%taux/tauy (Pa, on the C-grid u/v points
  averaged to the T point), forces%ustar (m/s), tau_mag (Pa).
- **Units:** MKS (all conversion factors are 1.0 in this run; applied
  explicitly via `US%kg_m2s_to_RZ_T`, `US%W_m2_to_QRZ_T`, etc.).
- **Cadence:** one row per traced cell per coupling step (24 steps/day);
  `dump_first_steps` limits to the first 12.

## MOM6's own diagnostic output (PLAN.md §4.3 item 2) — checked, cannot supply these boundaries

The C_JRA `diag_table` writes `mom6.h.*` files at **monthly** frequency; the
1-day reference run produced **no history files at all** (only restarts and
`ocean.stats.nc`). Even with the cadence fixed, the requested boundaries are
not in the diagnostic list as instantaneous per-cell values: `KPP_OBLdepth`
is available (as a time-mean), but per-layer KPP Kd/Kv/non-local transport,
drho/dT-drho/dS, and the forcing import→flux mapping are not
(`available_diags.000000` in the run dir confirms). **Hooks are required**;
they are additive (dump calls only, solution-neutral).

## Instrumentation (fork `stage1/EqESM-cesm-mom6`, branch `instrument/mom6-dumps`)

- Fork = git worktree of the pristine `stage1/cesm` (e59f251,
  cesm3.0-alphabranch), submodules pinned per `stage1/runs/README.md`
  (cime cime6.7.1, mom mi_260909 → MOM6 dev/ncar_260909 @301879336 →
  CVMix-src @fce4221, GSW-Fortran @29e64d6; cdeps 1.0.112, cmeps 1.1.73,
  FMS fi_260310, parallelio pio2_6_8, share 1.1.22, ccs_config 1.0.92).
  Verified `git submodule status` on the fork and on the pristine checkout.
- New module `MOM6/src/framework/mom6_dump_mod.F90` (auto-picked-up: the
  MOM6 Filepath list in `components/mom/cime_config/buildlib` includes
  `src/framework`). Namelist gate: MOM6 reads its own parameters from
  `MOM_override` (user_nl_mom) at runtime — `dump_mom6`, `dump_dir`,
  `dump_first_cells`, `dump_first_steps` are read with `get_param` in
  `dump_init`, so **no CIME namelist-definition change is needed** (unlike
  CICE, MOM6 has no buildnml hard-coded group list; MOM_input comes from
  `MOM_input.json` + user_nl_mom → MOM_override).
- Wiring (additive dump calls only): (1) `MOM_CVMix_KPP.F90` KPP_compute_BLD
  after line 1574 (OBL dump) and KPP_calculate at the end of the column loop
  (KPP dump + EOS dump of the KPP-computed densities; the EOS dump uses the
  same per-layer T/S/p as the model's `calculate_density` call at line 1437
  plus explicit `calculate_density_derivs_scalar` calls); (2)
  `mom_ocean_model_nuopc.F90` after `set_net_mass_forcing` (forcing dump).
  `dump_init`/`dump_finalize` wired into `ocean_model_init`/`ocean_model_end`
  (via `mom_cap_methods.F90` or `mom_ocean_model_nuopc.F90`).
- The three source files + submodule pointers are committed to the fork
  (local, unpushed): `instrument/mom6-dumps`.

## Validation

1. **Solution neutrality:** instrumented run's `mom6.r`, `mom6.r_stoch`,
   `cpl.r`, `datm.r`, `dice.r`, `drof.r` restarts + `ocean.stats.nc`
   bit-identical to the reference `C_JRA_t62` run's, modulo timestamps/
   titles.
2. **CSV hygiene:** no blank lines, arity = `# columns` count, row counts =
   cells × steps as selected.
3. **Law checks (hand, vs Fortran):**
   - EOS: ρ(T,S,p) and drho/dT, drho/dS reproduced by the Wright_FULL
     formula on all rows (Python reimplementation of
     `density_elem_Wright_full`), and cross-checked against
     `MOM_EOS_Wright_full.F90` constants.
   - KPP: kOBL consistency (OBL depth within the kOBL cell),
     Kt/Ks/Kv update rule (KPP_IS_ADDITIVE=False ⇒ K replaces interior
     values where non-zero), nonLocalTrans zero where
     surfBuoyFlux≥0, NLT shape profile (1-σ)²·(1+2σ) per MATCH_TECHNIQUE/
     NLT_SHAPE on the traced columns.
   - Forcing: fluxes%sw = sum of the four sw components; latent =
     q_flux·Lv − fprec·Lf − frunoff·Lf (sign conventions); salt_flux =
     kg_m2_s_conversion·IOB%salt_flux; u/v_flux rotation
     (cos_rot/sin_rot) reproduced from the mediator values.

## Proposed stage-2 target

**EOS (WRIGHT_FULL) first** — it is a closed-form pointwise law with 15
constants, no state, and no discretization: the natural first MOM6 stage-2
component (`mom6_eos.esm`). The dump **suffices to pin it to 1e-12
relative**: CSV values are written with `es26.17` (17 significant digits,
full double precision), T/S/p inputs and ρ/dρ outputs are exact doubles of
the Fortran computation, and the Wright_FULL formula is a fixed sequence of
~20 arithmetic operations (round-off ~1e-15 relative for a naive
transcription). The derivative columns (drho/dT, drho/dS) let the tests pin
both the law and its derivatives. KPP is the second target (per-column
CVMix law with OBL depth + profile functions); the surface forcing import
is the third (pure field mapping, no physics).
