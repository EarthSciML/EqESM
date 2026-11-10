# MOM6 lane — stage-1 instrumentation (W1)

Status: **run complete, dumps validated, all checks green** (2026-11-02).
Fork `stage1/EqESM-cesm-mom6` (branch `instrument/mom6-dumps`), case
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
     q_flux·Lv − fprec·Lf − frunoff·Lf − frunoff_glc·Lf (sign
     conventions); salt_flux_in = kg_m2_s_conversion·IOB%salt_flux and
     salt_flux = salt_flux_in + salinity-restoring flux
     (RESTORE_SALINITY=True; the mediator's import is 0 in this run);
     AGRID tau_mag/ustar/omega built from the imported u_flux/v_flux.

## Run results (2026-11-02)

- **Run 1** (`case.run` 11980255, COMPLETED 0:0): dumps were written but
  with two bugs found during inspection - (a) `dump_open` wrote to
  `trim(dump_dir)//trim(fname)` and DUMP_DIR has no trailing slash, so the
  files landed one directory up as `mom6mom6_*.csv`; (b) the coupling-step
  counter was incremented once per traced *cell* instead of once per
  coupling interval, so the trace covered 1 full + 1 partial coupling step
  and let the KPP dumps run 2 intervals past the intended budget. Fixed in
  MOM6 84f9c9a66 (dump_open path separator; counter advanced by the first
  traced cell of each interval, re-gated after increment; init
  confirmation messages via `MOM_error(NOTE)` so they land in the
  component log). Solution neutrality of run 1 was verified anyway
  (`compare_mom6_restarts.py`: mom6.r/mom6.r_stoch/cpl.r/datm.r/dice.r/
  drof.r restarts + ocean.stats.nc + Vertical_coordinate.nc all
  bit-identical to the reference, 0 differences).
- **Run 2** (11986853): rerun with the fixed module - correct filenames
  and step semantics, but the validation pass exposed an **uninitialized
  `LangEnhK`** in `MOM_CVMix_KPP.F90`: with LT_K_ENHANCEMENT=False the
  local remained a denormal (`1.636e-311`) that leaked into the kpp CSV.
  Fixed with an explicit `else; LangEnhK = 1.0` branch (fork commit; the
  original trailing `endif` after the Langmuir do-loop was removed so the
  loop stays inside the IF). Run-2 restarts also bit-identical to the
  reference.
- **Run 3** (11992759, COMPLETED 0:0, 7:43 wall; st_archive 11992761):
  the final, valid run. `validate_mom6_dumps.py` passes all checks
  (results below) and the run-3 restarts are bit-identical to the
  reference (all 6 restart files + ocean.stats.nc + Vertical_coordinate.nc,
  0 differences).
- **Traced domain:** PET 0 is the SW corner block of the 6x6
  decomposition; only 12 of its cells are ocean (the 77 southernmost grid
  rows are Antarctic land mask), so the traced set is all 12 ocean cells
  (dump_first_cells=32 requested, 12 found).
- **Trace pattern (CESM startup ocean-lag, mom_cap.F90:1835-1856):** the
  first coupling timestep is skipped, and the second call steps back one
  interval and advances by two - so for the 12 traced coupling steps the
  forcing dump fired 12 times (OS%Time tods 0, 7200..43200; the first is
  the lag call whose OS%Time is still 0) while the KPP dumps fired on all
  13 thermo steps covered (Time_end tods 3600..46800; cpl_step 1 appears
  twice - the first advance covers two thermo steps with the same
  forcing). Row counts: eos 12x75x13 = 11700; kpp_obl/kpp 12x13 = 156
  each; sfc_forcing 12x12 = 144.

## Validation results (run 3, `validate_mom6_dumps.py`, all green)

- **EOS law:** Wright_FULL ρ(T,S,p) and drho/dT, drho/dS bit-exact (max
  rel err 0.0) on all 11700 rows (Python mirror of the Fortran operation
  order).
- **Hydrostatic pressure:** p_{k+1} − p_k = g·H_to_RZ·h_k with a constant
  ratio 1.0143e4 = RHO_0·g = 1035·9.8 (stdev 1.2e-7 relative, round-off
  only) over all 11544 layer interfaces; p(k=1) == p_surf exactly on all
  156 columns. (The pressure column is the Boussinesq pRef at the top
  interface of each layer, accumulated in KPP_compute_BLD:
  pRef += g_Earth·H_to_RZ·h.)
- **KPP update:** with KPP_IS_ADDITIVE=False, Kt_out = Kd_heat where
  Kd_heat≠0 else Kt_in (same for Ks/Kv) on all 11856 interfaces.
- **NLT:** NLT_heat == NLT_salt everywhere; NLT = 0 at the surface
  interface (k=1) and everywhere surfBuoyFlux ≥ 0. All 156 traced columns
  are stably stratified (surfBuoyFlux ≥ 0), so the unstable-column
  NLT_SHAPE branch is not exercised by this trace.
- **OBL:** kbl = int(kOBL) in [1,ke], fractional part in [0,1), and
  −OBLdepth bracketed by the interface pair around kOBL, on all 156
  columns; z_inter identical between the OBL and KPP dumps for the same
  (cell, step).
- **Forcing:** latent (q_flux·Lv − fprec·Lf − frunoff·Lf − frunoff_glc·Lf)
  and the three latent diags, sw = sum of the 4 components, evap =
  q_flux, lprec/fprec/lrunoff/frunoff and seaice_melt(+heat) mapping,
  p_surf_full = p (both fluxes and forces), ice_fraction,
  tau_mag_o = GUST_CONST + |u_flux,v_flux|, ustar_o = sqrt(tau_mag_o/1035),
  omega_w2x = atan2(v_flux,u_flux): all 144/144. **salt_flux:**
  RESTORE_SALINITY=True ⇒ fluxes%salt_flux = salinity-restoring flux +
  coupler import, with the import alone recorded in fluxes%salt_flux_in
  (convert_IOB_to_fluxes, mom_surface_forcing_nuopc.F90:686-689); the
  mediator sends salt_flux = 0 in this run, so salt_flux_in == salt_flux
  == 0 on all 144 rows and salt_flux_o (4.90e-7..1.02e-6) is purely the
  restoring flux (SALT_RESTORE_FILE=state_restore_tx2_3_20230416.nc,
  ADJUST_NET_SRESTORE_TO_ZERO=True).

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

## Wave-2 results (2026-11-05)

Wave-2 boundary dumps (five new boundaries) were added to the fork
(`MOM6` submodule HEAD 5df49d288 on top of the wave-1 219c8ccee; fork HEAD
234c605) and validated on the final rerun of `C_JRA_t62_mom6_inst`
(case.run 12139764, st_archive 12139765, COMPLETED 2026-11-05 18:00:59):

- **mom6_sw_abs.csv** (156 rows = 12 cells × 13 steps; nsw=2): per-band
  surface split `sw_pen_band_n`, exponential band decay
  `pen_sw_rem(n,k) == pen_sw_rem(n,k-1)·exp(-opacity(n,k)·h(k))` above the
  min_SW_heat cutoff (3431 testable layers, all exact), net penetrating SW
  `netPen_k`, layer heating `dT_k`; total absorption identity vs
  `penSW_bnd` (surface in − bottom leftover) holds on all 156 rows.
- **mom6_frazil.csv** (312 rows = 2 calls × 156): heat conservation
  `frazil_out − frazil_in == −RHO_0·C_P·Σ h·(T_after − T_before)` on all
  312 rows; freezing-point law (LINEAR: 0.0 − 0.054·S − 7.75e-8·p) on all
  312; `T_after ≥ T_freeze` on all 23400 layer values.  No frazil heat
  forms in this window (early-summer Antarctic cells; sum of changes is 0)
  — the formation/reclaim branches need a winter window.
- **mom6_bflux.csv** (156 rows): mass/salt/heat update laws hold on all
  rows (heat: `Σ h·T_after − h·T_before == netHeat + Σ penSW_bnd −
  leftover`, max rel 1.6e-11; `nonpenSW` is inside `netHeat`; the earlier
  draft that added `nonpenSW` separately was wrong).  netMassIn ≡ 0
  (AGGREGATE_FW_FORCING); netMassOut −1.6e-5..2.0e-3 m.
- **mom6_diabatic_tri.csv** (156 rows): conservation and the tridiagonal
  law hold (max residual 1.2e-12 relative to an operand-magnitude scale —
  convective interfaces carry ent ~ 1e6 m, which makes the plain h·T
  relative scale meaningless); cross-file equality with bflux
  (h_after, T_after, S_after == tri h, T_before, S_before) on all 156 rows.
- **mom6_kpp_wide.csv** (8164 rows over 33 PEs): unstable (surfBuoyFlux<0)
  or warm (T_sfc > 20 °C) columns of the first 32 ocean cells of every PE,
  per-PE per-cell per-step rows merged at finalize; selection bits,
  Kt/Ks/Kv update rule, NLT structure, and z_inter/kOBL/OBLdepth geometry
  verified on all 8164 rows.  `step` is a per-PE ordinal over the PE's
  *qualifying* steps (the counter increments only when a traced column
  qualifies), so PET 0's Antarctic shelf cells — stable until late in the
  day — contribute rows at tods 64800..75600 with ordinals 1..4; `tod` is
  the model time.

Bugs found and fixed during wave-2 validation (all in the dump code, not
the model): wide-file merge line buffer truncation (rows are ~24 KB, not
~5 KB — 20000 → 40000 chars), bflux header band-10 name (`i1` overflow
wrote `*`; `i2.2` now zero-pads all band names), and the wave-1-style
`penSW_bnd_ 1` leading space (`i2` → `i2.2`).  Neutrality of the final run
re-verified: mom6.r/mom6.r_stoch/cpl.r/datm.r/dice.r/drof.r restarts +
ocean.stats.nc + Vertical_coordinate.nc all bit-identical to the reference
(8 files, 0 differences).  Wave-1 validator still green on the regenerated
CSVs (eos 11700, kpp 156, kpp_obl 156, sfc_forcing 144 rows).

Wave-2 stage-2 targets (see `stage2/README-mom6-w2.md`): shortwave
absorption, frazil formation, boundary flux application, tridiagonal
mixing, and widened KPP tests (unstable/warm columns).

## Wave-3 results (2026-11-09)

Wave-3 dynamics boundary dumps (continuity, Coriolis/momentum advection,
pressure gradient force) were added to the fork (`MOM6` submodule HEAD
4ba241193 = a1360523b dynamics dumps + c3519857f extended stencils +
4ba241193 layer-k fix; fork HEAD 8dd0e0e, branch `instrument/mom6-dumps`)
in case `C_JRA_t62_mom6_dyn_inst` (clone of the wave-2 case; same 12
traced cells of PET 0, 12 dynamics steps × 75 layers):

- **mom6_continuity.csv** (21600 rows = 12 cells × 12 steps × 75 layers ×
  2 calls; calls 1 and 3 — call 2 is a skipped duplicate in the
  STORE_CORIOLIS_ACCEL=True configuration): directionally-split PPM
  thickness advection (MOM_continuity_PPM.F90, FIRST_DIRECTION=x_first,
  MONOTONIC_CONTINUITY=False, SIMPLE_2ND_PPM=False, UPWIND_1ST=False,
  VOL_CFL=False).  Call 1 is the pure BT_cont layer-flux call (uh/vh obey
  the PPM flux law); call 3 is the BT-adjusted final call.  Per-cell
  columns: hin/h_out, h_inter (zonal-pass intermediate thickness),
  h_inter_jm4..jp3 (7-cell meridional stencil), uh_W/E, vh_S/N, u/v on
  the 4 faces, edge_ok, dt, x_first, h_min/h_min_recon
  (=Angstrom_H/2·Angstrom_H).
- **mom6_coradv.csv** (10800 rows = 12 × 12 × 75, call 2 only): h/u/v/uh/
  vh/Area/KE/KEx-KEy/CAu-CAv on the cell stencil (CORIOLIS_SCHEME=
  SADOURNY75_ENERGY, KE_SCHEME=KE_ARAKAWA, BOUND_CORIOLIS=True,
  no_slip=False), dvdx/dudy/rv/av/q at the 4 corners, Ih_q, grid metrics.
- **mom6_pgf.csv** (10800 rows): e_bot from the layer-k rho profile,
  e_sfc/e cross-layer law, pa_top at k=1 (constants for the reference run:
  GxRho_ref = 9.8·1035, I_Rho0 = 1/1035, Z_ref = 0.0, h_neglect = 1e-30,
  H_to_Z = 1).

**Dump bug found and fixed (the only model-side change of the wave):** the
coradv dump passed the full 3-D arrays `u, v, h, uh, vh, CAu, CAv` to 2-D
dummy arguments, so every written k-row carried the layer-1 values of
those arrays while the per-layer private locals (rv, KE, hArea, q) were
correct — the cross-layer laws for rv/hArea/KE/CAu-CAv failed on exactly
the aliased columns.  Fixed at the call site (MOM_CoriolisAdv.F90:1210) by
passing layer-k slices; committed as 4ba241193 (MOM6) → b2e4a9e (mom) →
8dd0e0e (fork root).  Rebuilt and reran: case.run 12282697, st_archive
12282698, COMPLETED 2026-11-09 10:30 (fresh start from 0001-01-01, restarts
removed).

**Not fixed (documented limitation):** the continuity dump's call-3
`hin`/`h_inter` columns alias the same array as `h` (the call passes
`continuity(u_inst, v_inst, h, h, ...)`), so call-3 hin is the final h.
h is bit-identical between calls 1 and 3 (verified), so the call-3
convergence identities are tested against the call-1 hin as base (bit
exact, 0 bad); the aliasing is a dump artifact, not a model issue.

**Neutrality re-confirmed after the fix:** the fixed run's 9 original
CSVs are bit/value-identical to the wave-2 reference run — mom6_kpp_wide
(203 MB, 33 PEs, 8164 rows) byte-identical to the Nov-5 final reference,
and bflux/diabatic_tri/frazil/sw_abs identical on all overlapping rows
(the `w2_test/` Nov-5 copies; the reference-run backups in
/scratch.local/ctessum/mom6_dumps_refrun/ + mom6_dumps_dyninst/ were lost
when /scratch.local was cleared on a job requeue — the surviving w2_test/
copies served as the diff baseline).  eos/kpp/kpp_obl/sfc_forcing are
written by unchanged dump code in the same bit-identical run.

**Validator fixes (validator-side transcription bugs, found while
revalidating the fixed data):** the coradv `dvdx` reconstruction used the
wrong v-stencil at the 'se'/'sw' corners (v_ne·dyCv_ne − v_s·dyCv_s and
v_se·dyCv_se − v_sw·dyCv_sw; correct is v_se·dyCv_se − v_s·dyCv_s and
v_s·dyCv_s − v_sw·dyCv_sw — the model's dvdx(I,J) = v(i+1,J)·dyCv(i+1,J) −
v(i,J)·dyCv(i,J) mapped to the dump's per-cell corner columns), and the
CAv_s raw used the CAv_n uh-stencil (uh_w+uh_wn / uh_e+uh_en; correct is
uh_ws+uh_w / uh_es+uh_e for the v-point (i, J−1)).  Both were wrong on
every row until the dump fix made the per-layer data correct; the fixed
validator also pins the dump's dvdx/dudy columns directly against the u/v
stencil (43200/43200).

**Final validation (`stage1/analysis/validate_mom6_dyn.py`, ALL CHECKS
PASSED on the fixed run's CSVs):** continuity — convergence identities
0/21600 bad, reconstruction pins 21600/21600, flux-law pins 7200/7200
(edge_ok cells (3,78),(4,78) only; PET 0's compute domain ends at local
j=84); coradv — rv 43200/43200, av 43200/43200, q 43200/43200,
dvdx/dudy 43200/43200, Ihq_ne 10800/10800, hArea 86400/86400, KE_c
10800/10800, CAu/CAv 43200/43200, grid 54000/54000; pgf — e_bot
32400/32400, e_sfc 144/144, cross-layer 21312/21312, pa_top(k=1)
144/144, constants 4500/4500.  Land masks: Area_h == 0 (coradv stencil)
extended with the h ≡ 0.001 land-thickness fill (40 land cells).  No
bugs found in the reference C++/Fortran model code during this wave
(bugs.md unchanged); all fixes were in the dump instrumentation and the
validator.

Wave-3 stage-2 target (see `stage2/README-mom6-dyn.md`): the continuity
PDE ∂h/∂t = −∇·(u h) with the directionally-split PPM reconstruction +
positive limiter as expression templates (PPM templates authored in
`stage2/mom6_templates.esm`), the C-grid flux-divergence rule (candidate
EarthSciDiscretizations PR), and the coradv/pgf laws pinned by informative
points from these dumps.
