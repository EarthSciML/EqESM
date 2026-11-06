# MOM6 stage 2 — wave 2: shortwave, frazil, boundary fluxes, tridiagonal mixing, widened KPP

Stage-2 implementation proposal for the five wave-2 MOM6 boundaries of the
MOM6 lane (`stage1/runs/MOM6_LANE.md`), as esm 1.3.0 files with `tests`
blocks pinned to the stage-1 wave-2 instrumentation dumps
(`stage1/dumps/mom6/mom6_sw_abs.csv`, `mom6_frazil.csv`, `mom6_bflux.csv`,
`mom6_diabatic_tri.csv`, `mom6_kpp_wide.csv`).

The dumps were produced by the instrumented fork
`stage1/EqESM-cesm-mom6` (branch `instrument/mom6-dumps`, MOM6
`dev/ncar_260909 @301879336` + the wave-2 dump commits), run
`C_JRA_t62_mom6_inst` (1 day, 0001-01-01; T62_t232, NK=75 z\* ALE, DT=1800,
DT_THERM=3600, 72 PEs; solution-neutral — restarts bit-identical to the
reference `C_JRA_t62` run).  Traced domain: the 12 ocean cells of PET 0
(the SW-corner block) × 13 thermo steps (tods 3600..46800; cpl_step 1
covers two thermo steps, the CESM startup ocean lag), plus — for the KPP
widening — the first `dump_first_cells=32` ocean cells of **every** PE,
restricted to columns that are unstable (surfBuoyFlux < 0) or warm
(T_sfc > DUMP_WIDE_WARM_T = 20 °C) at that step (8164 rows over 33 PEs:
3490 unstable-only, 1401 warm-only, 3273 both; T_sfc −1.86..29.96 °C,
surfBuoyFlux −4.5e-7..5.9e-7, OBLdepth 5..397 m).

Sources: `components/mom/MOM6/src/parameterizations/vertical/MOM_opacity.F90`
(`sumSWoverBands` :875-1016, `absorbRemainingSW` :602-869),
`MOM_diabatic_aux.F90` (`make_frazil` :120-240, `applyBoundaryFluxesInOut`
:464-961), `MOM_diabatic_driver.F90` (`diabatic_ALE` :1692-1735, the
`tracer_vertdiff_Eulerian` calls), `MOM_tracer_diabatic.F90`
(`tracer_vertdiff_Eulerian` :226-417), `MOM_CVMix_KPP.F90` (`KPP_calculate`
:1056), `MOM_TFreeze.F90` (linear freezing point), and
`src/framework/mom6_dump_mod.F90` (the dump writer).

## File map (proposed)

| File | Component | Boundary | Fortran source |
|---|---|---|---|
| `mom6_shortwave.esm` | `ShortwaveAbsorption` — per-band penetrating SW: surface split `sw_pen_band(n)`, exponential band decay through `opacity(n,k)`, layer heating `dT_k`, net penetrating SW on interfaces `netPen_k` | `mom6_sw_abs.csv` | `MOM_opacity.F90`: sumSWoverBands :875-1016, absorbRemainingSW :602-869 |
| `mom6_frazil.esm` | `FrazilFormation` — column frazil heat accounting: supercooling → T raised to T_freeze, heat accumulated in `frazil`; reclaim into warm layers; thin-layer exception | `mom6_frazil.csv` | `MOM_diabatic_aux.F90`: make_frazil :120-240 |
| `mom6_bflux.esm` | `BoundaryFluxes` — surface flux application: netMassIn/Out thickness change, netHeat/netSalt/nonpenSW/penSW per-layer updates, h/T/S before/after | `mom6_bflux.csv` | `MOM_diabatic_aux.F90`: applyBoundaryFluxesInOut :464-961 |
| `mom6_diabatic_tri.esm` | `TridiagonalMixing` — the implicit vertical diffusion solve `(h_k+ent_k+ent_{k+1})T'_k = h_k T_k + ent_k T'_{k-1} + ent_{k+1} T'_{k+1}` with `ent = dt·Kd/dz_int` | `mom6_diabatic_tri.csv` | `MOM_tracer_diabatic.F90`: tracer_vertdiff_Eulerian :226-417 (no-sink branch :384-413); driver :1692-1735 |
| `mom6_kpp_wide.esm` | extension of the wave-1 `mom6_kpp.esm` tests to unstable/warm columns (the wave-1 traced PET-0 cells are all stably stratified, so the NLT and OBL branches for unstable columns were untested) | `mom6_kpp_wide.csv` | `MOM_CVMix_KPP.F90`: KPP_calculate :1056 |

The wave-1 `mom6_templates.esm` gains a `linear_tfreeze` expression
template (`T_Fr = (TFREEZE_S0_P0 + DTFREEZE_DS·S) + DTFREEZE_DP·p`, with the
run's coefficients 0.0 / −0.054 / −7.75e-8) shared by the frazil component.

## Dump formats and identities (what the tests pin)

All values in physical units (all internal rescale factors are 1.0 in the
reference run: internal == MKS).  Common context columns on every row:
`cpl_step,ymd,tod,gi,gj,i,j,lon,lat`.  Rows are identified by
(cpl_step, tod, gi, gj) — cpl_step 1 appears twice (two thermo steps under
the first coupling interval) so `tod` is the unique key.

### mom6_sw_abs.csv — 156 rows (12 cells × 13 steps), nsw = 2 bands

Columns: `nsw`, `sw_pen_band_1..nsw` [W m-2] (the surface band split from
set_pen_shortwave), `opacity_n_k` [m-1], `pen_sw_rem_n_k` [degC m]
(integrated penetrating SW remaining below layer k, per band),
`dT_k` [degC] (T change caused by the SW in this call),
`netPen_01..76` [degC m s-1] (net penetrating SW at the interfaces from
sumSWoverBands; 0 when the buffered profile is not fresh — it is fresh on
every traced row).

Identities (all hold to round-off in the dump; verified by
`stage1/analysis/validate_mom6_w2.py`):

- surface three-way consistency (sw_abs × bflux × sw_abs):
  `sw_pen_band(n) == RHO_0·C_P·penSW_bnd(n)/dt` and
  `netPen_1 == Σ_n penSW_bnd(n)/dt` (RHO_0 = 1035 kg m-3, C_P = 3992 J kg-1
  K-1 from the run's MOM_input);
- exponential band decay for every layer above the min_SW_heat cutoff
  (pen_sw_rem(n,k-1) > 1e-6 degC m — below that the model clamps or zeroes
  SW_trans):
  `pen_sw_rem(n,k) == pen_sw_rem(n,k-1)·exp(−opacity(n,k)·h(k))` with h
  from bflux `h_after` (the h passed to absorbRemainingSW);
- total absorption: `Σ_k dT_k·h_k == Σ_n (penSW_bnd(n) − pen_sw_rem(n,75))`
  (absorbed = surface in minus bottom leftover; exact model arithmetic,
  including the thin-column `T_chg` redistribution path).

### mom6_frazil.csv — 312 rows (12 cells × 13 steps × 2 calls)

Columns: `call` (1 = pre-diabatic, 2 = post-diabatic), `p_surf` [Pa],
`frazil_in`/`frazil_out` [J m-2] (tv%frazil before/after the call),
`fraz_col` [J m-2] (column heat requirement accumulated by this call),
per layer `T_before_k`/`T_after_k` [degC], `S_k` [ppt], `p_k` [Pa]
(mid-layer pressure), `T_freeze_k` [degC] (calculate_TFreeze at (S_k,p_k),
re-executed by the dump with the same arguments), `h_k` [m].

Identities (verified in the dump):

- heat conservation per call:
  `frazil_out − frazil_in == −RHO_0·C_P·Σ_k h_k·(T_after_k − T_before_k)`;
- freezing-point law: `T_freeze_k == (0.0 + (−0.054)·S_k) + (−7.75e-8)·p_k`
  (TFREEZE_FORM=LINEAR; TFREEZE_S0_P0/DTFREEZE_DS at their defaults,
  DTFREEZE_DP from MOM_input);
- `T_after_k ≥ T_freeze_k − 1e-9` on every layer (make_frazil never leaves
  a layer below freezing; the thin-layer branch and the reclaim branch both
  end at or above T_freeze).

**Coverage note**: no frazil heat forms in this window (sum of
frazil_out − frazil_in over all 312 rows is 0 — the traced Antarctic cells
are in early summer).  The stage-2 tests pin the no-op branches (T_after ==
T_before, frazil unchanged) plus the freezing-point law on all layers; the
formation/reclaim branches need a winter window (or the widened trace).

### mom6_bflux.csv — 156 rows

Columns: `dt` [s], `netMassInOut`/`netMassIn`/`netMassOut` [m] (net mass
entering/leaving over the step), `netHeat` [degC m] (heat flux over the
step, excluding penetrative SW and the heat content of outgoing mass),
`netSalt` [ppt m], `nonpenSW` [degC m], `penSW_bnd_1..10` [degC m] (per-band
penetrating SW entering the water column; bands 3..10 are the
MAX_SW_BANDS=10 padding, zero when nsw=2), per layer `h_before_k`/`h_after_k`
[m], `T_before_k`/`T_after_k` [degC], `S_before_k`/`S_after_k` [ppt].

Identities (all hold to round-off in the dump):

- mass: `Σ_k (h_after − h_before) == netMassInOut`;
- salt: `Σ_k (h_after·S_after − h_before·S_before) == netSalt`;
- heat: `Σ_k (h_after·T_after − h_before·T_before) ==
  netHeat + Σ_n penSW_bnd(n) − leftover`, where `netHeat` already includes
  the non-penetrating SW (`nonpenSW` is the non-penetrating part of the SW,
  inside netHeat), and leftover is the bottom `pen_sw_rem` from
  mom6_sw_abs.csv (zero for the deep traced columns — the reference run's
  traced cells absorb all SW).  Max relative error in the dump: 1.6e-11.

The stage-2 component consumes the ShortwaveAbsorption component's outputs
(dT, leftover) — the two components share the same h/T/S state per step.

### mom6_diabatic_tri.csv — 156 rows

Columns: `dt` [s], per interface (76) `ent_t_k`/`ent_s_k` [m] (tridiagonal
entrainment coefficients dt·Kd/dz_int for T and S; zero at the surface and
bottom interfaces), `Kd_heat_k`/`Kd_salt_k` [m2 s-1] (post-KPP interface
diffusivities), per layer (75) `h_k` [m] (unchanged by the solve),
`T_before_k`/`T_after_k` [degC], `S_before_k`/`S_after_k` [ppt].

Identities (verified in the dump):

- conservation (no surface/bottom fluxes in this call):
  `Σ_k h_k·(T_after_k − T_before_k) == 0` and the same for S;
- tridiagonal residual per layer:
  `h_k·(T'_k − T_k) == ent_t_k·(T'_{k−1} − T'_k) + ent_t_{k+1}·(T'_{k+1} − T'_k)`
  with ent_1 = ent_76 = 0 (the defining equation of the solve).  Convective
  interfaces carry ent ~ 1e6 m, so the residual is scaled by the operand
  magnitudes (|lhs| + |e_in| + |e_out| + |ent·T'| terms): max residual
  1.2e-12 relative — round-off (a Python Thomas solve from the dumped
  h/ent/T_before reproduces T_after to 2e-9 absolute);
- cross-file: `bflux.h_after == tri.h`, `bflux.T_after == tri.T_before`,
  `bflux.S_after == tri.S_before` exactly (only diagnostics run between
  applyBoundaryFluxesInOut and the solve in diabatic_ALE).

### mom6_kpp_wide.csv — 8164 rows (33 PEs; per-PE per-cell per-step)

Columns: `pet` (ocean PE index), `step` (per-PE ordinal over the traced
thermo steps in which the PE had at least one qualifying traced column —
**not** a global step number; `tod` is the model time), then the
mom6_kpp.csv context (`ymd,tod,gi,gj,i,j,lon,lat`), `sel` (selection bits:
1 = unstable, surfBuoyFlux < 0; 2 = warm, T_sfc > 20 °C), `T_sfc`/`S_sfc`
[degC/ppt], `OBLdepth` [m], `kOBL`, `uStar` [m s-1], `surfBuoyFlux`,
`La_SL`, `LangEnhK`, then per interface (76): `z_inter_k` [m],
`Kd_heat_k`/`Kd_salt_k`/`Kv_k` [m2 s-1], `NLT_heat_k`/`NLT_salt_k`,
`Kt_in_k`/`Ks_in_k`/`Kv_in_k` (pre-KPP interior values), `Kt_out_k`/
`Ks_out_k`/`Kv_out_k` (final).

The per-PE files `mom6_kpp_wide.petNNNN.csv` are merged into
`mom6_kpp_wide.csv` by the root PE at finalize (`sync_PEs` barrier, PE
order, headers skipped, per-PE files deleted).

Identities (verified in the dump):

- selection bits: `(sel & 1) == (surfBuoyFlux < 0)`,
  `(sel & 2) == (T_sfc > 20.0)`, sel != 0;
- the wave-1 KPP update rule on every interface of every row:
  `Kt_out == Kd_heat` where Kd_heat ≠ 0, else `Kt_in` (and same for
  salt/vis) — this run's PET-0 trace is all stable columns, so the wide
  file is where the non-trivial branches (NLT shape, OBL depth, Kt/Ks/Kv
  updates on unstable columns) get exercised;
- NLT structure: NLT_heat == NLT_salt; NLT(1) == 0 (surface interface);
  NLT == 0 everywhere on stable rows;
- geometry: z_inter(1) == 0, interfaces below the surface negative and
  increasing; kOBL in [1,75]; z_inter(kOBL) ≥ −OBLdepth > z_inter(kOBL+1);
  uStar ≥ 0; OBLdepth > 0.

## Units

| Quantity | Symbol | Units |
|---|---|---|
| Potential temperature | `T`, `T_sfc`, `T_freeze` | degC |
| Salinity | `S`, `S_sfc` | ppt |
| Pressure | `p_k`, `p_surf` | Pa |
| Layer thickness | `h_k` | m |
| Interface height | `z_inter` | m |
| Diffusivities | `Kd_heat`, `Kd_salt`, `Kv`, `Kt_*`, `Ks_*`, `Kv_*` | m²/s |
| Penetrating SW (integrated) | `pen_sw_rem`, `penSW_bnd`, `netHeat`, `nonpenSW` | degC·m |
| Penetrating SW (surface split) | `sw_pen_band` | W/m² |
| Penetrating SW (rate) | `netPen` | degC·m/s |
| Opacity | `opacity` | 1/m |
| Frazil heat | `frazil_in`, `frazil_out`, `fraz_col` | J/m² |
| Mass fluxes | `netMassInOut`, `netMassIn`, `netMassOut` | m |
| Salt flux | `netSalt` | ppt·m |
| Entrainment coefficients | `ent_t`, `ent_s` | m |
| Boundary-layer depth | `OBLdepth` | m |
| Surface buoyancy flux | `surfBuoyFlux` | m²/s³ |
| Friction velocity | `uStar` | m/s |
| Langmuir numbers | `La_SL`, `LangEnhK` | 1 |

## Testing

Run with

    esm validate stage2/mom6_templates.esm stage2/mom6_shortwave.esm \
        stage2/mom6_frazil.esm stage2/mom6_bflux.esm \
        stage2/mom6_diabatic_tri.esm
    esm test stage2/mom6_shortwave.esm stage2/mom6_frazil.esm \
        stage2/mom6_bflux.esm stage2/mom6_diabatic_tri.esm

All `tests` blocks are `time_span {0,0}` instantaneous pins; every expected
value is a real dump row, identified by (tod, gi, gj) — for the wide file,
by (pet, tod, gi, gj).  The dump writes `es26.17` and the reference run is
r8, so the transcriptions are pinned to round-off (rel 1e-12).  The run
configuration is documented in the file metadata (nsw=2, OHLMANN_03
opacity bands; RHO_0 = 1035, C_P = 3992; TFREEZE_FORM=LINEAR with
TFREEZE_S0_P0 = 0, DTFREEZE_DS = −0.054, DTFREEZE_DP = −7.75e-8;
AGGREGATE_FW_FORCING=True, so netMassIn ≡ 0; all unit rescale factors 1.0).

### Coverage limitations (this data window)

- **sw_abs**: the deep layers of every traced profile fall below the
  min_SW_heat cutoff (pen_sw_rem < 1e-6 degC m), where the model clamps or
  zeroes SW_trans — the decay law is pinnable only above the cutoff (3431
  of 23400 layer values; the remaining layers pin at exactly 0 or at the
  clamped value); the thin-column redistribution path (h_heat <
  H_limit_fluxes) is not exercised (all traced columns are deep); only
  nsw=2 bands (OHLMANN_03) are covered.
- **frazil**: no frazil heat forms (early-summer Antarctic cells) — the
  formation/reclaim/thin-layer branches are unpinned; the freezing-point
  law and the T_after ≥ T_freeze invariant are pinned on all 312 rows.
- **bflux**: netMassIn ≡ 0 (AGGREGATE_FW_FORCING) so the step-A
  netMassIn branch is pinned at zero only; netMassOut is nonzero (up to
  2e-3 m) so the mass-removal heat accounting is exercised; all traced
  columns absorb all SW (leftover = 0), so the shallow-column leftover
  branch is unpinned.
- **tri**: no surface/bottom flux terms in this call (ent_1 = ent_76 = 0
  pinned); Kd_heat == Kd_salt on every row (KPP sets both), so the
  T/S-differing-diffusivity branch is unpinned; the `diabatic_diff_tendency_diag`
  and salt-adjust branches are not in the call path.
- **kpp_wide**: the traced window is one January day — the warm rows come
  from the tropics and the unstable rows from mid/high latitudes (33 of 36
  PEs have qualifying columns; PET 0's Antarctic shelf cells qualify only
  late in the day); the NLT shape law on unstable columns is exercised but
  the CVMix surface-layer matching (sigma_0) details are only pinned
  through the raw NLT values, not through a closed-form shape fit; the
  NLT_SHAPE=CVMix profile is config-dependent (MATCH_TECHNIQUE), so the
  shape check in `validate_mom6_w2.py` is structural (NLT_heat == NLT_salt,
  surface zero, stable zero) rather than a closed-form fit.
- The B-grid/C-grid wind_stagger alternatives, other opacity parameterizations
  (2-band vs 3-band vs 4-band), and other TFREEZE_FORM choices are not
  transcribed.

## Next steps

1. `stage2/mom6_shortwave.esm`, `mom6_frazil.esm`, `mom6_bflux.esm`,
   `mom6_diabatic_tri.esm` — the four new components above; the bflux
   component consumes the shortwave component (dT/leftover), and the tri
   component consumes the KPP component of wave 1 (`mom6_kpp.esm`) for its
   Kd/ent inputs.
2. Extend `mom6_kpp.esm` with `mom6_kpp_wide.esm` tests (unstable/warm
   columns) — the wide file is the only place the unstable-column branches
   are exercised.
3. A wider stage-1 window (a winter month for frazil formation, a shallow
   column for the SW leftover branch, a non-aggregate-FW configuration).
4. Mount into `eqesm.esm` with `earthsciintegrator` (stage 3) and compare
   against full-CESM runs.
