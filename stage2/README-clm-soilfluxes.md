# CLM soil fluxes — SoilFluxes

Stage-2 implementation of `SoilFluxesMod.F90:SoilFluxes` (ctsm5.4.063,
cesm3.0-alphabranch; checkout `stage1/cesm/components/clm`): the update of the
surface energy and moisture fluxes to the new soil temperature — the
previous-time-step ground temperature and T⁴ blends, the sensible/latent flux
corrections, the liquid/solid partition of ground evaporation, the
snow-evaporation moisture limit, the ground heat flux, the totals, and the
outgoing/net longwave radiation with the bare-ground skin temperature.
`clm_soilfluxes.esm` is the component file (model `SoilFluxes`, 32 equations),
with the shared laws imported by reference from the `clm_templates.esm`
library; the tests pin the stage-1 instrumentation dumps (`stage1/dumps/clm/`:
`clm_soilfluxes.csv`, `clm_bareground.csv`, `clm_soiltemp.csv`,
`clm_surfrad.csv`, `clm_run_imports.csv`) and, for the downward longwave
forcing of columns 1–5, `stage1/dumps/cdeps/datm_gswp3_exports.csv`.

The traced run is the D-class single-point column experiment of README-clm.md:
16 patches on 16 columns (bare, deep snowpacks), 12 steps of dtime = 1800 s.
The SoilFluxes dump traces 192 rows (16 patches × 12 steps) with p == c on
every row, so the model is written all-`patch`-shaped with the column-level
quantities of the Fortran replicated patch-wise (the stage-1 probe confirmed
the esm rejects equations mixing patch- and column-shaped operands).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `clm_templates.esm` | 4 new pointwise laws: `ground_temperature_blend` (t = frac_sno_fluxes·t_snow + (1 − frac_sno_fluxes − frac_h2osfc)·t_soil + frac_h2osfc·t_h2osfc, for t_grnd0 and lw_grnd), `stefan_boltzmann_emission` (emg·sb·T4), `flux_temperature_correction` (pre + tinc·deriv, the five flux corrections of :192-204) and `stefan_boltzmann_emission_dT` (4·emg·sb·T³, the dT-derivative whose product with tinc is the longwave increment) | `SoilFluxesMod.F90`:176-178, 344-346, 350, 471-472, 192-204; `shr_const_mod.F90`:23 |
| `clm_soilfluxes.esm` | Model `SoilFluxes` — the full SoilFluxes update, 32 equations (the top-layer water h2osoi_liq+h2osoi_ice is named once as the `h2osoi_top` intermediate, shared by the evaporation limit and the partition) | `SoilFluxesMod.F90`:38-518 |
| `stage1/analysis/clm_soilfluxes_derived.py` | Produces the expected values of the derived (not dump-column) assertions (lw_grnd, qflx_ev_snow/ev_soil, evaporation_limit, qflx_soliddew_to_top_layer, lwrad_emission/increment, eflx_lwrad_out, qflx_evap_tot) by evaluating the transcribed equations at the pinned dump inputs | — |

The `qflx_*` moisture fluxes are declared as mass fluxes in `kg/m^2/s` (1 mm/s
≡ 1 kg/m^2/s of water): the esm unit checker enforces dimensional consistency
of every equation, and only the mass-flux form makes the latent-heat products
qflx·htvp [J/kg] = W/m^2 and the flux corrections pre + tinc·cgrndl (cgrndl in
kg/m^2/s/K) consistent. The values are numerically identical to the Fortran's
mm/s.

## Laws

**Previous-time-step ground temperature and T⁴ (SoilFluxesMod.F90:175-183,
344-346).** With j = snl+1 the top layer,

```
t_grnd0 = frac_sno_fluxes·tssbef(c,snl+1) + (1 − frac_sno_fluxes − frac_h2osfc)·tssbef(c,1) + frac_h2osfc·t_h2osfc_bef
tinc    = t_grnd − t_grnd0
lw_grnd = frac_sno_fluxes·tssbef(c,snl+1)⁴ + (1 − frac_sno_fluxes − frac_h2osfc)·tssbef(c,1)⁴ + frac_h2osfc·t_h2osfc_bef⁴
```

(the Fortran's no-snow branch of t_grnd0, F90:180, is the same blend at
frac_sno_fluxes = 0; lw_grnd is computed unconditionally at :344-346). Both
blends are the `ground_temperature_blend` template. On every traced row
frac_sno_fluxes = 1, frac_h2osfc = 0 and snl = −9, so t_grnd0 = lw_grnd¹ᐟ⁴ =
tssbef(c,−8), and t_grnd0 equals the bareground t_grnd to round-off on all
192 rows.

**Flux corrections (F90:189-206).** The pre-correction fluxes from
BareGroundFluxes/CanopyFluxes are shifted to the present soil temperature, all
five by the same law pre + tinc·deriv (the `flux_temperature_correction`
template, applied with deriv = cgrnds for the sensible flux and cgrndl for the
four mass fluxes):

```
eflx_sh_grnd  += tinc·cgrnds
qflx_evap_soi += tinc·cgrndl
qflx_ev_snow  += tinc·cgrndl
qflx_ev_soil  += tinc·cgrndl
qflx_ev_h2osfc += tinc·cgrndl
```

(non-urban branch; the urban branch sets qflx_ev_soil = qflx_ev_h2osfc = 0 and
qflx_ev_snow = qflx_evap_soi, F90:197-200 — see gaps). The traced
clm_bareground.csv rows carry the pre values and clm_soilfluxes.csv the post
values; the pins reproduce the post values to round-off on all 192 rows.
cgrnds/cgrndl are inputs (computed by FrictionVelocity/BareGroundFluxes; the
dump carries identical values in both files). The Fortran comments cgrndl as
w/m2/k but applies it directly to the mass fluxes, so its effective units are
mm/s/K ≡ kg/m^2/s/K.

**Liquid/solid partition (F90:208-274).** The four partition fluxes
(qflx_liqevap_from_top_layer, qflx_solidevap_from_top_layer,
qflx_soliddew_to_top_layer, qflx_liqdew_to_top_layer) are computed from the
signed qflx_ev_snow: for evaporation (≥ 0) split by the top-layer
liquid/(liquid+ice) ratio, for dew (< 0) all solid when t_grnd < tfrz and all
liquid otherwise. Every traced row is dew at t_grnd ≈ 238-248 K < tfrz, so
qflx_soliddew_to_top_layer = |qflx_ev_snow| and the other three are 0 on all
192 rows. The evaporation branch (and the soil solidevap limit of F90:317-327)
is transcribed but unpinned (see gaps).

**Snow-evaporation limit (F90:276-297).** For snow layers (j < 1),
evaporation_limit = (h2osoi_ice(c,j) + h2osoi_liq(c,j))/(frac_sno_fluxes·dtime)
and qflx_ev_snow is capped at it, with the excess returned to the soil
evaporation and sensible heat (energy-conserving):

```
if qflx_ev_snow > evaporation_limit:
    qflx_ev_snow  = evaporation_limit
    qflx_evap_soi -= frac_sno_fluxes·(demand − limit)
    eflx_sh_grnd  += frac_sno_fluxes·(demand − limit)·htvp
```

The .esm transcribes this without a self-referencing conditional as
`qflx_ev_snow_excess = max(qflx_ev_snow_raw − evaporation_limit, 0)` and the
limited flux, the two corrections and the top-layer partition written linearly
in the excess — exact, because the limiting branch's partition values are the
general-branch values evaluated at the limited qflx_ev_snow (F90:292-293 vs
:238-243 coincide). The limit never fires on the traced rows (all dew; the
limit is ≈ 3e-3 mm/s against a −2e-7 mm/s dew rate), so the firing branch is
unpinned (see gaps).

**Ground heat flux (F90:341-371).**

```
eflx_soil_grnd = (1 − frac_sno_fluxes)·sabg_soil + frac_sno_fluxes·sabg_snow
               + dlrad + (1 − frac_veg_nosno)·emg·forc_lwrad
               − emg·sb·lw_grnd − 4·emg·sb·t_grnd0³·tinc
               − (eflx_sh_grnd + qflx_evap_soi·htvp)
```

with the two emission terms factored as the `stefan_boltzmann_emission`
template and the increment written as
`lwrad_increment = stefan_boltzmann_emission_dT(emg, sb, t_grnd0)·tinc` — the
dT-derivative of the emission law evaluated at the previous-time-step
temperature, times the temperature increment (F90:350,472; the urban
eflx_lwrad_del of :361 is the same linearization at t_grnd).

**Totals and history variables (F90:373-394).**

```
eflx_sh_tot   = eflx_sh_veg + eflx_sh_grnd + eflx_sh_stem
qflx_evap_tot = qflx_evap_veg + qflx_evap_soi
eflx_lh_tot   = hvap·qflx_evap_veg + htvp·qflx_evap_soi
qflx_evap_can = qflx_evap_veg − qflx_tran_veg
eflx_lh_vege  = qflx_evap_can·hvap      (reuses the :390 intermediate)
eflx_lh_vegt  = qflx_tran_veg·hvap
eflx_lh_grnd  = qflx_evap_soi·htvp
```

**Outgoing longwave and skin temperature (F90:451-495).**

```
eflx_lwrad_out = ulrad + (1 − frac_veg_nosno)·(1 − emg)·forc_lwrad
               + (1 − frac_veg_nosno)·emg·sb·lw_grnd + 4·emg·sb·t_grnd0³·tinc
eflx_lwrad_net = eflx_lwrad_out − forc_lwrad
if frac_veg_nosno = 0:  t_skin = (lw_grnd)¹ᐟ⁴
```

(the vegetation branch of t_skin is set in CanopyFluxes and enters as the
input t_skin_veg).

## Dump analysis

### The Stefan-Boltzmann constant: sb = 5.67e-8 exactly

CESM's `SHR_CONST_STEBOL = 5.67e-8` (shr_const_mod.F90:23) — the rounded
value, not the CODATA 5.670374419e-8 — and the dumps pin it: with sb = 5.67e-8
the eflx_soil_grnd and eflx_lwrad_net transcriptions reproduce
clm_soilfluxes.csv to 7.9e-15 / 4.6e-16 relative on all 60 columns-1-5 rows,
and the sb implied by the dump (from eflx_soil_grnd + eflx_lwrad_net, which
cancel everything except the emission terms) is 5.67000000000000048e-08 — the
5.67e-8 literal bit-exactly. This is a provenance note (the value follows the
Fortran), not a bug; it is recorded here because a stage-3 coupler or a future
component that assumes the CODATA value would break the pins.

### forc_lwrad: sourced from the DATM dump

The traced SoilFluxes run receives its downward longwave from the DATM
component (med passes Faxa_lwdn with merge 'copy'), and the CLM dump does not
carry forc_lwrad. The datm_gswp3_exports.csv Faxa_lwdn column supplies it,
with the step offset verified through forc_t: clm step s consumes the DATM
exports of datm step s−1 (clm forc_t equals datm Sa_tbot to round-off on all
192 rows; the datm dump's points 1-16 are the identical pole gridcells, so
point n = 1 stands for every traced column). This identity holds only for the
traced columns 1-5 — columns 6-16 have no forc_lwrad, and their
forc_lwrad-dependent outputs are unpinned (see gaps).

### The 192-row uniformity

Every traced row is a bare istsoil column (ptype = 0, ctype = 1) with
snl = −9, frac_sno_fluxes = 1, frac_h2osfc = 0, frac_veg_nosno = 0, emg =
0.97 (frac_sno_albedo = 1), htvp = 2,834,700 J/kg (the frozen top layer with
h2osoi_liq = 0, h2osoi_ice > 0: htvp = hsub = LATICE + LATVAP =
3.337e5 + 2.501e6, BiogeophysPreFluxCalcsMod.F90:354-358), dlrad = ulrad = 0
(BareGroundFluxesMod.F90:313-314) and qflx_ev_snow < 0 (dew). The inputs are
read from the dump rows as follows:

| Input | Source |
|---|---|
| tssbef_snow, tssbef_soil | clm_soiltemp.csv entry stage, j = −8 / j = 1 (t_soisno) |
| frac_sno_fluxes, snl | clm_soiltemp.csv entry stage (all rows 1.0 / −9) |
| t_grnd | clm_soilfluxes.csv (the post-solve input) |
| cgrnds, cgrndl | clm_soilfluxes.csv (= clm_bareground.csv bit-exactly) |
| eflx_sh_grnd_pre, qflx_evap_soi_pre, qflx_ev_snow_pre, qflx_ev_soil_pre | clm_bareground.csv |
| h2osoi_liq_top, h2osoi_ice_top | clm_soiltemp.csv exit stage, j = −8 |
| sabg_soil, sabg_snow | clm_surfrad.csv |
| emg, htvp | clm_run_imports.csv (preflux stage) |
| forc_lwrad | datm_gswp3_exports.csv Faxa_lwdn, datm step s−1, n = 1 (columns 1-5 only) |

The quantities that are not dump columns (qflx_ev_snow/qflx_ev_soil after
correction, the partition, evaporation_limit, qflx_ev_snow_excess, lw_grnd,
t_grnd0, tinc, lwrad_emission/lwrad_increment, eflx_lwrad_out, qflx_evap_tot)
are the f64 evaluations of the transcribed equations at these inputs,
cross-checked against every dump column they touch (all 192 rows, max rel err
7.9e-15). The .esm tests pin the subset that can fail (lw_grnd, qflx_ev_snow,
qflx_ev_soil, evaporation_limit, qflx_soliddew_to_top_layer,
lwrad_emission/lwrad_increment, eflx_lwrad_out, qflx_evap_tot) and label them
"derived (not a dump column)"; their expected values are produced by the
committed script `stage1/analysis/clm_soilfluxes_derived.py` (run it to
recompute the tuples: `python3 stage1/analysis/clm_soilfluxes_derived.py`).
The remaining derived quantities are not asserted because they cannot fail on
the traced rows: t_grnd0 is the fs = 1 blend (a copy of the overridden
tssbef_snow, equal to the dumped t_skin), tinc is the difference of the two
overrides t_grnd − tssbef_snow, qflx_ev_snow_excess and
qflx_liqevap/qflx_solidevap_from_top_layer/qflx_liqdew_to_top_layer are
identically zero on every traced row (all dew at t_grnd < tfrz, limit never
firing), and qflx_evap_can/eflx_lh_vege/eflx_lh_vegt vanish with the zero
vegetation inputs. Their branches are listed as coverage gaps below.

## Tests

All 4 tests use time_span {0,0} (instantaneous pins), rel 1e-12, and cite
their dump rows in the description:

| Test | Pins (assertions) | Dump rows |
|---|---|---|
| `clm_soilfluxes_step1_p1` | 17 assertions: the dump columns eflx_sh_grnd, qflx_evap_soi, eflx_soil_grnd, eflx_lwrad_net, t_skin, eflx_sh_tot, eflx_lh_tot, eflx_lh_grnd plus the derived (not a dump column) lw_grnd, qflx_ev_snow, qflx_ev_soil, evaporation_limit, qflx_soliddew_to_top_layer, lwrad_emission/increment, eflx_lwrad_out, qflx_evap_tot (clm_soilfluxes_derived.py) | step 1, patch/column 1; forc_lwrad from datm step 0 |
| `clm_soilfluxes_step12_p1` | 15 assertions, the same pins at a second longwave value (8 dump columns + 7 derived) | step 12, patch/column 1; forc_lwrad from datm step 11 |
| `clm_soilfluxes_step7_p8` | 13 assertions: the corrections/partition/totals at a positive tinc (+2.44e-1 K); eflx_soil_grnd/eflx_lwrad_* unpinned (no forc_lwrad for columns 6-16) | step 7, patch/column 8 |
| `clm_soilfluxes_step6_p16` | 13 assertions, the same no-LW pin at the coldest traced column (largest |tinc| of the pins) | step 6, patch/column 16 |

The step-1/step-12 column-1 pins fix sb = 5.67e-8 exactly (see above); the
step-7/step-8 and step-6/step-16 pins cover the correction/partition/totals
laws at columns without a longwave input.

## Coverage gaps

The dumps support exactly the pins above; the following branches are
transcribed or documented but not dump-pinned:

- **forc_lwrad for columns 6-16.** Only the traced columns 1-5 have a
  DATM longwave identity; eflx_soil_grnd, eflx_lwrad_out and eflx_lwrad_net
  are unpinned at columns 6-16 (the remaining equations are pinned there).
  A dump carrying forc_lwrad per column would close this.
- **qflx_ev_h2osfc and qflx_ev_h2osfc_pre.** The surface-water evaporation is
  not dumped (frac_h2osfc = 0 on all rows, so it is 0-weight anyway); the
  equation is transcribed with qflx_ev_h2osfc_pre as a default-0 input.
- **The snow-evaporation limit firing (F90:288-296).** Never fires on the
  traced rows (all dew; limit ≈ 3e-3 mm/s vs dew ≈ −2e-7 mm/s), so
  qflx_ev_snow_excess is identically zero and its assertion was dropped
  (2026-11-06 review). The transcription is exact by linearity; a traced step
  with evaporating snow would pin it. Requires h2osoi at the pre-SoilFluxes
  state, which the exit-stage soiltemp dump approximates (the traced rows'
  top-layer water is bit-identical from exit(step) to entry(step+1)).
- **The evaporation branch of the partition (F90:234-243) and the soil
  solidevap limit (F90:317-327, j = 1).** Every row is dew; the liquid-evap/
  sublimation split by the top-layer liquid/(liquid+ice) ratio is untested
  (qflx_liqevap/qflx_solidevap_from_top_layer are identically zero on all
  traced rows; their assertions were dropped 2026-11-06).
- **The liquid-dew branch (F90:248, t_grnd ≥ tfrz with qflx_ev_snow < 0).**
  All traced rows have t_grnd ≈ 238-248 K < tfrz, so qflx_liqdew_to_top_layer
  is identically zero; its assertion was dropped (2026-11-06 review). A
  thawing row (t_grnd > tfrz, dew) would pin it.
- **The urban branch (F90:197-200, 356-371, 490-495) and the `_u`/`_r`
  aliases.** The traced columns are all istsoil; the urban branch (with
  eflx_lwrad_del, the wasteheat/AC/traffic/ventilation terms and the urban
  t_skin of F90:506-510) is not transcribed.
- **frac_h2osfc > 0 and t_h2osfc_bef, and a mixed snow/soil blend
  (0 < frac_sno_fluxes < 1).** frac_h2osfc = 0 and frac_sno_fluxes = 1 on all
  rows, so the surface-water terms of the blends are 0-weight and the blends
  degenerate to t_grnd0 = tssbef_snow (a copy of the overridden input — the
  t_grnd0/tinc assertions were dropped for that reason, 2026-11-06 review; the
  blend law itself is still pinned through lw_grnd, and the dumped t_skin
  pins the fourth-root relation). The no-snow (snl = 0) branches of t_grnd0
  and t_skin and the vegetation branch (frac_veg_nosno = 1: t_skin_veg,
  eflx_sh_veg/stem, qflx_evap_veg/tran_veg and the derived
  qflx_evap_can/eflx_lh_vege/eflx_lh_vegt) are unpinned (the last three were
  previously asserted only at their all-zero endpoint, which cannot fail, so
  those assertions were dropped 2026-11-06).
- **errsoi (F90:399-446).** The soil energy-balance check needs the 27-layer
  sums over (t_soisno − tssbef)/fact, xmf and the h2osfc terms; not dumped
  (and its conservation content is already covered by the
  eflx_soil_grnd/eflx_lwrad_net pins for columns 1-5, whose sum is the
  forc_lwrad-independent energy-balance identity sabg_eff − eflx_sh_grnd −
  qflx_evap_soi·htvp, verified to 6.2e-16 on all 192 rows).
