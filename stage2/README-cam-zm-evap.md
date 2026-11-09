# CAM ZM q1q2 tendencies and rain evaporation (stage 2)

This lane implements the two Zhang–McFarlane deep-convection pieces of the
evaporation branch of the ZM chain (PLAN.md lane `cam_zm_evap`), pinned by the
stage-1 SCAM dumps `cam_zm.csv` / `cam_zm_scal.csv` (window steps 18–23 and
40–48, all columns g = 1, jt = 15, maxg = 32 = pver):

- `stage2/cam_zm_q1q2.esm` — the q1q2_pjr temperature/moisture tendencies
  (`components/cam/src/atmos_phys/schemes/zhang_mcfarlane/zm_convr.F90:2584-2622`,
  called at F90:718-724).
- `stage2/cam_zm_evap.esm` — the zm_conv_evap_run rain evaporation step
  (`zhang_mcfarlane/zm_conv_evap.F90:118-236`, called from
  `physics/cam/zm_conv_intr.F90:671-679`).
- `stage2/cam_templates.esm` — appended (append-only): `zm_mass_flux_excess`,
  `zm_evap_flux_step`, `zm_evap_snow_frac`.

Levels are the CAM convention (lev 1 = model top, lev 32 = surface), exactly as
dumped; the models declare the level set `lev` over the file-level
metaparameter `PVER = 32`, so dumped arrays enter the tests unchanged.

## Laws

**q1q2** (all per level k, with ktm = jt and kbm = maxg of the single-column
gathered loop): the net condensation `emc = -cme` (F90:2586-2587, cme = cu − evp
of cldprp), and over the convection layer jt ≤ k < pver

```
dqdt = emc + ( mu(k+1)(qu−qhat)(k+1) − mu(k)(qu−qhat)(k)
             + md(k+1)(qd−qhat)(k+1) − md(k)(qd−qhat)(k) ) / dp      (F90:2596-2601)
dsdt = −rl/cp·emc + ( mu(k+1)(su−shat)(k+1) − mu(k)(su−shat)(k)
                   + md(k+1)(sd−shat)(k+1) − md(k)(sd−shat)(k) ) / dp (F90:2589-2594)
dl   = du(k)·qcde(k+1)                                                 (F90:2603-2604)
```

with `dp = 0.01·pdel` (F90:717, mb). At k = maxg the tendencies are replaced by
the subcloud form `(1/dsubcld)(−mu(x−xhat) − md(y−xhat))` (F90:2611-2619) and
for k > maxg copied down from k−1 (F90:2620-2622). The four flux terms are four
applications of the shared template `zm_mass_flux_excess`; the esm model
evaluates the laws pointwise per level with a level faq, the k+1 neighbor
gather clamped to `min(k+1, PVER)`, and the piecewise mask `ifelse(k >= maxg,
subcloud, ifelse(k >= jt, flux form, 0))`.

**Evap** (per level k, top down): the evaporation is the min of the Sundqvist
rate and three upper bounds (F90:158-178):

```
evpprec = min( min( min( max(0, (qs−q)/deltat),               flxprec_lo·gravit/pdel ),
                    (prec_gen−evpvint)·gravit/pdel ),                        evpprec_raw )
evpprec_raw = ke·(1−cldfrc)·max(1−q/qs, 0)·sqrt(flxprec_lo)                   (kemask = ke, F90:160-164)
```

then the snow chain with `old_snow = .true.` (F90:108): the melt split
`flxsntm = 0 / snowmlt = flxsnow_lo·gravit/pdel` where t1 > tmelt, else
`flxsntm = flxsnow_lo / snowmlt = 0` (F90:128-134); `evpsnow = evpprec·work1`
with work1 = min(max(0, flxsntm/flxprec_lo), 1) when flxprec_lo > 0 (F90:186-192);
`ntprprd = rprd − evpprec` (F90:198); `ntsnprd = rprd·work2 − evpsnow − snowmlt`
with work2 = max(fsnow_conv, work1b), zeroed when snowmlt > 0 (F90:204-215);
and the flux steps with the roundoff guard
`flxprec_hi = max(flxprec_lo + ntprprd·pdel/gravit, 0)` (F90:219-224), applied
once for the rain flux and once (with ntsnprd) for the snow flux via the shared
template `zm_evap_flux_step`. The surface rates are the bottom-interface fluxes
converted by density_fresh_water = 1000 kg/m³ (F90:234-236). The qs comes from
the estbl table lookup + `svp_to_qsat` (see below). All `ifelse`/`min`/`max`
nestings reproduce the Fortran's expression order so the evaluation is
bit-faithful where the inputs are bit-faithful.

## The precipitation-flux integral and the dependency cycle

The precipitation flux is the vertical integral of the net production
`flxprec(k+1) = flxprec(k) + ntprprd(k)·pdel(k)/gravit`, accumulated top down
(F90:118-120, 219-224). The coupling is *two-way*: the evaporation at a level
depends on the flux entering it (the falling-flux limit F90:173 and the
Sundqvist sqrt factor F90:164), and the flux entering a level is the integral
of production-minus-evaporation above it. The esm CLI cannot evaluate this
top-down recurrence; the coupled formulation is rejected with the exact error

```
✗ /models/<m>: Observed dependency cycle: evpprec -> flxprec_hi -> ntprprd -> evpprec.
  Each of these observed variables is defined in terms of the next, so no
  evaluation order satisfies every definition (esm-spec §4.9.6). Break the
  cycle by splitting one observed into a pre-value and a post-value.
```

(reproduced with `scratch/probe_cycle.esm`; the task's per-level-step fallback
applies). `cam_zm_evap.esm` therefore pins the per-level step: the accumulated
boundary flux `flxprec_lo`/`flxsnow_lo` enters as a dumped input at every level
(the recurrence boundary), and the Fortran's running evaporation accumulator
`evpvint` is supplied as a **derived input** (it is not a dump column):

```
evpvint(k) = Σ_{j<k} (rprd(j) − ntprprd(j))·pdel(j)/gravit
```

Each increment equals `evpprec(j)·pdel(j)/gravit` of F90:195 exactly because
`ntprprd = rprd − evpprec` (F90:198); the derivation is in
`scratch/extract_cam_zm_tests.py` and reproduces the Fortran accumulator
bit-exactly. The `zm_evap_flux_step` template documents the per-level step as
the integral of the production over the layer mass, so the vertical coupling is
preserved as data in the tests and can be re-imposed as an ESD integral in
stage 3.

## The estbl saturation table

`qs` is `estblf(t1)` (wv_saturation.F90:191-254, the 251-knot Goff–Gratch
mixed-phase table at T = 127.16 + (i−1) K, with the estblf clamps to
[127.16, 375.16] K) followed by `svp_to_qsat` (eps·es/(p − omeps·es),
wv_saturation.F90:264-274, reused by reference from `cam_templates.esm`).
The estbl table is duplicated as a `function_tables` entry in
`cam_zm_evap.esm` because the CLI rejects importing a function table by
reference (function tables are data, not template libraries —
`template_import_not_library`). The copy is extracted verbatim from the
committed `cam_wv_saturation.esm` table (script-verified identical, axis and
251 values), i.e. the same data the Fortran interpolates, and the evap tests
pin the qs dependence through `tend_q` (see Coverage).

## Notes on the dumped du units

The dump description labels `eu/du` as 1/s and the validator's mu-continuity
check (validate_cam_zm_dump.py) confirms it: `mu(k) − mu(k+1) =
(eu(k) − du(k))·dp(k)` with dp in mb holds with eu/du in their final 1/s
units. `cam_zm_q1q2.esm` declares `du` as 1/s accordingly, so
`dlf = du·qcde(k+1)` is dimensionally a mixing-ratio tendency (kg/kg/s),
consistent with the rliq identity `rliq = Σ dlf·dpp/gravit` of
zm_convr.F90:760-766. The mu/md fluxes are declared in hPa/s (mb = hPa).

## Tests

Assertion values come from the dump rows only (cited per test:
`cam_zm.csv`/`cam_zm_scal.csv` step and lev); `evpvint` is the derived input
above, `prec_gen = prec_convr·1000` is computed in-model (F90:111). Tolerances
are rel 1e-12; the identically-zero guard pins (flxprec_hi at the guard and
bottom levels, flxsnow_hi at the melting level, prec_evap/snow_evap of the
evanesced columns) carry abs 1e-24 in addition.

| File | Column | Pins |
|---|---|---|
| cam_zm_q1q2 | steps 18, 40, 48 | dqdt/dsdt at levs 15 (jt), 17, 25, 32 (maxg); dlf at levs 15-16 |
| cam_zm_evap | steps 18, 22, 40, 48 | ntprprd/flxprec_hi across the raining column, the guard level 31 (steps 18/22), ntsnprd/tend_s_snwprd/tend_s_snwevmlt incl. the melting level 22, tend_q/tend_s at the active-evaporation and Sundqvist-binding levels, prec_evap/snow_evap |

The step-18 and step-22 columns evaporate their precipitation completely
(prec_evap = 0; the max roundoff guard of F90:222-224 fires at lev 31 on
steps 18, 21, 22, 23 — residual fluxes of −1.2e-19 to −3.0e-19 kg/m²/s — so the
dumped flxprec_hi(31) is exactly 0). Steps 40/48 are the active columns
(prec_evap = 4.5e-9 and 1.0e-8 m/s).

## Coverage

**Exercised by the dump window** (480 levels over 15 columns): both melt-branch
sides — the melting level 22 of *every* column (t1 = 274.3–274.7 K > tmelt with
flxsnow_lo > 0: snowmlt > 0, work2 zeroed at F90:212, ntsnprd = −snowmlt,
flxsnow_hi = 0; step 43 also melts at lev 23, t1 = 280.2 K) and the non-melting
levels 17–21 (flxsntm = flxsnow, evpsnow = evpprec·work1 with work1 = 1, the
snow share of the evaporated rain); the Sundqvist raw limit (binds on 401 of
480 levels, 157 of them non-degenerately with flxprec_lo > 0 — steps 20–23
levs 23–26, steps 40–48 levs 17–32), the falling-flux limit (45 levels), the
remaining-precipitation limit (34 levels), the roundoff max guard, the q1q2
subcloud branch (lev 32 = maxg), the dl detrainment (levs 15–16), and the
surface conversion.

**Not exercised** (documented gaps): the deltat saturation limit
`max(0, (qs−q)/deltat)` (0 of 480 levels; its branch is still in the min
chain), the k > maxg copy-down of q1q2 (maxg = pver on every dumped column),
the `new_snow` branch of zm_conv_evap (dead code: old_snow = .true. is
hardcoded at F90:108), intermediate clamp values of work1/work1b (the dumped
snow fractions are 0 or 1 only), and a negative evpprec (the old_snow path has
no max(0, ·) on evpprec itself; it never goes negative in this window).

## Verification pipeline

- `scratch/extract_cam_zm_tests.py` — reproduces the full Fortran chains
  (q1q2_pjr and zm_conv_evap_run with the estbl of cam_wv_saturation.esm) in
  Python against the dump: q1q2 to ~1e-15 rel (dl bit-exact), evap bit-exact
  at step 18 and ~1e-15 elsewhere; writes the test tuples and the derived
  evpvint.
- `scratch/inject_cam_zm_tests.py` — materializes the dump-extracted tuples
  into the hand-authored tests blocks (the sanctioned extraction use; the
  equations, structure, descriptions and pin choices are hand-authored, and
  the script is idempotent).
- Probes: `scratch/probe_tbl3.esm` (lev-shaped table_lookup + svp_to_qsat,
  4/4 at rel 1e-12 vs cam_qsat.csv), `scratch/probe_units.esm` (sqrt and
  index-at-PVER), `scratch/probe_cycle.esm` (the exact cycle rejection above).
- Suite: 1970P/0F/0E (1868 baseline + 102 new assertions from this lane).

## bugs.md

The zm_conv_evap `ke_lnd`/`landfrac` dead-argument finding is recorded as
bugs.md entry 33 (kemask = ke unconditionally, identical in the vanilla CESM
checkout — not an instrumentation artifact).
