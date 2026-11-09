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
  `zm_flux`, `zm_flux_divergence`, `zm_subcloud_tendency`, `zm_evap_flux_step`,
  `zm_evap_snow_frac`, `estbl_es0`, `estbl_tmax`.

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
for k > maxg copied down from k−1 (F90:2620-2622). The two tendencies share
their structure and are factored into three templates of `cam_templates.esm`
(one law, one template — no suffix-repeated blocks): the layer-centered plume
flux `F = zm_flux = mu(xu−xhat) + md(xd−xhat)` (two applications of
`zm_mass_flux_excess`), the pressure-coordinate flux divergence
`zm_flux_divergence = (F[min(k+1, PVER)] − F[k])/dp[k]` (the k+1 gather clamped
to the surface level), and the subcloud tendency
`zm_subcloud_tendency = −F[k]/dsubcld`; each is applied once for the q species
(fd_q, sub_dqdt) and once for the s species (fd_s, sub_dsdt), and the
piecewise mask `ifelse(k >= maxg, subcloud, ifelse(k >= jt, flux form, 0))`
selects between them. The `D(F, lev)` discretization rule of
EarthSciDiscretizations (`face_flux_D_lev_supplied_faces`) was not used: it
lowers to `(F[k+1] − F[k])/dz[k]` on the geometric-height `column_nonuniform_1d`
grid, while the ZM divergence is per pressure thickness `dp` (hPa) on a
grid-less CAM level set, so the template carries the pressure-coordinate
discretization explicitly.

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
template `zm_evap_flux_step`. The surface rate is the bottom-interface rain
flux converted by density_fresh_water = 1000 kg/m³ (F90:234-236); `snow_evap`
is the same conversion of the snow flux and is 0 in this window (see Tests).
The qs comes from the estbl table lookup + `svp_to_qsat` (see below). All
`ifelse`/`min`/`max` nestings reproduce the Fortran's expression order so the
evaluation is bit-faithful where the inputs are bit-faithful.

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
applies). The ESD integral route was probed as well
(`scratch/probe_cumulative.esm`, the `integral_lev_cumulative_to_top` rule of
`column_nonuniform_1d`, with the triangular recurrence f = 1 − b,
b = c·sqrt(integral(f, lev..z_top) − f·dz/2)): the rule lowers, but the
coupled multi-variable recurrence is rejected the same way, with the exact
error

```
✗ /models/Probe: Observed dependency cycle: b -> f -> b. Each of these observed
  variables is defined in terms of the next, so no evaluation order satisfies
  every definition (esm-spec §4.9.6). Break the cycle by splitting one observed
  into a pre-value and a post-value.
```

so the vertical coupling cannot be written as an ESD integral in this CLI
either (esm-spec §4.9.6, the known multi-variable-recurrence limit).
`cam_zm_evap.esm` therefore pins the per-level step: the accumulated boundary
flux `flxprec_lo`/`flxsnow_lo` enters as a dumped input at every level (the
recurrence boundary), and the Fortran's running evaporation accumulator
`evpvint` is supplied as a **derived input** (it is not a dump column):

```
evpvint(k) = Σ_{j<k} (rprd(j) − ntprprd(j))·pdel(j)/gravit
```

Each increment equals `evpprec(j)·pdel(j)/gravit` of F90:195 exactly because
`ntprprd = rprd − evpprec` (F90:198); the derivation is in
`stage1/analysis/extract_cam_zm_tests.py` and reproduces the Fortran
accumulator bit-exactly. The `zm_evap_flux_step` template documents the
per-level step as the integral of the production over the layer mass, so the
vertical coupling is preserved as data in the tests and can be re-imposed as
an ESD integral in stage 3.

## The estbl saturation table

`qs` is `estblf(t1)` (wv_saturation.F90:191-254, the 251-knot Goff–Gratch
mixed-phase table at T = 127.16 + (i−1) K, with the estblf clamps to
[127.16, 375.16] K) followed by `svp_to_qsat` (eps·es/(p − omeps·es),
wv_saturation.F90:264-274, reused by reference from `cam_templates.esm`).
The lookup law is single-sourced: one application of the `estbl_es0` template
of `cam_templates.esm` in each of `cam_wv_saturation.esm` (WvSaturation,
bound to T) and `cam_zm_evap.esm` (ZMevap, bound to the layer temperature t1),
with the 375.16 K clamp as the named zero-parameter fragment `estbl_tmax`
(wv_saturation.F90:96). Both attempts to share the *table data* itself failed
in the current CLI (esm 0.7.0), with these exact errors:

- (a) table + lookup law moved into the library `cam_templates.esm`, consumed
  via `expression_template_imports` — validation of the importing document:

```
✗ : Failed to load ESM file: Structural validation error:
  [table_lookup_unknown_table] `table_lookup` references table `estbl`, which
  the document's `function_tables` block does not declare
  (at /models/WvSaturation/equations/0/rhs)
```

- (b) `cam_zm_evap.esm` composing WvSaturation by a `subsystems` ref — the
  mounted child's own lookup is checked against the *mounting* document's
  tables:

```
✗ : Failed to load ESM file: Structural validation error:
  [table_lookup_unknown_table] `table_lookup` references table `ft`, which the
  document's `function_tables` block does not declare
  (at /models/Parent/subsystems/Child/equations/0/rhs)
```

(probe with a minimal child table at `scratch/probe_subsys.esm` /
`scratch/probe_child.esm`.) `table_lookup` resolves only against the
declaring document's `function_tables` block and there is no by-reference
table mechanism, so the 251 values remain a per-document data copy: the
identical `function_tables.estbl` entry in `cam_wv_saturation.esm` and
`cam_zm_evap.esm` (data-identical, script-verified — the copy is extracted
verbatim by `stage1/analysis/inject_cam_zm_tests.py`), i.e. the same data the
Fortran interpolates. The evap tests pin the qs dependence through `tend_q`
(see Coverage).

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
bottom levels, flxsnow_hi at the melting level, the evanesced-column
prec_evap) carry abs 1e-24 in addition. The near-zero-region pins are either
kept with an explicit failure-mode explanation or dropped: `ntsnprd = 0`
below the snow region is kept because its zero is a three-term composite
(rprd·work2 − evpsnow − snowmlt) that fails if any of the branch states
(fsnow_conv = 0, flxsntm = 0, snowmlt = 0, flxsnow_hi = 0) misbehaves — the
test descriptions spell this out per column; `snow_evap = 0` is dropped as a
pin because it is the same zero region one step downstream of the pinned
`flxsnow_hi(22) = 0` melt-vanishing pin (documented in the descriptions).

| File | Column | Pins |
|---|---|---|
| cam_zm_q1q2 | steps 18, 40, 48 | dqdt/dsdt at levs 15 (jt), 17, 25, 32 (maxg); dlf at levs 15-16 |
| cam_zm_evap | steps 18, 22, 40, 48 | ntprprd/flxprec_hi across the raining column, the guard level 31 (steps 18/22), ntsnprd/tend_s_snwprd/tend_s_snwevmlt incl. the melting level 22 and the zero below the snow region (lev 25, lev 24 on step 22), tend_q/tend_s at the active-evaporation and Sundqvist-binding levels, prec_evap |

The step-18 and step-22 columns evaporate their precipitation completely
(prec_evap = 0; the max roundoff guard of F90:222-224 fires at lev 31 on
steps 18, 21, 22, 23 — residual fluxes of −1.2e-19 to −3.0e-19 kg/m²/s — so the
dumped flxprec_hi(31) is exactly 0). Steps 40/48 are the active columns
(prec_evap = 4.5e-9 and 1.0e-8 m/s). Every column melts its snow at lev 22, so
no snow reaches the surface in the whole window and `snow_evap` is unpinned
everywhere (see Coverage for the melt-region pins that do cover the
mechanism).

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
surface rain conversion. The zero `ntsnprd` below the snow region is exercised
at levs 24–25 of every column (see Tests for the failure modes).

**Not exercised** (documented gaps): the deltat saturation limit
`max(0, (qs−q)/deltat)` (0 of 480 levels; its branch is still in the min
chain), the k > maxg copy-down of q1q2 (maxg = pver on every dumped column),
the `new_snow` branch of zm_conv_evap (dead code: old_snow = .true. is
hardcoded at F90:108), intermediate clamp values of work1/work1b (the dumped
snow fractions are 0 or 1 only), a negative evpprec (the old_snow path has
no max(0, ·) on evpprec itself; it never goes negative in this window), and
`snow_evap` (the melting level 22 of every column melts all snow, so the
surface snow flux is identically zero in the window — the mechanism is pinned
by `flxsnow_hi(22) = 0` instead, and the `snow_evap = 0` assertion was dropped
as a redundant downstream zero).

## Verification pipeline

- `stage1/analysis/extract_cam_zm_tests.py` — reproduces the full Fortran
  chains (q1q2_pjr and zm_conv_evap_run with the estbl of
  cam_wv_saturation.esm) in Python against the dump, with a hard copy check
  at the assertion tolerance (rel 1e-12, abs 1e-24 for identically-zero dump
  values — never `==`): q1q2 to ~1e-15 rel (dl bit-exact), evap to ≤1.4e-14
  rel over all checked columns and prec_evap; writes the test tuples and the
  derived evpvint. The script fails loudly if a dump or law change would
  break the committed pins.
- `stage1/analysis/inject_cam_zm_tests.py` — materializes the dump-extracted
  tuples into the hand-authored tests blocks (the sanctioned extraction use;
  the equations, structure, descriptions and pin choices are hand-authored,
  and the script is idempotent — re-running it reproduces the committed tests
  verbatim, including the estbl table copy and the binder counts in the
  descriptions). Both scripts are cited in the test descriptions.
- Probes: `scratch/probe_tbl3.esm` (lev-shaped table_lookup + svp_to_qsat,
  4/4 at rel 1e-12 vs cam_qsat.csv), `scratch/probe_units.esm` (sqrt and
  index-at-PVER), `scratch/probe_cycle.esm` (the exact faq cycle rejection
  above), `scratch/probe_cumulative.esm` (the exact integral-rule cycle
  rejection above), `scratch/probe_subsys.esm`/`scratch/probe_child.esm`
  (the table-import failure above).
- Suite: baseline after the main merge 3175P/0F/0E; final suite on this
  branch 3171P/0F/0E (the four dropped `snow_evap = 0` pins) — see the commit
  messages for the per-file counts; this lane's files contribute 163
  assertions: cam_zm_q1q2 30, cam_zm_evap 68, cam_wv_saturation 65.

## bugs.md

The zm_conv_evap `ke_lnd`/`landfrac` dead-argument finding is recorded as
bugs.md entry 33 (kemask = ke unconditionally, identical in the vanilla CESM
checkout — not an instrumentation artifact); after the main merge it was
re-numbered `## TBD. ZM: ...` (main's #33–#37 kept, the manager numbers this
one).
