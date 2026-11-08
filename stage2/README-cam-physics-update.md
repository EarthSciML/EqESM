# CAM6 `physics_update` state law — stage-2 `.esm` implementation

`cam_physics_update.esm` implements the single state-update step that CAM6
applies after **every** parameterization, `physics_update` of
`components/cam/src/physics/cam/physics_types.F90:207-495`, pinned by the
stage-1 instrumented dumps of the SCAM run `C_JRA_t62_scam_inst` (all 48
steps of the dump window, `dt = 1800 s`).  The instrumented fork is the
cam6_4_212 checkout (`stage1/EqESM-cesm-cam`, branch
`instrument/cam-dumps`); dump documentation and validation: `stage1/runs/
CAM_LANE.md` (tphysbc/tphysac dump section), `stage1/analysis/
validate_cam_phys_dump.py`.

This is the *law only*: it consumes the `ptend` produced by a call site and
the pre-call state, and produces the post-call state.  It does not contain
any parameterization physics, so it is the first fully-assembled, fully
coupled piece of the CAM stage-2 chain — every other CAM component in this
repo feeds or consumes it.

## The law, per call

Given the pre-call state `(u, v, T, s, zm, zi, phis, pmid, pdel)` and the
constituents `q`, plus the call's `ptend` (`ptend_s/u/v`, `ptend_q`, flags
`ls/lu/lv/lq`, `top_level/bot_level`) and the timestep `dt`:

1. **Window**: `win(k) = 1` for `NLEV+1-bot_level <= k <= NLEV+1-top_level`
   (the Fortran's `do k = top_level..bot_level` in TOA-first indexing,
   F90:289; model indices are surface-first, so the window is reversed).
   Every active site in this run uses the full column (1, 32), so `win = 1`
   everywhere in the tests.
2. **u, v**: `u_after = u + lu*win*ptend_u*dt`, `v_after = v + lv*win*ptend_v*dt`
   (F90:287-301).
3. **Constituents** (F90:313-334, the `constituent_update` template applied
   over the `[lev, cnst]` index set): the raw update
   `qraw = q + win*ptend_q*dt` is computed **once** as a named intermediate
   (the Fortran writes it in both clamp branches), then, where `lq` is set,
   either the qneg3 floor `qclamp = max(qraw, qmin)` over the whole column,
   or the number-concentration clamp
   `min(max(qraw, 1e-12), 1e10)` inside the window for
   NUMICE/NUMLIQ/NUMRAI/NUMSNO (`num_clamp`, F90:320-329); `lq = 0`
   short-circuits to `q` unchanged.  The `qmin` table is the chemistry.F90
   one, carried as the `[cnst]` parameter default: Q = 1e-12, the chem
   gases and aerosol-mass species = 1e-36, num_a1..a4 = 1e-5, and
   CLDLIQ/CLDICE/RAINQM/SNOWQM = 0.
4. **CLDLIQ/CLDICE zeroing** (F90:353-375, the `state_cnst_min_nz_mask`
   template): at the `convect_deep`/`zm_conv_tend` sites only
   (`cld_min_nz` = 1 at the mass constituents CLDLIQ and CLDICE), wherever
   the post-clamp `qclamp(., CLDLIQ) < 1e-36` the (CLDLIQ, NUMLIQ) pair is
   zeroed together, and likewise (CLDICE, NUMICE) — expressed with the
   `[cnst]`-shaped pair masks `pair_cldliq`/`pair_cldice` selecting the
   members of each pair in the single `q_after` equation.  This run:
   CLDLIQ zeroed on ~29 of 32 levels, CLDICE on ~8.
5. **Temperature** (F90:415-421): `T_after = T + ls*win*ptend_s*dt/cpairv`,
   with `cpairv = cpair` because `dry_air_species = ''` (air_composition.
   F90:354).
6. **Geopotential** (F90:423-449, geopotential_temp.F90:133-216): when
   `derive_new_geopotential` is true — `ptend%ls` set, or any
   thermodynamically active water species (this run all five: Q, CLDLIQ,
   CLDICE, RAINQM, SNOWQM — the `wet_mask` subset) has a set `lq` — the
   heights are recomputed hydrostatically from the surface up:
   `f_tv = rog*tv` with `rog = rairv/gravit`,
   `tv = T_after*tvfac` (the `virtual_temperature_factor` template folds
   the wet-to-dry `qfac` and the dry-mixing-ratio sum; its `q_Q` input is
   the Q column of `q_after`, selected with a constant index over `cnst`),
   `q_wet_sum = sum over the wet_mask constituents of q_after` (one masked
   `faq` reduction over the `cnst` axis, not a hand-listed sum),
   `zm_after = integral(f_tv, z_sfc -> lev)` and
   `zi_after = integral(f_tv, z_sfc -> lev_nodes)` — cumulative integrals
   with measure `dz = pdel/pmid` (the SE non-lagrangian branch; the
   lagrangian `piln` branch is not used by the SE dycore).  Otherwise
   `zm_after = zm`, `zi_after = zi` (identity branch, not asserted — see
   the copy rule below).
7. **Static energy** (F90:445-448): `s_after = T_after*cpairv + gravit*zm_after + phis`
   in the window, unchanged outside it.
8. **Tend accumulation**: `tend%dudt += ptend_u`, `tend%dvdt += ptend_v`,
   `tend%dtdt += ptend_s/cpairv` in the window when the flag is set
   (F90:287-301, 417-419).  The per-call increment is the `dudt`/`dvdt`/
   `dtdt` equation of `PhysicsUpdate`; the accumulation over the call
   sequence is the **`TendAccumulation`** model in the same file: one
   masked `faq` sum over the `site` axis (the consecutive call sites of
   the phase), `dudt = sum over sites of lu*win*ptend_u` — the linear
   accumulation written as a sum over the call sequence.  `tend%flx_net`
   is *not* touched by `physics_update` (accumulated in
   physpkg.F90:1527/3048 outside the law).

All of these are gated per constituent/site by the dumped flags, so
zero-tendency sites (energy_fixer, aoa_tracers, co2_cycle, rayleigh, gw,
qbo, lunar, iondrag, …) need no pinning: their `physics_update` calls are
identity transformations, which the chaining check covers for free.

## Axis conventions

- **`lev` is surface-first**: index 1 = the surface layer = dump lev 32,
  the orientation of the EarthSciDiscretizations `column_nonuniform_1d`
  grid whose cumulative-from-bottom integral rules lower the hydrostatic
  sums (the CAM dump is TOA-first: dump lev 1, `pmid ≈ 277 Pa`,
  `zm(1) ≈ 38070 m`).  Every `[lev]` test array is the dump's lev-1..32
  column reversed; `[lev_nodes]` (zi, zi_after) is reversed the same way.
- **`cnst` is the dump's pcnst order** (34 constituents, 1-based): Q,
  CLDLIQ, CLDICE, NUMLIQ, NUMICE, RAINQM, SNOWQM, NUMRAI, NUMSNO, then the
  aerosol/chem species (bc_a1, bc_a4, DMS, dst_a1..3, H2O2, H2SO4,
  ncl_a1..3, num_a1..4, pom_a1, pom_a4, SO2, so4_a1..3, soa_a1..2, SOAE,
  SOAG).  Assertion coords name constituents by position (`{"cnst": j}`);
  the equations touch the axis only through the `[cnst]`-shaped parameter
  tables and the two constant column selections CLDLIQ = 2, CLDICE = 3.
- **`site` is the tphysbc call order** (15 sites: energy_fixer, dadadj,
  convect_deep, convect_shallow, ssatcontrail, clubb#1..3, microp#1..3,
  aero_calcsize, aero_wetdep, convect_deep_2, radiation), used only by
  the `TendAccumulation` model.

The geopotential integrals are written as `integral` (not aggregates) and
lowered by the imported rules
`column_nonuniform_1d/rules/integral_lev_cumulative_from_bottom.esm` and
`integral_lev_nodes_cumulative_from_bottom.esm` (exclusive prefix + half
layer for `lev`, exclusive prefix for `lev_nodes` — the exact Fortran forms
`zm = zi + rog*tv*hkk`, `zi(kint) = zi(kint-1) + rog*tv*hkl`, with the
`dz = pdel/pmid` measure bound to the rules' `dz` free name and `z_sfc = 0`
to their `z_sfc`).  The imported rules constrain the integrand to a bare
lev-shaped variable reference, so the integrand `f_tv = rog*tv` is a named
variable (the `ifelse(derive, integral(...), zm)` selection is applied
around the whole integral).

## Reuse

- `constituent_update` (qneg3/num-clamp step, now consuming the
  once-computed `qraw` named intermediate), `state_cnst_min_nz_mask` (the
  1e-36 zeroing masks, applied to the CLDLIQ/CLDICE columns) and
  `virtual_temperature_factor` (the wet-to-dry tv factor) are applied
  **by reference** from `cam_templates.esm`, as is the
  `column_nonuniform_1d` integral-rule pair imported from
  EarthSciDiscretizations.
- Constants (`rair`, `gravit`, `cpair`, `rh2o`, `z_sfc`) are parameters
  with the physconst values; `rairv = cpairv = constant` because
  `dry_air_species = ''` and there are no subcolumns (air_composition.
  F90:354), documented on the variables.

## Tests (169 assertions, rel 1e-12, `time_span {0,0}`)

Eight tests pin the law at six call sites plus the tend accumulation over
the whole tphysbc phase, all from `stage1/dumps/cam/cam_phys_tend.csv`,
`cam_phys_tend_q.csv`, `cam_tend_total.csv`:

| test | assertions | pins |
|---|---|---|
| `convect_deep_step18` | 36 | the full law (ls/lu/lv + 9 dumped constituents) at proc zm_conv_tend, incl. the CLDLIQ/CLDICE pair zeroing, the number clamp and the geopotential; step 18, a ZM-active step |
| `convect_deep_step40` | 37 | the same at step 40 (steps 1-17, 24-39, 41-48 are ZM-inactive with all-zero ptends and are not pinned) |
| `clubb_macmic1_step1` | 33 | the first CLUBB substep (macmic_it 1 of 3, ptend pre-scaled by 1/3), incl. the chem gases |
| `microp_macmic1_step1` | 18 | the first MG2 substep (ls only), incl. the aerosol floors |
| `aero_calcsize_step1` | 8 | the q-only site: the num_a1 1e-5 floor (9.588e-6 → 1e-5) and the 1e-36 mass floors |
| `radiation_step1` | 15 | ls-only site, the last tphysbc call (T/zm/s chained to aoa_tracers) |
| `vdiff_step1` | 13 | lu/lv/ls + all 34 constituents at the single active level (dump lev 32) |
| `tend_accumulation_bc_step18` (TendAccumulation) | 9 | dudt/dvdt/dtdt over the 15 consecutive tphysbc calls as a multi-site sum vs the cam_tend_total.csv phase=bc row |

Suite status: `./esm test stage2/` runs **2831P/0F/0E** — the 169
physics_update assertions (rel 1e-12) on top of the 2662-assertion stage-2
baseline (was 3234P with the old 572-assertion generator-built file).

**Pinning rules** (the review fixes): every assertion is pinned at a few
informative levels (the active region, one deep level) — never every
level — and **only where the law changes the value**: the extractor
(`stage1/analysis/extract_cam_physics_update_pins.py`) drops any assertion
whose expected value equals the overridden input it derives from
(T_after vs T, u_after vs u, zm_after vs zm, s_after vs s, q_after vs q).
That removes the old unfailable pins: the all-zero convect_deep steps 1/2/24,
the zero-`ptend` levels of every profile, the per-call tend increments
(dudt = ptend_u with the summed ptend fed in is a copy), the derive=0
`zm_after = zm` copy at aero_calcsize, and the Python-derived vdiff
expectations (vdiff is now chained to the `gw` site's dump rows — the
intervening calls rayleigh/aero_drydep/carma_tend leave the state
unchanged, and gw is the first later site that dumps all 34 constituents).

**Chaining**: every expected value is chained to a later call site's dump
rows (`state_before(next site) == state_after(this site)`): convect_deep →
convect_shallow, clubb#1 → microp#1, microp#1 → clubb#2 / microp#2 (the
aerosols, which clubb#2 does not touch), aero_calcsize → aero_wetdep,
radiation → aoa_tracers, vdiff → gw.  Independently, the committed
extractor re-verifies the full 48-step × 26-site chain for every level and
constituent (worst relative error 0 for the column state, 4.3e-16 for the
geopotential zm/s) and the bit-identity of the phase=bc row with the
15-site sum, and fails loudly on any mismatch — the dump columns are
mutually consistent with the law.

## Coverage gaps and notes

- **ac-phase `dtdt` SE scaling** (physpkg.F90:2026-2030): between the bc
  and ac phases the SE dycore scales the temperature tendency by
  `cpairv/cp_or_cv_dycore`, and `cp_or_cv_dycore` is not dumped.  The ac
  `dtdt` row of cam_tend_total.csv therefore does not equal the raw bc+ac
  `ptend_s/cpairv` sum (scaling ≈ 0.985-1.0, site-dependent).  The old
  `tend_accumulation_ac_step1` test (dudt/dvdt only) was dropped in the
  review rework: its inputs were the phase sums, making dudt = ptend_u a
  copy assertion; the dtdt scaling gap remains open.  A future dump of
  `cp_or_cv_dycore` would close it.
- **zi not dumped**: the identity branch `zi_after = zi` is covered only
  structurally; zi itself is a coverage gap (the derive=1 branch recomputes
  zi from the integral, which is fully pinned).
- **Derived-input labelling**: `radiation_step1`'s five water-species q
  inputs are chained from the running state (water species are not dumped
  at radiation; the chain propagation was verified against every dump) and
  the test description labels them as derived inputs, not dump columns.
  No expected value in the suite is derived from the Python
  reimplementation anymore (vdiff is chained to gw).
- **sslt_rebin_adv mutation**: `sslt_rebin_adv` mutates the ncl_/num_
  species between convect_shallow and ssatcontrail (rebinning), so the
  chain has a documented discontinuity there (ncl_a1/a2/a3, num_a1/a2/a3);
  the chaining checks skip that pair and the law itself is unaffected.
- **Dump 'E'-dropping quirk**: gfortran `es26.17` writes 3-digit exponents
  without the 'E' (`2.39688771960016664-252`) in cam_phys_tend_q.csv (310
  rows); the parsers insert the 'E' before the sign.  Dump quirk, not a
  CAM bug — recorded here, not in `bugs.md`.
- **qmin provenance**: the qmin table is the chemistry.F90 one
  (chemistry.F90 qmin array), read off for the 34 dumped constituents and
  carried as the `qmin` parameter default; Q = 1e-12 never engages in this
  run (Q stays ≫ 1e-12), the num_a1 floor engages at aero_calcsize, and
  the num clamps engage at convect_deep (NUMICE 9.9999999e-13 → 1e-12).
- **Window orientation**: the `win` window is transcribed in the model's
  surface-first indexing (`NLEV+1-bot_level <= k <= NLEV+1-top_level`);
  all dumped windows are the full column (1, 32), so the orientation is
  untestable with this dump but faithful to F90:289.
- **macmic pre-scaling**: the CLUBB/MG2 ptends are dumped already carrying
  the 1/3 macmic pre-scaling (physpkg.F90:2680/2723/2886), so the model
  consumes them as-is; the `clubb_macmic1_step1`/`microp_macmic1_step1`
  descriptions state this.
- **Untranscribed branches**: the WACCM-X H/H2 caps and
  `cam_thermo_dry_air_update` (F90:459-488) and the lagrangian `piln`
  geopotential branch are not exercised by this run and are not
  transcribed.

## Regeneration

`stage1/analysis/extract_cam_physics_update_pins.py` re-verifies the whole
law against every dump row and emits the numeric `parameter_overrides` +
`assertions` of the tests (8 tests, 169 assertions) as JSON fragments; the
fragments are pasted into the hand-authored tests blocks of
`cam_physics_update.esm` — the model body (variables, equations, index
sets, templates) is hand-authored in the `.esm` file itself, and the
script touches numbers only.  The `.esm` file is the artifact; running the
script reproduces the exact committed tuples (and fails loudly if the
dumps ever disagree with the law).
