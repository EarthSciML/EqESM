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
constituents `q(c)`, plus the call's `ptend` (`ptend_s/u/v`, `ptend_q(c)`,
flags `ls/lu/lv/lq(c)`, `top_level/bot_level`) and the timestep `dt`:

1. **Window**: `win(k) = 1` for `top_level <= k <= bot_level` (1..32 at
   every active site in this run); everything below is inside the window,
   and the update is applied only there (F90:288-329, 416-420, 445-448).
2. **u, v**: `u_after = u + lu*win*ptend_u*dt`, `v_after = v + lv*win*ptend_v*dt`
   (F90:287-301).
3. **Constituents** (F90:313-334, the `constituent_update` template): each
   `q(c)` is updated by the qneg3 floor `q_after = max(q + lq*win*ptend_q*dt, qmin)`
   — except the number-concentration species NUMICE/NUMLIQ/NUMRAI/NUMSNO,
   which are **clamped** instead of floored: `min(max(q + lq*win*ptend_q*dt, 1e-12), 1e10)`
   (they are not qneg3'd; qmin is not consulted).  The `qmin` table is the
   chemistry.F90 one: Q = 1e-12, the chem gases and aerosol-mass species =
   1e-36, num_a1..a4 = 1e-5, and CLDLIQ/CLDICE/RAINQM/SNOWQM = 0 (qmin = 0
   is a no-op floor; those two water species are instead handled by step 5).
4. **CLDLIQ/CLDICE zeroing** (F90:353-375, the `state_cnst_min_nz_mask`
   template): at the `convect_deep`/`zm_conv_tend` sites only
   (`cld_min_nz` = 1), wherever the post-clamp `qclamp_CLDLIQ < 1e-36` the
   pair (CLDLIQ, NUMLIQ) is zeroed together, and likewise (CLDICE, NUMICE)
   with the 1e-36 threshold.  This run: CLDLIQ zeroed on ~29 of 32 levels,
   CLDICE on ~8; the lev-21 CLDLIQ value 2.0e-28 is *not* zeroed (it
   exceeds the threshold).
5. **Temperature** (F90:415-421): `T_after = T + ls*win*ptend_s*dt/cpairv`,
   with `cpairv = cpair` because `dry_air_species = ''` (air_composition.
   F90:354).
6. **Geopotential** (F90:423-449, geopotential_temp.F90:133-216): when
   `derive_new_geopotential` is true — `ptend%ls` set, or any
   thermodynamically active water species (this run all five: Q, CLDLIQ,
   CLDICE, RAINQM, SNOWQM) has a set `lq` — the heights are recomputed
   hydrostatically from the surface up:
   `f_tv = rog*tv` with `rog = rairv/gravit`,
   `tv = T_after*tvfac` (the `virtual_temperature_factor` template folds
   the wet-to-dry `qfac` and the dry-mixing-ratio sum),
   `zm_after = integral(f_tv, z_sfc -> lev)` and
   `zi_after = integral(f_tv, z_sfc -> lev_nodes)` — cumulative integrals
   with measure `dz = pdel/pmid` (the SE non-lagrangian branch; the
   lagrangian `piln` branch is not used by the SE dycore).  Otherwise
   `zm_after = zm`, `zi_after = zi` (identity branch, exercised at
   aero_calcsize).
7. **Static energy** (F90:445-448): `s_after = T_after*cpairv + gravit*zm_after + phis`
   in the window, unchanged outside it.
8. **Tend accumulation**: `tend%dudt += ptend_u`, `tend%dvdt += ptend_v`,
   `tend%dtdt += ptend_s/cpairv` in the window when the flag is set
   (F90:287-301, 417-419).  `tend%flx_net` is *not* touched by
   `physics_update` (accumulated in physpkg.F90:1527/3048 outside the law).

All four of these steps are gated per constituent/site by the dumped flags,
so zero-tendency sites (energy_fixer, aoa_tracers, co2_cycle, rayleigh,
gw, qbo, lunar, iondrag, …) need no pinning: their `physics_update` calls
are identity transformations, which the chaining check covers for free.

## Axis convention: `lev` is surface-first

The model's `lev` axis is **surface-first**: index 1 = the surface layer =
dump lev 32.  This is the orientation of the EarthSciDiscretizations
`column_nonuniform_1d` grid, whose cumulative-from-bottom integral rules
lower the hydrostatic sums (the CAM dump is TOA-first: dump lev 1,
`pmid ≈ 277 Pa`, `zm(1) ≈ 38070 m`).  Every `[lev]` test array is
therefore the dump's lev-1..32 column reversed; `[lev_nodes]` (zi, zi_after)
is reversed the same way, with model interface 1 = the surface (0) and
interface 33 = the top.  The convention is asserted in the model metadata
and documented on every shape-bearing variable.

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

- `constituent_update` (qneg3/num-clamp step), `state_cnst_min_nz_mask`
  (the 1e-36 zeroing masks) and `virtual_temperature_factor` (the
  wet-to-dry tv factor) are applied **by reference** from
  `cam_templates.esm` (appended, append-only file), as is the
  `column_nonuniform_1d` integral-rule pair imported from
  EarthSciDiscretizations.
- Constants (`rair`, `gravit`, `cpair`, `rh2o`, `z_sfc`) are parameters
  with the physconst values; `rairv = cpairv = constant` because
  `dry_air_species = ''` and there are no subcolumns (air_composition.
  F90:354), documented on the variables.

## Tests (572 assertions, rel 1e-12, `time_span {0,0}`)

Eleven tests pin the law at six call sites across the step and the tend
accumulation over both physics phases, all from `stage1/dumps/cam/
cam_phys_tend.csv`, `cam_phys_tend_q.csv`, `cam_tend_total.csv`:

| test | assertions | pins |
|---|---|---|
| `convect_deep_step1/2/24/48` | 57 each | the full law (ls/lu/lv + all 34 constituents) at proc zm_conv_tend, incl. the CLDLIQ/CLDICE zeroing; steps 1/2/24/48 span the dump window |
| `clubb_macmic1_step1` | 73 | the first CLUBB substep (macmic_it 1 of 3, ptend pre-scaled by 1/3) |
| `microp_macmic1_step1` | 72 | the first MG2 substep (ls only) |
| `aero_calcsize_step1` | 61 | the q-only site (ls = lu = lv = 0): the derive=0 geopotential identity branch and the num_a1 1e-5 floor (9.588e-6 → 1e-5) |
| `radiation_step1` | 17 | ls-only site, the last tphysbc call |
| `vdiff_step1` | 96 | lu/lv/ls + all 34 constituents, the last tphysac call |
| `tend_accumulation_bc_step1` | 15 | dudt/dvdt/dtdt over the 21 tphysbc calls (15 unique sites; CLUBB and microphysics each called 3× with their 1/3-scaled ptends) |
| `tend_accumulation_ac_step1` | 10 | dudt/dvdt over all 26 call sites vs the cam_tend_total.csv phase=ac row (dtdt excluded — see the SE-scaling gap below) |

Site sequence (26 physics_update calls per step): energy_fixer, dadadj,
convect_deep, convect_shallow, ssatcontrail, clubb#1..3 / microp#1..3
interleaved, aero_calcsize, aero_wetdep, convect_deep_2, radiation (all
tphysbc) | aoa_tracers, co2_cycle, chem, vdiff, rayleigh, aero_drydep,
carma_tend, gw, qbo, lunar, iondrag (tphysac).

**Chaining**: every test's expected values are chained to the *next* call
site's dump rows (`state_before(next site) == state_after(this site)`), as
labelled consistency checks; the forward-scan `chain_q` handles
constituents absent from the next site's lq (chem gases chain
clubb#1 → clubb#2, aerosols microp#1 → microp#2).  Independently,
`scratch/verify_chain.py` verified the full 48-step × 26-site chain for
every level and constituent at ≤ 4e-16 relative (geopotential zm/s
included) — the dump columns are mutually consistent with the law.

## Coverage gaps and notes

- **ac-phase `dtdt` SE scaling** (physpkg.F90:2026-2030): between the bc
  and ac phases the SE dycore scales the temperature tendency by
  `cpairv/cp_or_cv_dycore`, and `cp_or_cv_dycore` is not dumped.  The ac
  `dtdt` row of cam_tend_total.csv therefore does not equal the raw bc+ac
  `ptend_s/cpairv` sum (scaling ≈ 0.985-1.0, site-dependent), so
  `tend_accumulation_ac_step1` pins only dudt/dvdt (exact, 0 error) and
  documents the gap; `dtdt` ac accumulation is transcribed as the unscaled
  law and the scaling is *not* guessed.  A future dump of
  `cp_or_cv_dycore` would close it.
- **sslt_rebin_adv mutation**: `sslt_rebin_adv` mutates the ncl_/num_
  species between convect_shallow and ssatcontrail (rebinning), so the
  chain has a documented discontinuity there (ncl_a1/a2/a3, num_a1/a2/a3);
  the chaining checks skip that pair and the law itself is unaffected.
- **Dump 'E'-dropping quirk**: gfortran `es26.17` writes 3-digit exponents
  without the 'E' (`2.39688771960016664-252`) in cam_phys_tend_q.csv (310
  rows); the parsers insert the 'E' before the sign.  Dump quirk, not a
  CAM bug — recorded here, not in `bugs.md`.
- **qmin provenance**: the qmin table is the chemistry.F90 one
  (chemistry.F90 qmin array), read off for the 34 dumped constituents;
  Q = 1e-12 never engages in this run (Q stays ≫ 1e-12), the num_a1 floor
  engages at aero_calcsize, and the num clamps engage at convect_deep
  (NUMICE 9.9999999e-13 → 1e-12).
- **zi not dumped**: the identity branch `zi_after = zi` is pinned only
  structurally (aero_calcsize pins `zm_after = zm`); zi itself is a
  coverage gap (the derive=1 branch recomputes zi from the integral, which
  is fully pinned).
- **Derived-input labelling**: `vdiff_step1` expectations are *derived*
  (the law's output; there is no next-site dump within the step, vdiff is
  the last call) and `radiation_step1` q-inputs are chained from the
  running state (water species are not dumped at radiation; the chain
  propagation was verified against every dump).  Both tests are labelled
  as such in their descriptions.
- **macmic pre-scaling**: the CLUBB/MG2 ptends are dumped already carrying
  the 1/3 macmic pre-scaling (physpkg.F90:2680/2723/2886), so the model
  consumes them as-is; the `clubb_macmic1_step1`/`microp_macmic1_step1`
  descriptions state this.
- **Untranscribed branches**: the WACCM-X H/H2 caps and
  `cam_thermo_dry_air_update` (F90:459-488) and the lagrangian `piln`
  geopotential branch are not exercised by this run and are not
  transcribed.

## Regeneration

`scratch/gen_pins.py` extracts the test tuples from the dumps into
`scratch/pins/tests.json` (11 tests, 572 assertions; `scratch/verify_chain.py`
is the full-law cross-check).  The model body is hand-authored in
`scratch/build_cam_physics_update.py` (variables + equations), which
assembles `cam_physics_update.esm` = body + pins; the `.esm` file itself is
the artifact — the build script only pastes the two halves together.
