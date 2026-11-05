# CAM parallel lane (atmosphere, milestone M2 — CAM6 physics in SCAM)

Scope: **CAM6 physics** (cam6_4_212) as exercised by the single-column
FSCAMARM97 compset on the ne3np4 grid. Deliverable: per-process instantaneous
tendencies for every `physics_update` call in `tphysbc`/`tphysac`, cam_in/cam_out
surface exchange per step, the `wv_saturation` qsat law, and accumulated
tendency totals — as CSV dumps from an instrumented run, verified
solution-neutral against a pristine run, and mapped to Fortran `file:line`
so stage-2 `.esm` components can be authored and tested against them.

## Ownership / collision guards

- **Ours:** `stage2/cam_*.esm`, fork `stage1/EqESM-cesm-cam` (branch
  `instrument/cam-dumps`), `stage1/dumps/cam/`, `stage1/runs/cam_instrument.json`,
  `stage1/analysis/compare_cam_hist.py`, `stage1/analysis/validate_cam_dumps.py`,
  this file.
- **Do not touch:** pristine `stage1/cesm` (read-only reference), the CICE
  (`EqESM-cesm-cice`) / CTSM (`EqESM-cesm-ctsm`) forks, `PROGRESS.md`/`PLAN.md`/
  `jobs.json`.
- **Fork chain (all local, unpushed):** cam submodule `4f1ed524` →
  `bcfd15b3` (line-length/cnst_name fixes) → superproject `95d771a` (cam
  pointer) on `instrument/cam-dumps`. The fork is a `git worktree` of
  `stage1/cesm`; all 17 top-level submodules plus the nested ones needed by
  this build are checked out at the pinned commits (nested submodules that
  the pristine checkout had initialized: `cdeps/{fox,share/genf90}`,
  `cice/src` + `src/icepack`, `clm/src/fates`, `cam-sima` (sparse: 2 files),
  `cam/src/atmos_phys/schemes/{pumas,mmm,rrtmgp,ext,tempo}`, `cam/src/hemco/HEMCO`
  (+docs), `geoschem_src/KPP/standalone`, `cime/CIME/non_py/cprnc`).

## Run configuration

| | reference | instrumented |
|---|---|---|
| case | `FSCAMARM97_ne3` | `FSCAMARM97_ne3_cam_inst` |
| compset / grid | `FSCAMARM97` = 2000_CAM60%SCAMARM97_CLM50%SP_CICE%PRES_DOCN%DOM_SROF_SGLC_SWAV_SESP; `ne3_ne3_mg37` | same |
| source | pristine `stage1/cesm` | fork `stage1/EqESM-cesm-cam` (cam @ bcfd15b3) |
| CIME scripts | pristine `stage1/cesm/cime` (cime6.7.1) | same (case CIMEROOT = pristine) |
| pes | NTASKS=1 (SCAM: "single column mode must be run on 1 pe") | same |
| run | start 1997-06-18 23:29:45 (START_TOD=84585), STOP_OPTION=nhours STOP_N=24 (48 physics steps of 1800 s), REST_OPTION=never, CLM_FORCE_COLDSTART=on | same |
| physics | deep=ZM, shallow/macrop=CLUBB_SGS, cld_macmic_num_steps=3, MG microphysics (PUMAS), rrtmg radiation (SW/LW every 2 steps), trop_mam4 chemistry ACTIVE, use_gw_* all false, scm_relaxation on, iopfile ARM97_4scam.nc (2089 hourly times, first 1997-06-18 00Z), pver=32 | same |
| namelist deltas | — | `user_nl_cam`: `dump_cam=.true., dump_dir=.../stage1/dumps/cam/, dump_first_steps=48` |

History output: `cam.h0i` at 30-min cadence (48 records), plus `cam.h0a`.

### Environment incident (2026-11-02): `share/domains` missing from inputdata

Both initial submissions (jobs 11978227/11978230) crashed with SIGSEGV in
`esm_get_single_column_attributes` (cmeps `esm.F90:1245+`) because
`stage1/inputdata/share/domains/domain.lnd.ne3np4_gx3v7.230718.nc` had
disappeared from the shared inputdata tree (the whole `share/domains/`
directory is gone; `share/meshes/` intact). The file is not present on the
CESM inputdata SVN server (only fv/ne16-ne240 domain files are listed in
`share/domains/`). It was reconstructed from the grid's own mesh
(`share/meshes/ne3np4_ESMFmesh_c230714_cdf5.nc`, same July-2023 grid
generation): unstructured domain file (ni=488, nj=1) with `xc/yc` =
mesh element centers, `area` = mesh element areas.

Two versions were needed:

- **v1 (`mask=1, frac=1`, land column):** the driver's nearest-neighbor
  search picks element 348 (96.708°W, 36.520°N) for
  PTS_LAT=36.6/PTS_LON=262.5 and the mesh lookup matches it exactly
  (≤1e-4) — the same column the reference archive shows. But the archived
  reference logs (`drv.log.11827138`) show the *original* file produced
  `scol_lndfrac=0.0`/`scol_ocnfrac=1.0` at the column (the archived run
  was an all-ocean single column), and with v1's land column CLM got an
  active land point and its `check_for_errors`
  (`components/clm/src/cpl/utils/lnd_import_export_utils.F90:107-130`)
  aborted both runs at step ~17: CAM's exported longwave-down
  (`flwds`) collapses to 0.0 from step ~5 onward in this configuration
  (radiation steps 1/3 give 385.5/92.9 W/m², then 0), and CLM rejects
  `flwds ≤ 0`. The reference survived only because its land fraction was
  zero (no land gridcells → check loop empty). Recorded in `bugs.md` #27.
- **v2 (`mask=0, frac=0`, ocean column, in use):** reproduces the
  reference configuration (`scol_lndfrac=0`, `scol_ocnfrac=1`) while
  keeping the correct element-348 column coordinates; both runs complete
  the full 48 steps. The v1 land-column run's 16 complete steps of dumps
  are discarded (they are also from the pre-fix dump module, see below).

Verified: the archived reference `cam.h0i` has lat/lon = 36.51977632 /
263.29179607 (element 348), area 0.03817521 sr, so v1 and v2 select the
same column as the reference regardless of mask/frac. The original file
was not byte-restorable (not on the SVN server); flagged for the inputdata
owner. The driver's unchecked-error handling behind all of this is logged
as a CESM bug (`bugs.md` #27).

### Dump-module bug found and fixed (fork `d6230637`)

The first instrumented run exposed a bug in `cam_dump_mod.F90`: for call
sites that pass an *empty* `physics_ptend` (disabled schemes such as
`lunar`/`iondrag`; `active=0`), `ptend%name`, `ptend%top_level` and
`ptend%bot_level` are uninitialized, and the dump wrote the garbage bytes
(embedded NULs, non-ASCII) into `cam_phys_tend.csv` (512 rows per
inactive site per 16 steps; the q file was unaffected). Fixed by writing
`'none'`/`0`/`0` when the ptend has no active fields
(`cam_dump_mod.F90:344-356`). Rerun required for clean CSVs.

## Instrumentation (cam_dump_mod.F90, cam submodule @ bcfd15b3)

- `src/physics/cam/cam_dump_mod.F90` — dump module (namelist `cam_dump_nml`:
  `dump_cam` [logical, default .false.], `dump_dir`, `dump_first_steps`;
  read in `cam_dump_readnl` from `runtime_opts.F90` `read_namelist`, broadcast
  masterproc→all). Gated by `dump_enabled()` = master PET ∧ init ∧
  1 ≤ nstep ≤ dump_first_steps; all writes inside one `!$omp critical`.
  When `dump_cam=.false.` the module is a no-op ⇒ default builds are
  bit-identical to pristine CAM.
- `src/physics/cam/physpkg.F90` — 34 dump call sites (see table).
- `src/control/runtime_opts.F90` — `use cam_dump_mod, only: cam_dump_readnl` +
  call in `read_namelist`.
- `bld/namelist_files/namelist_definition.xml` — `cam_dump_nml` group
  (build-namelist writes all groups into `atm_in` automatically).
- Nested submodule gotchas hit during the build (recorded for future lanes):
  mkSrcfiles dedups sources by basename and is only regenerated when
  Filepath/Srcfiles/Depends are missing, so initializing a nested submodule
  *after* a failed build leaves a stale `Depends`; always wipe the component
  obj dir (`bld/atm`, `bld/ice`, `bld/gnu/.../clm`) after touching submodules.

## Call sequence (physpkg.F90, cam6_4_212 + instrumentation)

`tphysac` = lines 1352–2076 (after dynamics), `tphysbc` = 2078–3087
(before dynamics). Order within a step: **tphysbc** sites first, then
**tphysac** sites. clubb/microp run inside the macro–micro co-substep
loop (`cld_macmic_num_steps=3`), so they appear 3× per step
(`macmic_it=1..3`, ptend scale = 1/3), interleaved (clubb_i, microp_i).

**Sites gated OFF in this configuration** (dump calls exist in the source
but the enclosing conditional is false — confirmed against the dumps):

- `macrop` (physpkg.F90:2717): `macrop_scheme='CLUBB_SGS'` ⇒
  `use_macrop=.false.` (CLUBB does the macro); not called.
- `microp_sc` (2892): `use_subcol_microp=.false.`; not called.
- `carma_wetdep` (3042) and `carma_emis` (1561): no CARMA
  (`use_carma=.false.`); not called.
- `nudging` (1963): gated by `Nudge_Model .and. Nudge_ON`, not by
  `scm_relaxation`; both false. The SCM relaxation is applied in the
  **dynamics** (`apply_iop_forcing.F90`), not in physics, so there is no
  nudging physics site in this case.

Effective per-step site sequence (26 site-calls, 22 unique):

| # | site label | dump call | physics_update | process (file:line) |
|---|---|---|---|---|
| bc | (surface bc_in) | `cam_dump_surface('bc_in')` @2325 | — | cam_in at entry to tphysbc |
| bc | (qsat bc) | `cam_dump_qsat('bc')` @2326 | — | qsat(T,pmid) at tphysbc entry |
| 1 | energy_fixer | 2399 | 2400 | energy_fixer (physics/cam/energy_fixer.F90) |
| 2 | dadadj | 2466 | 2467 | dadadj (dry adiabatic adjustment) |
| 3 | convect_deep | 2506 | 2507 | ZM deep convection (convect_deep.ZM) |
| 4 | convect_shallow | 2558 | 2559 | CLUBB_SGS shallow (via clubb) |
| 5 | ssatcontrail | 2670 | 2671 | ssatcontrail |
| 6×3 | clubb | 2765 | 2767 | CLUBB_SGS (clubb_tend_cam) |
| 7×3 | microp | 2931 | 2933 | MG/PUMAS micro (state) |
| 8 | aero_calcsize | 3003 | 3004 | modal_aero_calcsize |
| 9 | aero_wetdep | 3024 | 3025 | wetdep (modal_aero_wetdep) |
| 10 | convect_deep_2 | 3049 | 3050 | second deep-convection call (deep convection for aerosols) |
| 11 | radiation | 3098 | 3099 | radiation (rrtmg SW/LW) |
| bc | (surface bc_out) | `cam_dump_surface('bc_out')` @3119 | — | cam_out after cam_export |
| bc | (tend_total bc) | `cam_dump_tend_total('bc')` @3132 | — | tend%dtdt/dudt/dvdt at end of tphysbc |
| ac | (qsat ac) | `cam_dump_qsat('ac')` @1491 | — | qsat(T,pmid) at tphysac entry |
| 12 | aoa_tracers | 1592 | 1593 | aoa_tracers |
| 13 | co2_cycle | 1608 | 1609 | co2_cycle |
| 14 | chem | 1644 | 1645 | CAM-chem (trop_mam4, superfast) |
| 15 | vdiff | 1715 | 1716 | vertical diffusion (vdiff) |
| 16 | rayleigh | 1757 | 1758 | rayleigh friction |
| 17 | aero_drydep | 1788 | 1789 | dry deposition (aero_drydep) |
| 18 | carma_tend | 1810 | 1811 | CARMA tendencies (carma_timestep_tend) |
| 19 | gw | 1845 | 1846 | gravity waves (disabled: use_gw_* all false) |
| 20 | qbo | 1876 | 1877 | QBO (disabled) |
| 21 | lunar | 1895 | 1896 | lunar tides (disabled) |
| 22 | iondrag | 1930 | 1931 | ion drag (disabled) |

Not called in this config: macrop (2717), microp_sc (2892), carma_wetdep
(3042), carma_emis (1561), nudging (1963) — see the gating notes above.
| ac | (surface ac_out) | `cam_dump_surface('ac_out')` @2095 | — | cam_out at end of tphysac |
| ac | (tend_total ac) | `cam_dump_tend_total('ac')` @2096 | — | full-step accumulated tend |

Notes:
- Site 4 (convect_shallow) has a **dead-code duplicate** `physics_update` at
  (old) line 2571: `ptend` is deallocated right after the first update at
  2559, so the second call is a no-op on a deallocated ptend. Documented
  quirk, not instrumented, not a bug (the state update happens once).
- Snapshot-only `physics_update` calls at 2606 (rk_stratiform_tend snapshot),
  2644, 2828 (state_snap copies) are for `cam_take_snapshot` diagnostics and
  are not instrumented (they never run in this configuration).
- `ptend%name` is `'none'` when a site returns an empty ptend (e.g. some
  disabled schemes); `active=1` only if ls/lu/lv/any-lq is set.
- Radiation runs every 2nd step (radtlw/radtsw=2): on the alternating step
  the `radiation` rows still appear (site called with an empty/zero ptend).

## Dump files (stage1/dumps/cam/, gitignored)

| file | row key | content |
|---|---|---|
| `cam_phys_tend.csv` | (step, site, lev) | state BEFORE each physics_update: T,u,v,s,pmid,pdel,ps,phis,zm,q_Q + ptend_s/u/v/q_Q (+0 outside [top_level,bot_level] or when flag unset), dt, scale, macmic_it/num, lat/lon, active, nlq |
| `cam_phys_tend_q.csv` | (step, site, lev, cnst) | q_before + ptend_q for every constituent in ptend%lq |
| `cam_surface_exchange.csv` | (step, phase) | bc_in: cam_in (asdir..cflx_Q); bc_out: cam_out after cam_export; ac_out: cam_out at end of tphysac |
| `cam_surface_cflx.csv` | (step, phase, cnst) | cam_in%cflx per constituent |
| `cam_qsat.csv` | (step, phase, lev) | qsat law: T, pmid → es, qs, gam (wv_saturation qsat_2D, Goff-Gratch mixed phase) |
| `cam_tend_total.csv` | (step, phase, lev) | tend%dtdt/dudt/dvdt at end of tphysbc (bc) and tphysac (ac) |
| `cam_vdiff_prof.csv` | (step, lev) | vdiff (site 15) solver context after tendencies_run: solve-grid mid/ifc_lo/ifc_hi (Coords1D(state%pint)), interface rhoi/tint/dpidz_sq, kvh/kvm/kvq + cgh/cgs at lo/hi interfaces, tau_damp_rate, dtk, post-solve u1/v1/s1, ptend_s/u/v, pdel/pdeldry, dragblj, tauresx/y, ksrftms, clubbtop, taux/tauy/shflux, dse_top |
| `cam_vdiff_q.csv` | (step, lev, cnst) | vdiff tracer context: q0 (pre-solve), q1 (post-solve q_tmp), ptend_q, cflx, do_diffusion flag — one row per constituent (pcnst = 34) |

Row counts for the 48-step run (verified): `cam_phys_tend` 48 × 26 × 32 =
39936 rows (26 = 20 single sites + clubb×3 + microp×3); `cam_phys_tend_q`
48 × Σ(site lq count) × 32; `cam_qsat`/`cam_tend_total` 48 × 2 × 32 =
3072 each; `cam_surface_exchange` 48 × 3 = 144; `cam_surface_cflx`
48 × pcnst = 1632 (pcnst = 34); `cam_vdiff_prof` 48 × 32 = 1536;
`cam_vdiff_q` 48 × 32 × 34 = 52224.

All CSVs: `# description =` / `# columns =` headers naming columns + units +
Fortran file:line, one row per sample, no blank lines, no trailing commas,
es26.17 reals.

## Validation

- [x] **Solution neutrality (VERIFIED 2026-11-02):** instrumented h0i
  bit-identical to pristine over all 48 records and all variables, modulo
  the `case` global attribute and the `time_written` char variable
  (write timestamps) — `compare_cam_hist.py` reports exactly those two
  differences. Verified twice: (pristine 11999519 vs inst 11999518) and
  (pristine 11999519 vs inst 12001684, the final CSV format). Re-verified
  2026-11-03 after adding the vdiff solver-context dump (run 12048430):
  bit-identical again (same two differences only).
- [x] **Structural:** `python3 stage1/analysis/validate_cam_dumps.py` —
  ALL STRUCTURAL CHECKS PASSED. Row counts: cam_phys_tend 39936 =
  48 × 26 × 32; cam_phys_tend_q 511488 = 48 × Σ(calls×lq×32); qsat and
  tend_total 3072 each; surface_exchange 144 = 48 × 3; surface_cflx
  1632 = 48 × 34 (pcnst=34); per-step site order OK (26 site-calls,
  22 unique, clubb/microp interleaved ×3); macmic_it=1..3, scale=1/3 OK;
  q-block nlq counts consistent with the ptend `nlq` column.
  vdiff dumps (added 2026-11-03): cam_vdiff_prof 1536 = 48 × 32,
  cam_vdiff_q 52224 = 48 × 32 × 34.
- [x] **qsat law:** es=estblf(T) (Goff-Gratch mixed-phase table, linear
  interpolation), qs=epsilo·es/(p−omeps·es), es=min(es,p), gam=(hltalt/cpair)·dqsdt
  with pcf transition polynomial — 3072 rows checked, max|es err| =
  1.819e-12 Pa, max|qs err| = 1.110e-16, max|gam err| = 4.263e-14, bad=0
  (round-off level).
- [x] **Surface:** bc_out cam_out == ac_out cam_out on all 48 steps
  (nothing in tphysac rewrites cam_out between cam_export and tphysac
  exit); bc_in rows present for all steps; cflx_Q cross-checked against
  the cam_surface_cflx Q row.
- [x] **Precip:** 15 of 48 steps have nonzero precip at bc_out
  (precc max 1.01e-8 m/s at step 48, precl max 6.5e-14 m/s at step 47);
  the IOP observed-precip forcing (from IOP hour 23.5) drives the
  precipitating regime as expected.

## Stage-2 ranking (proposed pinning order)

Ranked by (a) standalone-ness (pure state→tendency laws, no pbuf),
(b) importance for the SCAM ARM97 test (precip/clouds), (c) tractability
of the Fortran chain:

1. **wv_saturation qsat law** (wv_saturation.F90 qsat_2D + wv_sat_methods.F90)
   — pure table+formula law, already fully pinned by `cam_qsat.csv`; first
   target for a hand-authored `.esm` component with tests from the dumps.
2. **energy_fixer** (energy_fixer.F90) — small, closed-form, applies to
   whole column; `cam_phys_tend.csv` site 1.
3. **rayleigh friction** (rayleigh_friction.F90) — tiny analytic law.
4. **convect_deep (ZM)** — the dominant precip source in this case; large
   chain (convect_deep.F90) but high value.
5. **microp (MG/PUMAS)** — the largest and most coupled; needs
   cam_phys_tend_q (lq constituents) plus pbuf context; defer.
6. **CLUBB_SGS (macrop/clubb/shallow)** — state-dependent closure, needs
   pbuf; defer behind ZM.

Rationale detail in the ranking section of cam_instrument.json.

## Bugs / quirks logged

- physpkg.F90: duplicate no-op `physics_update` for shallow convection
  (dead code, see above) — quirk, not a bug.
- `cam_abortutils` (src/utils/cam_abortutils.F90) does not export
  `check_allocate` while CAM-SIMA's copy does; the cam-sima sparse checkout
  (2 files) is what keeps the build from compiling orbital_data.F90
  (src/utils/cam-sima/...). Not a CESM bug — build layout note.
- Environment: `stage1/inputdata/share/domains/` was found missing
  (2026-11-02); domain file reconstructed from the ne3np4 mesh (see above).
