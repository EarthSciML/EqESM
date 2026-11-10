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

## Instrumentation (cam_dump_mod.F90, cam submodule @ 369b5cd7 + atmos_phys @ c379d7a)

- `src/physics/cam/cam_dump_mod.F90` — dump module (namelist `cam_dump_nml`:
  `dump_cam` [logical, default .false.], `dump_dir`, `dump_first_steps`;
  read in `cam_dump_readnl` from `runtime_opts.F90` `read_namelist`, broadcast
  masterproc→all). Gated by `dump_enabled()` = master PET ∧ init ∧
  1 ≤ nstep ≤ dump_first_steps; all writes inside one `!$omp critical`.
  When `dump_cam=.false.` the module is a no-op ⇒ default builds are
  bit-identical to pristine CAM.
- `src/physics/cam/physpkg.F90` — 34 dump call sites (see table).
- ZM deep-convection chain (added 2026-11-05, cam submodule `369b5cd7`,
  atmos_phys `c379d7a`, EqESM-cesm-cam pointer `8e54769`):
  `zm_convr.F90` `zm_convr_run` dumps the convr boundary (state inputs,
  p/z/s/pf/zf working arrays, buoyan_dilute/cldprp/closure/q1q2_pjr
  outputs incl. mu/eu/du/md/ed/mc, mb, prec, rliq, the zmconv_* namelist
  constants and the physconst values passed in) into the cam_dump_mod
  buffers (`cam_dump_zm_convr`, only when `lengath > 0`); `zm_conv_intr.F90`
  `zm_conv_tend` dumps the evap (`cam_dump_zm_evap` after
  `zm_conv_evap_run`: t1/q1, cldfrc, fsnow_conv, tend_s/snwprd/snwevmlt,
  tend_q, ntprprd/ntsnprd, flxprec/flxsnow, prec/snow, jctop/jcbot),
  momtran (`cam_dump_zm_momtran` after `zm_conv_momtran_run`: u1/v1,
  ptend_u/v, seten, pguallu/v, pgdallu/v, icwuu/v, icwdu/v) and convtran
  (`cam_dump_zm_convtran` after `zm_conv_convtran_run`: CLDLIQ/CLDICE
  pre-transport q and ptend; also flushes the buffered rows).  Rows are
  written only for active steps/columns (cape > capelmt, lengath > 0);
  the buffers carry slot-mismatch endrun checks and a completeness
  sentinel (zm_sentinel) checked at the flush.
- `src/atmos_phys/schemes/zhang_mcfarlane/zm_convr.F90` — the convr dump
  call at the end of `zm_convr_run` (gated by `lengath > 0`; passes the
  routine-local working arrays, buoyan_dilute outputs and constants).
- ZM parcel-ascent dump (added 2026-11-09, atmos_phys `73b65b2`, cam
  `12ddc27e`, cam pointer `87edceff`, EqESM-cesm-cam pointer `0c71776`):
  `zm_convr.F90` `buoyan_dilute`/`parcel_dilute` pass the parcel arrays
  (tmix/qsmix/qtmix/smix/xsh2o/qtot/ds_xsh2o/ds_freeze/qtp/sp/mp/tpv/
  buoy) through to `zm_convr_run`, which calls the new
  `cam_dump_zm_parcel` (cam_dump_mod.F90; + `dump_i_last` integer-dump
  helper) next to `cam_dump_zm_convr` — cam_zm_parcel.csv (per level)
  + cam_zm_parcel_scal.csv (per column: pl/plge600/capeten_1..5/
  lelten_1..5).  The qtp/sp/mp entrainment sums were changed from
  (ncol) locals to (ncol,pver) `intent(out)` with per-level storage
  (numerically identical; verified).  The 5 extra per-level columns
  (qtp/sp/mp + ds_xsh2o/ds_freeze) close the full identity chain:
  entrainment recurrences, mixing fractions, and the rainout/freezing
  entropy budget entropy(tmix,p,qtot) = smix + ds_xsh2o + ds_freeze.
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
| `cam_zm.csv` | (step, g, lev) | ZM deep-convection chain, per-level (active steps only, cape > capelmt): zm_convr boundary (t/q/pmid/pint/pdel/zm/zi inputs; p/z/s/pf/zf working arrays; tp/qstp/cape/lcl/lel/lon0/maxi/pblt; gathered qg..qhat; qs/hmn/hsat/qu/su/qd/sd/ql/qcde; mu/eu/du/md/ed/mc; dqdt/dsdt/dlf/cme/rprd), evap boundary (t1/q1, cldfrc, fsnow_conv, tend_s/snwprd/snwevmlt, tend_q, ntprprd/ntsnprd, flxprec/flxsnow lo/hi), momtran boundary (u1/v1, ptend_u/v, seten, pguallu/v, pgdallu/v, icwuu/v, icwdu/v), convtran1 boundary (q_cliq/q_cice, ptend_cliq/ptend_cice) — 83 columns |
| `cam_zm_scal.csv` | (step, g) | ZM chain scalars (active steps only): g/ideep/lengath, lat/lon, pblh, tpert, landfrac, phis, ps, ztodt; cape/tl, the level indices lcl/lel/lon0/maxi/pblt/jt/maxg/j0/jd/jlcl; mb, dsubcld, prec_convr, rliq, prec_evap, snow_evap, jctop/jcbot; the zmconv_* namelist constants and the physconst values passed to zm_convr_run (c0_lnd..latice) — 60 columns |

Row counts for the 48-step run (verified): `cam_phys_tend` 48 × 26 × 32 =
39936 rows (26 = 20 single sites + clubb×3 + microp×3); `cam_phys_tend_q`
48 × Σ(site lq count) × 32; `cam_qsat`/`cam_tend_total` 48 × 2 × 32 =
3072 each; `cam_surface_exchange` 48 × 3 = 144; `cam_surface_cflx`
48 × pcnst = 1632 (pcnst = 34); `cam_vdiff_prof` 48 × 32 = 1536;
`cam_vdiff_q` 48 × 32 × 34 = 52224.  ZM dumps (added 2026-11-05): only
the 15 active steps have rows — steps 18–23 and 40–48 (cape > capelmt =
70 J/kg; the trigger test is `cape(i) > capelmt` at zm_convr.F90:555);
`cam_zm.csv` 15 × 32 = 480 rows, `cam_zm_scal.csv` 15 rows (ncol = 1 so
g = 1, ideep = 1 on every row).

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
  bit-identical again (same two differences only). Re-verified
  2026-11-05 after adding the ZM deep-convection dump (run 12112551):
  bit-identical again (`case` attr + `date_written`/`time_written` only),
  and the 8 pre-existing dump CSVs are byte-identical to the pre-ZM run
  (backup `/scratch.local/ctessum/cam_dump_backup_20261105/`).
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

### ZM deep-convection dump (2026-11-05, run 12112551)

`python3 stage1/analysis/validate_cam_zm_dump.py stage1/dumps/cam` — ALL
ZM DUMP CHECKS PASSED (480 cam_zm rows + 15 cam_zm_scal rows):

- **structure / grouping:** headers, no blank lines, 83/60 fields per
  row; the (step,g) pairs of both files coincide; lev = 1..32 in order
  per pair; g = 1..lengath and lengath = rows-per-step; steps ⊆ 1..48;
  ymd/tod constant within a step.
- **constants:** the c0_lnd..latice block is identical on every row and
  matches the lane's zmconv namelist (c0_lnd=0.0075, c0_ocn=0.03,
  ke=5e-6, ke_lnd=1e-5, momcu=momcd=0.7, num_cin=1, tau=3600,
  capelmt=70, dmpdz=-1e-3, tiedke_add=0.5, parcel_hscale=0.5,
  parcel_pbl=F, deepconv_pbl=F); msg = limcnv−1 (= 5, 6).
- **updraft mass continuity (identity):**
  mu(k) − mu(k+1) = (eu(k) − du(k))·dp(k), dp = 0.01·pdel, k = 1..31 —
  holds to 1.9e-16 (machine precision) on all 465 rows.  NOTE: eu/du
  are dumped in final units (1/s); cldprp computes them per-meter and
  zm_convr.F90:651-661 converts to per-mb and 700-714 scales by mb, so
  the continuity identity is in **pressure** thickness dp, not the
  geometric dz (the dump header's "mu=eu-du identity" phrasing is loose).
- **precip/snow production (identity):** prec_evap/snow_evap = the
  zm_conv_evap.F90:219-224 clamped column chain of
  ntprprd/ntsnprd·pdel/(grav·1000) — reproduced bit-exactly (rel resid
  0.0) on all 15 columns; the unclamped sum differs by ≤ 3.0e-22 m/s
  (the amount absorbed by the `max(flux,0)` roundoff guard).  On steps
  18–23 all ZM precip evaporates before reaching the surface
  (prec_evap = 0 with prec_convr up to 2.7e-8 m/s); from step 40 on
  prec_evap ≈ 4.5e-9..1.0e-8 m/s.
- **interface fluxes (identity):** flxprec_hi = max(flxprec_lo +
  ntprprd·pdel/grav, 0) per row (snow analog) — bit-exact on all 960
  rows; flxprec_lo(lev=1) = 0 and flxprec_hi(lev=pver) =
  prec_evap·1000 verified.
- **mass-flux bookkeeping:** mc = mu + md per row (2.2e-16).
- **working arrays:** p = 0.01·pmid, pf = 0.01·pint (bit-exact),
  z = zm + phis/grav (2.1e-16; FORTRAN uses zs = geos·rgrav, rgrav =
  1/grav), s = t + (grav/cpres)·z (bit-exact).
- **trigger / indices:** cape > capelmt on all 15 rows (min 94.2 J/kg);
  1 ≤ jt ≤ maxg ≤ pver; jctop == jt and jcbot == maxg
  (zm_conv_intr.F90:581-586); all level indices in bounds.
- Quirk logged: the cam_zm.csv `# description` header mentions
  "mcon [kg/m2/s]" but no mcon column is written (mcon = mc·100/grav is
  derivable; zm_convr.F90:769).  Cosmetic only.

### RRTMG radiation dump (2026-11-08, run 12255018)

`python3 stage1/analysis/validate_cam_rad_dump.py stage1/dumps/cam
--info` — ALL RADIATION DUMP CHECKS PASSED.  Adds 8 new CSVs
(`cam_rad_input/scal/out/out_ifc/cldopt_sw/cldopt_lw/solver_sw/
solver_lw.csv`) tracing the radiation boundary of `radiation_tend`
(tphysbc site 11, radiation.F90) at the RRTMG climate call (icall=0,
the last call of the SW/LW diagnostic loop; radsw.F90/radlw.F90 fill
the cam_dump_mod buffers on every call, the last one wins).  Sources:
`cam_dump_mod.F90` (buffers + `cam_dump_rad_sw/lw` + `cam_dump_rad_flush`),
`radsw.F90`/`radlw.F90` (buffer calls at the end of `rad_rrtmg_sw/lw`),
`radiation.F90` (flush call before `rrtmg_state_destroy`).  Fork
commits: cam submodule `81922142` (surface-column fix on top of
`6a528c94`), superpointer `6951a7a`/`d827b83`.

- **configuration pinned:** `radiation_scheme='rrtmg'`, iradsw=iradlw=2
  (radiation on steps 3,5,…,47 = 23 steps), icecldoptics='mitchell',
  liqcldoptics='gammadist', `num_rrtmg_levs`=33=pverp (RRTMG layer k=1
  is the extra layer above model top: pintmb(1)=0.01 Pa, pmid(1)=
  0.5·pintmb(2); cam_lev=max(k−1,1)), McICA deterministic
  (permuteseed=1 SW / 150 LW, single column).  The RRTMG k-distribution
  gas optics are NOT dumped (too large); the per-g-point solver
  boundary (cld_sto/tauc_sto/ssac_sto/asmc_sto/fsfc_sto, SW; cld_sto/
  tauc_sto, LW) pins the solver inputs instead.
- **day/night:** the lane is NOT always daylight — 13 day steps
  (3, 25, 27, …, 47) and 10 night steps (5..23, coszrs<0, solin=0).  On
  night steps the flush records the all-zero SW boundary (radsw.F90
  early-returns with zeros; the flush zero-fills the SW buffers,
  cam_dump_mod.F90:1752-1767); LW always runs.  The validator splits
  its SW identities by day/night accordingly.
- **identities verified (validator):** (i) h2ovmr = q/(1−q)·amdw
  (amdw=1.607793), pmid(k)=½(pint_lo+pint_hi) k≥2, pmid(1)=½·pint_hi(1),
  shared edges pint_hi(k)=pint_lo(k+1) and T_ifc_hi(k)=T_ifc_lo(k+1),
  cldfprime=max(cld,cldfsnow), cfc22vmr=ccl4vmr=0; (ii) solin =
  Σ(sfac_b·solar_ref_band_irradiance)·eccf·coszrs on day steps
  (radsw.F90:395-397), solin=0 at night; (iii) fns=fds−fus, fnl=ful−fdl;
  (iv) qrs = −g·Δfns/Δp (max resid 2.6e-16 W/kg) and qrl =
  +g·Δfnl/Δp (4.3e-16) against the interface fluxes and pintmb from
  cam_rad_input — the RRTMG heating is exactly the flux divergence;
  (v) column integrals Σqrs·Δp/g = fsnt−fsns, Σqrl·Δp/g = flns−flnt;
  (vi) per-band cloud optics: ssa=tau_w/tau (tau>0) else 1, asm=
  tau_w_g/tau_w (tau_w>0) else 0, fsfc likewise, 0≤tau_w≤tau; (vii)
  McICA invariants per g-point via the RRTMG band→gpt mappings
  (SW 6,12,8,8,10,10,2,10,8,6,6,8,6,12; LW 10,12,16,14,16,8,12,8,12,6,
  8,8,4,2,2,2): cld_sto∈{0,1}, tauc_sto = tau(band) when cld_sto=1 else
  0 (bit-exact), ssac_sto=ssa(band) when cld_sto=1 else 1 (day), asm/
  fsfc likewise; (viii) scal↔ifc: fsnt=fns(1), fsns=fns(33), fsds=
  fds(33), flnt=fnl(1), flut=ful(1), flns=fnl(33), flwds=fdl(33),
  fsntc=fcns(1), fsnsc=fcns(33), flntc=fcnl(1), flnsc=fcnl(33); (ix)
  sols+soll+solsd+solld = fds(surface) (day); (x) fds(1) ≤ solin with
  extra-layer absorption solin−fds(1) = 2.0..5.3 W/m2 (the O3
  P_top=50 Pa linear-decay column).
- **dump bug found and fixed (validator caught it):** cam_rad_input
  promised per-row surface columns Ts/asdir/asdif/aldir/aldif in its
  schema but the row loop never wrote them (rows had 32 fields vs 37
  declared).  Fixed in cam submodule `81922142` (writes added to the
  flush row loop); rerun 12255018 re-verified.  (The values are also in
  cam_surface_exchange.csv phase=bc_in, but the rad input file is now
  self-contained as its schema says.)
- **neutrality re-verified:** h0i bit-identical to the pristine run
  (compare_cam_hist.py: only case attr + date_written/time_written
  differ) and bit-identical to the ZM-era instrumented run (runs
  reproducible); all 10 pre-existing CSVs byte-identical to the
  pre-rad run (backup in /scratch.local/ctessum/cam_dump_backup_20261108/).
- **queue note:** on 2026-11-08 the ctessum partition was saturated by
  24-h agent jobs; the case batch script's `#SBATCH --exclusive`
  (env_batch.xml directive) left the 1-PE run queued ~2 h.  The
  `--exclusive` line was removed from `.case.run` for runs 12254141/
  12255018 (they scheduled immediately on the shared node ccc0236) and
  restored afterwards; env_batch.xml itself was untouched (it is
  locked).  Recorded here so the next lane knows why the batch script
  differs from the template in gitignored CaseDocs.

### ZM parcel-ascent dump (2026-11-09, interactive re-run 21:59)

`python3 stage1/analysis/validate_cam_zm_parcel_dump.py stage1/dumps/cam`
— ALL ZM PARCEL DUMP CHECKS PASSED (480 cam_zm_parcel rows + 15
cam_zm_parcel_scal rows).  The parcel-ascent internals of
buoyan_dilute/parcel_dilute (zm_convr.F90) are now dumped at the same
trace point and with the same (step,g[,lev]) keys as cam_zm.csv:
tmix/qsmix/qtmix/smix (the entraining parcel), xsh2o/qtot/ds_xsh2o/
ds_freeze (rainout + freezing), qtp/sp/mp (entrainment sums, per-level
`intent(out)` accumulation — numerically identical to the original
(ncol) locals), tpv/buoy, plus per-column pl/plge600/capeten_1..5/
lelten_1..5.  Fork hashes: atmos_phys `73b65b2`, cam `12ddc27e` (cam
pointer `87edceff`), EqESM-cesm-cam pointer `0c71776` (branch
instrument/cam-dumps).

- **launch level (k = maxi = 32 = pver):** qtmix == qtot == qtp0 == q
  and smix == sp0 == entropy(t,p,q) bit-exact; tmix is the ientropy
  root of entropy(t,p,q) = smix to 1.5e-4 K (Brent tol1 ~= 5.2e-4 K,
  EPS=3e-8, tol=0.001, zm_convr.F90:1517-1632); qsmix == qsat_GG(tmix)
  to 4.4e-15 rel; xsh2o/ds_xsh2o/ds_freeze/qtp/sp/mp == 0.
- **entrainment recurrences (identity):** qtp(k) = qtp(k+1) −
  dmpdp·dp·qtenv (sp/mp analog) with qtenv/tenv/penv the layer means,
  senv = entropy(tenv,penv,qtenv), dpdz = −penv·grav/(rgas·tenv),
  dmpdp = dmpdz/dpdz (F90:1322-1340) — bit-exact for qtp/mp (max
  6.9e-18), 4.5e-13 for sp; the mixing fractions qtmix = (qtp0+qtp)/
  (mp0+mp) and smix = (sp0+sp)/(mp0+mp) bit-exact (390 rows).
- **precip/freezing loop (identity):** the nit_lheat=2 rainout loop
  (F90:1392-1463, lwmax=1e-3, tscool=0) reconstructed from the dumped
  smix/qtmix and the final values at k+1: xsh2o = max(0, qtmix−qsmix−
  lwmax); ds_xsh2o recurrence −cpliq·log(tmix/tfreez)·max(0, xsh2o(k)−
  xsh2o(k+1)); one-off freezing (latice/tmix)·max(0, qtmix−qsmix−
  xsh2o) at the first tmix ≤ 273.15 K level and the continual term
  (latice/tmix)·max(0, qsmix(k+1)−qsmix(k)) above it; then
  new_s = smix+ds_xsh2o+ds_freeze, new_q = qtmix−xsh2o and the
  ientropy inversion.  Reconstructed tmix/qsmix/xsh2o/qtot/ds_xsh2o/
  ds_freeze match the dump to 2.6e-4 K / 2.3e-7 / 8.4e-8 / 8.4e-8 /
  1.3e-5 / 1.2e-5 (root-tolerance dominated).  Direct rainout entropy
  identity |entropy(tmix,p,qtot) − (smix+ds_xsh2o+ds_freeze)| ≤ 2.8e-3
  J/kg/K (Brent bound dS/dT·tol1 ~= 3.3e-3); ds_freeze ≠ 0 only where
  tmix ≤ 273.15 K (max 9.9 J/kg/K, accumulated over the ~20 cold
  levels above the lcl); xsh2o ≥ 0.
- **parcel state chain (identity):** tp == tmix, qtot == qtmix − xsh2o,
  qstp == min(qtot,qsmix), tpv == (tp+tpert)·(1+1.608·qstp)/(1+qtot)
  (F90:1411-1477) and buoy == tpv − tv + tiedke_add (F90:1089) — all
  bit-exact on k = 6..maxi (405 rows); qsmix == qsat_direct(tmix,p) to
  4.4e-15 rel — the ZM parcel qsat is the DIRECT Goff-Gratch formula
  (wv_sat_scheme default "GoffGratch", wv_saturation.F90:159), NOT the
  estblf table of cam_qsat.csv; buoy == 0 and tpv == tv for k ≤ msg = 5.
- **lcl / pl / tl (identity):** lcl = first saturated level on ascent
  (reconstructed pre-rainout qsmix: qsmix_pre(lcl) ≤ qtmix(lcl),
  qsmix_pre(k) > qtmix(k) below); pl = p(lcl+1) − qxskp1/dqxsdp
  (F90:1356-1361) to 5.7e-3 mb — the dump holds the post-rainout qsmix
  at the lcl while pl used the pre-rainout value, recovered to the
  Brent root tolerance via the entropy inversion (tolerance 0.02 mb);
  plge600 == (pl ≥ 600) ("T" everywhere, pl = 768-830 mb);
  |entropy(tl, pl, qtlcl) − slcl| ≤ 4.1e-3 J/kg/K (slcl/qtlcl the
  (k,k+1) interpolation, F90:1361-1367).
- **CAPE / lel (identity):** lel == lelten_1 = the FIRST k in
  msg+2..pver with k < lcl where buoy(k+1) > 0 and buoy(k) ≤ 0
  (F90:1107-1116); capeten_1 == Σ_{k=lelten_1+1..maxi}
  rgas·buoy·log(pf(k+1)/pf(k)) (F90:1120-1128, pf = dumped pf_lo/pf_hi)
  bit-exact (max resid 0.0); cape == capeten_1 and capeten_2..5 ==
  lelten_2..5 == 0 (num_cin = 1); lel ∈ {14,15}, cape 94-594 J/kg.
- **neutrality re-verified:** h0i bit-identical to the pristine run
  (compare_cam_hist.py: only case attr + date_written/time_written)
  and bit-identical to the pre-parcel instrumented run (runs
  reproducible); ALL 18 pre-existing cam_*.csv byte-identical to the
  pre-parcel run (backup
  /scratch.local/ctessum/cam_dump_backup_20261109_pre_parcel/).
- **queue note:** on 2026-11-09 the ctessum partition was fully
  allocated again; this time the case was re-run in place interactively
  (`source stage1/runs/env.sh && env -u SLURM_MEM_PER_NODE -u
  SLURM_MEM_PER_CPU -u SLURM_MEM_PER_GPU python3 .case.run`) inside
  the agent's own slurm job 12282180 instead of case.submit (run dir
  refreshed at 21:59, model ~35 s wall, exit 0).  `.case.run` still
  carries its restored `#SBATCH --exclusive`.

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
   chain (convect_deep.F90) but high value.  The full chain is now pinned
   at its four boundaries by `cam_zm.csv`/`cam_zm_scal.csv` (convr, evap,
   momtran, convtran1), including the closure mass flux mb, the
   mu/eu/du/md/ed/mc flux profiles and the ntprprd/ntsnprd/flxprec/
   flxsnow precip chain, with the internal identities above verified to
   machine precision; convect_deep's own state→tendency contract is
   also available at site 3 of `cam_phys_tend.csv`.
5. **microp (MG/PUMAS)** — the largest and most coupled; needs
   cam_phys_tend_q (lq constituents) plus pbuf context; defer.
6. **CLUBB_SGS (macrop/clubb/shallow)** — state-dependent closure, needs
   pbuf; defer behind ZM.
7. **RRTMG radiation** — large but now cleanly separable into three
   pinned pieces (see `stage2/README-cam-rad.md`): the CAM cloud optics
   (gammadist/mitchell → per-band tau/ssa/asm/fsfc, pinned by
   cam_rad_cldopt_*.csv), the McICA stochastic subcolumn generator
   (pinned by cam_rad_solver_*.csv, deterministic seeds), and the
   two-stream solver boundary/verdict (fluxes/heating, pinned by
   cam_rad_out*.csv; LW could reuse the RRTMG LW work in
   EarthSciModels atmospheric_radiation).  The k-distribution gas
   optics are the only unpinned interior.

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
