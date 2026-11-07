# CAM6 ZM deep convection — proposed stage-2 `.esm` decomposition

Stage-2 plan for the Zhang–McFarlane deep-convection chain
(`convect_deep.ZM`, tphysbc site 3), to be implemented as `.esm`
components pinned by the stage-1 ZM dumps `cam_zm.csv` /
`cam_zm_scal.csv` (15 active steps of the 48-step SCAM run: 18–23,
40–48; rows only where cape > capelmt).  Sources: the instrumented
cam6_4_212 fork (`stage1/EqESM-cesm-cam`, branches
`instrument/cam-dumps` + `instrument/cam-zm-dump`):
`components/cam/src/atmos_phys/schemes/zhang_mcfarlane/zm_convr.F90`
(the chain driver + buoyan_dilute + cldprp + closure + q1q2_pjr),
`components/cam/src/physics/cam/zm_conv_evap.F90`,
`zm_conv_momtran.F90`, `zm_conv_convtran.F90`, and the wrapper
`zm_conv_intr.F90` (`zm_conv_tend`).  Boundary/column documentation:
`stage1/runs/CAM_LANE.md` (ZM dump section), `stage1/runs/
cam_instrument.json`; dump validation: `stage1/analysis/
validate_cam_zm_dump.py` (ALL ZM DUMP CHECKS PASSED 2026-11-05).

## The chain and its data flow

```
state (t,q,pmid,pint,pdel,zm,zi,phis,pblh,tpert,landfrac,ztodt)
  |  zm_convr_run (zm_convr.F90): trigger, gather, working arrays
  |    buoyan_dilute  -> tp,qstp,tl,cape,lcl,lel,lon0,maxi,pblt
  |    gather         -> qg,tg,pg,zg,sg,zfg,tpg,qstpg,shat,qhat  (i=ideep)
  |    cldprp         -> mu,eu,du,md,ed,mc,qu,su,qd,sd,qs,hmn,hsat,
  |                      ql,qcde,cmeg,cug,rprdg,evpg,jt,j0,jd,jlcl
  |    (rescales eu/du/ed/cug/cmeg/rprdg/evpg by dz/dp: 1/m -> 1/mb)
  |    closure        -> mb (closure mass flux), dsubcld
  |    (limits mb by mumax and no_deep_pbl; scales all fluxes by mb)
  |    q1q2_pjr       -> dqdt, dsdt, dlg; prec, rliq, mcon
  |  physics_update (state + ptend_s/q)          -> state1 (evap inputs)
  |    zm_conv_evap_run -> prec_evap, snow_evap, ntprprd, ntsnprd,
  |                        flxprec/flxsnow (interfaces), tend_s/snwprd/
  |                        snwevmlt, tend_q (evapcdp), cldfrc, fsnow_conv
  |    zm_conv_momtran_run -> ptend_u/v, seten, pguallu/v, pgdallu/v,
  |                           icwuu/v, icwdu/v
  |    zm_conv_convtran_run (convtran1) -> CLDLIQ/CLDICE transport
  v
state' + precip fluxes
```

`cam_zm.csv` carries one row per (step, g, lev) with the convr, evap,
momtran and convtran1 boundary values; `cam_zm_scal.csv` one row per
(step, g) with the scalars (indices, mb, prec_evap/snow_evap, jctop/
jcbot = jt/maxg, and the zmconv_* + physconst constants).  The pbuf
fields ZM_JT/ZM_MAXG and the CLDTOP/CLDBOT of convect_deep equal the
dumped jt/maxg (zm_conv_intr.F90:581-586), so the chain is pinned
without pbuf reads except for `cld` (pbuf CLD, dumped as `cldfrc`).

## Proposed components (each its own `.esm`, tests from the dumps)

| component | Fortran source | pinned by | internal identities (validator checks) |
|---|---|---|---|
| `cam_zm_parcel` | buoyan_dilute (zm_convr.F90:890-1490) | cape, tl, tp, qstp, lcl, lel, lon0, maxi, pblt (scal + prof) | — |
| `cam_zm_cldprp` | cldprp (zm_convr.F90:1690-2320) | mu/eu/du/md/ed/mc, qu/su/qd/sd, qs/hmn/hsat, ql/qcde, jt/j0/jd/jlcl | mu(k) − mu(k+1) = (eu(k) − du(k))·dp(k), dp = 0.01·pdel (check 4); mc = mu + md (check 7) |
| `cam_zm_closure` | closure (zm_convr.F90:2340-2450) | mb, dsubcld, cape | — (mb enters every flux ×mb scaling) |
| `cam_zm_q1q2` | q1q2_pjr (zm_convr.F90:2530-2630) | dqdt, dsdt, dlf, cme, rprd | — |
| `cam_zm_evap` | zm_conv_evap_run (zm_conv_evap.F90) | prec_evap, snow_evap, ntprprd, ntsnprd, flxprec/flxsnow lo/hi, tend_s/snwprd/snwevmlt, tend_q | prec_evap/snow_evap = clamped column chain of ntprprd/ntsnprd·pdel/(grav·1000) (check 5); flxprec_hi = max(flxprec_lo + ntprprd·pdel/grav, 0) (check 6) |
| `cam_zm_momtran` | zm_conv_momtran_run (zm_conv_momtran.F90) | u1/v1, ptend_u/v, seten, pguallu/v, pgdallu/v, icwuu/v, icwdu/v | — |
| `cam_zm_convtran` | zm_conv_convtran_run, convtran1 (zm_conv_convtran.F90) | q_cliq/q_cice, ptend_cliq/ptend_cice | — (other convtran1 constituents obey the same law with the same fluxes; their combined tendencies are at site convect_deep of cam_phys_tend_q.csv) |
| `cam_zm_driver` | zm_convr_run skeleton + zm_conv_intr.zm_conv_tend | the full chain: state → state' + prec/rliq | the working-array laws p = 0.01·pmid, pf = 0.01·pint, z = zm + phis/grav, s = t + (grav/cpres)·z (check 8); trigger/indices jctop = jt, jcbot = maxg (check 9) |

## Unit and scaling conventions that must be reproduced exactly

- Working pressures are **mb**: p = 0.01·pmid, pf = 0.01·pint, dp =
  0.01·pdel; geometric heights z/zf/zfg are absolute (surface height
  included, z = zm + phis/grav, zs = geos·rgrav with rgrav = 1/grav);
  s = t + (grav/cpres)·z.
- cldprp computes mu normalized by the closure mass flux (mu(jb) = 1)
  and eu/du per **meter** (continuity: (eu − du)·dz = mu(k) − mu(k+1),
  dz = zf_lo − zf_hi); zm_convr.F90:651-661 converts eu/du/ed/cug/cmeg/
  rprdg/evpg to per-mb by ×dz/dp, closure returns mb, and 700-714
  scales mu/md/mc/du/eu/ed/cmeg/rprdg/cug/evpg by mb.  The dumped
  mu/md/mc are **mb/s** and eu/du/ed are **1/s**, and the continuity
  identity in dumped units is mu(k) − mu(k+1) = (eu(k) − du(k))·dp(k).
- The evap chain accumulates interface fluxes as flxprec(k+1) =
  max(flxprec(k) + ntprprd(k)·pdel(k)/gravit, 0) (the max is a
  roundoff guard; it fires on 5 of the 15 dumped columns by up to
  3.0e-22 m/s) and divides by 1000 only at the surface.
- mcon (kg/m²/s) = mc·100/grav is not a dumped column (header prose
  only); rliq = Σ dlf·pdel/grav/1000 over the column.

## Test plan (mirroring the validator identities)

All tests `time_span {0,0}`, rel 1e-12 (as in the other CAM stage-2
files), citing dump file, step, g, lev:

1. **mu continuity** (pins cldprp + the two rescalings): for k = 1..31,
   mu(k) − mu(k+1) − (eu(k) − du(k))·0.01·pdel(k) = 0.  This is the
   strongest pin: it exercises the 1/m→1/mb conversion and the ×mb
   scaling together with the flux-profile computation.
2. **evap chain** (pins zm_conv_evap_run): prec_evap and snow_evap
   against the clamped column chain of the dumped ntprprd/ntsnprd, and
   per-level flxprec_hi = max(flxprec_lo + ntprprd·pdel/grav, 0).
3. **closure + trigger**: cape > capelmt and mb > 0 on the dumped rows;
   mb from the closure formula (mumax limit and no_deep_pbl branches
   are inactive in this lane: lparcel_pbl = no_deep_pbl = F).
4. **driver plumbing**: p/pf/z/s working-array laws; jctop/jcbot =
   jt/maxg; constants block equal to the namelist values.
5. Compiled-Fortran reference pins for the branches the dump window
   does not exercise (no_deep_pbl = T, mb = 0 paths), in the
   `scratch/` driver style of README-cam.md.

Data volume: 15 columns × 32 levels per identity; the dump window
spans cape 94–594 J/kg, prec_convr 4.3e-9..7.4e-8 m/s, prec_evap
0..1.0e-8 m/s (all precip evaporates on steps 18–23).

## Coverage gaps / notes

- **buoyan_dilute** is pinned only through its outputs; the parcel
  ascent loop internals (capeten levels, qsmix, parcel_pbl branches)
  are not separately observable.  lparcel_pbl = F in this lane, so the
  parcel_pbl branch is unpinned by data.
- **cldprp's** downdraft branch and the `mu(k) < 0.02` cutoff are
  exercised (jt/maxg span 15..32) but the branch structure is not
  directly observable; the continuity identity holds through all of it.
- **convtran1** is pinned only via CLDLIQ/CLDICE (q_cliq/q_cice,
  ptend_cliq/ptend_cice); the other transported constituents (Q, CLDICE
  sources, …) obey the same flux law and their combined tendencies are
  at site convect_deep of cam_phys_tend_q.csv — a cross-check for the
  assembled `cam_zm_driver`.
- The `cldfrc` input to the evap is pbuf CLD from the previous
  processes; it is dumped, so `cam_zm_evap` needs no pbuf read.
- fsnow_conv (the evap snow-fraction input) comes from
  `cloud_fraction_fice_run` (zm_conv_intr.F90:669); dumped as an input
  to the evap component.
