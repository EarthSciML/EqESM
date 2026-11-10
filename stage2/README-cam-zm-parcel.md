# CAM6 ZM deep convection — stage-2 `.esm`: parcel/CAPE + closure

Stage-2 implementation of the two lowest blocks of the Zhang–McFarlane
deep-convection chain (see `README-cam-zm.md` for the full-chain
decomposition): **buoyan_dilute** (the parcel/CAPE calculation) and
**closure** (the cloud-base mass flux).  Three files:

| file | contents |
|---|---|
| `stage2/cam_zm_templates.esm` | shared expression-template library (17 templates): the ZM qsat (zm_qsat_hPa), the moist static energies hmn/hsat, the virtual temperatures tv/tpv, the Raymond–Blyth entropy, the closure per-level laws (thetavp/thetavm, dqsdtp, dtpdt, dtmdt/dqmdt at jt and mid, dboydt above/below the LCL), each with its lel/lcl/jt/mx region guard inside the template body |
| `stage2/cam_zm_parcel.esm` | model `Parcel`: buoyan_dilute (zm_convr.F90:890-1490) + the parcel_dilute launch/entropy laws (F90:1414-1415, 1476, 1489-1513) + the cldprp hmn/hsat laws (F90:1829-1830) |
| `stage2/cam_zm_closure.esm` | model `Closure`: closure (zm_convr.F90:2318-2506) + the driver's dsubcld (F90:597-603) |

Both models import the templates by reference
(`expression_template_imports`), so every shared law has a single
definition.  The CAPE and dadt integrals lower through the ESD rule
`integral_lev_whole` (`${ESD_ROOT}/grids/column_nonuniform_1d/
rules/integral_lev_whole.esm`, binding NLEV = 32) with the consumer's
free-name layer measure `dz`.  In `Parcel`, `dz` is the true d(ln p)
layer measure log(pf_hi/pf_lo) of the CAPE integral; in `Closure`, `dz`
is the true height layer measure zf_lo − zf_hi of the dadt integral
(`dsubcld` is not an integral at all: it is the telescoped closed form
pf_hi(pver) − pf_lo(mx) of the driver's dp sum, see below).

Sources: `stage1/dumps/cam/cam_zm.csv` (480 profile rows),
`cam_zm_scal.csv` (15 scalar rows, steps 18–23 and 40–48, g = 1),
`cam_zm_parcel.csv` (480 rows: the parcel-ascent internals
tmix/qsmix/qtmix/smix/xsh2o/qtot/ds_xsh2o/ds_freeze/qtp/sp/mp/tpv/buoy)
and `cam_zm_parcel_scal.csv` (15 rows: pl/plge600/capeten_1..5/
lelten_1..5) of the instrumented cam6_4_212 fork
`components/cam/src/atmos_phys/schemes/zhang_mcfarlane/zm_convr.F90`
(line numbers cited per law below).

## Model `Parcel` (buoyan_dilute + entropy)

Variables/equations: `qst`, `hmn`, `hsat`, `tv`, `tpv`, `buoy`,
`integrand`, `dz`, `cape` (an `integral`), `ent`, `qstp_launch`.

- `qst` = zm_qsat_hPa(t, p): the qsat_hPa wrapper (F90:2634-2649), its
  qsat_water core (wv_sat_methods.F90:227-250), reusing
  `goff_gratch_water`/`svp_to_qsat` of `cam_templates.esm` by reference
  (bit-identical to the qsat pinned by `cam_wv_saturation.esm`).
- `hmn`/`hsat` = cpres·t + grav·z + rl·q / rl·qst (F90:1829-1830; the
  same law is hmn_lev of the parcel_pbl branch, F90:979, and of the hmax
  search, F90:1018).
- `tv` = t·(1+1.608·q)/(1+q) (F90:1005-1006, 1088); `tpv` =
  (tp+tpert)·(1+1.608·qstp)/(1+qtot) for msg+1 ≤ k ≤ maxi, else the
  environmental tv (F90:1415 launch, 1476 ascent; the tpv = tv / buoy = 0
  initializations of F90:1004-1007 for k ≤ msg — the mask lives in the
  `zm_parcel_virtual_temperature` template body, bound to the `msg`
  parameter); `buoy` = tpv − tv + tiedke_add for msg+1 ≤ k ≤ maxi, else 0
  (F90:1088-1089).
- `integrand` = rgas·buoy and `dz` = log(pf_hi/pf_lo) (d(ln p)); `cape`
  = max(0, integral) over lel < k ≤ maxi (F90:1124, 1146), lowered by
  `integral_lev_whole` (the measure dz is dimensionless, so the integral
  is a plain sum of the Fortran log increments).  The lel < k bound makes
  the msg+1 lower bound of buoy redundant in the integral (lel ≥ 14).
- `ent` = the Raymond–Blyth moist entropy (F90:1489-1513): (cpres +
  qtot·cpliq)·log(TK/tfreez) − rgas·log((p−e)/pref) + L·qv/TK −
  qv·rh2o·log(qv/qst), with L = rl − (cpliq−cpwv)·(TK−tfreez),
  qv = min(qtot, qst), e = qv·p/(eps1+qv), pref = 1000 hPa.
- `qstp_launch` = q(maxi): carries the parcel_dilute launch branch
  qstp = q (F90:1414) as the bare index q(maxi) — a definition, not
  pinned by any test (its Fortran value is a dumped parameter input).

The logarithms are taken of **dimensionless ratios** (TK/tfreez,
(p−e)/pref, qv/qst, t/1K, eb/1hPa) as the dimensional checker requires;
the Fortran takes them of dimensioned values.  Working pressures are
**hPa** (the mb of the Fortran); hPa↔Pa conversion happens only inside
zm_qsat_hPa via the declared 100 Pa/hPa constant.

### What is a parameter here (and why)

- `maxi` (the launch level) — the argmax of hmn over k ≥ pblt
  (F90:1018-1020); argmax is not evaluable by the CLI interpreter
  (probed), so maxi is a parameter pinned by the dump column.
- `tp`, `qstp` and `qtot` = new_q — outputs of the **ientropy ascent
  solve** (F90:1517-1632).  The solve is an iterative Brent root-find,
  not expressible as closed-form laws, so its outputs are dumped and
  carried as parameters (cam_zm_parcel.csv columns tp/qstp/qtot), which
  pins `tpv`, `buoy` and `cape` against the dumped values.  qtot
  defaults to the step-18 profile (the step-18 test runs on the
  defaults); the dumped qtot is 0 below msg+1 = 6, where the parcel
  loops never run.  The other solve outputs (tl, lcl) are not dumped
  and are unused by `Parcel`.
- `lel` — with num_cin = 1 on all dumped rows, the candidate selection
  (F90:1136-1140) collapses to lel = lelten(1), so lel is a parameter
  (the dumped column).  The collapse is verified against the dumped
  capeten/lelten columns of cam_zm_parcel_scal.csv: cape =
  capeten(1) and lel = lelten(1) on every active row.
- `landfrac` is not read: lparcel_pbl = F in this lane (the parcel_pbl
  branch, F90:960-1030, is unused — no hmn_lev/maxi coupling needed).
- The `pl ≥ 600` gate (F90:1079-1080) is omitted: pl is dumped
  (768–830 mb on the active rows, cam_zm_parcel_scal.csv) and plge600 =
  T on all of them, so the gate never fires.

## Model `Closure` (closure + dsubcld)

Variables/equations: the dp conversion, `dsubcld` (the closed form
pf_hi(pver) − pf_lo(mx)), the launch-level extractions t_mx/q_mx/p_mx,
`eb`, `dtbdt`, `dqbdt`, `debdt`, `dtldt`,
`dtmdt_jt`/`dtmdt_mid`/`dtmdt`, `dqmdt_jt`/`dqmdt_mid`/`dqmdt`,
`thetavp_above`/`thetavp_below`/`thetavp`, `thetavm`, `dqsdtp`,
`dtpdt`, `dboydt_above`/`dboydt_below`/`dboydt`, `dz` (the dadt
measure), `dadt` (an `integral`), `mb`.

- `dp` = 0.01·pdel (the README-cam-zm working identity; pdel dumped in
  Pa, dp in hPa).
- `dsubcld` = Σ dp over k ≥ mx (F90:597-603), written in its telescoped
  **closed form** `dsubcld = pf_hi(pver) − pf_lo(mx)`, the
  interface-pressure difference over the subcloud layer: dp(k) =
  pf_hi(k) − pf_lo(k) with pf_hi(k) = pf_lo(k+1) on the interface grid,
  so the sum collapses up to round-off and needs no integral rule or measure.
  With mx = 32 = pver it is the single layer thickness dp(32), pinned
  to the dumped column (worst rel 8.3e-15 on the 15 rows;
  14.026758 hPa on step 18, the pf_hi(32) − pf_lo(32) difference of the
  defaults).
- Launch-level block (F90:2401-2410): eb = p(mx)·q(mx)/(eps1+q(mx)),
  dtbdt/dqbdt from the mu/md fluxes at mx, debdt =
  eps1·p(mx)/(eps1+q(mx))²·dqbdt, dtldt =
  −2840·(3.5/t(mx)·dtbdt − debdt/eb)/(3.5·log(t(mx)) − log(eb) −
  4.805)² (the Fortran's dimensioned logs are transcribed as the
  dimensionless ratios t(mx)/1K and eb/1hPa).
- Cumulus heating/drying (F90:2423-2446): the jt interface terms
  dtmdt_jt/dqmdt_jt and the interior terms dtmdt_mid/dqmdt_mid with
  beta = 0 (F90:2433), which collapses the condensate weights
  (beta·ql(k) + (1−beta)·ql(k+1)) to ql(k+1); dtmdt = dtmdt_jt +
  dtmdt_mid and dqmdt = dqmdt_jt + dqmdt_mid (the Fortran loops write
  disjoint level sets, k = jt and jt < k < mx).
- Virtual potential temperatures (F90:2457, 2462, 2474, 2480): thetavp
  = tp·(1000/p)^(rd/cp)·(1 + 1.608·qstp − q(mx)) above the LCL and
  tp·(1000/p)^(rd/cp)·(1 + 0.608·q(mx)) below; thetavm =
  t·(1000/p)^(rd/cp)·(1 + 0.608·q) in both branches (the same law, so
  one unknown, zeroed outside lel ≤ k < mx).
- dqsdtp = qstp·(1 + qstp/eps1)·eps1·rl/(rd·tp²) (F90:2458); dtpdt =
  tp/(1 + rl/cp·(dqsdtp − qstp/tp))·(dtbdt/t(mx) + rl/cp·(dqbdt/tl −
  q(mx)/tl²·dtldt)) (F90:2460-2462).
- dboydt above the LCL (F90:2464-2468) and below (F90:2482-2485);
  `dz` = zf_lo − zf_hi (the Fortran layer weight of the summand,
  F90:2498), and `dadt` = `integral` of dboydt over that height measure
  on the active levels lel ≤ k < mx (F90:2490-2499; dboydt is zero
  outside that range, so the whole-column `integral_lev_whole` sum
  equals the Fortran sum).  The measure is real: the integrand carries
  only the buoyancy-change rate dboydt [m/(s²·hPa)] and the layer
  thickness lives in dz [m], so dadt comes out in m²/(s²·hPa).
- `mb` = max((capelmt − cape)/(tau·dadt), 0), zero when dadt = 0
  (F90:2500-2503: dltaa = −(cape − capelmt), mb = max(dltaa/tau/dadt,
  0)).

The closure consumes the **NORMALIZED** fluxes: mu/md/mc are
dimensionless fractions of mb and du is per hPa (1/hPa) — the form the
closure receives after the zm_convr_run mb scaling (F90:648-714); the
dumped mu/md/mc/du are the mb-scaled values, and the tests feed
dump/mb (the mb pin; du is read at jt < k < mx = 16..31, where du(16)
is nonzero on every row and du(17..31) = 0).  The lel/lcl/jt/mx
region masks of the branch laws live inside the `zm_*` template bodies
(zeroing the inactive levels — the `surface_layer_source` pattern), not
in the consumer equations.  Omitted with justification: `zg` and `qs`
(closure arguments never read by its body — grep-verified); the msg
loop lower bound (loops start at msg+1 = 6, F90:2401 etc., inactive
because every branch lower bound is lel ≥ 14 > 6); the mumax flux
limit (F90:675-690, inactive on all 15 rows — verified, mumax-limited
mb equals the dumped mb); lel/lcl/jt/mx enter as dumped parameters
(they are computed by other blocks of the chain).

## Tests

All tests `time_span {0,0}`, rel 1e-12, citing dump file, step, g and
lev (280 assertions):

- **Parcel — 240 assertions** (15 steps × 16): hmn and hsat at lev
  {1, 15, 25, 32} against the dumped columns (cam_zm.csv rows
  (step, 1, k)); qst at lev {25, 32} against the dumped `qs` column of
  cam_zm.csv (same rows; also cross-pinned by hsat through hsat =
  cpres·t + grav·z + rl·qst); tpv at lev {25, 32} and buoy at lev
  {15, 18, 25} against the dumped cam_zm_parcel.csv columns; cape
  against the dumped cape of cam_zm_scal.csv (== capeten_1 of
  cam_zm_parcel_scal.csv under the num_cin = 1 collapse, verified).
  Every test overrides the 9 profile arrays T/q/z/p/tp/qstp/qtot/
  pf_hi/pf_lo and the scalars maxi/lel (the step-18 values equal the
  model defaults, so the default run *is* the step-18 state).
- **Parcel — 10 ent/tv assertions** (5 states × 2): ent at (18,32),
  (18,25), (18,15), (40,26), (48,32) and tv at the same states, against
  the entropy (F90:1489-1513) and virtual-temperature (F90:1088)
  replicas (DERIVED — not dump columns).
- **Closure — 30 assertions** (15 steps × 2): dsubcld (worst rel
  8.3e-15 — the pf_lo/pf_hi es26.17 dump round-trip of the interface
  difference) and mb (worst rel 2.5e-16) against the dumped columns,
  from the NORMALIZED fluxes dump/mb.  The step-18 test runs on the
  model defaults (every step-18 profile equals a default, so the default
  run *is* the step-18 state); the other steps override all 26
  profile/scalar arrays.

Data volume: 15 steps × 32 levels per profile array; the dump window
spans cape
94–594 J/kg, lel 14/15, lcl 25/26/27, jt 15, mx 32.

## Numeric tuple provenance

`stage1/analysis/zm_parcel_closure_tuples.py` extracts the numeric
arrays of the tests blocks from the dumps (per the stage-2 workflow:
scripts may extract numeric test tuples into an existing hand-authored
tests block; all equations and structure are hand-authored in the .esm
files).  `scratch/splice_zm_tests.py` regenerates only the `tests`
blocks of the two .esm files from the tuples, in the exact
json.dumps(indent=1) format the region was originally written in, so
`git diff` shows only test-block changes.  The tuples script also
*verifies* every law replica against the dump before writing the
tuples: worst rel |mb| 2.5e-16, |dsubcld| 0.0, |dsubcld closed form|
8.3e-15, |hmn|, |hsat|, |tpv|, |buoy|, |cape| 0.0; the pl reconstruction (from the dumped pf_hi/pf_lo
interface pressures) is exact to worst |pl| 5.7e-3 mb and the ientropy
launch root reproduces the dumped tmix to worst |tmix − root| 1.5e-4 K
(run with `ZM_DUMP_DIR`/`ZM_TUPLE_DIR` env overrides; default output
`scratch/zm_tuples/tuples.json`, gitignored).

## Coverage gaps / notes

- **The parcel ascent solve is not expressible**: the ientropy Brent
  solve (F90:1517-1632) fixes tp, qstp, tl, lcl, pl and the internal
  tmix/qsmix/xsh2o.  `tp`, `qstp`, `qtot` are dumped and carried as
  parameters (pinning tpv/buoy/cape); `tl` and `lcl` are not dumped and
  are unused by `Parcel`.  The solve's internals tmix/qsmix/xsh2o are
  dumped and verified against the stage-1 replicas (launch root
  |tmix − root| ≤ 1.5e-4 K, the ascent recurrences bit-exact), but the
  iterative solve itself remains unexpressed.  The per-level
  closed-form laws are transcribed and pinned separately.
- **argmax (maxi) and the Brent root-find are not evaluable** by the
  CLI interpreter (probed); `maxi` is a parameter (dump column; the
  qstp = q(mx) launch law is carried via qstp_launch, a definition).
- pl ≥ 600 gate (F90:1079-1080) omitted (pl dumped 768–830 mb, plge600
  = T on all rows); landfrac unused (lparcel_pbl = F); mumax limit
  (F90:675-690) and the msg loop bound omitted (inactive on all 15
  rows); zg/qs omitted (unread); beta = 0 (F90:2433) collapses the ql
  weights to ql(k+1).
- The Fortran's `dtldt` (F90:2407-2410) and entropy (F90:1489-1513)
  take logs of dimensioned values; the .esm laws use the numerically
  identical dimensionless ratios (t(mx)/1K, eb/1hPa, TK/tfreez,
  (p−e)/pref, qv/qst).
- Transcription note: binary `-` in the esm operator set is
  subtraction, so negations are written `0 − x`, never `x − 0` (the
  latter evaluates to +x; this bit the dtldt and entropy rgas terms
  during development and is regression-covered by the tests above).
- The closure tests feed level values 1–13 zeroed (the closure reads
  only lel..mx = 14/15..32 on this lane; thetavm is zeroed outside
  [lel, mx), everything else below lel is unused).  The step-18 test
  omits the overrides entirely and exercises the zeroed defaults, which
  is why the region guards live inside the template bodies rather than
  mask-multiplied in the consumer: a mask multiply would evaluate the
  (1000/p)^(rd/cp) law at p = 0 (levels 1–13) and turn mb into NaN,
  while the ifelse guard short-circuits the untaken branch.
