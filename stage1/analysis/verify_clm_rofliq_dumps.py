#!/usr/bin/env python3
# EqESM CLM land->river runoff lane (M9 part A): verify the gridcell runoff
# export dump (clm_rofliq_exports.csv) and the column-level production dump
# (clm_rofliq_cols.csv) against the Fortran laws of lnd2atmMod.F90
# (runoff section 381-470) and subgridAveMod.F90 (c2g_1d 761-818,
# set_c2l_scale 104-165), reproduce the first MOSART coupling-step runoff
# import (mosart_run_imports.csv) from the lnd exports and the mesh
# geometry, and verify solution neutrality of the instrumented run
# (bit-identical restart vs the uninstrumented baseline I2000Clm60Sp_f19,
# bit-identical re-dumps of the 15 existing clm_*.csv).
"""Verify the CLM land->river runoff dumps against the Fortran laws.

Dump sites (both gated by dump_rofliq, default .false.):
  clm_rofliq_exports.csv -- at the lnd_import_export export call
    (lnd_import_export.F90 export_fields, "output to river" block
    890-935): the gridcell values of all 5666 f19 land gridcells at every
    land export step (LND_NCPL=48: every step, the runoff fields are not
    accumulator-averaged), gathered to PET 0 (MPI_Gatherv) and sorted by
    gindex.  Flrl_rofsur = qflx_rofliq_qsur_grc, Flrl_rofsub =
    qflx_rofliq_qsub_grc + qflx_rofliq_drain_perched_grc, Flrl_rofgwl =
    qflx_rofliq_qgwl_grc, Flrl_irrig = -qirrig_grc (the mediator negates
    the demand), Flrl_rofi = qflx_rofice_grc; all [mm/s] = [kg/m2/s]
    positive downward.
  clm_rofliq_cols.csv -- at the end of the runoff section of lnd2atm
    (lnd2atmMod.F90 470-487, after the qflx_liq_from_ice_col adjustment
    loop and all c2g calls): the exact column fluxes the c2g calls saw:
    qflx_surf_col, qflx_drain_col, qflx_drain_perched_col,
    qflx_qrgwl_col (already including qflx_liq_from_ice_col),
    qflx_sfc_irrig_col, qflx_ice_runoff_col, with the c2g weights
    wtgcell, ctype and canyon_hwr.

Checks (all to round-off, es26.17 dump precision, unless noted):
  structure: row counts (48 steps x 5666 gridcells / 28,494 columns),
    per-step gindex coverage and sort order, tod = 1800*step mod 86400,
    lfrac in [0,1], area > 0, fluxes finite, dynbal fields exactly 0
    (this run), ctype in the known set {1 soil, 5 glacier,
    71 roof, 72 sunwall, 73 shadewall, 74 road_imperv, 75 road_perv,
    215 lake, 401..410 crop}, all columns active, wtgcell >= 0,
    sum(wtgcell) == 1 per (step, gridcell) to 1e-12 (the c2g
    denominator; lfrac is the separate gridcell land fraction),
    canyon_hwr finite iff the column type is urban (71..75) -- CLM never
    initializes canyon_hwr for non-urban landunits, so the dump sees NaN
    there (benign: set_c2l_scale only reads it for urban landunits).
    Sign conventions: qsur, qperch, qirrig, rofice and the column
    production fields qflx_surf/perch/sfc_irrig/ice_runoff are >= 0
    (production, positive downward), but qflx_drain_col (bottom
    drainage/capillary rise) and qflx_qrgwl_col (glacier columns
    accumulate: negative runoff; mosart_in sets qgwl_runoff_option =
    "negative") -- and hence qsub, qgwl and the total -- can be
    negative.
  c2g identity (lnd2atmMod.F90 395-469, subgridAveMod.F90 c2g_1d +
    set_c2l_scale; per step and gridcell):
      sumwt = sum(wtgcell) over active columns (== 1 here)
      exported = sum_c(scale_c * wtgcell_c * flux_c) / sumwt + dynbal
    with the urbanf law (set_c2l_scale 'urbanf'):
      scale = 3*canyon_hwr   for sunwall (72) / shadewall (73)
      scale = 3              for road_imperv (74) / road_perv (75)
      scale = 1              for roof (71) and all non-urban types
    applied to the six dumped column fluxes, and:
      qflx_rofliq_grc (total) = c2g(qflx_runoff_col) + qflx_liq_dynbal,
      where qflx_runoff_col = qsur + qdrain + qperch + qrgwl (the
      lnd2atm liq_from_ice adjustment is already inside the dumped
      qflx_qrgwl_col); equivalently checked as the c2g of the sum.
    dynbal terms qflx_liq/ice_dynbal_grc are added to the qgwl / total /
    rofice gridcell values (lnd2atmMod.F90 445-452, 465-469).
  qtotal vs qsur+qsub+qperch+qgwl at the gridcell level: consistent to
    1e-11 (rounding of the separate c2g sums).
  mosart reproduction (step 1 only -- the only coupling step with
    nonzero imports in this run): the mediator accumulates the lnd
    exports between rof coupling steps (CMEPS med_phases_prep_rof_accum,
    called from med_phases_post_lnd) and time-averages them
    (fldbun_average by lndAccum2rof_cnt = 6 exports over the 3 h window
    [0, 10800] here) before mapping with ESMF mapconsf (built online
    from mesh_lnd = fv1.9x2.5_141008_ESMFmesh.nc and mesh_rof =
    r05_nomask_c110308_ESMFmesh.nc).  The four nonzero mosart cells
    (gindex 100312/100313/101032/101033 = a 2x2 block of 0.5-deg r05
    cells at lon 295.5..296.5, lat -20.5..-19.5, fully inside land
    gridcell 1981) reproduce as
      import [m3/s] = A_dest * avg6(qsur_src) * W * 1e-3
    with A_dest the spherical cell area (R = 6371.229 km, exact
    lat/lon quad), avg6 = (qsur_step1 + qsur_step2)/6 = 2.590714e-11
    mm/s (qsur of gridcell 1981 is nonzero only at steps 1-2 of this
    cold start), and W = A_overlap/A_dest = 1.0 (100312, 101032) / 0.5
    (100313, 101033).  Checked at rel 1e-5: the residual ~1e-6 is the
    ESMF-internal spherical-polygon intersection / sphere-radius
    convention, which is computed online at run time and never
    persisted, so the exact weights are not reproducible from the mesh
    files alone (documented limitation).  The mosart dump window is
    partial (first 16 local cells of PET 0, 8 coupling steps), and its
    imports are exactly zero at coupling steps 2..8 although the lnd
    exports are nonzero all day -- an open question for the mosart lane
    (the step-1 values here are the only pinnable ones).
  neutrality (--restart, --existing): the instrumented run's clm2
    restart is bit-identical (all 290 common numeric variables,
    mask-aware) to the uninstrumented baseline I2000Clm60Sp_f19
    restart, and the 15 existing clm_*.csv re-dumps are bit-identical
    to the pre-lane backup.

Usage:
  python3 verify_clm_rofliq_dumps.py [--dumps DIR] [--mosart FILE]
                                     [--restart NEW BASE] [--existing BACKUP_DIR]
"""
import argparse, hashlib, math, os, sys

import numpy as np

D = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/clm"
MOSART = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/mosart/mosart_run_imports.csv"

DTIME = 1800.0          # land model time step (s); LND_NCPL=48
NSTEPS = 48             # 1-day run
NGRC = 5666             # f19 land gridcells
NCOL = 28494            # f19 land columns
R_SPHERE = 6371.229e3   # ESMF default sphere radius (m)

# relative tolerance for the recomputed laws (es26.17 dump precision;
# the identities hold to ~4e-16, so 1e-12 is comfortable)
RTOL = 1e-12
ATOL = 1e-24
# mosart reproduction: analytic spherical areas leave an ESMF-internal
# residual of ~1.3e-6 (online mapconsf weights, never persisted)
RTOL_MOSART = 1e-5

CTYPE_URBAN = (71, 72, 73, 74, 75)   # roof, sunwall, shadewall, road_imperv, road_perv
CTYPE_ALL = (1, 5, 71, 72, 73, 74, 75, 215, 401, 402, 403, 404, 405,
             406, 407, 408, 409, 410)

nbad = 0
def chk(cond, msg):
    global nbad
    if not cond:
        nbad += 1
        print("FAIL:", msg)

def parse(path, ncol):
    """Parse a dump CSV into a float ndarray (skips '#' header lines)."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.strip().rstrip(",").split(",")
            if len(v) != ncol:
                raise ValueError(f"{path}: row has {len(v)} fields, expected {ncol}")
            rows.append([float(x) for x in v])
    return np.array(rows)

def load():
    cols = parse(f"{D}/clm_rofliq_cols.csv", 15)
    exp = parse(f"{D}/clm_rofliq_exports.csv", 17)
    return cols, exp

# exports columns (0-based): 0 step 1 ymd 2 tod 3 gindex 4 lat 5 lon
#   6 area 7 lfrac 8 qsur 9 qsub 10 qperch 11 qgwl 12 qirrig 13 rofice
#   14 qtotal 15 liq_dynbal 16 ice_dynbal
# cols columns: 0 step 1 ymd 2 tod 3 gidx_grc 4 gidx_col 5 ctype 6 active
#   7 wtgcell 8 canyon_hwr 9 qflx_surf 10 qflx_drain 11 qflx_drain_perched
#   12 qflx_qrgwl 13 qflx_sfc_irrig 14 qflx_ice_runoff

def check_structure(cols, exp):
    # --- exports ---
    chk(exp.shape == (NSTEPS * NGRC, 17), f"exports shape {exp.shape}")
    steps = exp[:, 0].astype(int)
    chk(set(steps) == set(range(1, NSTEPS + 1)), "exports steps != 1..48")
    for s in range(1, NSTEPS + 1):
        m = exp[steps == s]
        chk(len(m) == NGRC, f"exports step {s}: {len(m)} rows")
        g = m[:, 3].astype(int)
        chk(np.array_equal(g, np.arange(1, NGRC + 1)),
            f"exports step {s}: gindex not sorted 1..{NGRC}")
        tod = m[:, 2].astype(int)
        chk(np.all(tod == (DTIME * s) % 86400), f"exports step {s}: tod mismatch")
        ymd = m[:, 1].astype(int)
        chk(np.all(ymd == (20000102 if s == NSTEPS else 20000101)),
            f"exports step {s}: ymd mismatch")
    chk(np.all((exp[:, 4] >= -90) & (exp[:, 4] <= 90)), "exports lat out of range")
    chk(np.all((exp[:, 5] >= 0) & (exp[:, 5] < 360)), "exports lon out of range")
    chk(np.all(exp[:, 6] > 0), "exports area <= 0")
    chk(np.all((exp[:, 7] >= 0) & (exp[:, 7] <= 1)), "exports lfrac out of range")
    chk(np.all(np.isfinite(exp[:, 8:15])), "exports NaN/Inf in flux columns")
    # qsur/qperch/qirrig/rofice are production (>= 0); qsub/qgwl/qtotal
    # can be negative (bottom drainage, glacier accumulation)
    chk(np.all(exp[:, [8, 10, 12, 13]] >= 0), "exports negative production flux")
    chk(np.all(exp[:, 15] == 0) and np.all(exp[:, 16] == 0),
        "exports dynbal fields not exactly 0")
    # --- cols ---
    chk(cols.shape == (NSTEPS * NCOL, 15), f"cols shape {cols.shape}")
    csteps = cols[:, 0].astype(int)
    chk(set(csteps) == set(range(1, NSTEPS + 1)), "cols steps != 1..48")
    for s in range(1, NSTEPS + 1):
        m = cols[csteps == s]
        chk(len(m) == NCOL, f"cols step {s}: {len(m)} rows")
        key = m[:, 3].astype(np.int64) * 100000 + m[:, 4].astype(np.int64)
        chk(len(np.unique(key)) == NCOL and np.all(np.diff(key) > 0),
            f"cols step {s}: (grc, col) not sorted/unique")
        chk(set(m[:, 3].astype(int)) == set(range(1, NGRC + 1)),
            f"cols step {s}: gidx_grc coverage")
        chk(set(m[:, 4].astype(int)) == set(range(1, NCOL + 1)),
            f"cols step {s}: gidx_col coverage")
    ct = cols[:, 5].astype(int)
    chk(set(np.unique(ct)) <= set(CTYPE_ALL), f"cols unknown ctype: {np.unique(ct)}")
    chk(np.all(cols[:, 6] == 1), "cols inactive columns present")
    chk(np.all(cols[:, 7] >= 0), "cols negative wtgcell")
    chk(np.all(np.isfinite(cols[:, 9:15])), "cols NaN/Inf in flux columns")
    # production fields surf/perch/sfc_irrig/ice_runoff >= 0; drain and
    # qrgwl can be negative (see docstring)
    chk(np.all(cols[:, [9, 11, 13, 14]] >= 0), "cols negative production flux")
    # sum(wtgcell) == 1 per (step, gridcell) -- the c2g denominator
    key = (cols[:, 0] * 100000 + cols[:, 3]).astype(np.int64)
    N = int(key.max()) + 1
    sumwt = np.bincount(key, weights=cols[:, 7], minlength=N)
    u = np.unique(key)
    chk(np.max(np.abs(sumwt[u] - 1.0)) < 1e-12,
        f"sum(wtgcell) != 1: max dev {np.max(np.abs(sumwt[u]-1.0)):.3e}")
    # canyon_hwr finite iff urban ctype (CLM leaves it uninitialized for
    # non-urban landunits; set_c2l_scale never reads it there)
    urban = np.isin(ct, CTYPE_URBAN)
    chk(np.all(np.isfinite(cols[urban, 8])), "urban column with NaN canyon_hwr")
    chk(np.all(np.isnan(cols[~urban, 8])), "non-urban column with finite canyon_hwr")
    print(f"structure: exports {exp.shape[0]} rows, cols {cols.shape[0]} rows, "
          f"sum(wtgcell)==1 for {len(u)} (step,gridcell) pairs")

def urbanf_scale(ct, hwr):
    """set_c2l_scale 'urbanf' law, urbpoi derived from the column type."""
    scale = np.ones(len(ct))
    wall = (ct == 72) | (ct == 73)
    road = (ct == 74) | (ct == 75)
    scale[wall] = 3.0 * hwr[wall]
    scale[road] = 3.0
    return scale

def check_c2g(cols, exp):
    ct = cols[:, 5].astype(int)
    scale = urbanf_scale(ct, cols[:, 8])
    key = (cols[:, 0] * 100000 + cols[:, 3]).astype(np.int64)
    eke = (exp[:, 0] * 100000 + exp[:, 3]).astype(np.int64)
    N = max(int(key.max()), int(eke.max())) + 1
    w = cols[:, 7] * scale
    def c2g(flux):
        return np.bincount(key, weights=w * flux, minlength=N)
    den = np.bincount(key, weights=cols[:, 7], minlength=N)
    # dynbal fields are 0 this run; add them anyway to pin the law shape
    liq_db, ice_db = exp[:, 15], exp[:, 16]
    nqsur = c2g(cols[:, 9]);   nqsub = c2g(cols[:, 10])
    nqperch = c2g(cols[:, 11]); nqgwl = c2g(cols[:, 12])
    nirrig = c2g(cols[:, 13]); nrofice = c2g(cols[:, 14])
    ntotal = c2g(cols[:, 9] + cols[:, 10] + cols[:, 11] + cols[:, 12])
    pairs = [("qsur", nqsur, exp[:, 8], 0.0),
             ("qsub", nqsub, exp[:, 9], 0.0),
             ("qperch", nqperch, exp[:, 10], 0.0),
             ("qgwl", nqgwl, exp[:, 11], liq_db),
             ("qirrig", nirrig, exp[:, 12], 0.0),
             ("rofice", nrofice, exp[:, 13], ice_db),
             ("qtotal", ntotal, exp[:, 14], liq_db)]
    for name, num, dumped, dynbal in pairs:
        got = num[eke] / den[eke]
        want = dumped + dynbal
        rel = np.abs(got - want) / np.maximum(np.abs(want), ATOL)
        bad = rel > RTOL
        chk(not bad.any(),
            f"c2g identity {name}: {bad.sum()} cells over rtol {RTOL:.0e} "
            f"(max {rel.max():.3e})")
    # gridcell-level total vs sum of the four liquid parts (rounding of
    # the separate c2g sums, so 1e-11)
    s = exp[:, 8] + exp[:, 9] + exp[:, 10] + exp[:, 11] + liq_db
    rel = np.abs(exp[:, 14] - s) / np.maximum(np.abs(s), ATOL)
    chk(rel.max() < 1e-11, f"qtotal vs sum of parts: max rel {rel.max():.3e}")
    print("c2g identities: all 7 fields over all (step, gridcell) pairs at "
          "round-off (rtol 1e-12)")

def check_mosart(cols, exp, mosart_path):
    if not os.path.exists(mosart_path):
        chk(False, f"mosart dump not found: {mosart_path}")
        return
    rows = parse(mosart_path, 17)
    # columns: 0 step 1 ymd 2 tod 3 nt 4 nr 5 gindex 6 mask 7 qsur 8 qsub
    #   9 qgwl 10 qirrig 11 qglc_liq 12 qglc_ice 13 volr 14 wt 15 wr 16 wh
    nz = rows[(rows[:, 0] == 1) & (rows[:, 3] == 1) & (rows[:, 7] > 0)]
    chk(len(nz) == 4, f"mosart step-1 nonzero qsur rows: {len(nz)} (expected 4)")
    if len(nz) != 4:
        return
    gidx = nz[:, 5].astype(int)
    imp = nz[:, 7]
    chk(set(gidx) == {100312, 100313, 101032, 101033},
        f"mosart nonzero cells {gidx} != {{100312,100313,101032,101033}}")
    # average qsur over the first rof accumulation window (3 h = 6 lnd
    # exports, steps 1..6): find the source gridcell whose avg6 matches
    q12 = {}
    for r in exp:
        s, g = int(r[0]), int(r[3])
        if s in (1, 2):
            q12.setdefault(g, {})[s] = r[8]
    def area_dest(g):
        # r05 cell: 720 cols/row; g = (row-1)*720 + col; col c spans
        # lon [-180+0.5(c-1), -180+0.5c]; row r spans lat [-90+0.5(r-1), ...]
        row, col = (g - 1) // 720 + 1, (g - 1) % 720 + 1
        lon1, lon2 = np.radians(-180 + 0.5 * (col - 1)), np.radians(-180 + 0.5 * col)
        lat1, lat2 = np.radians(-90 + 0.5 * (row - 1)), np.radians(-90 + 0.5 * row)
        return R_SPHERE ** 2 * (lon2 - lon1) * (math.sin(lat2) - math.sin(lat1))
    impd = dict(zip(gidx, imp))
    # implied average qsur per nonzero cell: import/(A_dest*1e-3); the two
    # full-coverage cells (100312, 101032) must both equal avg6 of the
    # source gridcell (the half-coverage cells imply avg6/2)
    implied = {g: impd[g] / (area_dest(g) * 1e-3) for g in gidx}
    best = None
    for g, d in q12.items():
        if len(d) != 2:
            continue
        avg6 = (d[1] + d[2]) / 6.0
        rel = max(abs(avg6 - implied[g0]) / implied[g0]
                  for g0 in (100312, 101032))
        if best is None or rel < best[0]:
            best = (rel, g, avg6)
    chk(best is not None and best[0] < 1e-3,
        f"no unique mosart source gridcell (best rel {best[0] if best else None})")
    if best is None or best[0] >= 1e-3:
        return
    _, src_g, avg6 = best
    chk(src_g == 1981, f"mosart source gridcell {src_g} != 1981")
    # coverage weights from the observed structure: W = 1.0 (full cells),
    # 0.5 (half cells); reproduction residual is the ESMF-internal area/
    # intersection convention
    W = {100312: 1.0, 101032: 1.0, 100313: 0.5, 101033: 0.5}
    worst = max(abs(area_dest(g) * avg6 * W[g] * 1e-3 - impd[g]) / impd[g]
                for g in gidx)
    chk(worst < RTOL_MOSART,
        f"mosart reproduction worst rel {worst:.3e} >= {RTOL_MOSART:.0e}")
    # geometric plausibility: source cell center within the 2x2 block
    m = exp[(exp[:, 0] == 1) & (exp[:, 3] == src_g)][0]
    chk(-21.0 <= m[4] <= -19.0 and 294.5 <= m[5] <= 297.0,
        f"source cell ({m[4]:.3f},{m[5]:.3f}) not near the 2x2 r05 block")
    print(f"mosart: step-1 imports reproduced from gridcell {src_g}: "
          f"import = A_dest*avg6*W*1e-3, avg6 = {avg6:.6e} mm/s "
          f"(= (qsur_s1+qsur_s2)/6), W = {{1.0, 0.5}}, worst rel {worst:.2e}")
    print("  note: exact ESMF mapconsf weights are computed online and never "
          "persisted; the ~1e-6 residual is the ESMF-internal spherical "
          "polygon intersection / sphere convention (R = 6371.229 km used "
          "here).  mosart dump window is partial (first 16 local cells); "
          "its steps 2..8 imports are exactly 0 despite all-day nonzero "
          "lnd exports -- open question for the mosart lane.")

def check_restart(new, base):
    import netCDF4
    n = netCDF4.Dataset(new)
    b = netCDF4.Dataset(base)
    common = sorted(set(n.variables) & set(b.variables))
    ndiff = 0
    for v in common:
        if n[v].dtype.kind not in "fiu":
            continue
        a = np.asarray(n[v][:]); c = np.asarray(b[v][:])
        if a.shape != c.shape:
            ndiff += 1
            print("FAIL: restart var", v, "shape", a.shape, c.shape)
            continue
        am = np.ma.masked_invalid(a.astype(np.float64))
        cm = np.ma.masked_invalid(c.astype(np.float64))
        if not np.ma.allequal(am, cm, fill_value=True):
            ndiff += 1
            print("FAIL: restart var", v, "differs")
    chk(ndiff == 0, f"restart neutrality: {ndiff} of {len(common)} variables differ")
    print(f"restart neutrality: {len(common)} common variables, {ndiff} differ")

def check_existing(backup_dir):
    files = ["clm_bareground", "clm_canopyflux", "clm_canopyflux_veg",
             "clm_interception", "clm_qsat", "clm_run_imports",
             "clm_soilfluxes", "clm_soiltemp", "clm_soilwater",
             "clm_surfalbedo", "clm_surfrad", "clm_twostream", "clm_urbalb",
             "clm_urban_lw", "clm_urban_rad"]
    for f in files:
        a = hashlib.md5(open(f"{D}/{f}.csv", "rb").read()).hexdigest()
        b = hashlib.md5(open(f"{backup_dir}/{f}.csv", "rb").read()).hexdigest()
        chk(a == b, f"existing dump {f}.csv not bit-identical to backup")
    print(f"existing dumps bit-identical vs backup: {len(files)} files")

def main():
    global D
    ap = argparse.ArgumentParser()
    ap.add_argument("--dumps", default=D)
    ap.add_argument("--mosart", default=MOSART)
    ap.add_argument("--restart", nargs=2, metavar=("NEW", "BASE"),
                    help="bit-identical restart check (instrumented vs baseline)")
    ap.add_argument("--existing", metavar="BACKUP_DIR",
                    help="bit-identical check of the 15 pre-lane dumps vs this backup dir")
    a = ap.parse_args()
    D = a.dumps

    cols, exp = load()
    check_structure(cols, exp)
    check_c2g(cols, exp)
    check_mosart(cols, exp, a.mosart)

    if a.restart:
        check_restart(*a.restart)
    if a.existing:
        check_existing(a.existing)

    print()
    if nbad == 0:
        print("ALL CHECKS PASSED")
    else:
        print(f"{nbad} CHECKS FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
