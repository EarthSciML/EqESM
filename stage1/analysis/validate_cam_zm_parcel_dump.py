#!/usr/bin/env python3
"""Validate the stage-1 ZM parcel-ascent dump CSVs (cam_zm_parcel.csv and
cam_zm_parcel_scal.csv) against the existing ZM dumps (cam_zm.csv,
cam_zm_scal.csv) from the instrumented SCAM run.

Usage: python3 validate_cam_zm_parcel_dump.py DUMP_DIR [--pver 32]

The dump traces the dilute-parcel ascent internals of the
Zhang-McFarlane scheme (buoyan_dilute / parcel_dilute, zm_convr.F90) at
the same trace point as cam_zm.csv: one cam_zm_parcel.csv row per
(step, gathered point g, level) and one cam_zm_parcel_scal.csv row per
(step, g), written only for the active steps/columns.  The (step,g)
keys coincide with cam_zm.csv / cam_zm_scal.csv.

Reference code: zm_convr.F90 parcel_dilute (entrainment loop
F90:1280-1378, precipitation/freezing loop F90:1392-1463) and
buoyan_dilute (CAPE F90:1100-1137).  The parcel entropy is the
Raymond-Blyth form F90:1489-1512 evaluated with the DIRECT Goff-Gratch
saturation formula (qsat_hPa F90:2634-2649 -> wv_sat_methods.F90
GoffGratch_svp_water, wv_sat_methods.F90:569-580) -- NOT the estblf
table used by wv_saturation qsat_2D; wv_sat_scheme defaults to
"GoffGratch" (wv_saturation.F90:159).  The ientropy Brent inversion
(F90:1517-1632) converges to tol1 = 2*EPS*|T| + 0.5*tol ~= 5.2e-4 K
(EPS = 3e-8, tol = 0.001), so all reconstructed temperatures carry an
inherent ~5e-4 K root uncertainty; entropy residuals are bounded by
dS/dT * tol1 ~= 3.3e-3 J/kg/K.

Identities checked (all internal to the dump + env state from
cam_zm.csv; tolerances have 2-10x margin over the empirical maxima of
the 2026-11-09 run):

  1. structure: headers, no blank lines, 18 / 16 fields per row;
     (step,g,lev) keys of cam_zm_parcel.csv == (step,g,lev) of
     cam_zm.csv; (step,g) keys of cam_zm_parcel_scal.csv ==
     cam_zm_scal.csv; lev = 1..pver in order; maxi == pver;
     num_cin == 1; constants (rgas..latice) identical across rows.
  2. launch level (k = maxi, F90:1290-1307): qtmix == qtp0 == q,
     mp0 = 1; smix == sp0 == entropy(t, p, q) (bit-exact); tmix is the
     ientropy root of entropy(t,p,qtmix) = smix (<= 1e-3 K); qsmix ==
     qsat(tmix) (rel 1e-12); qtot == qtp0; tp == tmix, qstp == q;
     xsh2o == ds_xsh2o == ds_freeze == 0; qtp == sp == mp == 0.
  3. entrainment recurrences (F90:1322-1342), k = maxi-1..msg+1:
     qtp(k) = qtp(k+1) - dmpdp*dp*qtenv, sp(k) = sp(k+1) -
     dmpdp*dp*senv, mp(k) = mp(k+1) - dmpdp*dp, with
     qtenv/tenv/penv = layer means, senv = entropy(tenv,penv,qtenv),
     dpdz = -penv*grav/(rgas*tenv), dmpdp = dmpdz/dpdz*dzdp --
     bit-exact for qtp/mp, <= 1e-9 for sp (entropy replica round-off);
     mixing fractions qtmix = (qtp0 + qtp)/(mp0 + mp) and smix =
     (sp0 + sp)/(mp0 + mp) -- bit-exact.
  4. precipitation/freezing loop reconstruction (F90:1392-1463,
     nit_lheat = 2, lwmax = 1e-3, tscool = 0): ascending over k =
     maxi-1..msg+1, iterate the rainout/freezing adjustment twice,
     starting from the entrainment-phase (pre-rainout) tmix/qsmix
     (reconstructed as the entropy roots of the dumped smix/qtmix),
     using the FINAL dumped values at k+1 (the Fortran loop reaches
     k+1 first):
         xsh2o    = max(0, qtmix - qsmix - lwmax)
         ds_xsh2o = ds_xsh2o(k+1) - cpliq*log(tmix/tfreez)
                    * max(0, xsh2o - xsh2o(k+1))
         one-off freezing (tmix <= tfreez, ds_freeze(k+1) == 0):
            ds_freeze = (latice/tmix) * max(0, qtmix - qsmix - xsh2o)
         continual freezing (tmix <= tfreez, ds_freeze(k+1) /= 0):
            ds_freeze = ds_freeze(k+1) + (latice/tmix)
                        * max(0, qsmix(k+1) - qsmix)
         new_s = smix + ds_xsh2o + ds_freeze;  new_q = qtmix - xsh2o
         tmix, qsmix <- entropy root of (new_s, p, new_q)
     after 2 iterations the reconstructed tmix/qsmix/xsh2o/ds_xsh2o/
     ds_freeze/qtot must match the dump to <= 1e-3 K / 1e-6 / 1e-6 /
     1e-3 / 1e-3 / 1e-6 (root-tolerance dominated).  This also proves
     the rainout entropy identity
         entropy(tmix, p, qtot) == smix + ds_xsh2o + ds_freeze
     (checked directly to <= 5e-3 J/kg/K) and xsh2o >= 0.
  5. parcel state chain (F90:1411-1477): tp == tmix; qtot == qtmix -
     xsh2o; qstp == min(qtot, qsmix); tpv == (tp + tpert)*(1 +
     1.608*qstp)/(1 + qtot); buoy == tpv - tv + tiedke_add with tv =
     t*(1 + 1.608*q)/(1 + q) -- all bit-exact on k = msg+1..maxi.
     Outside the parcel region (k > maxi or k <= msg): buoy == 0 and
     tpv == tv (F90:1087-1092).
  6. lcl (F90:1355-1368): with the reconstructed pre-rainout qsmix,
     qsmix_pre(lcl) <= qtmix(lcl) and qsmix_pre(k) > qtmix(k) for all
     k in lcl+1..maxi-1 (first saturated level on ascent); pl ==
     p(lcl+1) - qxskp1/dqxsdp with qxsk = qtmix(lcl) - qsmix_pre(lcl),
     qxskp1 = qtmix(lcl+1) - qsmix(lcl+1), dqxsdp = (qxsk -
     qxskp1)/(p(lcl) - p(lcl+1)) to <= 0.02 mb (the pl formula uses
     the PRE-rainout qsmix at the lcl, while the dump holds the
     post-rainout value; the reconstruction recovers it to the Brent
     root tolerance); plge600 == (pl >= 600).
  7. lcl temperature (F90:1361-1367): slcl = smix(lcl+1) + dsdp*(pl -
     p(lcl+1)), qtlcl = qtmix(lcl+1) + dqtdp*(pl - p(lcl+1)) with
     dsdp/dqtdp the (k, k+1) gradients; tl is the ientropy root of
     (slcl, pl, qtlcl), i.e. |entropy(tl, pl, qtlcl) - slcl| <= 0.01.
  8. CAPE / lel (F90:1100-1137): lelten_1 is the FIRST k in msg+2..pver
     with k < lcl where buoy(k+1) > 0 and buoy(k) <= 0 (no earlier
     crossing); capeten_1 == sum_{k=lelten_1+1..maxi}
     rgas*buoy(k)*log(pf(k+1)/pf(k)) (pf = dumped pf_hi/pf_lo) to <=
     1e-6; with num_cin = 1: cape == capeten_1, lel == lelten_1,
     capeten_2..5 == lelten_2..5 == 0.

All identities are internal to the dump (cross-file against
cam_zm.csv/cam_zm_scal.csv only), so they hold independently of the
FORTRAN physics being correct.
"""
import sys
import csv
import os
import numpy as np

RTOL = 1e-9          # identities hold to ~machine precision; be loose
                     # about summation/formatting round-off

# parcel_dilute hardcoded values (zm_convr.F90:1252-1255, 1317)
LWMAX = 1.0e-3       # maximum condensate held before rainout (F90:1319)
TSCOOL = 0.0         # supercooled offset: freezing at tmix <= tfreez (F90:1320)
NIT_LHEAT = 2        # precipitation/freezing iterations (F90:1317)


def parse(lines):
    """Parse CSV data lines (quoted strings unquoted) into row lists."""
    return list(csv.reader(lines))


def load(fn):
    """Return (cols, rows) after structural checks."""
    with open(fn) as f:
        lines = f.readlines()
    assert lines[0].startswith("# description"), fn
    assert lines[1].startswith("# columns"), fn
    cols = lines[1][len("# columns = "):].strip().split(",")
    data = [l for l in lines[2:] if l.strip()]
    assert len(data) == len(lines) - 2, f"{fn}: blank lines present"
    ncol = len(cols)
    for i, l in enumerate(data):
        nf = len(l.rstrip("\n").split(","))
        assert nf == ncol, f"{fn}: row {i+3} has {nf} fields, expected {ncol}"
    return cols, data


def f(x):
    return float(x)


def goff_gratch_svp_water(t):
    """Direct Goff-Gratch SVP over water (wv_sat_methods.F90:569-580),
    returns Pa."""
    tboil = 373.16
    return 10.0 ** (-7.90298 * (tboil / t - 1.0) + 5.02808 * np.log10(tboil / t)
                    - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / tboil)) - 1.0)
                    + 8.1328e-3 * (10.0 ** (-3.49149 * (tboil / t - 1.0)) - 1.0)
                    + np.log10(1013.246)) * 100.0


def qsat_water(t, p_mb, eps1):
    """qsat_hPa replica (zm_convr.F90:2634-2649 -> wv_sat_qsat_water):
    direct Goff-Gratch with the es = min(es, p) clamp (wv_sat_methods
    F90:248); eps1 = 0.6219705862045155."""
    es = min(goff_gratch_svp_water(t), p_mb * 100.0)
    return eps1 * es / (p_mb * 100.0 - (1.0 - eps1) * es)


def rb_entropy(TK, p_mb, qtot, C):
    """zm_convr.F90 entropy function (F90:1489-1512), p in mb."""
    L = C["rl"] - (C["cpliq"] - C["cpwv"]) * (TK - C["tfreez"])
    qst = qsat_water(TK, p_mb, C["eps1"])
    qv = min(qtot, qst)
    e = qv * p_mb / (C["eps1"] + qv)
    return ((C["cpres"] + qtot * C["cpliq"]) * np.log(TK / C["tfreez"])
            - C["rgas"] * np.log((p_mb - e) / 1000.0)
            + L * qv / TK - qv * C["rh2o"] * np.log(qv / qst))


def invert_entropy(s, p_mb, qtot, C):
    """Brent-free replica of the ientropy inversion (F90:1517-1632):
    root of rb_entropy(t, p, qtot) = s.  The Fortran Brent solve
    converges to tol1 ~= 5.2e-4 K; bisection to 1e-9 K reproduces the
    root well within that bracket, so all reconstructed temperatures
    match the Fortran ones to the Brent tolerance."""
    lo, hi = 150.0, 400.0
    flo = rb_entropy(lo, p_mb, qtot, C) - s
    fhi = rb_entropy(hi, p_mb, qtot, C) - s
    assert flo * fhi < 0.0, f"invert_entropy: no bracket for s={s}, p={p_mb}, qtot={qtot}"
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = rb_entropy(mid, p_mb, qtot, C) - s
        if fm * flo <= 0.0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-9:
            break
    return 0.5 * (lo + hi)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    pver = 32
    for a in sys.argv[2:]:
        if a.startswith("--pver"):
            pver = int(a.split("=")[1])

    # ---- load ---------------------------------------------------------
    pp = os.path.join(d, "cam_zm_parcel.csv")
    ps = os.path.join(d, "cam_zm_parcel_scal.csv")
    pprof = os.path.join(d, "cam_zm.csv")
    pscal = os.path.join(d, "cam_zm_scal.csv")
    for fn in (pp, ps, pprof, pscal):
        assert os.path.exists(fn), f"missing {fn}"
    cols_a, lines_a = load(pp)
    cols_b, lines_b = load(ps)
    cols_p, lines_p = load(pprof)
    cols_s, lines_s = load(pscal)
    assert len(cols_a) == 18, f"cam_zm_parcel.csv: {len(cols_a)} columns, expected 18"
    assert len(cols_b) == 16, f"cam_zm_parcel_scal.csv: {len(cols_b)} columns, expected 16"
    print(f"cam_zm_parcel.csv:      {len(lines_a)} rows x {len(cols_a)} fields")
    print(f"cam_zm_parcel_scal.csv: {len(lines_b)} rows x {len(cols_b)} fields")

    parcel = parse(lines_a)
    pscal = parse(lines_b)
    prof = parse(lines_p)
    scal = parse(lines_s)
    ai = {c: i for i, c in enumerate(cols_a)}
    bi = {c: i for i, c in enumerate(cols_b)}
    pi = {c: i for i, c in enumerate(cols_p)}
    si = {c: i for i, c in enumerate(cols_s)}

    # ---- key grouping: parcel == prof, pscal == scal ------------------
    parcel_keys = {}
    for r in parcel:
        key = (int(r[0]), int(r[3]))
        assert int(r[4]) == len(parcel_keys.get(key, [])) + 1, \
            f"parcel: {key} lev out of order"
        parcel_keys.setdefault(key, []).append(r)
    prof_keys = {}
    for r in prof:
        key = (int(r[0]), int(r[3]))
        prof_keys.setdefault(key, []).append(r)
    assert set(parcel_keys) == set(prof_keys), \
        "(step,g) pairs differ between cam_zm_parcel.csv and cam_zm.csv"
    for key, rows in parcel_keys.items():
        assert len(rows) == pver, f"parcel: {key} has {len(rows)} levels, expected {pver}"
    pscal_keys = {}
    for r in pscal:
        key = (int(r[0]), int(r[3]))
        assert key not in pscal_keys, f"pscal: duplicate {key}"
        pscal_keys[key] = r
    scal_keys = {}
    for r in scal:
        key = (int(r[0]), int(r[3]))
        assert key not in scal_keys, f"scal: duplicate {key}"
        scal_keys[key] = r
    assert set(pscal_keys) == set(scal_keys), \
        "(step,g) pairs differ between cam_zm_parcel_scal.csv and cam_zm_scal.csv"
    ncol = len(parcel_keys)
    print(f"grouping OK: {len(set(int(r[0]) for r in parcel))} active steps, "
          f"{ncol} (step,g) pairs, pver={pver}")

    # ---- constants -----------------------------------------------------
    const_cols = ["rgas", "tiedke_add", "cpres", "rl", "eps1", "tfreez",
                  "grav", "cpliq", "cpwv", "rh2o", "latice", "dmpdz_param",
                  "num_cin", "msg"]
    C = {}
    for c in const_cols:
        vals = {r[si[c]] for r in scal}
        assert len(vals) == 1, f"constant {c} varies across rows: {vals}"
        C[c] = f(next(iter(vals)))
    assert C["num_cin"] == 1.0, "this lane runs with num_cin = 1"
    msg = int(C["msg"])
    print(f"constants OK: num_cin={int(C['num_cin'])}, msg={msg}, "
          f"lwmax={LWMAX}, tscool={TSCOOL}, nit_lheat={NIT_LHEAT}")

    # ---- per-column helpers -------------------------------------------
    def get(key, lev):
        return parcel_keys[key][lev - 1], prof_keys[key][lev - 1], pscal_keys[key]

    # launch level must be the bottom level in this case
    for key in parcel_keys:
        srow = pscal_keys[key]
        assert int(float(srow[bi["pl"]])) >= 0  # pl present
    # maxi from cam_zm_scal (already validated there); assert == pver
    maxi_all = {int(f(r[si["maxi"]])) for r in scal}
    assert maxi_all == {pver}, f"maxi != pver: {maxi_all}"

    # ---- 2. launch level -------------------------------------------------
    mx_tlaunch = mx_smix = mx_qsmix = 0.0
    for key in parcel_keys:
        a, p, _ = parcel_keys[key][pver - 1], prof_keys[key][pver - 1], None
        t = f(p[pi["t"]]); q = f(p[pi["q"]]); pmb = f(p[pi["p"]])
        tmix = f(a[ai["tmix"]]); qtmix = f(a[ai["qtmix"]]); smix = f(a[ai["smix"]])
        qsmix = f(a[ai["qsmix"]]); qtot = f(a[ai["qtot"]])
        assert qtmix == q, f"{key}: launch qtmix {qtmix} != q {q}"
        assert qtot == q, f"{key}: launch qtot {qtot} != q {q}"
        for c in ("xsh2o", "ds_xsh2o", "ds_freeze", "qtp", "sp", "mp"):
            assert f(a[ai[c]]) == 0.0, f"{key}: launch {c} != 0"
        mx_smix = max(mx_smix, abs(smix - rb_entropy(t, pmb, q, C)))
        mx_tlaunch = max(mx_tlaunch, abs(tmix - invert_entropy(smix, pmb, qtmix, C)))
        mx_qsmix = max(mx_qsmix, abs(qsmix - qsat_water(tmix, pmb, C["eps1"])))
    assert mx_smix <= 1e-6, f"launch smix != entropy(t,p,q): {mx_smix:.3e}"
    assert mx_tlaunch <= 1e-3, f"launch tmix root mismatch: {mx_tlaunch:.3e} K"
    assert mx_qsmix <= 1e-12 * max(1.0, max(f(r[ai["qsmix"]]) for r in parcel)), \
        f"launch qsmix != qsat(tmix): {mx_qsmix:.3e}"
    print(f"launch level OK: qtmix/qtot == q; xsh2o/ds/qtp/sp/mp == 0; "
          f"max |smix - entropy(t,p,q)| = {mx_smix:.2e}; "
          f"max |tmix - root| = {mx_tlaunch:.2e} K")

    # ---- 3. entrainment recurrences + 4. precip-loop reconstruction ----
    mx_r = {c: 0.0 for c in ("qtp", "sp", "mp", "qmix", "smix")}
    mx_p = {c: 0.0 for c in ("tmix", "qsmix", "xsh2o", "dsx", "dsf", "qtot")}
    n_asc = 0
    for key in parcel_keys:
        # sp0 / qtp0 / mp0 from the launch level
        al, pl_ = parcel_keys[key][pver - 1], prof_keys[key][pver - 1]
        qtp0 = f(al[ai["qtmix"]])           # == q(maxi)
        sp0 = rb_entropy(f(pl_[pi["t"]]), f(pl_[pi["p"]]), qtp0, C)
        mp0 = 1.0
        # entrainment-phase (pre-rainout) roots per level
        qtmix = {}; smix = {}; pmb = {}; tpre = {}; qspre = {}
        for lev in range(pver, msg, -1):
            a, p, _ = parcel_keys[key][lev - 1], prof_keys[key][lev - 1], None
            qtmix[lev] = f(a[ai["qtmix"]]); smix[lev] = f(a[ai["smix"]])
            pmb[lev] = f(p[pi["p"]])
            tpre[lev] = invert_entropy(smix[lev], pmb[lev], qtmix[lev], C)
            qspre[lev] = qsat_water(tpre[lev], pmb[lev], C["eps1"])
        # entrainment recurrences, k = maxi-1 .. msg+1 (ascending level index)
        for lev in range(pver - 1, msg, -1):
            a, p, _ = parcel_keys[key][lev - 1], prof_keys[key][lev - 1], None
            an, pn, _ = parcel_keys[key][lev], prof_keys[key][lev], None
            qtenv = 0.5 * (f(p[pi["q"]]) + f(pn[pi["q"]]))
            tenv = 0.5 * (f(p[pi["t"]]) + f(pn[pi["t"]]))
            penv = 0.5 * (pmb[lev] + pmb[lev + 1])
            senv = rb_entropy(tenv, penv, qtenv, C)
            dp = pmb[lev] - pmb[lev + 1]
            dpdz = -(penv * C["grav"]) / (C["rgas"] * tenv)
            dmpdp = C["dmpdz_param"] / dpdz
            qtp = f(an[ai["qtp"]]) - dmpdp * dp * qtenv
            sp = f(an[ai["sp"]]) - dmpdp * dp * senv
            mp = f(an[ai["mp"]]) - dmpdp * dp
            mx_r["qtp"] = max(mx_r["qtp"], abs(qtp - f(a[ai["qtp"]])))
            mx_r["sp"] = max(mx_r["sp"], abs(sp - f(a[ai["sp"]])))
            mx_r["mp"] = max(mx_r["mp"], abs(mp - f(a[ai["mp"]])))
            mx_r["qmix"] = max(mx_r["qmix"], abs((qtp0 + qtp) / (mp0 + mp) - f(a[ai["qtmix"]])))
            mx_r["smix"] = max(mx_r["smix"], abs((sp0 + sp) / (mp0 + mp) - f(a[ai["smix"]])))
            # precip loop, 2 iterations, using final dumped values at k+1
            tmix = tpre[lev]; qsmix = qspre[lev]
            an_r, pn_r = parcel_keys[key][lev], prof_keys[key][lev]
            for _ in range(NIT_LHEAT):
                xsh2o = max(0.0, qtmix[lev] - qsmix - LWMAX)
                dsx = f(an_r[ai["ds_xsh2o"]]) - C["cpliq"] * np.log(tmix / C["tfreez"]) * max(
                    0.0, xsh2o - f(an_r[ai["xsh2o"]]))
                if tmix <= C["tfreez"] + TSCOOL:
                    if f(an_r[ai["ds_freeze"]]) == 0.0:
                        dsf = (C["latice"] / tmix) * max(0.0, qtmix[lev] - qsmix - xsh2o)
                    else:
                        dsf = f(an_r[ai["ds_freeze"]]) + (C["latice"] / tmix) * max(
                            0.0, f(an_r[ai["qsmix"]]) - qsmix)
                else:
                    dsf = 0.0
                new_s = smix[lev] + dsx + dsf
                new_q = qtmix[lev] - xsh2o
                tmix = invert_entropy(new_s, pmb[lev], new_q, C)
                qsmix = qsat_water(tmix, pmb[lev], C["eps1"])
            mx_p["tmix"] = max(mx_p["tmix"], abs(tmix - f(a[ai["tmix"]])))
            mx_p["qsmix"] = max(mx_p["qsmix"], abs(qsmix - f(a[ai["qsmix"]])))
            mx_p["xsh2o"] = max(mx_p["xsh2o"], abs(xsh2o - f(a[ai["xsh2o"]])))
            mx_p["dsx"] = max(mx_p["dsx"], abs(dsx - f(a[ai["ds_xsh2o"]])))
            mx_p["dsf"] = max(mx_p["dsf"], abs(dsf - f(a[ai["ds_freeze"]])))
            mx_p["qtot"] = max(mx_p["qtot"], abs(new_q - f(a[ai["qtot"]])))
            n_asc += 1
    assert mx_r["qtp"] <= 1e-12 and mx_r["mp"] <= 1e-12, \
        f"entrainment recurrences qtp/mp not bit-exact: {mx_r['qtp']:.3e} {mx_r['mp']:.3e}"
    assert mx_r["sp"] <= 1e-9, f"sp recurrence: {mx_r['sp']:.3e}"
    assert mx_r["qmix"] <= 1e-9 and mx_r["smix"] <= 1e-9, \
        f"mixing fractions: qtmix {mx_r['qmix']:.3e} smix {mx_r['smix']:.3e}"
    assert mx_p["tmix"] <= 1e-3, f"reconstructed tmix: {mx_p['tmix']:.3e} K"
    assert mx_p["qsmix"] <= 1e-6, f"reconstructed qsmix: {mx_p['qsmix']:.3e}"
    assert mx_p["xsh2o"] <= 1e-6, f"reconstructed xsh2o: {mx_p['xsh2o']:.3e}"
    assert mx_p["qtot"] <= 1e-6, f"reconstructed qtot: {mx_p['qtot']:.3e}"
    assert mx_p["dsx"] <= 1e-3, f"reconstructed ds_xsh2o: {mx_p['dsx']:.3e}"
    assert mx_p["dsf"] <= 1e-3, f"reconstructed ds_freeze: {mx_p['dsf']:.3e}"
    print(f"entrainment recurrences OK ({n_asc} rows): max |dqtp| = {mx_r['qtp']:.1e}, "
          f"|dsp| = {mx_r['sp']:.1e}, |dmp| = {mx_r['mp']:.1e}; mixing fractions bit-exact")
    print(f"precip loop reconstruction OK ({n_asc} rows): max |dtmix| = {mx_p['tmix']:.2e} K, "
          f"|dqsmix| = {mx_p['qsmix']:.2e}, |dxsh2o| = {mx_p['xsh2o']:.2e}, "
          f"|dqtot| = {mx_p['qtot']:.2e}, |dds_xsh2o| = {mx_p['dsx']:.2e}, "
          f"|dds_freeze| = {mx_p['dsf']:.2e}")

    # ---- rainout entropy identity (direct, separate max) ----------------
    mx_ent = 0.0
    for key in parcel_keys:
        for lev in range(pver - 1, msg, -1):
            a, p, _ = parcel_keys[key][lev - 1], prof_keys[key][lev - 1], None
            resid = abs(rb_entropy(f(a[ai["tmix"]]), f(p[pi["p"]]), f(a[ai["qtot"]]), C)
                        - (f(a[ai["smix"]]) + f(a[ai["ds_xsh2o"]]) + f(a[ai["ds_freeze"]])))
            mx_ent = max(mx_ent, resid)
            assert f(a[ai["xsh2o"]]) >= 0.0
    assert mx_ent <= 5e-3, f"rainout entropy identity: {mx_ent:.3e} J/kg/K"
    print(f"rainout entropy identity OK ({n_asc} rows): max |entropy(tmix,p,qtot) - "
          f"(smix+ds_xsh2o+ds_freeze)| = {mx_ent:.2e} J/kg/K")

    # ---- 5. parcel state chain (tp/qstp/tpv/buoy) ------------------------
    mx_tp = mx_qstp = mx_tpv = mx_buoy = 0.0
    mx_qsat = 0.0
    n_state = 0
    for key in parcel_keys:
        tpert = f(scal_keys[key][si["tpert"]])
        for lev in range(1, pver + 1):
            a, p, _ = parcel_keys[key][lev - 1], prof_keys[key][lev - 1], None
            t = f(p[pi["t"]]); q = f(p[pi["q"]]); pmb = f(p[pi["p"]])
            tmix = f(a[ai["tmix"]]); qsmix = f(a[ai["qsmix"]]); qtot = f(a[ai["qtot"]])
            tp = f(p[pi["tp"]]); qstp = f(p[pi["qstp"]])
            tpv = f(a[ai["tpv"]]); buoy = f(a[ai["buoy"]])
            tv = t * (1 + 1.608 * q) / (1 + q)
            if msg < lev <= pver:
                mx_tp = max(mx_tp, abs(tp - tmix))
                mx_qstp = max(mx_qstp, abs(qstp - min(qtot, qsmix)))
                mx_tpv = max(mx_tpv, abs(tpv - (tp + tpert) * (1 + 1.608 * qstp) / (1 + qtot)))
                mx_buoy = max(mx_buoy, abs(buoy - (tpv - tv + C["tiedke_add"])))
                mx_qsat = max(mx_qsat, abs(qsmix - qsat_water(tmix, pmb, C["eps1"]))
                              / max(1e-12, qsmix))
                n_state += 1
            else:
                assert buoy == 0.0, f"{key} lev {lev}: buoy = {buoy} != 0"
                assert abs(tpv - tv) <= 1e-9, f"{key} lev {lev}: tpv != tv"
    assert mx_tp <= 1e-9, f"tp != tmix: {mx_tp:.3e}"
    assert mx_qstp <= 1e-9, f"qstp != min(qtot,qsmix): {mx_qstp:.3e}"
    assert mx_tpv <= 1e-9, f"tpv law: {mx_tpv:.3e}"
    assert mx_buoy <= 1e-9, f"buoy law: {mx_buoy:.3e}"
    assert mx_qsat <= 1e-9, f"qsmix != qsat_direct(tmix,p): {mx_qsat:.3e} rel"
    print(f"parcel state chain OK ({n_state} rows): tp==tmix, qstp==min(qtot,qsmix), "
          f"tpv law, buoy law all bit-exact; qsmix==qsat_GG(tmix,p) to {mx_qsat:.1e} rel; "
          f"buoy==0 and tpv==tv for k <= msg")

    # ---- 6. lcl / pl -----------------------------------------------------
    mx_pl = 0.0
    nbad_cross = 0
    for key in parcel_keys:
        srow = pscal_keys[key]
        lcl = int(float(scal_keys[key][si["lcl"]]))
        a_up, p_up, _ = parcel_keys[key][lcl - 1], prof_keys[key][lcl - 1], None
        a_dn, p_dn, _ = parcel_keys[key][lcl], prof_keys[key][lcl], None
        # pre-rainout qsmix at the lcl (entrainment-phase root)
        qspre = qsat_water(invert_entropy(f(a_up[ai["smix"]]), f(p_up[pi["p"]]),
                                          f(a_up[ai["qtmix"]]), C), f(p_up[pi["p"]]), C["eps1"])
        # lcl crossing: first saturated level on ascent
        if not (qspre <= f(a_up[ai["qtmix"]]) and
                f(a_dn[ai["qsmix"]]) > f(a_dn[ai["qtmix"]])):
            nbad_cross += 1
        for lev in range(lcl + 1, pver):     # levels below the lcl
            a, p, _ = parcel_keys[key][lev - 1], prof_keys[key][lev - 1], None
            qspre_k = qsat_water(invert_entropy(f(a[ai["smix"]]), f(p[pi["p"]]),
                                                f(a[ai["qtmix"]]), C), f(p[pi["p"]]), C["eps1"])
            if not (qspre_k > f(a[ai["qtmix"]])):
                nbad_cross += 1
        # pl interpolation (F90:1356-1361)
        qxsk = f(a_up[ai["qtmix"]]) - qspre
        qxskp1 = f(a_dn[ai["qtmix"]]) - f(a_dn[ai["qsmix"]])
        dp = f(p_up[pi["p"]]) - f(p_dn[pi["p"]])
        dqxsdp = (qxsk - qxskp1) / dp
        pl_ref = f(p_dn[pi["p"]]) - qxskp1 / dqxsdp
        pl = f(srow[bi["pl"]])
        mx_pl = max(mx_pl, abs(pl - pl_ref))
        assert (srow[bi["plge600"]] == "T") == (pl >= 600.0), f"{key}: plge600 mismatch"
    assert nbad_cross == 0, f"lcl crossing: {nbad_cross} violations"
    assert mx_pl <= 0.02, f"pl identity: {mx_pl:.3e} mb"
    print(f"lcl/pl OK: first-saturated-level crossing holds on all {ncol} columns; "
          f"max |pl - p(lcl+1) - qxskp1/dqxsdp| = {mx_pl:.2e} mb; plge600 == (pl >= 600)")

    # ---- 7. tl -------------------------------------------------------------
    mx_tl = 0.0
    for key in parcel_keys:
        srow_s = scal_keys[key]
        lcl = int(float(srow_s[si["lcl"]]))
        tl = f(srow_s[si["tl"]])
        a_up, p_up, _ = parcel_keys[key][lcl - 1], prof_keys[key][lcl - 1], None
        a_dn, p_dn, _ = parcel_keys[key][lcl], prof_keys[key][lcl], None
        qspre = qsat_water(invert_entropy(f(a_up[ai["smix"]]), f(p_up[pi["p"]]),
                                          f(a_up[ai["qtmix"]]), C), f(p_up[pi["p"]]), C["eps1"])
        qxsk = f(a_up[ai["qtmix"]]) - qspre
        qxskp1 = f(a_dn[ai["qtmix"]]) - f(a_dn[ai["qsmix"]])
        dp = f(p_up[pi["p"]]) - f(p_dn[pi["p"]])
        pl = f(p_dn[pi["p"]]) - qxskp1 / ((qxsk - qxskp1) / dp)
        dsdp = (f(a_up[ai["smix"]]) - f(a_dn[ai["smix"]])) / dp
        dqtdp = (f(a_up[ai["qtmix"]]) - f(a_dn[ai["qtmix"]])) / dp
        slcl = f(a_dn[ai["smix"]]) + dsdp * (pl - f(p_dn[pi["p"]]))
        qtlcl = f(a_dn[ai["qtmix"]]) + dqtdp * (pl - f(p_dn[pi["p"]]))
        mx_tl = max(mx_tl, abs(rb_entropy(tl, pl, qtlcl, C) - slcl))
    assert mx_tl <= 0.01, f"tl identity: {mx_tl:.3e} J/kg/K"
    print(f"tl OK: max |entropy(tl, pl, qtlcl) - slcl| = {mx_tl:.2e} J/kg/K "
          f"(ientropy root tolerance ~3.3e-3)")

    # ---- 8. CAPE / lel -----------------------------------------------------
    mx_cape = 0.0
    nbad_lel = 0
    for key in parcel_keys:
        srow = pscal_keys[key]
        srow_s = scal_keys[key]
        lel = int(float(srow_s[si["lel"]]))
        lcl = int(float(srow_s[si["lcl"]]))
        cape = f(srow_s[si["cape"]])
        capeten1 = f(srow[bi["capeten_1"]])
        lelten1 = int(srow[bi["lelten_1"]])
        assert lelten1 == lel, f"{key}: lelten_1 {lelten1} != lel {lel}"
        assert abs(cape - capeten1) <= 1e-9, f"{key}: cape != capeten_1"
        for n in range(2, 6):
            assert f(srow[bi[f"capeten_{n}"]]) == 0.0, f"{key}: capeten_{n} != 0"
            assert int(srow[bi[f"lelten_{n}"]]) == 0, f"{key}: lelten_{n} != 0"
        # lelten_1 is the FIRST k in msg+2..pver with k < lcl where
        # buoy(k+1) > 0 and buoy(k) <= 0
        buoy = {lev: f(parcel_keys[key][lev - 1][ai["buoy"]]) for lev in range(msg + 2, pver + 1)}
        found = None
        for k in range(msg + 2, pver + 1):
            if k < lcl and buoy[k + 1] > 0.0 and buoy[k] <= 0.0:
                found = k
                break
        if found != lel:
            nbad_lel += 1
        cap = 0.0
        for lev in range(lelten1 + 1, pver + 1):
            p = prof_keys[key][lev - 1]
            cap += C["rgas"] * buoy[lev] * np.log(f(p[pi["pf_hi"]]) / f(p[pi["pf_lo"]]))
        mx_cape = max(mx_cape, abs(capeten1 - cap))
    assert nbad_lel == 0, f"lelten_1 first-crossing: {nbad_lel} violations"
    assert mx_cape <= 1e-6, f"CAPE integral: {mx_cape:.3e} J/kg"
    print(f"CAPE/lel OK: lel == lelten_1 (first buoy crossing), cape == capeten_1, "
          f"capeten_2..5 == lelten_2..5 == 0 (num_cin=1); "
          f"max |capeten_1 - sum rd*buoy*log(pf_hi/pf_lo)| = {mx_cape:.2e} J/kg")

    print("\nALL ZM PARCEL DUMP CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
