#!/usr/bin/env python3
"""Validate the stage-1 ZM deep-convection dump CSVs (cam_zm.csv and
cam_zm_scal.csv) from the instrumented SCAM run.

Usage: python3 validate_cam_zm_dump.py DUMP_DIR [--pver 32] [--nsteps 48]

The dump traces the full Zhang-McFarlane deep-convection chain of
convect_deep.ZM (tphysbc site 3, zm_conv_intr.zm_conv_tend) at its four
boundaries: zm_convr (zm_convr.F90 zm_convr_run), zm_conv_evap,
zm_conv_momtran and zm_conv_convtran (zm_conv_intr.F90).  Rows are
written only for active steps/columns (cape > capelmt, lengath > 0):
one cam_zm.csv row per (step, gathered point g, level), one
cam_zm_scal.csv row per (step, g).  Context columns step,ymd,tod are
prepended to every row by dump_ctx and lead the `# columns` headers,
so cam_zm.csv has 83 fields and cam_zm_scal.csv 60.

Checks:
  1. structure: files exist, `# description`/`# columns` headers, no
     blank lines, consistent field counts (83 / 60);
  2. row grouping: the (step,g) pairs of the two files coincide; every
     pair has exactly pver rows with lev = 1..pver in order; within a
     step g runs 1..lengath in order and lengath equals the row count
     of that step; steps in 1..nsteps; ymd/tod constant within a step;
  3. the ZM-constants block of cam_zm_scal.csv (c0_lnd..latice) is
     identical on every row and matches the zmconv namelist values of
     this lane; msg == limcnv-1; lparcel_pbl / no_deep_pbl == F;
  4. updraft mass continuity (cldprp, zm_convr.F90:2026-2107, with the
     1/m -> 1/mb conversion at zm_convr.F90:651-661 and the closure
     scaling by mb at zm_convr.F90:700-714):
         mu(k) - mu(k+1) = (eu(k) - du(k)) * dp(k),   dp = 0.01*pdel
     for k = 1..pver-1.  NOTE: eu/du are dumped in their final 1/s
     units (per-mb rates scaled by mb), so the continuity identity is
     in pressure thickness dp, not geometric depth dz;
  5. precipitation / snow production sums (zm_conv_evap.F90:219-224,
     240-241, flxprec(1) = 0): the FORTRAN accumulates the interface
     fluxes as flxprec(k+1) = max(flxprec(k) + ntprprd(k)*pdel(k)/
     gravit, 0) (the max is a roundoff guard) and converts to m/s only
     at the surface, so:
         prec_evap  = clamped column chain of ntprprd*pdel/(grav*1000)
         snow_evap  = clamped column chain of ntsnprd*pdel/(grav*1000)
  6. interface flux increments (zm_conv_evap.F90:219-224), per row:
         flxprec_hi = max(flxprec_lo + ntprprd*pdel/grav, 0)   (snow
     analog) plus the boundary values flxprec_lo(lev=1) = 0 and
     flxprec_hi(lev=pver) = prec_evap*1000 (snow analog);
  7. mass-flux bookkeeping, per row: mc = mu + md (zm_convr.F90:2311,
     702-704).  (The dump `# description` also mentions mcon [kg/m2/s]
     but no mcon column is written; mcon = mc*100/grav, zm_convr.F90:
     769.)
  8. working arrays (zm_convr.F90:498-518): p = 0.01*pmid,
     pf_lo/pf_hi = 0.01*pint_lo/pint_hi, z = zm + phis/grav
     (zs = geos*rgrav, rgrav = 1/grav), s = t + (grav/cpres)*z;
  9. trigger / index sanity: cape > capelmt on every row; 1 <= jt <=
     maxg <= pver; jctop == jt and jcbot == maxg (zm_conv_intr.F90:
     581-586); other level indices within [0, pver].

All identities are internal to the dump (no cross-file references), so
they hold independently of the FORTRAN physics being correct.
"""
import sys
import csv
import os
import numpy as np

RTOL = 1e-9          # identities hold to ~machine precision; be loose
                     # about summation/formatting round-off


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


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    pver = 32
    nsteps = 48
    for a in sys.argv[2:]:
        if a.startswith("--pver"):
            pver = int(a.split("=")[1])
        elif a.startswith("--nsteps"):
            nsteps = int(a.split("=")[1])

    # ---- load ---------------------------------------------------------
    pprof = os.path.join(d, "cam_zm.csv")
    pscal = os.path.join(d, "cam_zm_scal.csv")
    assert os.path.exists(pprof), f"missing {pprof}"
    assert os.path.exists(pscal), f"missing {pscal}"
    cols_p, lines_p = load(pprof)
    cols_s, lines_s = load(pscal)
    assert len(cols_p) == 83, f"cam_zm.csv: {len(cols_p)} columns, expected 83"
    assert len(cols_s) == 60, f"cam_zm_scal.csv: {len(cols_s)} columns, expected 60"
    print(f"cam_zm.csv:      {len(lines_p)} rows x {len(cols_p)} fields")
    print(f"cam_zm_scal.csv: {len(lines_s)} rows x {len(cols_s)} fields")

    prof = parse(lines_p)          # step,ymd,tod + 80 named columns
    scal = parse(lines_s)
    ci = {c: i for i, c in enumerate(cols_s)}   # index into the full row
    pi = {c: i for i, c in enumerate(cols_p)}

    # ---- row grouping --------------------------------------------------
    # scal: group by step; check g == 1..lengath and lengath == nrows
    scal_steps = {}
    for r in scal:
        st, g, leng = int(r[0]), int(r[3]), int(r[5])
        assert 1 <= st <= nsteps, f"scal: step {st} out of 1..{nsteps}"
        assert g == len(scal_steps.get(st, [])) + 1, \
            f"scal: step {st} g {g} out of order"
        assert leng == int(r[5]) and leng >= 1, f"scal: step {st} bad lengath"
        scal_steps.setdefault(st, []).append(r)
    for st, rows in scal_steps.items():
        assert len(rows) == int(rows[0][5]), \
            f"scal: step {st} has {len(rows)} rows, lengath={rows[0][5]}"

    # prof: (step,g) pairs must match scal; lev 1..pver in order per pair
    prof_pairs = {}
    for r in prof:
        key = (int(r[0]), int(r[3]))
        assert int(r[4]) == len(prof_pairs.get(key, [])) + 1, \
            f"prof: {key} lev out of order"
        prof_pairs.setdefault(key, []).append(r)
    scal_pairs = set()
    for st, rows in scal_steps.items():
        for r in rows:
            scal_pairs.add((st, int(r[3])))
    assert set(prof_pairs) == scal_pairs, \
        "(step,g) pairs differ between cam_zm.csv and cam_zm_scal.csv"
    for key, rows in prof_pairs.items():
        assert len(rows) == pver, f"prof: {key} has {len(rows)} levels, expected {pver}"
    ncol_pairs = len(scal_pairs)
    print(f"grouping OK: {len(scal_steps)} active steps, "
          f"{ncol_pairs} (step,g) pairs, pver={pver}")

    # ymd/tod constant within each step
    for st, rows in scal_steps.items():
        ymd = {r[1] for r in rows}
        tod = {r[2] for r in rows}
        assert len(ymd) == 1 and len(tod) == 1, \
            f"scal: step {st} ymd/tod not constant"
    print("ymd/tod constant within each step")

    # ---- constants block ----------------------------------------------
    const_cols = ["c0_lnd", "c0_ocn", "ke", "ke_lnd", "momcu", "momcd",
                  "num_cin", "tau", "capelmt", "dmpdz_param", "tiedke_add",
                  "parcel_hscale", "lparcel_pbl", "no_deep_pbl", "limcnv",
                  "msg", "cpres", "rl", "rgas", "eps1", "tfreez", "grav",
                  "cpliq", "cpwv", "rh2o", "latice"]
    consts = {}
    for c in const_cols:
        vals = {r[ci[c]] for r in scal}
        assert len(vals) == 1, f"constant {c} varies across rows: {vals}"
        v = next(iter(vals))
        consts[c] = v if c in ("lparcel_pbl", "no_deep_pbl") else f(v)
    # namelist values of this lane (zm_convr_run readnl, FSCAMARM97_ne3
    # user_nl_cam): zmconv_c0_lnd/ocn, zmconv_ke/ke_lnd, zmconv_momcu/momcd,
    # zmconv_num_cin, zmconv_tau, zmconv_capelmt, zmconv_dmpdz,
    # zmconv_tiedtke_add, zmconv_parcel_hscale, zmconv_parcel_pbl,
    # phys_deepconv_pbl
    expected = {
        "c0_lnd": 0.0075, "c0_ocn": 0.03, "ke": 5e-6, "ke_lnd": 1e-5,
        "momcu": 0.7, "momcd": 0.7, "num_cin": 1, "tau": 3600.0,
        "capelmt": 70.0, "dmpdz_param": -1e-3, "tiedke_add": 0.5,
        "parcel_hscale": 0.5,
    }
    for c, v in expected.items():
        assert abs(consts[c] - v) <= RTOL * max(abs(v), 1.0), \
            f"{c} = {consts[c]}, expected {v}"
    assert consts["lparcel_pbl"] == "F" and consts["no_deep_pbl"] == "F", \
        "lparcel_pbl / no_deep_pbl should be F"
    assert consts["msg"] == consts["limcnv"] - 1, \
        f"msg={consts['msg']} != limcnv-1={consts['limcnv']-1}"
    print(f"constants OK: namelist values match; msg={consts['msg']:.0f}, "
          f"limcnv={consts['limcnv']:.0f}, grav={consts['grav']:.5f}")

    grav = consts["grav"]
    cpres = consts["cpres"]
    capelmt = consts["capelmt"]

    # helper: rows of one (step,g) profile, keyed by lev
    def profile(key):
        return {int(r[4]): r for r in prof_pairs[key]}

    # ---- 4. updraft mass continuity ------------------------------------
    nbad = 0
    maxres = 0.0
    for key, rows in prof_pairs.items():
        pr = profile(key)
        for k in range(1, pver):                 # skip lev=pver (no k+1)
            r = pr[k]
            rn = pr[k + 1]
            mu_k, eu_k, du_k = f(r[pi["mu"]]), f(r[pi["eu"]]), f(r[pi["du"]])
            mu_n = f(rn[pi["mu"]])
            dp = 0.01 * f(r[pi["pdel"]])
            lhs = mu_k - mu_n
            rhs = (eu_k - du_k) * dp
            scale = max(abs(mu_k), abs(mu_n), abs(rhs))
            res = abs(lhs - rhs)
            maxres = max(maxres, res / scale if scale > 0 else (1.0 if res else 0.0))
            if res > RTOL * scale:
                nbad += 1
                if nbad <= 5:
                    print(f"  mu continuity {key} lev {k}: {lhs:.6e} vs {rhs:.6e}")
    assert nbad == 0, f"mu continuity: {nbad} bad rows"
    print(f"mu continuity OK: mu(k)-mu(k+1) = (eu-du)*0.01*pdel "
          f"({ncol_pairs*(pver-1)} rows, max rel resid {maxres:.2e})")

    # ---- 5. precip/snow production sums --------------------------------
    # The FORTRAN accumulates flxprec(k+1) = max(flxprec(k) +
    # ntprprd*pdel/gravit, 0) level by level and divides by 1000 only at
    # the end (zm_conv_evap.F90:219-224, 240-241); the max() is a
    # roundoff guard that makes the plain sum differ from the final flux
    # by at most a few ulps when the flux dips negative.  Reproduce the
    # chain (same order of operations -> bit exact), and report the
    # unclamped integral as a diagnostic.
    nbad = 0
    maxpre = maxsno = 0.0
    maxunc = 0.0
    for key in sorted(scal_pairs):
        pr = profile(key)
        r = scal_steps[key[0]][key[1] - 1]
        prec_evap = f(r[ci["prec_evap"]])
        snow_evap = f(r[ci["snow_evap"]])
        for name, got, prod, acc in (("prec_evap", prec_evap, "ntprprd", "prec"),
                                     ("snow_evap", snow_evap, "ntsnprd", "snow")):
            flx = 0.0
            unc = 0.0
            for k in range(1, pver + 1):
                term = f(pr[k][pi[prod]]) * f(pr[k][pi["pdel"]]) / grav
                flx = max(flx + term, 0.0)          # the FORTRAN chain
                unc += term                          # unclamped integral
            w = flx / 1000.0
            scale = max(abs(got), abs(w))
            res = abs(got - w) / scale if scale > 0 else (1.0 if got != w else 0.0)
            if acc == "prec":
                maxpre = max(maxpre, res)
            else:
                maxsno = max(maxsno, res)
            if res > RTOL:
                nbad += 1
                if nbad <= 5:
                    print(f"  {name} {key}: {got:.6e} vs chain {w:.6e}")
            if abs(unc / 1000.0 - got) > maxunc:
                maxunc = abs(unc / 1000.0 - got)
    assert nbad == 0, f"precip/snow sums: {nbad} bad columns"
    print(f"precip sums OK: prec_evap/snow_evap = clamped column chain of "
          f"ntprprd/ntsnprd*pdel/(grav*1000) "
          f"({ncol_pairs} cols, max rel resid prec {maxpre:.2e}, "
          f"snow {maxsno:.2e}, max |unclamped - dump| {maxunc:.2e})")

    # ---- 6. interface flux increments ----------------------------------
    # Same law as check 5, per level: flx_hi = max(flx_lo + prod*pdel/grav, 0).
    nbad = 0
    maxfi = 0.0
    nzero = nbot = 0
    for key, rows in prof_pairs.items():
        pr = profile(key)
        r = scal_steps[key[0]][key[1] - 1]
        prec_evap = f(r[ci["prec_evap"]])
        snow_evap = f(r[ci["snow_evap"]])
        for k in range(1, pver + 1):
            rr = pr[k]
            for prod, lo, hi in (("ntprprd", "flxprec_lo", "flxprec_hi"),
                                 ("ntsnprd", "flxsnow_lo", "flxsnow_hi")):
                fl_lo = f(rr[pi[lo]])
                fl_hi = f(rr[pi[hi]])
                prod_pdel = f(rr[pi[prod]]) * f(rr[pi["pdel"]]) / grav
                want = max(fl_lo + prod_pdel, 0.0)   # the FORTRAN chain
                scale = max(abs(fl_lo), abs(fl_hi), abs(prod_pdel))
                res = abs(fl_hi - want) / scale if scale > 0 else (1.0 if fl_hi else 0.0)
                maxfi = max(maxfi, res)
                if res > RTOL:
                    nbad += 1
            if k == 1:
                if f(pr[1][pi["flxprec_lo"]]) == 0.0:
                    nzero += 1
            if k == pver:
                if abs(f(pr[pver][pi["flxprec_hi"]]) - prec_evap * 1000.0) <= \
                        RTOL * max(abs(prec_evap * 1000.0), 1.0):
                    nbot += 1
                if abs(f(pr[pver][pi["flxsnow_hi"]]) - snow_evap * 1000.0) <= \
                        RTOL * max(abs(snow_evap * 1000.0), 1.0):
                    nbot += 1
    assert nbad == 0, f"flux increments: {nbad} bad rows"
    assert nzero == ncol_pairs, f"flxprec_lo(lev=1) != 0 in {ncol_pairs - nzero} columns"
    assert nbot == 2 * ncol_pairs, "flx_hi(lev=pver) != prec/snow*1000"
    print(f"flux increments OK: flx_hi-flx_lo = ntprprd/ntsnprd*pdel/grav "
          f"({ncol_pairs*pver*2} rows, max rel resid {maxfi:.2e}); "
          f"top flux 0 and bottom flux = prec/snow*1000 checked")

    # ---- 7. mass-flux bookkeeping --------------------------------------
    # mc = mu + md per row (zm_convr.F90:2311, then all three scaled by
    # mb at 702-704).  (mcon is mentioned in the dump `# description`
    # prose but is not a dumped column; it is mc*100/grav.)
    nbad = 0
    maxmc = 0.0
    for key, rows in prof_pairs.items():
        for r in rows:
            mu, md, mc = (f(r[pi[c]]) for c in ("mu", "md", "mc"))
            scale = max(abs(mu), abs(md), abs(mc))
            res = abs(mc - (mu + md)) / scale if scale > 0 else (1.0 if mc else 0.0)
            maxmc = max(maxmc, res)
            if res > RTOL:
                nbad += 1
    assert nbad == 0, f"mass-flux bookkeeping: {nbad} bad rows"
    print(f"mass-flux bookkeeping OK: mc = mu+md ({len(prof)} rows, "
          f"max rel resid {maxmc:.2e})")

    # ---- 8. working arrays p/z/s ----------------------------------------
    nbad = 0
    worst = {}
    for key, rows in prof_pairs.items():
        r = scal_steps[key[0]][key[1] - 1]
        phis = f(r[ci["phis"]])
        for rr in rows:
            t = f(rr[pi["t"]])
            pmid = f(rr[pi["pmid"]])
            p = f(rr[pi["p"]])
            for a, b, name in ((p, 0.01 * pmid, "p"),
                               (f(rr[pi["pf_lo"]]), 0.01 * f(rr[pi["pint_lo"]]), "pf_lo"),
                               (f(rr[pi["pf_hi"]]), 0.01 * f(rr[pi["pint_hi"]]), "pf_hi")):
                scale = max(abs(a), abs(b))
                res = abs(a - b) / scale if scale > 0 else (1.0 if a != b else 0.0)
                worst[name] = max(worst.get(name, 0.0), res)
                if res > RTOL:
                    nbad += 1
            z = f(rr[pi["z"]])
            want_z = f(rr[pi["zm"]]) + phis / grav
            scale = max(abs(z), abs(want_z))
            res = abs(z - want_z) / scale if scale > 0 else (1.0 if z else 0.0)
            worst["z"] = max(worst.get("z", 0.0), res)
            if res > RTOL:
                nbad += 1
            s = f(rr[pi["s"]])
            want_s = t + (grav / cpres) * z
            scale = max(abs(s), abs(want_s))
            res = abs(s - want_s) / scale if scale > 0 else (1.0 if s else 0.0)
            worst["s"] = max(worst.get("s", 0.0), res)
            if res > RTOL:
                nbad += 1
    assert nbad == 0, f"working arrays: {nbad} bad rows"
    print(f"working arrays OK: p=0.01*pmid, pf=0.01*pint, z=zm+phis/grav, "
          f"s=t+(grav/cpres)*z ({len(prof)} rows; max rel resid "
          + ", ".join(f"{k} {v:.1e}" for k, v in worst.items()) + ")")

    # ---- 9. trigger / index sanity -------------------------------------
    nbad = 0
    mincape = 1e300
    for key in sorted(scal_pairs):
        r = scal_steps[key[0]][key[1] - 1]
        cape = f(r[ci["cape"]])
        mincape = min(mincape, cape)
        if not (cape > capelmt):
            nbad += 1
        jt = int(f(r[ci["jt"]]))
        maxg = int(f(r[ci["maxg"]]))
        if not (1 <= jt <= maxg <= pver):
            nbad += 1
        if int(f(r[ci["jctop"]])) != jt or int(f(r[ci["jcbot"]])) != maxg:
            nbad += 1
        for c in ("lcl", "lel", "lon0", "maxi", "pblt", "jlcl", "j0", "jd"):
            v = int(f(r[ci[c]]))
            if not (0 <= v <= pver):
                nbad += 1
    assert nbad == 0, f"index sanity: {nbad} bad rows"
    print(f"trigger/index sanity OK: cape > capelmt on all {ncol_pairs} rows "
          f"(min cape {mincape:.3g} J/kg); 1<=jt<=maxg<=pver; "
          f"jctop==jt, jcbot==maxg; level indices in bounds")

    print("\nALL ZM DUMP CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
