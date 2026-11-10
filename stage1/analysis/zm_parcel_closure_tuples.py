#!/usr/bin/env python3
"""Extract the numeric test tuples for cam_zm_parcel.esm and cam_zm_closure.esm
from the stage-1 ZM dumps, and verify every law transcription against the dump.

The .esm models (stage2/cam_zm_parcel.esm, stage2/cam_zm_closure.esm,
stage2/cam_zm_templates.esm) are hand-authored; the numeric arrays of their
`tests` blocks are generated here, per the stage-2 workflow (scripts may
extract numeric test tuples from the instrumented-model dumps into an
existing hand-authored tests block).  Every law replica below is the
zm_convr.F90 law in double precision; the verification section prints the
worst relative error against the dumped values and exits nonzero on any
mismatch beyond 1e-12.

Output: writes the tuples used by the test blocks as JSON under the output
directory (default scratch/zm_tuples of the repo where the script runs):
  parcel/   per-step full t/q/z/p/tp/qstp/qtot profiles (dump columns t
            [K], q [kg/kg], z [m], p [mb] = the hPa working pressure,
            tp/qstp the parcel temperature and saturation specific humidity
            of the ientropy solve; qtot from cam_zm_parcel.csv, the parcel
            total water new_q, 0 below msg+1 = 6), the launch level maxi,
            and the
            expected values at the pinned levels: hmn/hsat at lev
            {1,15,25,32} (cam_zm.csv dump columns), qst at lev {25,32}
            (DERIVED: the qsat_hPa replica of zm_convr.F90:2634-2649,
            cross-pinned by hsat), tpv at lev {25,32} and buoy at lev
            {15,18,25} (cam_zm_parcel.csv dump columns), and cape (the
            dumped cam_zm_scal.csv cape, == capeten_1 of
            cam_zm_parcel_scal.csv under num_cin = 1); plus the five
            (step, lev) entropy/virtual-temperature states (DERIVED: the
            entropy replica of F90:1489-1513 and tv = t (1+1.608 q)/(1+q),
            F90:1088).  qstp_launch = q(maxi) (parcel_dilute, F90:1414) is
            NOT emitted: it is a definition (index(q, maxi)), unfailable as
            a pin, and its Fortran value qstp(maxi) is a dumped parameter
            input, not a model output.

  closure/  per-step scalars (cape, tl, lel, lcl, jt, mx, mb, dsubcld from
            cam_zm_scal.csv) and the 18 closure input profiles at levels
            14..32 (zero outside; the closure reads only k >= lel = 14..15
            through mx = 32 in this lane): qg/tg/pg/sg/tpg/qstpg/qu/su/qd/sd/
            qhat/shat/ql/pdel/zf_lo/zf_hi (dump columns; zg and qs are NOT
            emitted -- the closure body never reads them) and the normalized
            fluxes mu/md/mc/du (dump value / mb, the form the closure
            receives: mu/md/mc dimensionless, du 1/hPa; du is read at
            jt < k < mx = 16..31, where du(16) is nonzero on every row and
            du(17..31) = 0).  pdel is dumped in Pa and the model converts
            with dp = 0.01 hPa/Pa x pdel (the README-cam-zm working identity
            dp = 0.01*pdel).

Sources: stage1/dumps/cam/cam_zm.csv and cam_zm_scal.csv (480 profile rows,
15 scalar rows, steps 18-23 and 40-48), stage1/dumps/cam/cam_zm_parcel.csv
(480 rows: the parcel-ascent internals tmix/qsmix/qtmix/smix/xsh2o/qtot/
ds_xsh2o/ds_freeze/qtp/sp/mp/tpv/buoy) and cam_zm_parcel_scal.csv (15 rows:
pl/plge600/capeten_1..5/lelten_1..5), and the instrumented cam6_4_212 fork
components/cam/src/atmos_phys/schemes/zhang_mcfarlane/zm_convr.F90.
"""
import csv
import json
import math
import os
import sys

DUMP = os.environ.get(
    "ZM_DUMP_DIR",
    "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cam")
OUT = os.environ.get("ZM_TUPLE_DIR", os.path.join("scratch", "zm_tuples"))


def load(fn):
    lines = open(fn).readlines()
    cols = lines[1][len("# columns = "):].strip().split(",")
    rows = [r for r in csv.reader(lines[2:]) if r]
    return cols, rows


def invert_entropy(s, p_mb, qtot, rb_entropy):
    """Brent-free replica of the ientropy inversion (F90:1517-1632): root of
    rb_entropy(t, p_mb, qtot) = s.  The Fortran Brent solve converges to
    tol1 ~= 5.2e-4 K; bisection to 1e-9 K reproduces the root well within
    that bracket, so all reconstructed temperatures match the Fortran ones
    to the Brent tolerance."""
    lo, hi = 150.0, 400.0
    flo = rb_entropy(lo, p_mb, qtot) - s
    fhi = rb_entropy(hi, p_mb, qtot) - s
    assert flo * fhi < 0.0, f"invert_entropy: no bracket for s={s}, p={p_mb}, qtot={qtot}"
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = rb_entropy(mid, p_mb, qtot) - s
        if fm * flo <= 0.0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-9:
            break
    return 0.5 * (lo + hi)


def main():
    cols_p, rows_p = load(f"{DUMP}/cam_zm.csv")
    cols_s, rows_s = load(f"{DUMP}/cam_zm_scal.csv")
    cols_a, rows_a = load(f"{DUMP}/cam_zm_parcel.csv")
    cols_b, rows_b = load(f"{DUMP}/cam_zm_parcel_scal.csv")
    pi = {c: i for i, c in enumerate(cols_p)}
    si = {c: i for i, c in enumerate(cols_s)}
    ai = {c: i for i, c in enumerate(cols_a)}
    bi = {c: i for i, c in enumerate(cols_b)}

    C = {}
    for c in ("cpres", "rl", "rgas", "eps1", "tfreez", "grav", "cpliq", "cpwv",
              "rh2o", "latice", "capelmt", "tiedke_add", "tau"):
        C[c] = float(rows_s[0][si[c]])

    prof = {}
    for r in rows_p:
        prof.setdefault((int(r[0]), int(r[3])), {})[int(r[4])] = r
    scal = {(int(r[0]), int(r[3])): r for r in rows_s}
    pdump = {}
    for r in rows_a:
        pdump.setdefault((int(r[0]), int(r[3])), {})[int(r[4])] = r
    pscal = {(int(r[0]), int(r[3])): r for r in rows_b}
    steps = sorted({int(r[0]) for r in rows_s})

    # ---------------- law replicas (zm_convr.F90, double precision) --------

    def goff_gratch_water(t):
        tboil = 373.16
        return 10.0 ** (-7.90298 * (tboil / t - 1.0) + 5.02808 * math.log10(tboil / t)
                        - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / tboil)) - 1.0)
                        + 8.1328e-3 * (10.0 ** (-3.49149 * (tboil / t - 1.0)) - 1.0)
                        + math.log10(1013.246)) * 100.0

    def qsat_hPa(t, p_mb):
        """zm_convr.F90:2634-2649 (qsat_water at p*100 Pa, es back in hPa)."""
        es_pa = goff_gratch_water(t)
        p_pa = p_mb * 100.0
        es_pa = min(es_pa, p_pa)
        if (p_pa - es_pa) <= 0.0:
            qs = 1.0
        else:
            qs = 0.6219705862045155 * es_pa / (p_pa - 0.3780294137954845 * es_pa)
        return qs

    def rb_entropy(TK, p_mb, qtot):
        """zm_convr.F90:1489-1513."""
        L = C["rl"] - (C["cpliq"] - C["cpwv"]) * (TK - C["tfreez"])
        qst = qsat_hPa(TK, p_mb)
        qv = min(qtot, qst)
        e = qv * p_mb / (C["eps1"] + qv)
        pref = 1000.0
        return ((C["cpres"] + qtot * C["cpliq"]) * math.log(TK / C["tfreez"])
                - C["rgas"] * math.log((p_mb - e) / pref)
                + L * qv / TK - qv * C["rh2o"] * math.log(qv / qst))

    def env_tv(t, q):
        return t * (1.0 + 1.608 * q) / (1.0 + q)

    def closure(pr, sr):
        """Replicate closure (zm_convr.F90:2318-2506) + dsubcld (597-603) at
        the dumped state, on the NORMALIZED (per unit mb) fluxes."""
        def A(k, name):
            return float(pr[k][pi[name]])

        mx = int(sr[si["maxg"]])
        jt = int(sr[si["jt"]])
        lel = int(sr[si["lel"]])
        lcl = int(sr[si["lcl"]])
        cape = float(sr[si["cape"]])
        tl = float(sr[si["tl"]])
        pver = 32
        mb_dump = float(sr[si["mb"]])

        dsubcld = 0.0
        for k in range(6, pver + 1):
            if k >= mx:
                dsubcld += 0.01 * A(k, "pdel")

        def flx(name, k):
            return A(k, name) / mb_dump

        q_mx = A(mx, "qg"); t_mx = A(mx, "tg"); p_mx = A(mx, "pg")
        eb = p_mx * q_mx / (C["eps1"] + q_mx)
        dtbdt = (1.0 / dsubcld) * (flx("mu", mx) * (A(mx, "shat") - A(mx, "su"))
                                   + flx("md", mx) * (A(mx, "shat") - A(mx, "sd")))
        dqbdt = (1.0 / dsubcld) * (flx("mu", mx) * (A(mx, "qhat") - A(mx, "qu"))
                                   + flx("md", mx) * (A(mx, "qhat") - A(mx, "qd")))
        debdt = C["eps1"] * p_mx / (C["eps1"] + q_mx) ** 2 * dqbdt
        dtldt = (-2840.0 * (3.5 / t_mx * dtbdt - debdt / eb)
                 / (3.5 * math.log(t_mx) - math.log(eb) - 4.805) ** 2)

        dtmdt = [0.0] * (pver + 1)
        dqmdt = [0.0] * (pver + 1)
        for k in range(6, pver):
            if k == jt:
                dtmdt[k] = (flx("mu", k + 1) * (A(k + 1, "su") - A(k + 1, "shat")
                            - C["rl"] / C["cpres"] * A(k + 1, "ql"))
                            + flx("md", k + 1) * (A(k + 1, "sd") - A(k + 1, "shat"))) \
                    / (0.01 * A(k, "pdel"))
                dqmdt[k] = (flx("mu", k + 1) * (A(k + 1, "qu") - A(k + 1, "qhat")
                            + A(k + 1, "ql"))
                            + flx("md", k + 1) * (A(k + 1, "qd") - A(k + 1, "qhat"))) \
                    / (0.01 * A(k, "pdel"))
            if k > jt and k < mx:
                dp = 0.01 * A(k, "pdel")
                dtmdt[k] = ((flx("mc", k) * (A(k, "shat") - A(k, "sg"))
                             + flx("mc", k + 1) * (A(k, "sg") - A(k + 1, "shat"))) / dp
                            - C["rl"] / C["cpres"] * flx("du", k) * A(k + 1, "ql"))
                dqmdt[k] = ((flx("mu", k + 1) * (A(k + 1, "qu") - A(k + 1, "qhat")
                             + C["cpres"] / C["rl"] * (A(k + 1, "su") - A(k, "sg")))
                             - flx("mu", k) * (A(k, "qu") - A(k, "qhat")
                             + C["cpres"] / C["rl"] * (A(k, "su") - A(k, "sg")))
                             + flx("md", k + 1) * (A(k + 1, "qd") - A(k + 1, "qhat")
                             + C["cpres"] / C["rl"] * (A(k + 1, "sd") - A(k, "sg")))
                             - flx("md", k) * (A(k, "qd") - A(k, "qhat")
                             + C["cpres"] / C["rl"] * (A(k, "sd") - A(k, "sg")))) / dp
                            + flx("du", k) * A(k + 1, "ql"))

        dboydt = [0.0] * (pver + 1)
        for k in range(6, pver + 1):
            if lel <= k <= lcl:
                thetavp = A(k, "tpg") * (1000.0 / A(k, "pg")) ** (C["rgas"] / C["cpres"]) \
                    * (1.0 + 1.608 * A(k, "qstpg") - q_mx)
                thetavm = A(k, "tg") * (1000.0 / A(k, "pg")) ** (C["rgas"] / C["cpres"]) \
                    * (1.0 + 0.608 * A(k, "qg"))
                dqsdtp = A(k, "qstpg") * (1.0 + A(k, "qstpg") / C["eps1"]) * C["eps1"] \
                    * C["rl"] / (C["rgas"] * A(k, "tpg") ** 2)
                dtpdt = A(k, "tpg") / (1.0 + C["rl"] / C["cpres"]
                                       * (dqsdtp - A(k, "qstpg") / A(k, "tpg"))) \
                    * (dtbdt / t_mx + C["rl"] / C["cpres"]
                       * (dqbdt / tl - q_mx / tl ** 2 * dtldt))
                dboydt[k] = ((dtpdt / A(k, "tpg")
                              + 1.0 / (1.0 + 1.608 * A(k, "qstpg") - q_mx)
                              * (1.608 * dqsdtp * dtpdt - dqbdt))
                             - (dtmdt[k] / A(k, "tg")
                                + 0.608 / (1.0 + 0.608 * A(k, "qg")) * dqmdt[k])) \
                    * C["grav"] * thetavp / thetavm
            elif lcl < k < mx:
                thetavp = A(k, "tpg") * (1000.0 / A(k, "pg")) ** (C["rgas"] / C["cpres"]) \
                    * (1.0 + 0.608 * q_mx)
                thetavm = A(k, "tg") * (1000.0 / A(k, "pg")) ** (C["rgas"] / C["cpres"]) \
                    * (1.0 + 0.608 * A(k, "qg"))
                dboydt[k] = (dtbdt / t_mx + 0.608 / (1.0 + 0.608 * q_mx) * dqbdt
                             - dtmdt[k] / A(k, "tg")
                             - 0.608 / (1.0 + 0.608 * A(k, "qg")) * dqmdt[k]) \
                    * C["grav"] * thetavp / thetavm

        dadt = 0.0
        for k in range(lel, mx):
            dadt += dboydt[k] * (A(k, "zf_lo") - A(k, "zf_hi"))
        dltaa = -1.0 * (cape - C["capelmt"])
        mb = 0.0 if dadt == 0.0 else max(dltaa / C["tau"] / dadt, 0.0)
        return dict(dsubcld=dsubcld, mb=mb)

    # ---------------- verification (all laws vs the dump) ------------------

    pver = 32
    worst_mb = worst_ds = worst_hmn = worst_hsat = 0.0
    worst_tpv = worst_buoy = worst_cape = worst_pl = worst_tlaunch = 0.0
    for st in steps:
        pr = prof[(st, 1)]
        sr = scal[(st, 1)]
        res = closure(pr, sr)
        mb_dump = float(sr[si["mb"]])
        ds_dump = float(sr[si["dsubcld"]])
        worst_mb = max(worst_mb, abs(res["mb"] - mb_dump) / max(abs(res["mb"]), abs(mb_dump)))
        worst_ds = max(worst_ds, abs(res["dsubcld"] - ds_dump) / max(abs(ds_dump), 1.0))
        mx = int(sr[si["maxi"]])
        for k in (1, 15, 25, 32):
            r = pr[k]
            hmn = C["cpres"] * float(r[pi["t"]]) + C["grav"] * float(r[pi["z"]]) \
                + C["rl"] * float(r[pi["q"]])
            qst = qsat_hPa(float(r[pi["t"]]), float(r[pi["p"]]))
            hsat = C["cpres"] * float(r[pi["t"]]) + C["grav"] * float(r[pi["z"]]) \
                + C["rl"] * qst
            worst_hmn = max(worst_hmn, abs(float(r[pi["hmn"]]) - hmn)
                            / max(abs(float(r[pi["hmn"]])), abs(hmn)))
            worst_hsat = max(worst_hsat, abs(float(r[pi["hsat"]]) - hsat)
                             / max(abs(float(r[pi["hsat"]])), abs(hsat)))
        assert float(pr[mx][pi["qstp"]]) == float(pr[mx][pi["q"]]), \
            f"step {st}: qstp(mx) != q(mx)"
        # ---- parcel-ascent laws (cam_zm_parcel.csv / cam_zm_parcel_scal.csv)
        pa = pdump[(st, 1)]
        brow = pscal[(st, 1)]
        tpert = float(sr[si["tpert"]])
        # launch level (F90:1347-1353): qtmix = qtot = q, tmix = ientropy root
        assert float(pa[mx][ai["qtmix"]]) == float(pr[mx][pi["q"]]), \
            f"step {st}: qtmix(mx) != q(mx)"
        assert float(pa[mx][ai["qtot"]]) == float(pr[mx][pi["q"]]), \
            f"step {st}: qtot(mx) != q(mx)"
        worst_tlaunch = max(worst_tlaunch, abs(
            float(pa[mx][ai["tmix"]]) - invert_entropy(
                float(pa[mx][ai["smix"]]), float(pr[mx][pi["p"]]),
                float(pa[mx][ai["qtmix"]]), rb_entropy)))
        # tpv/buoy laws at all levels (F90:1415 launch, 1476 ascent, 1088-1089;
        # tpv = tv and buoy = 0 below msg+1 = 6 and above mx)
        for k in range(1, pver + 1):
            a = pa[k]
            tpv = float(a[ai["tpv"]]); buoy = float(a[ai["buoy"]])
            tv = env_tv(float(pr[k][pi["t"]]), float(pr[k][pi["q"]]))
            if 6 <= k <= mx:
                tpv_ref = (float(pr[k][pi["tp"]]) + tpert) \
                    * (1.0 + 1.608 * float(pr[k][pi["qstp"]])) \
                    / (1.0 + float(a[ai["qtot"]]))
                worst_tpv = max(worst_tpv, abs(tpv - tpv_ref))
                worst_buoy = max(worst_buoy, abs(buoy - (tpv - tv + C["tiedke_add"])))
            else:
                assert buoy == 0.0, f"step {st} lev {k}: buoy = {buoy} != 0"
                worst_tpv = max(worst_tpv, abs(tpv - tv))
        # CAPE integral over the dumped buoy (F90:1169-1178) and the
        # num_cin = 1 collapse (F90:1180-1194)
        lelten1 = int(brow[bi["lelten_1"]])
        cap = 0.0
        for k in range(lelten1 + 1, mx + 1):
            cap += C["rgas"] * float(pa[k][ai["buoy"]]) \
                * math.log(float(pr[k][pi["pf_hi"]]) / float(pr[k][pi["pf_lo"]]))
        worst_cape = max(worst_cape, abs(cap - float(brow[bi["capeten_1"]])))
        assert float(sr[si["cape"]]) == float(brow[bi["capeten_1"]]), \
            f"step {st}: cape != capeten_1"
        assert int(sr[si["lel"]]) == lelten1, f"step {st}: lel != lelten_1"
        for n in range(2, 6):
            assert float(brow[bi[f"capeten_{n}"]]) == 0.0, \
                f"step {st}: capeten_{n} != 0"
            assert int(brow[bi[f"lelten_{n}"]]) == 0, \
                f"step {st}: lelten_{n} != 0"
        # pl: qtmix-qsmix zero-crossing interpolation (F90:1356-1361) with the
        # PRE-rainout qsmix at the lcl (entrainment-phase ientropy root; the
        # dumped qsmix at the lcl is post-loop, though xsh2o(lcl) = 0)
        lcl = int(sr[si["lcl"]])
        a_up, p_up = pa[lcl], pr[lcl]
        a_dn, p_dn = pa[lcl + 1], pr[lcl + 1]
        qspre = qsat_hPa(invert_entropy(float(a_up[ai["smix"]]),
                                        float(p_up[pi["p"]]),
                                        float(a_up[ai["qtmix"]]), rb_entropy),
                         float(p_up[pi["p"]]))
        qxsk = float(a_up[ai["qtmix"]]) - qspre
        qxskp1 = float(a_dn[ai["qtmix"]]) - float(a_dn[ai["qsmix"]])
        dqxsdp = (qxsk - qxskp1) / (float(p_up[pi["p"]]) - float(p_dn[pi["p"]]))
        pl_ref = float(p_dn[pi["p"]]) - qxskp1 / dqxsdp
        pl_dump = float(brow[bi["pl"]])
        worst_pl = max(worst_pl, abs(pl_dump - pl_ref))
        assert (brow[bi["plge600"]] == "T") == (pl_dump >= 600.0), \
            f"step {st}: plge600 mismatch"
    print(f"verification: worst rel |mb| {worst_mb:.2e}, |dsubcld| {worst_ds:.2e}, "
          f"|hmn| {worst_hmn:.2e}, |hsat| {worst_hsat:.2e}")
    print(f"verification: worst |tpv| {worst_tpv:.2e}, |buoy| {worst_buoy:.2e}, "
          f"|cape| {worst_cape:.2e}, |pl| {worst_pl:.2e} mb, "
          f"|tmix-launch root| {worst_tlaunch:.2e} K")
    if max(worst_mb, worst_ds, worst_hmn, worst_hsat) > 1e-12:
        print("FAIL: verification mismatch above 1e-12")
        sys.exit(1)
    if max(worst_tpv, worst_buoy) > 1e-9 or worst_cape > 1e-6 \
            or worst_pl > 0.02 or worst_tlaunch > 1e-3:
        print("FAIL: parcel-law verification mismatch (tpv/buoy 1e-9, cape 1e-6, "
              "pl 0.02 mb, launch root 1e-3)")
        sys.exit(1)

    # ---------------- tuple extraction -------------------------------------

    pver = 32
    parcel = {}
    closure_tuples = {}
    for st in steps:
        pr = prof[(st, 1)]
        sr = scal[(st, 1)]
        mx = int(sr[si["maxi"]])
        arr = {c: [float(pr[k][pi[c]]) for k in range(1, pver + 1)]
               for c in ("t", "q", "z", "p", "tp", "qstp", "pf_hi", "pf_lo")}
        pa = pdump[(st, 1)]
        arr["qtot"] = [float(pa[k][ai["qtot"]]) for k in range(1, pver + 1)]
        asserts = {}
        for k in (1, 15, 25, 32):
            asserts[f"hmn{k}"] = float(pr[k][pi["hmn"]])
            asserts[f"hsat{k}"] = float(pr[k][pi["hsat"]])
        for k in (25, 32):
            asserts[f"qst{k}"] = qsat_hPa(float(pr[k][pi["t"]]), float(pr[k][pi["p"]]))
        for k in (25, 32):
            asserts[f"tpv{k}"] = float(pa[k][ai["tpv"]])
        for k in (15, 18, 25):
            asserts[f"buoy{k}"] = float(pa[k][ai["buoy"]])
        asserts["cape"] = float(sr[si["cape"]])
        parcel[str(st)] = dict(overrides=arr, maxi=mx,
                               lel=int(sr[si["lel"]]), asserts=asserts)

        mb_dump = float(sr[si["mb"]])
        flux_vars = ("mu", "md", "mc", "du")
        prof_vars = ("qg", "tg", "pg", "sg", "tpg", "qstpg", "qu", "su",
                     "qd", "sd", "qhat", "shat", "ql", "pdel", "zf_lo", "zf_hi")
        norm = {}
        for v in prof_vars:
            norm[v] = [0.0] * 13 + [float(pr[k][pi[v]]) for k in range(14, pver + 1)]
        for v in flux_vars:
            norm[v] = [0.0] * 13 + [float(pr[k][pi[v]]) / mb_dump
                                    for k in range(14, pver + 1)]
        closure_tuples[str(st)] = dict(
            scalars=dict(cape=float(sr[si["cape"]]), tl=float(sr[si["tl"]]),
                         lel=int(sr[si["lel"]]), lcl=int(sr[si["lcl"]]),
                         jt=int(sr[si["jt"]]), mx=int(sr[si["maxg"]])),
            norm=norm,
            asserts=dict(dsubcld=float(sr[si["dsubcld"]]), mb=mb_dump))

    ent_tv_states = [(18, 32), (18, 25), (18, 15), (40, 26), (48, 32)]
    ent_tv = []
    for st, k in ent_tv_states:
        r = prof[(st, 1)][k]
        t = float(r[pi["t"]]); p = float(r[pi["p"]]); q = float(r[pi["q"]])
        ent_tv.append(dict(step=st, lev=k, t=t, p=p, q=q,
                           ent=rb_entropy(t, p, q), tv=env_tv(t, q)))

    constants = {c: C[c] for c in ("cpres", "rl", "rgas", "eps1", "tfreez",
                                   "grav", "cpliq", "cpwv", "rh2o",
                                   "capelmt", "tiedke_add", "tau")}
    out = dict(parcel=dict(per_step=parcel, ent_tv=ent_tv),
               closure=dict(per_step=closure_tuples, constants=constants))
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/tuples.json", "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {OUT}/tuples.json")
    print(f"parcel: {len(parcel)} steps x 16 asserts + {len(ent_tv)} entropy/tv states; "
          f"closure: {len(closure_tuples)} steps x (dsubcld, mb)")


if __name__ == "__main__":
    main()
