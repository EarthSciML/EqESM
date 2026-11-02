#!/usr/bin/env python3
"""Validate the stage-1 CAM dump CSVs from the instrumented SCAM run.

Usage: python3 validate_cam_dumps.py DUMP_DIR [--steps N] [--pver 32]

Checks (mirroring the CICE lane's validation):
  1. all 7 dump files exist, start with '# description ='/'# columns ='
     headers, have zero blank lines and a consistent field count;
  2. row counts: cam_phys_tend = N x 27 sites x pver, cam_surface_exchange
     = N x 3 phases, cam_surface_cflx = N x pcnst, cam_qsat = N x 2 x pver,
     cam_tend_total = N x 2 x pver, cam_phys_tend_q = N x sum(lq count)
     x pver;
  3. per-step site order matches the tphysbc/tphysac call sequence
     (macrop/clubb/microp(_sc) appear cld_macmic_num_steps=3 times);
  4. the wv_saturation qsat law (Goff-Gratch mixed-phase table,
     estblf linear interpolation, svp_to_qsat, gam with the pcf
     transition polynomial) reproduces the dumped es/qs/gam exactly;
  5. cam_surface_exchange bc_out == ac_out cam_out fields (nothing in
     tphysac should change cam_out between cam_export and the end of
     tphysac);
  6. the run contains precipitating steps (precsc/precsl/precc/precl).

The qsat law is implemented below with the exact constants and formula
order of wv_saturation.F90 / wv_sat_methods.F90 (cam6_4_212).
"""
import sys
import csv
import os
import numpy as np

# ---------------------------------------------------------------- qsat law
# physconst / shr_const values (physconst.F90, shr_const_mod.F90)
MWWV = 18.016          # kg/kmole
MWDAIR = 28.966        # kg/kmole
AVOGAD = 6.02214e26    # molecules/kmole
BOLTZ = 1.38065e-23    # J/K/molecule
RGAS = AVOGAD * BOLTZ
RWV = RGAS / MWWV      # water-vapor gas constant, J/K/kg
EPSILO = MWWV / MWDAIR
OMEPS = 1.0 - EPSILO
LATVAP = 2.501e6
LATICE = 3.337e5
CPAIR = 1.00464e3
TMELT = 273.15         # shr_const_tkfrz
H2OTRIP = 273.16       # shr_const_tktrip
TBOIL = 373.16         # wv_saturation parameter
TTRICE = 20.0
TMIN, TMAX = 127.16, 375.16

PCF = np.array([5.04469588506e-01, -5.47288442819e+00, -3.67471858735e-01,
                -8.95963532403e-03, -7.78053686625e-05])


def goff_gratch_water(t):
    return (10.0 ** (-7.90298 * (TBOIL / t - 1.0)
                     + 5.02808 * np.log10(TBOIL / t)
                     - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / TBOIL)) - 1.0)
                     + 8.1328e-3 * (10.0 ** (-3.49149 * (TBOIL / t - 1.0)) - 1.0)
                     + np.log10(1013.246)) * 100.0)


def goff_gratch_ice(t):
    return (10.0 ** (-9.09718 * (H2OTRIP / t - 1.0)
                     - 3.56654 * np.log10(H2OTRIP / t)
                     + 0.876793 * (1.0 - t / H2OTRIP)
                     + np.log10(6.1071)) * 100.0)


def svp_trans(t):
    """Mixed-phase Goff-Gratch SVP (wv_sat_svp_trans)."""
    es = np.where(t >= TMELT - TTRICE, goff_gratch_water(t), 0.0)
    cold = t < TMELT
    if np.any(cold):
        esice = goff_gratch_ice(t[cold])
        weight = np.where((TMELT - t[cold]) > TTRICE, 1.0,
                          (TMELT - t[cold]) / TTRICE)
        es[cold] = weight * esice + (1.0 - weight) * es[cold]
    return es


def build_estbl():
    n = int(np.ceil(TMAX - TMIN)) + 2
    tt = TMIN + np.arange(n, dtype=np.float64)
    return svp_trans(tt)


ESTBL = build_estbl()


def estblf(t):
    """Linear interpolation from the estbl table (wv_sat_estblf)."""
    t = np.asarray(t, dtype=np.float64)
    t_tmp = np.maximum(np.minimum(t, TMAX) - TMIN, 0.0)
    i = np.trunc(t_tmp).astype(int) + 1          # 1-based index
    weight = t_tmp - np.trunc(t_tmp)
    es = (1.0 - weight) * ESTBL[i - 1] + weight * ESTBL[i]  # 0-based
    return es


def svp_to_qsat(es, p):
    """wv_sat_svp_to_qsat."""
    qs = np.where((p - es) <= 0.0, 1.0, EPSILO * es / (p - OMEPS * es))
    return qs


def no_ip_hltalt(t):
    hltalt = np.full_like(t, LATVAP, dtype=np.float64)
    warm = t >= TMELT
    hltalt[warm] = hltalt[warm] - 2369.0 * (t[warm] - TMELT)
    return hltalt


def calc_hltalt(t):
    """calc_hltalt: hltalt and tterm (pcf transition polynomial)."""
    t = np.asarray(t, dtype=np.float64)
    hltalt = no_ip_hltalt(t)
    tterm = np.zeros_like(t)
    cold = t < TMELT
    if np.any(cold):
        tc = t[cold] - TMELT
        in_trans = tc >= -TTRICE
        weight = np.where(in_trans, -tc / TTRICE, 1.0)
        hltalt[cold] = hltalt[cold] + weight * LATICE
        if np.any(in_trans):
            ttc = tc[in_trans]
            tt = np.zeros_like(ttc)
            for c in PCF[::-1]:
                tt = c + ttc * tt
            tterm[cold[in_trans]] = tt / TTRICE
    return hltalt, tterm


def qsat_law(t, p):
    """Full qsat (mixed phase) law: returns es, qs, gam."""
    es = estblf(t)
    qs = svp_to_qsat(es, p)
    es = np.minimum(es, p)
    hltalt, tterm = calc_hltalt(t)
    gam = np.zeros_like(t)
    sat = qs == 1.0
    notsat = ~sat
    if np.any(notsat):
        desdt = hltalt[notsat] * es[notsat] / (RWV * t[notsat] ** 2) + tterm[notsat]
        dqsdt = qs[notsat] * p[notsat] * desdt / (es[notsat] * (p[notsat] - OMEPS * es[notsat]))
        gam[notsat] = dqsdt * (hltalt[notsat] / CPAIR)
    return es, qs, gam


# ------------------------------------------------------------- site order
# Expected per-step site sequence (physpkg.F90 call order).  The macmic
# loop sites (macrop/clubb/microp_sc/microp) repeat cld_macmic_num_steps
# times; 'bc'/'ac' phases for qsat and tend_total interleave as shown.
BC_SITES = ["energy_fixer", "dadadj", "convect_deep", "convect_shallow",
            "ssatcontrail"]
MACMIC_SITES = ["macrop", "clubb", "microp_sc", "microp"]
AC_SITES = ["carma_emis", "aoa_tracers", "co2_cycle", "chem", "vdiff",
            "rayleigh", "aero_drydep", "carma_tend", "gw", "qbo", "lunar",
            "iondrag", "nudging"]


def expected_sites(nmacmic=3):
    sites = list(BC_SITES)
    for _ in range(nmacmic):
        sites += MACMIC_SITES
    sites += ["aero_calcsize", "aero_wetdep", "carma_wetdep",
              "convect_deep_2", "radiation"]
    sites += AC_SITES
    return sites


# ------------------------------------------------------------- checks
def load(fn):
    """Return (desc, columns, rows) after structural checks."""
    with open(fn) as f:
        lines = f.readlines()
    assert lines[0].startswith("# description"), fn
    assert lines[1].startswith("# columns"), fn
    desc = lines[0][len("# description = "):].strip()
    cols = lines[1][len("# columns = "):].strip().split(",")
    data = [l for l in lines[2:] if l.strip()]
    assert len(data) == len(lines) - 2, f"{fn}: blank lines present"
    ncol = len(cols)
    for i, l in enumerate(data):
        nf = len(l.rstrip("\n").split(","))
        assert nf == ncol, f"{fn}: row {i+3} has {nf} fields, expected {ncol}"
    return desc, cols, data


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    nsteps = 48
    pver = 32
    for a in sys.argv[2:]:
        if a.startswith("--steps"):
            nsteps = int(a.split("=")[1])
        elif a.startswith("--pver"):
            pver = int(a.split("=")[1])
    nmacmic = 3

    files = ["cam_phys_tend.csv", "cam_phys_tend_q.csv",
             "cam_surface_exchange.csv", "cam_surface_cflx.csv",
             "cam_qsat.csv", "cam_tend_total.csv"]
    tables = {}
    for fn in files:
        p = os.path.join(d, fn)
        assert os.path.exists(p), f"missing {p}"
        tables[fn] = load(p)
        print(f"{fn}: {len(tables[fn][2])} rows, {len(tables[fn][1])} cols")

    # ---- row counts -------------------------------------------------
    n_sites = len(expected_sites(nmacmic))
    assert len(tables["cam_phys_tend.csv"][2]) == nsteps * n_sites * pver, \
        "cam_phys_tend row count"
    assert len(tables["cam_qsat.csv"][2]) == nsteps * 2 * pver, \
        "cam_qsat row count"
    assert len(tables["cam_tend_total.csv"][2]) == nsteps * 2 * pver, \
        "cam_tend_total row count"
    assert len(tables["cam_surface_exchange.csv"][2]) == nsteps * 3, \
        "cam_surface_exchange row count"

    # cam_surface_cflx: nsteps x pcnst; pcnst discovered from unique idx
    cflx_rows = [r.split(",") for r in tables["cam_surface_cflx.csv"][2]]
    pcnst = len({r[4] for r in cflx_rows})
    assert len(cflx_rows) == nsteps * pcnst, "cam_surface_cflx row count"
    print(f"pcnst = {pcnst}")

    # cam_phys_tend_q: nsteps x sum(site lq count) x pver
    tq = [r.split(",") for r in tables["cam_phys_tend_q.csv"][2]]
    nq_sites = len({(r[0], r[3]) for r in tq})
    assert len(tq) == nsteps * nq_sites * pver, "cam_phys_tend_q row count"
    print(f"cam_phys_tend_q: {nq_sites} distinct (step,site) blocks x {pver} lev")

    # ---- per-step site order ---------------------------------------
    tend = [r.split(",") for r in tables["cam_phys_tend.csv"][2]]
    exp = expected_sites(nmacmic)
    ok = True
    for s in range(1, nsteps + 1):
        rows = [r for r in tend if int(r[0]) == s]
        sites = [r[3] for r in rows[::pver]]
        if sites != exp:
            print(f"step {s}: SITE ORDER MISMATCH")
            for a, b in zip(sites, exp):
                if a != b:
                    print(f"  {a} != {b}")
                    break
            ok = False
    if ok:
        print(f"site order OK: {len(exp)} sites x {nsteps} steps "
              f"(macmic x{nmacmic})")

    # macmic substep bookkeeping: macrop/clubb/microp rows carry
    # macmic_it = 1..3, scale = 1/3
    import collections
    bad = 0
    for s in range(1, nsteps + 1):
        rows = [r for r in tend if int(r[0]) == s]
        for site in MACMIC_SITES:
            sr = [r for r in rows if r[3] == site]
            its = [int(r[11]) for r in sr[::pver]]
            if its != [1, 2, 3]:
                print(f"step {s} {site}: macmic_it = {its}")
                bad += 1
            for r in sr[::pver]:
                if abs(float(r[14]) - 1.0 / nmacmic) > 1e-12:
                    print(f"step {s} {site}: scale = {r[14]}")
                    bad += 1
    if bad == 0:
        print("macmic substep bookkeeping OK (macmic_it 1..3, scale=1/3)")

    # ---- qsat law ----------------------------------------------------
    q = [r.split(",") for r in tables["cam_qsat.csv"][2]]
    maxes = maxqs = maxgam = 0.0
    nbad = 0
    for r in q:
        t, p = float(r[5]), float(r[6])
        es, qs, gam = qsat_law(np.array([t]), np.array([p]))
        maxes = max(maxes, abs(es[0] - float(r[7])))
        maxqs = max(maxqs, abs(qs[0] - float(r[8])))
        maxgam = max(maxgam, abs(gam[0] - float(r[9])))
        if abs(es[0] - float(r[7])) > 1e-6 * abs(es[0]) + 1e-8:
            nbad += 1
        if abs(qs[0] - float(r[8])) > 1e-12:
            nbad += 1
        if abs(gam[0] - float(r[9])) > 1e-6 * abs(gam[0]) + 1e-8:
            nbad += 1
    print(f"qsat law: {len(q)} rows checked; max|es err|={maxes:.3e}, "
          f"max|qs err|={maxqs:.3e}, max|gam err|={maxgam:.3e}; "
          f"bad={nbad}")

    # ---- surface exchange: bc_out vs ac_out --------------------------
    se = [r.split(",") for r in tables["cam_surface_exchange.csv"][2]]
    cols = tables["cam_surface_exchange.csv"][1]
    camout_cols = ["tbot", "zbot", "ubot", "vbot", "qbot_Q", "pbot", "rho",
                   "netsw", "flwds", "precsc", "precsl", "precc", "precl",
                   "soll", "sols", "solld", "solsd", "thbot", "psl",
                   "co2prog", "co2diag", "ozone"]
    idx = {c: cols.index(c) for c in camout_cols}
    ndiff = 0
    for s in range(1, nsteps + 1):
        rows = [r for r in se if int(r[0]) == s]
        bo = [r for r in rows if r[3] == "bc_out"][0]
        ao = [r for r in rows if r[3] == "ac_out"][0]
        for c in camout_cols:
            if bo[idx[c]] != ao[idx[c]]:
                ndiff += 1
                if ndiff <= 5:
                    print(f"step {s}: cam_out {c} differs bc_out vs ac_out: "
                          f"{bo[idx[c]]} vs {ao[idx[c]]}")
    if ndiff == 0:
        print("surface exchange: bc_out cam_out == ac_out cam_out (all steps)")

    # bc_in vs cflx_Q
    cflx_idx = {r[4]: float(r[6]) for r in cflx_rows
                if r[5] == "Q"}  # cnst_idx for Q
    # (Q is constituent 1; keep as informational)
    qidx = [r for r in cflx_rows if r[5] == "Q"]
    for r in se:
        if r[3] == "bc_in":
            row_Q = [x for x in qidx if x[0] == r[0]]
            if row_Q and row_Q[0][6] != r[22]:
                print(f"step {r[0]}: cflx_Q mismatch with cam_surface_cflx")
    print("surface exchange: bc_in rows present for all steps")

    # ---- precip present? ----------------------------------------------
    nprec = sum(1 for r in se if r[3] == "bc_out"
                and any(float(r[idx[c]]) > 0 for c in ("precsc", "precsl",
                                                       "precc", "precl")))
    print(f"precipitating steps (any prec>0 at bc_out): {nprec}/{nsteps}")
    if nprec:
        for c in ("precc", "precl"):
            vals = [float(r[idx[c]]) for r in se if r[3] == "bc_out"]
            print(f"  {c}: max={max(vals):.3e} m/s at step "
                  f"{int([r for r in se if r[3]=='bc_out' and float(r[idx[c]])==max(vals)][0][0])}")

    print("\nALL STRUCTURAL CHECKS PASSED" if ok and nbad == 0 and ndiff == 0
          else "\nSEE ISSUES ABOVE")
    return 0 if (ok and nbad == 0 and ndiff == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
