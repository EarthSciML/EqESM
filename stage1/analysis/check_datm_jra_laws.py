#!/usr/bin/env python3
# EqESM DATM lane: Fortran-law verification of the JRA55 datamode stage-1 dump CSVs in
# stage1/dumps/cdeps/ against datm_datamode_jra_mod.F90 and dshr_strdata_mod.F90.
"""Verify every row of the JRA55 DATM dump CSVs against the Fortran laws.
Checks (all to round-off / tiny epsilon):
  streams layer (cdeps_strdata_streams_jra.csv):
    coszen : value = lb*flb (flb = cosz/tavCosz, fub = 0); value == 0 when flb<=0 (night)
    linear/nearest/lower/upper : value = flb*lb + fub*ub
  datamode layer (datm_jra_streams.csv + datm_jra_exports.csv), per row:
    Sa_pslv = Sa_pbot = strm_pslv; Sa_tbot = Sa_ptem = strm_tbot
    Sa_u = strm_u; Sa_v = strm_v; Sa_shum = strm_shum (no qsat law in the JRA datamode:
      the humidity passes through unchanged)
    Faxa_swdn = strm_swdn; Faxa_lwdn = strm_lwdn
    Sa_z = 10 m; Sa_u10m = Sa_u; Sa_v10m = Sa_v
    Sa_dens = pbot/(rdair*tbot*(1 + 0.608*shum)), rdair = 287.0423113650487
    precip: separates_rain_snow=0 -> step split of strm_prec by tbot < TKFRZ (all snow
      below, all rain at/above); rainc = snowc = 0; separates_rain_snow=1 -> Faxa_rainl
      = strm_prrn, Faxa_snowl = strm_prsn (prrn/prsn spval 1e30 when absent)
    radiation: swvdr/ndr/vdf/ndf = swdn*0.28/0.31/0.24/0.17
    avg_alb = 0.069 - 0.011*cos(2*lat*degtorad) (LY08), Faxa_swnet = swdn*(1 - avg_alb)
  cross-checks between the layers:
    every strm_* value equals the streams-layer interpolated value of the mapped
    stream field at the same (step, point): prec->Faxa_prec, lwdn->Faxa_lwdn,
    swdn->Faxa_swdn, q_10->Sa_shum, slp->Sa_pslv, t_10->Sa_tbot, u_10->Sa_u, v_10->Sa_v
    (datm.streams.xml stream_variables of the C_JRA_t62 case)
"""
import math
import numpy as np

D = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cdeps"
TKFRZ = 273.15
RDAIR = 287.0423113650487  # SHR_CONST_RDAIR = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ/SHR_CONST_MWDAIR (shr_const_mod.F90:24-29)
SPVAL = 1e30
DEGTORAD = math.pi / 180.0
# The JRA datamode writes the moisture coefficient and the LY08 albedo
# coefficients as UNSUFFIXED literals (0.608, 0.069, 0.011 in
# datm_datamode_jra_mod.F90:299,324), i.e. default-real = SINGLE precision;
# the laws therefore use the r4-rounded values (see bugs.md).
QCOEF = float(np.float32(0.608))
A0 = float(np.float32(0.069))
A1 = float(np.float32(0.011))
# stream field (netcdf var) -> datm_jra_streams.csv strm_* column (datm.streams.xml)
FIELD2STRM = {"prec": "strm_prec", "lwdn": "strm_lwdn", "swdn": "strm_swdn",
              "q_10": "strm_shum", "slp": "strm_pslv", "t_10": "strm_tbot",
              "u_10": "strm_u", "v_10": "strm_v"}


def load(path):
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    hdr = [l for l in open(path) if l.startswith("# columns")][0]
    cols = [c.strip() for c in hdr.split("=", 1)[1].strip().split(",")]
    rows = []
    for l in lines:
        v = l.strip().split(",")
        rows.append(dict(zip(cols, v)))
    return rows


nbad = 0
def chk(cond, msg):
    global nbad
    if not cond:
        nbad += 1
        print("FAIL:", msg)

def near(a, b, rtol=1e-12):
    return abs(a - b) <= rtol * max(1.0, abs(b))

# ---------- streams layer ----------
rows = load(f"{D}/cdeps_strdata_streams_jra.csv")
print(f"streams rows: {len(rows)}")
bykey = {}
for r in rows:
    lb, ub, flb, fub, val = (float(r[c]) for c in "lb ub flb fub value".split())
    algo = r["algo"].strip()
    if algo == "coszen":
        chk(fub == 0.0, f"coszen fub!=0 {r['ns']} {r['field']} {r['n']}")
        if flb > 0.0:
            chk(near(val, lb*flb), f"coszen value {r['field']} {r['n']}")
        else:
            chk(val == 0.0, f"coszen night value!=0 {r['field']} {r['n']} flb={flb}")
    else:
        expect = flb*lb + fub*ub
        chk(near(val, expect), f"{algo} value {r['ns']} {r['field']} {r['n']}")
    if algo == "linear":
        chk(abs(flb + fub - 1.0) < 1e-8, f"linear factors {r['field']} {r['n']}")
    if algo == "nearest":
        chk((flb == 1.0 and fub == 0.0) or (flb == 0.0 and fub == 1.0),
            f"nearest factors {r['field']} {r['n']}")
    if algo == "lower":
        chk(flb == 1.0 and fub == 0.0, f"lower factors {r['field']} {r['n']}")
    if r["field"] in FIELD2STRM:
        bykey[(r["step"], r["n"], FIELD2STRM[r["field"]])] = val

# ---------- datamode layer ----------
sin_ = load(f"{D}/datm_jra_streams.csv")
ex_ = load(f"{D}/datm_jra_exports.csv")
print(f"datamode rows: streams {len(sin_)}, exports {len(ex_)}")
assert len(sin_) == len(ex_), "streams/exports row-count mismatch"
for s, e in zip(sin_, ex_):
    key = (s["step"], s["n"])
    assert (e["step"], e["n"]) == key, "row order mismatch"
    n = s["n"]
    tbot = float(s["strm_tbot"]); pslv = float(s["strm_pslv"])
    u = float(s["strm_u"]); v = float(s["strm_v"]); shum = float(s["strm_shum"])
    prec = float(s["strm_prec"]); prrn = float(s["strm_prrn"]); prsn = float(s["strm_prsn"])
    swdn = float(s["strm_swdn"]); lwdn = float(s["strm_lwdn"])
    sep = int(s["separates_rain_snow"])
    lat = float(s["lat"])
    assert sep in (0, 1)
    # stream copies into the export state
    chk(near(float(e["Sa_pslv"]), pslv), f"Sa_pslv {n}")
    chk(near(float(e["Sa_pbot"]), pslv), f"Sa_pbot {n}")
    chk(near(float(e["Sa_tbot"]), tbot), f"Sa_tbot {n}")
    chk(near(float(e["Sa_ptem"]), tbot), f"Sa_ptem {n}")
    chk(near(float(e["Sa_u"]), u), f"Sa_u {n}")
    chk(near(float(e["Sa_v"]), v), f"Sa_v {n}")
    chk(near(float(e["Sa_shum"]), shum), f"Sa_shum {n}")
    chk(near(float(e["Faxa_swdn"]), swdn), f"Faxa_swdn {n}")
    chk(near(float(e["Faxa_lwdn"]), lwdn), f"Faxa_lwdn {n}")
    # constants and wind
    chk(float(e["Sa_z"]) == 10.0, f"Sa_z {n}")
    chk(near(float(e["Sa_u10m"]), u), f"Sa_u10m {n}")
    chk(near(float(e["Sa_v10m"]), v), f"Sa_v10m {n}")
    # density (the 0.608 coefficient is an unsuffixed single-precision
    # literal in the Fortran, hence QCOEF = r4(0.608))
    dens = pslv / (RDAIR * tbot * (1.0 + QCOEF*shum))
    chk(near(float(e["Sa_dens"]), dens), f"Sa_dens {n}")
    # precip split
    chk(float(e["Faxa_rainc"]) == 0.0, f"Faxa_rainc {n}")
    chk(float(e["Faxa_snowc"]) == 0.0, f"Faxa_snowc {n}")
    if sep == 1:
        chk(near(float(e["Faxa_rainl"]), prrn), f"Faxa_rainl(sep) {n}")
        chk(near(float(e["Faxa_snowl"]), prsn), f"Faxa_snowl(sep) {n}")
    else:
        chk(prrn == SPVAL and prsn == SPVAL, f"prrn/prsn not spval {n}")
        if tbot < TKFRZ:
            rainl, snowl = 0.0, prec
        else:
            rainl, snowl = prec, 0.0
        chk(near(float(e["Faxa_rainl"]), rainl), f"Faxa_rainl {n}")
        chk(near(float(e["Faxa_snowl"]), snowl), f"Faxa_snowl {n}")
    # shortwave split
    chk(near(float(e["Faxa_swvdr"]), 0.28*swdn), f"Faxa_swvdr {n}")
    chk(near(float(e["Faxa_swndr"]), 0.31*swdn), f"Faxa_swndr {n}")
    chk(near(float(e["Faxa_swvdf"]), 0.24*swdn), f"Faxa_swvdf {n}")
    chk(near(float(e["Faxa_swndf"]), 0.17*swdn), f"Faxa_swndf {n}")
    # LY08 albedo and net shortwave (the 0.069/0.011 coefficients are
    # unsuffixed single-precision literals in the Fortran, hence A0/A1)
    alb = A0 - A1*math.cos(2.0*lat*DEGTORAD)
    chk(near(float(e["avg_alb"]), alb), f"avg_alb {n}")
    chk(near(float(e["Faxa_swnet"]), swdn*(1.0 - alb)), f"Faxa_swnet {n}")
    # cross-check against the streams layer
    for col in ("strm_tbot", "strm_pslv", "strm_u", "strm_v", "strm_shum",
                "strm_prec", "strm_lwdn", "strm_swdn"):
        if (s["step"], n, col) in bykey:
            chk(near(float(s[col]), bykey[(s["step"], n, col)]), f"cross {col} {n}")

print(f"\nTOTAL FAILURES: {nbad}")
