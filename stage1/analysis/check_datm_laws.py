#!/usr/bin/env python3
# EqESM DATM lane: Fortran-law verification of the stage-1 dump CSVs in stage1/dumps/cdeps/ against datm_datamode_clmncep_mod.F90, dshr_strdata_mod.F90, shr_precip_mod.F90.
"""Verify every row of the DATM dump CSVs against the Fortran laws.
Checks (all to round-off / tiny epsilon):
  streams layer (cdeps_strdata_streams.csv):
    coszen : value = lb*flb (flb = cosz/tavCosz, fub = 0); value == 0 when flb<=0 (night)
    linear/nearest/lower/upper : value = flb*lb + fub*ub
  datamode layer (datm_gswp3_streams.csv + datm_gswp3_exports.csv), per row:
    u = v = wind/sqrt(2)
    tbot: if tbotmax < 50 add 273.15; clamp >= 180; ptem == tbot
    pbot: 0 -> 101325; pslv == pbot
    shum = min(qsat(esat(tbot)), strm_shum)
    dens from vp formula
    rainc0=0.1*precn, rainl0=0.9*precn; frac=clamp((tbot-273.15)*0.5,0,1)
    rainc=rainc0*frac, rainl=rainl0*frac, snowc=rainc0*(1-frac), snowl=rainl0*(1-frac)
    swndr=0.5*swdn; ratio_nir clamp poly; swndf=(1-ratio)*0.5*swdn; same vis
    swnet: 0 if anidrmax<1e-8 or anidrmax>0.9*1e30; else albedo-weighted; else plain sum
    Sa_topo = strm_topo or 1e30; Sa_z = 30.0 when strm_z absent (spval)
"""
import numpy as np, math

D = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cdeps"
TKFRZ = 273.15
PSTD = 101325.0
RDAIR = 287.0423113650487
SPVAL = 1e30

def esat(tK, tKbot):
    t = min(50.0, max(-50.0, tK - TKFRZ))
    if tKbot < TKFRZ:
        b = [6.109177956, 5.034698970e-1, 1.886013408e-2, 4.176223716e-4,
             5.824720280e-6, 4.838803174e-8, 1.838826904e-10]
        p = b[0] + t*(b[1] + t*(b[2] + t*(b[3] + t*(b[4] + t*(b[5] + t*b[6])))))
    else:
        a = [6.107799961, 4.436518521e-1, 1.428945805e-2, 2.650648471e-4,
             3.031240396e-6, 2.034080948e-8, 6.136820929e-11]
        p = a[0] + t*(a[1] + t*(a[2] + t*(a[3] + t*(a[4] + t*(a[5] + t*a[6])))))
    return 100.0 * p

def ratio_nir(s):
    return min(0.99, max(0.01, 0.29548 + 0.00504*s - 1.4957e-5*s*s + 1.4881e-8*s*s*s))

def ratio_vis(s):
    return min(0.99, max(0.01, 0.17639 + 0.00380*s - 9.0039e-6*s*s + 8.1351e-9*s*s*s))

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

# ---------- streams layer ----------
rows = load(f"{D}/cdeps_strdata_streams.csv")
print(f"streams rows: {len(rows)}")
for r in rows:
    lb, ub, flb, fub, val = (float(r[c]) for c in "lb ub flb fub value".split())
    algo = r["algo"].strip()
    if algo == "coszen":
        chk(fub == 0.0, f"coszen fub!=0 {r['ns']} {r['field']} {r['n']}")
        if flb > 0.0:
            chk(abs(val - lb*flb) <= 1e-10*max(1.0, abs(val)), f"coszen value {r['field']} {r['n']}")
        else:
            chk(val == 0.0, f"coszen night value!=0 {r['field']} {r['n']} flb={flb}")
    else:
        expect = flb*lb + fub*ub
        chk(abs(val - expect) <= 1e-10*max(1.0, abs(val)), f"{algo} value {r['ns']} {r['field']} {r['n']}")
    # factor sanity: linear factors sum to 1
    if algo == "linear":
        chk(abs(flb + fub - 1.0) < 1e-8, f"linear factors {r['field']} {r['n']}")
    if algo == "nearest":
        chk((flb == 1.0 and fub == 0.0) or (flb == 0.0 and fub == 1.0), f"nearest factors {r['field']} {r['n']}")
    if algo == "lower":
        chk(flb == 1.0 and fub == 0.0, f"lower factors {r['field']} {r['n']}")

# ---------- datamode layer ----------
sin_ = load(f"{D}/datm_gswp3_streams.csv")
ex_ = load(f"{D}/datm_gswp3_exports.csv")
print(f"datamode rows: streams {len(sin_)}, exports {len(ex_)}")
assert len(sin_) == len(ex_)
for s, e in zip(sin_, ex_):
    key = (s["step"], s["n"])
    assert (e["step"], e["n"]) == key, "row order mismatch"
    tbotmax = float(s["tbotmax"])
    strm_tbot = float(s["strm_tbot"])
    strm_pbot = float(s["strm_pbot"])
    wind = float(s["strm_wind"])
    strm_shum = float(s["strm_shum"])
    strm_lwdn = float(s["strm_lwdn"])
    strm_swdn = float(s["strm_swdn"])
    precn = float(s["strm_precn"])
    strm_topo = float(s["strm_topo"])
    strm_z = float(s["strm_z"])
    anidrmax = float(e["anidrmax"])
    n = s["n"]
    # tbot
    tbot = strm_tbot + (TKFRZ if tbotmax < 50.0 else 0.0)
    tbot = max(180.0, tbot)
    chk(abs(tbot - float(e["Sa_tbot"])) <= 1e-9*abs(tbot), f"Sa_tbot {n}")
    chk(float(e["Sa_ptem"]) == tbot, f"Sa_ptem {n}")
    # wind split
    u = wind/math.sqrt(2.0)
    chk(abs(u - float(e["Sa_u"])) <= 1e-9*abs(u), f"Sa_u {n}")
    chk(float(e["Sa_v"]) == float(e["Sa_u"]), f"Sa_v {n}")
    # pressure
    pbot = PSTD if strm_pbot == 0.0 else strm_pbot
    chk(float(e["Sa_pbot"]) == pbot, f"Sa_pbot {n}")
    chk(float(e["Sa_pslv"]) == pbot, f"Sa_pslv {n}")
    # topo/z
    chk(float(e["Sa_topo"]) == strm_topo, f"Sa_topo {n}")
    chk(float(e["Sa_z"]) == (30.0 if strm_z >= SPVAL else strm_z), f"Sa_z {n}")
    # shum
    es = esat(tbot, tbot)
    qsat = 0.622*es/(pbot - 0.378*es)
    shum = min(qsat, strm_shum)
    chk(abs(shum - float(e["Sa_shum"])) <= 1e-12*max(1.0, abs(shum)), f"Sa_shum {n}")
    # density
    vp = shum*pbot/(0.622 + 0.378*shum)
    dens = (pbot - 0.378*vp)/(tbot*RDAIR)
    chk(abs(dens - float(e["Sa_dens"])) <= 1e-9*abs(dens), f"Sa_dens {n}")
    # precip split
    rainc0, rainl0 = 0.1*precn, 0.9*precn
    frac = min(1.0, max(0.0, (tbot - TKFRZ)*0.5))
    chk(abs(max(0.0, rainc0*frac) - float(e["Faxa_rainc"])) <= 1e-12, f"Faxa_rainc {n}")
    chk(abs(max(0.0, rainl0*frac) - float(e["Faxa_rainl"])) <= 1e-12, f"Faxa_rainl {n}")
    chk(abs(max(0.0, rainc0*(1.0-frac)) - float(e["Faxa_snowc"])) <= 1e-12, f"Faxa_snowc {n}")
    chk(abs(max(0.0, rainl0*(1.0-frac)) - float(e["Faxa_snowl"])) <= 1e-12, f"Faxa_snowl {n}")
    # shortwave split
    swndr = 0.5*strm_swdn
    rn = ratio_nir(swndr)
    chk(abs(rn*swndr - float(e["Faxa_swndr"])) <= 1e-9*abs(swndr), f"Faxa_swndr {n}")
    chk(abs((1.0-rn)*swndr - float(e["Faxa_swndf"])) <= 1e-9*abs(swndr), f"Faxa_swndf {n}")
    rv = ratio_vis(swndr)
    chk(abs(rv*swndr - float(e["Faxa_swvdr"])) <= 1e-9*abs(swndr), f"Faxa_swvdr {n}")
    chk(abs((1.0-rv)*swndr - float(e["Faxa_swvdf"])) <= 1e-9*abs(swndr), f"Faxa_swvdf {n}")
    # swnet
    sx = [float(e[c]) for c in ("Sx_anidr", "Sx_avsdr", "Sx_anidf", "Sx_avsdf")]
    if anidrmax < 1e-8 or anidrmax > 0.9*SPVAL:
        swnet = 0.0
    elif all(x < SPVAL for x in sx):
        swnet = ((1.0-sx[0])*float(e["Faxa_swndr"]) + (1.0-sx[1])*float(e["Faxa_swvdr"]) +
                 (1.0-sx[2])*float(e["Faxa_swndf"]) + (1.0-sx[3])*float(e["Faxa_swvdf"]))
    else:
        swnet = (float(e["Faxa_swndr"]) + float(e["Faxa_swvdr"]) +
                 float(e["Faxa_swndf"]) + float(e["Faxa_swvdf"]))
    chk(abs(swnet - float(e["Faxa_swnet"])) <= 1e-6*max(1.0, abs(swnet)), f"Faxa_swnet {n}")
    # swdn/lwdn pass-through
    chk(float(e["Faxa_swdn"]) == strm_swdn, f"Faxa_swdn {n}")
    chk(float(e["Faxa_lwdn"]) == strm_lwdn, f"Faxa_lwdn {n}")

print(f"\nTOTAL FAILURES: {nbad}")
