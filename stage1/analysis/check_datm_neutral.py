#!/usr/bin/env python3
# EqESM DATM lane: restart-file solution-neutrality check (instrumented vs reference I2000Clm60Sp_f19 run); run with `source stage1/runs/env.sh` loaded python (netCDF4).
"""Compare restart files between the instrumented and reference runs.
Bit-identical data is expected (instrumentation is additive).
"""
import netCDF4, numpy as np, glob, os

ref = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/cime/cases/I2000Clm60Sp_f19/run"
inst = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/cime/cases/I2000Clm60Sp_cdeps_inst/run"

ok = True
for comp in ["clm2.r", "cpl.r", "datm.r", "dglc.r", "mosart.r"]:
    rp = sorted(glob.glob(os.path.join(ref, f"*{comp}.*.nc")))[-1]
    ip = sorted(glob.glob(os.path.join(inst, f"*{comp}.*.nc")))[-1]
    rd = netCDF4.Dataset(rp); id_ = netCDF4.Dataset(ip)
    diffs = []
    for v in rd.variables:
        if v not in id_.variables:
            diffs.append(f"  missing in inst: {v}"); continue
        rv, iv = rd.variables[v], id_.variables[v]
        if rv.dtype.kind in "SU" or rv.dtype == "S1":
            continue
        a, b = rv[:], iv[:]
        if not np.ma.allequal(a, b):
            a2 = np.ma.filled(a, np.nan).astype(float)
            b2 = np.ma.filled(b, np.nan).astype(float)
            d = np.abs(a2 - b2)
            d = d[np.isfinite(d)]
            diffs.append(f"  DIFF {v}: max|d| = {d.max() if d.size else 0:.6e}")
    status = "OK" if not diffs else "DIFFS"
    if diffs: ok = False
    print(f"{comp:10s} {os.path.basename(rp)} vs {os.path.basename(ip)}: {status}")
    for d in diffs: print(d)
    rd.close(); id_.close()
print("NEUTRAL" if ok else "NOT NEUTRAL")
