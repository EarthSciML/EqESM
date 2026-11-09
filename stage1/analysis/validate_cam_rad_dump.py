#!/usr/bin/env python3
"""Validate the CAM radiation (RRTMG) CSV dumps of the stage-1 instrumented run.

Reads the 8 cam_rad_*.csv files produced by cam_dump_mod.F90
(cam_dump_rad_flush, radiation.F90 radiation_tend, climate icall=0 call) and
checks them against the real identities of the CAM/RRTMG radiation code:

  * cam_rad_input.csv     (step, RRTMG layer k=1..33)
      - cam_lev = max(k-1,1); k=1 is the extra layer above model top with
        pint_lo = 0.01 Pa (rrtmg_state.F90:141-144)
      - pmid(k) = 0.5*(pint_lo+pint_hi) for k>=2 (rrtmg_state.F90:126),
        pmid(1) = 0.5*pint_hi(1) (rrtmg_state.F90:140)
      - shared edges: pint_hi(k) == pint_lo(k+1), T_ifc_hi(k) == T_ifc_lo(k+1)
      - h2ovmr = (q/(1-q)) * amdw, amdw = 1.607793 (rrtmg_state.F90:47,198)
      - cfc22vmr = ccl4vmr = 0 in this configuration (rrtmg_state.F90:222-223)
      - cldfprime = max(cld, cldfsnow) (radiation.F90:956-976)
      - T_ifc from surface upward flux / linear interpolation in ln p
        (rrtmg_state.F90:110-120)
  * cam_rad_scal.csv      (step)
      - num_rrtmg_levs = 33
      - solin = sum(sfac_b1..b14 * solar_ref_band_irradiance) * eccf * coszrs
        (radsw.F90:395-397; solar_ref from radconstants.F90:34-40)
      - cross-file: fsnt = fns(ifc=1), fsntc = fcns(ifc=1), flnt = fnl(ifc=1),
        flntc = fcnl(ifc=1), flut = ful(ifc=1), flutc = fsul(ifc=1),
        fsns = fns(ifc=33), fsnsc = fcns(ifc=33), fsds = fds(ifc=33),
        fsdsc = fdsc(ifc=33), flns = fnl(ifc=33), flnsc = fcnl(ifc=33),
        flwds = fdl(ifc=33), fldsc = fsdl(ifc=33)
      - fds(ifc=1) <= solin (absorption in the extra layer above model top)
      - fsntoa = solin - fus(ifc=1) to within the extra-layer absorption
        (info only)
  * cam_rad_out.csv       (step, lev=1..32)
      - qrs(k) = -g * (fns(k+1)-fns(k)) / (pintmb(k+2)-pintmb(k+1))
        (flux divergence; radsw.F90:573-574, dps=1/86400)
      - qrsc(k) same with fcns; qrl/qrlc with fnl/fcnl (LW net is
        upward-positive: qrl = +g*(fnl(k+1)-fnl(k))/dp, radlw.F90:244-245)
      - column integral: sum_k qrs(k)*dp_k/g = fsnt - fsns (SW)
  * cam_rad_out_ifc.csv   (step, ifc=1..33)
      - fns = fds - fus (radsw.F90:565-567), fnl = ful - fdl (radlw.F90:237-238)
        (the clear-sky up/down components fusc/fdsc/fsul/fsdl are not dumped;
        fcns/fcnl are cross-checked against the scal file at ifc 1 and 33)
      - sols+soll+solsd+solld = fds(ifc=33) (surface total down; radsw.F90:564-568)
      - monotone interfaces: fds decreasing downward, fdl increasing
  * cam_rad_cldopt_sw.csv (step, band=1..14, lev)
      - ssa = tau_w/tau if tau>0 else 1; asm = tau_w_g/tau_w if tau_w>0 else 0;
        fsfc = tau_w_f/tau_w if tau_w>0 else 0 (radsw.F90:437-475)
      - 0 <= tau_w <= tau, tau_w_g <= tau_w, tau_w_f <= tau_w
      - tau = 0 when cld = 0
  * cam_rad_cldopt_lw.csv (step, band=1..16, lev)
      - tau_lw >= 0; tau_lw = 0 when cld = 0
  * cam_rad_solver_sw.csv (step, gpt=1..112, lev)
      - cld_sto in {0,1} (mcica_subcol_gen_sw.f90:404-414)
      - tauc_sto = tau(band(gpt)) when cld_sto=1 else 0 (line 439-444)
      - ssac_sto = ssa(band(gpt)) when cld_sto=1 else 1
      - asmc_sto = asm(band(gpt)) when cld_sto=1 else 0
      - fsfc_sto = fsfc(band(gpt)) when cld_sto=1 else 0
      - band(gpt) from the RRTMG SW k-distribution (rrtmg_sw_init.f90:266-280,
        ngptsw=112: gpt counts 6,12,8,8,10,10,2,10,8,6,6,8,6,12)
  * cam_rad_solver_lw.csv (step, gpt=1..140, lev)
      - cld_sto in {0,1}; tauc_sto = tau_lw(band(gpt)) when cld_sto=1 else 0
        (mcica_subcol_gen_lw.f90:408-416; ngptlw=140: counts
        10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2)

Also checks run structure: 23 radiation steps (3,5,...,47, iradsw=iradlw=2),
row counts per step, and that all 8 files cover the same steps.

Usage: python3 validate_cam_rad_dump.py [--dumps-dir DIR] [--info]
Exit code 0 = all hard checks passed (info checks reported separately).
"""
import argparse
import sys

import numpy as np

# --- physical constants used by the dump identity checks -------------------
G = 9.80616          # gravitational acceleration (m/s2), physconst.F90
AMDW = 1.607793      # molecular weight of dry air / water vapor (rrtmg_state.F90:47)
SOLAR_REF = np.array([12.11, 20.3600000000001, 23.73, 22.43, 55.63, 102.93,
                      24.29, 345.74, 218.19, 347.20, 129.49, 50.15, 3.08,
                      12.89])  # W/m2, radconstants.F90:34-40 (sum = 1368.22)
# gpt -> band (1-based CAM dump band) mappings from the RRTMG k-distributions
SW_GPT_COUNTS = [6, 12, 8, 8, 10, 10, 2, 10, 8, 6, 6, 8, 6, 12]   # ngptsw=112
LW_GPT_COUNTS = [10, 12, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2]  # 140
SW_GPT_BAND = np.repeat(np.arange(1, 15), SW_GPT_COUNTS)
LW_GPT_BAND = np.repeat(np.arange(1, 17), LW_GPT_COUNTS)

PVER = 32            # CAM levels (PLEV=32)
PVERP = 33           # CAM interfaces
NUM_RRTMG_LEVS = 33  # pverp in this configuration
RAD_STEPS = list(range(3, 49, 2))  # iradsw=iradlw=2, radiation on odd steps 3..47


class Checker:
    def __init__(self, info=False):
        self.failures = []
        self.info_notes = []
        self.info = info

    def check(self, name, ok, detail=""):
        if not ok:
            self.failures.append(f"{name}: {detail}")

    def info_note(self, msg):
        if self.info:
            self.info_notes.append(msg)

    def close(self, name, a, b, rtol=1e-9, atol=0.0, what=""):
        """Return max |a-b| (0 if arrays are empty); record a failure if
        the mismatch exceeds rtol*|b| + atol."""
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        if a.size == 0:
            return 0.0
        d = np.abs(a - b)
        scale = np.maximum(np.abs(b), 1e-300)
        bad = d > rtol * scale + atol
        nbad = int(bad.sum())
        if nbad:
            i = np.argmax(d)
            self.check(f"{name}", False,
                       f"{what} max |a-b|={d[i]:.6e} at {i} "
                       f"(a={a[i]:.17e}, b={b[i]:.17e}), {nbad}/{d.size} bad")
        return float(d.max()) if d.size else 0.0


def read_csv(path):
    """Return (columns, data) where data is a dict column-name -> np.array."""
    cols = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# description ="):
                continue
            if line.startswith("# columns ="):
                cols = line[len("# columns ="):].strip().split(",")
                continue
            rows.append([x.strip() for x in line.split(",")])
    if cols is None:
        raise ValueError(f"{path}: no '# columns' header")
    data = {c: np.array([float(r[i]) for r in rows]) for i, c in enumerate(cols)}
    return cols, data


def by_step(data, step_col="step"):
    """Return {step: row-index array} for one step per unique step value."""
    steps = np.unique(data[step_col]).astype(int)
    return {int(s): np.where(data[step_col] == s)[0] for s in steps}, steps


def step_rows(rows_of, steps):
    return np.concatenate([rows_of[s] for s in steps])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dumps-dir", default="stage1/dumps/cam",
                    help="directory holding the cam_rad_*.csv files")
    ap.add_argument("--info", action="store_true",
                    help="also print info-level soft checks")
    args = ap.parse_args()

    chk = Checker(info=args.info)

    files = ["cam_rad_input.csv", "cam_rad_scal.csv", "cam_rad_out.csv",
             "cam_rad_out_ifc.csv", "cam_rad_cldopt_sw.csv",
             "cam_rad_cldopt_lw.csv", "cam_rad_solver_sw.csv",
             "cam_rad_solver_lw.csv"]
    D = {}
    for fn in files:
        try:
            cols, data = read_csv(f"{args.dumps_dir}/{fn}")
        except FileNotFoundError:
            print(f"FAIL: missing {args.dumps_dir}/{fn}")
            sys.exit(1)
        D[fn] = (cols, data)
        steps_uniq, steps = by_step(data)
        chk.check(f"{fn} steps", set(steps) == set(RAD_STEPS),
                  f"got {sorted(steps)} expected {RAD_STEPS}")

    # ---- common structure --------------------------------------------------
    exp_rows = {"cam_rad_input.csv": NUM_RRTMG_LEVS,
                "cam_rad_scal.csv": 1,
                "cam_rad_out.csv": PVER,
                "cam_rad_out_ifc.csv": PVERP,
                "cam_rad_cldopt_sw.csv": 14 * PVER,
                "cam_rad_cldopt_lw.csv": 16 * PVER,
                "cam_rad_solver_sw.csv": 112 * PVER,
                "cam_rad_solver_lw.csv": 140 * PVER}
    for fn, per in exp_rows.items():
        _, data = D[fn]
        steps_uniq, steps = by_step(data)
        got = {s: int(len(np.where(data["step"] == s)[0])) for s in steps}
        chk.check(f"{fn} rows/step", all(n == per for n in got.values()),
                  f"expected {per} rows per step, got {got}")
        # tod monotonic across steps (repeats within each step's rows)
        tod = data["tod"]
        _, first = np.unique(data["step"], return_index=True)
        chk.check(f"{fn} tod monotone", np.all(np.diff(tod[first]) > 0),
                  f"tods {tod[first]}")

    # ---- cam_rad_input.csv -------------------------------------------------
    _, inp = D["cam_rad_input.csv"]
    steps_inp, steps = by_step(inp)
    rows = step_rows(steps_inp, RAD_STEPS)
    k = inp["k"][rows]
    chk.check("input k range", k.min() == 1 and k.max() == NUM_RRTMG_LEVS,
              f"k in {k.min()}..{k.max()}")
    chk.check("input cam_lev", np.all(inp["cam_lev"][rows] ==
                                      np.maximum(k - 1, 1)),
              "cam_lev != max(k-1,1)")
    plo, phi = inp["pint_lo"][rows], inp["pint_hi"][rows]
    pmid = inp["pmid"][rows]
    # shared edges within each step
    for s in RAD_STEPS:
        r = steps_inp[s]
        chk.check(f"input pint shared edges (step {s})",
                  np.allclose(inp["pint_hi"][r][:-1], inp["pint_lo"][r][1:],
                              rtol=0, atol=1e-9),
                  "pint_hi(k) != pint_lo(k+1)")
        chk.check(f"input T_ifc shared edges (step {s})",
                  np.allclose(inp["T_ifc_hi"][r][:-1], inp["T_ifc_lo"][r][1:],
                              rtol=0, atol=1e-9),
                  "T_ifc_hi(k) != T_ifc_lo(k+1)")
    chk.check("input top pressure", np.allclose(plo[k == 1], 0.01, rtol=0,
                                                atol=1e-9),
              f"pint_lo(k=1) = {np.unique(plo[k==1])} (expect 0.01 Pa)")
    # pmid relations (k>=2: arithmetic mean of edges; k=1: half of top edge)
    m2 = (k >= 2)
    chk.close("input pmid mean k>=2", pmid[m2],
              0.5 * (plo[m2] + phi[m2]), rtol=1e-12, atol=1e-12,
              what="pmid vs (pint_lo+pint_hi)/2")
    chk.close("input pmid k=1", pmid[k == 1], 0.5 * phi[k == 1],
              rtol=1e-12, atol=1e-12, what="pmid(1) vs pint_hi(1)/2")
    # monotone pressure increasing with k (k=1 is the top layer; pmid only
    # for k>=2: the extra-layer midpoint pmid(1) = 0.5*pint_hi(1) sits
    # above pmid(2))
    for s in RAD_STEPS:
        r = steps_inp[s]
        chk.check(f"input pressure monotone (step {s})",
                  np.all(np.diff(inp["pint_hi"][r]) > 0) and
                  np.all(np.diff(inp["pmid"][r][1:]) > 0),
                  "pressure not monotone increasing with k")
    # h2ovmr identity
    q = inp["q_H2O"][rows]
    h2o = inp["h2ovmr"][rows]
    chk.close("input h2ovmr", h2o, (q / (1.0 - q)) * AMDW,
              rtol=1e-12, atol=1e-15, what="h2ovmr vs (q/(1-q))*amdw")
    chk.check("input cfc22vmr==0", np.all(inp["cfc22vmr"][rows] == 0.0),
              "cfc22vmr nonzero")
    chk.check("input ccl4vmr==0", np.all(inp["ccl4vmr"][rows] == 0.0),
              "ccl4vmr nonzero")
    # cldfprime = max(cld, cldfsnow)
    cld, cfs, cfp = inp["cld"][rows], inp["cldfsnow"][rows], inp["cldfprime"][rows]
    chk.close("input cldfprime", cfp, np.maximum(cld, cfs),
              rtol=1e-12, atol=1e-15, what="cldfprime vs max(cld,cldfsnow)")
    chk.check("input cld in [0,1]", np.all(cld >= -1e-12) and np.all(cld <= 1 + 1e-12),
              f"cld range {cld.min():.3e}..{cld.max():.3e}")
    chk.check("input cldfsnow in [0,1]",
              np.all(cfs >= -1e-12) and np.all(cfs <= 1 + 1e-12),
              f"cldfsnow range {cfs.min():.3e}..{cfs.max():.3e}")
    # cloud water / optics inputs
    for col in ("iclwp", "iciwp", "icswp", "dei", "lamc", "pgam", "des"):
        v = inp[col][rows]
        chk.check(f"input {col} >= 0", np.all(v >= -1e-12),
                  f"min {v.min():.3e}")
    chk.check("input Ts range", np.all((inp["Ts"][rows] > 200) &
                                       (inp["Ts"][rows] < 400)),
              f"Ts in {inp['Ts'][rows].min():.2f}..{inp['Ts'][rows].max():.2f}")
    for col in ("asdir", "asdif", "aldir", "aldif"):
        v = inp[col][rows]
        chk.check(f"input {col} in [0,1]", np.all(v >= -1e-9) and np.all(v <= 1 + 1e-9),
                  f"range {v.min():.3e}..{v.max():.3e}")
    # gas vmr sanity ranges.  In this SCAM case the greenhouse gas vmrs are
    # prescribed by the boundary data (scamMod.F90 readiopdata: "using column
    # value of co2vmr from boundary data"), so only shape/positivity checks.
    chk.check("input o3vmr > 0", np.all(inp["o3vmr"][rows] > 0),
              "o3vmr <= 0 somewhere")
    chk.check("input o2vmr ~0.2095", np.allclose(inp["o2vmr"][rows], 0.2095,
                                                 rtol=1e-3, atol=1e-5),
              f"o2vmr in {inp['o2vmr'][rows].min():.6f}..{inp['o2vmr'][rows].max():.6f}")
    chk.check("input gas vmrs positive",
              np.all(inp["co2vmr"][rows] > 0) and np.all(inp["ch4vmr"][rows] >= 0)
              and np.all(inp["n2ovmr"][rows] >= 0),
              f"co2vmr {inp['co2vmr'][rows].min():.3e}, "
              f"ch4vmr {inp['ch4vmr'][rows].min():.3e}, "
              f"n2ovmr {inp['n2ovmr'][rows].min():.3e}")

    # ---- cam_rad_scal.csv --------------------------------------------------
    _, sc = D["cam_rad_scal.csv"]
    steps_sc, _ = by_step(sc)
    rows = step_rows(steps_sc, RAD_STEPS)
    chk.check("scal num_rrtmg_levs", np.all(sc["num_rrtmg_levs"][rows] ==
                                            NUM_RRTMG_LEVS),
              f"{np.unique(sc['num_rrtmg_levs'][rows])}")
    cosz = sc["coszrs"][rows]
    chk.check("scal coszrs range", np.all(np.abs(cosz) <= 1.5),
              f"coszrs in {cosz.min():.6e}..{cosz.max():.6e}")
    day_steps = [s for s in RAD_STEPS if cosz[steps_sc[s][0]] > 0]
    night_steps = [s for s in RAD_STEPS if cosz[steps_sc[s][0]] <= 0]
    chk.info_note(f"day steps ({len(day_steps)}): {day_steps}; "
                  f"night steps ({len(night_steps)}): {night_steps}")
    eccf = sc["eccf"][rows]
    chk.check("scal eccf range", np.all((eccf > 0.95) & (eccf < 1.05)),
              f"eccf {eccf.min():.6f}..{eccf.max():.6f}")
    sfac = np.column_stack([sc[f"sfac_b{i}"][rows] for i in range(1, 15)])
    solin = sc["solin"][rows]
    for s in day_steps:
        r = steps_sc[s]
        chk.close(f"scal solin identity (step {s})", sc["solin"][r],
                  (sfac[r] @ SOLAR_REF) * eccf[r] * cosz[r],
                  rtol=1e-12, atol=1e-9,
                  what="solin vs sum(sfac*solar_ref)*eccf*coszrs")
    for s in night_steps:
        r = steps_sc[s]
        chk.check(f"scal solin=0 at night (step {s})",
                  np.all(sc["solin"][r] == 0.0),
                  f"solin = {sc['solin'][r]}")
        # night SW boundary is all zero (radsw.F90 early return; the dump
        # flush zero-fills the SW buffers, cam_dump_mod.F90:1752-1767)
        for col in ("fsnt", "fsns", "fsds", "fsntoa", "fsutoa", "fsntc",
                    "fsntoac", "fsnsc", "fsdsc", "fsnirtoa", "fsnrtoac",
                    "fsnrtoaq", "fsn200", "fsn200c", "sols", "soll",
                    "solsd", "solld"):
            chk.check(f"scal {col}=0 at night (step {s})",
                      np.all(sc[col][r] == 0.0), f"{sc[col][r]}")

    # ---- cam_rad_out_ifc.csv -----------------------------------------------
    _, oi = D["cam_rad_out_ifc.csv"]
    steps_oi, _ = by_step(oi)
    rows = step_rows(steps_oi, RAD_STEPS)
    fus, fds, fns, fcns = oi["fus"][rows], oi["fds"][rows], oi["fns"][rows], oi["fcns"][rows]
    ful, fdl, fnl, fcnl = oi["ful"][rows], oi["fdl"][rows], oi["fnl"][rows], oi["fcnl"][rows]
    chk.close("ifc fns=fds-fus", fns, fds - fus, rtol=1e-12, atol=1e-9,
              what="fns vs fds-fus")
    chk.close("ifc fnl=ful-fdl", fnl, ful - fdl, rtol=1e-12, atol=1e-9,
              what="fnl vs ful-fdl")
    chk.check("ifc fds >= fus >= 0", np.all(fds >= -1e-9) and
              np.all(fds - fus >= -1e-9) and np.all(fus >= -1e-9),
              f"fds {fds.min():.4e}..{fds.max():.4e}, fus {fus.min():.4e}..{fus.max():.4e}")
    chk.check("ifc fdl >= 0, ful >= 0", np.all(fdl >= -1e-9) and
              np.all(ful >= -1e-9),
              f"fdl {fdl.min():.4e}..{fdl.max():.4e}, ful {ful.min():.4e}..{ful.max():.4e}")
    # monotone: fds decreases downward, fdl increases downward (per step)
    for s in RAD_STEPS:
        r = steps_oi[s]
        chk.check(f"ifc fds monotone (step {s})",
                  np.all(np.diff(oi["fds"][r]) <= 1e-6),
                  "fds not monotone decreasing")
        chk.check(f"ifc fdl monotone (step {s})",
                  np.all(np.diff(oi["fdl"][r]) >= -1e-6),
                  "fdl not monotone increasing")
    # surface spectral components sum to total down (day steps only;
    # sols..solld are scal-file columns, compared with the surface fds)
    day_surf = np.concatenate([steps_oi[s][oi["ifc"][steps_oi[s]] == PVERP]
                               for s in day_steps])
    day_scal = np.concatenate([steps_sc[s] for s in day_steps])
    chk.close("ifc sols+soll+solsd+solld", sc["sols"][day_scal] +
              sc["soll"][day_scal] + sc["solsd"][day_scal] +
              sc["solld"][day_scal], oi["fds"][day_surf],
              rtol=1e-6, atol=1e-4, what="sols+soll+solsd+solld vs fds(surface)")

    # ---- cross-file: scal vs ifc -------------------------------------------
    for s in RAD_STEPS:
        ri = steps_oi[s]
        rs = steps_sc[s]
        one = lambda v: v[rs][0]
        def ifc(col, k):
            return oi[col][ri][oi["ifc"][ri] == k][0]
        chk.close(f"scal/ifc fsnt (step {s})", one(sc["fsnt"]), ifc("fns", 1),
                  rtol=1e-12, atol=1e-9, what="fsnt vs fns(1)")
        chk.close(f"scal/ifc fsntc (step {s})", one(sc["fsntc"]), ifc("fcns", 1),
                  rtol=1e-12, atol=1e-9, what="fsntc vs fcns(1)")
        chk.close(f"scal/ifc fsns (step {s})", one(sc["fsns"]), ifc("fns", PVERP),
                  rtol=1e-12, atol=1e-9, what="fsns vs fns(pverp)")
        chk.close(f"scal/ifc fsds (step {s})", one(sc["fsds"]), ifc("fds", PVERP),
                  rtol=1e-12, atol=1e-9, what="fsds vs fds(pverp)")
        chk.close(f"scal/ifc flnt (step {s})", one(sc["flnt"]), ifc("fnl", 1),
                  rtol=1e-12, atol=1e-9, what="flnt vs fnl(1)")
        chk.close(f"scal/ifc flntc (step {s})", one(sc["flntc"]), ifc("fcnl", 1),
                  rtol=1e-12, atol=1e-9, what="flntc vs fcnl(1)")
        chk.close(f"scal/ifc flut (step {s})", one(sc["flut"]), ifc("ful", 1),
                  rtol=1e-12, atol=1e-9, what="flut vs ful(1)")
        chk.close(f"scal/ifc flns (step {s})", one(sc["flns"]), ifc("fnl", PVERP),
                  rtol=1e-12, atol=1e-9, what="flns vs fnl(pverp)")
        chk.close(f"scal/ifc flwds (step {s})", one(sc["flwds"]), ifc("fdl", PVERP),
                  rtol=1e-12, atol=1e-9, what="flwds vs fdl(pverp)")
        chk.close(f"scal/ifc flnsc (step {s})", one(sc["flnsc"]), ifc("fcnl", PVERP),
                  rtol=1e-12, atol=1e-9, what="flnsc vs fcnl(pverp)")
        # fds at top <= solin (extra-layer absorption)
        chk.check(f"scal fds(1)<=solin (step {s})",
                  ifc("fds", 1) <= one(sc["solin"]) + 1e-9,
                  f"fds(1)={ifc('fds',1):.6e} solin={one(sc['solin']):.6e}")
        if s in day_steps:
            chk.info_note(
                f"step {s}: extra-layer absorption solin-fds(1) = "
                f"{one(sc['solin'])-ifc('fds',1):.4f} W/m2, "
                f"TOA net solin-fus(1) = {one(sc['solin'])-ifc('fus',1):.4f} "
                f"(fsntoa={one(sc['fsntoa']):.4f})")

    # ---- cam_rad_out.csv ---------------------------------------------------
    _, out = D["cam_rad_out.csv"]
    steps_out, _ = by_step(out)
    rows = step_rows(steps_out, RAD_STEPS)
    lev = out["lev"][rows]
    qrs, qrl, qrsc, qrlc = (out["qrs"][rows], out["qrl"][rows],
                            out["qrsc"][rows], out["qrlc"][rows])
    chk.check("out lev range", lev.min() == 1 and lev.max() == PVER,
              f"lev {lev.min()}..{lev.max()}")
    # build pintmb per step from cam_rad_input (0-based: pintmb[0]=pint_lo(1),
    # pintmb[i] = pint_hi(row k=i) for i=1..33)
    qrs_max = qrl_max = 0.0
    for s in RAD_STEPS:
        r_in = steps_inp[s]
        r_oi = steps_oi[s]
        r_out = steps_out[s]
        pintmb = np.concatenate([[inp["pint_lo"][r_in][0]],
                                 inp["pint_hi"][r_in]])
        fns_s = oi["fns"][r_oi]
        fcns_s = oi["fcns"][r_oi]
        fnl_s = oi["fnl"][r_oi]
        fcnl_s = oi["fcnl"][r_oi]
        qo = out["qrs"][r_out]
        dp = pintmb[1:] - pintmb[:-1]           # 33 layer thicknesses
        dplev = dp[1:]                          # CAM level thicknesses (32)
        div_sw = -G * (fns_s[1:] - fns_s[:-1]) / dplev
        div_sw_c = -G * (fcns_s[1:] - fcns_s[:-1]) / dplev
        div_lw = G * (fnl_s[1:] - fnl_s[:-1]) / dplev
        div_lw_c = G * (fcnl_s[1:] - fcnl_s[:-1]) / dplev
        qrs_max = max(qrs_max, chk.close(f"out qrs flux div (step {s})",
                                         qo, div_sw, rtol=1e-6, atol=1e-6,
                                         what="qrs vs -g*(fns(k+1)-fns(k))/dp"))
        chk.close(f"out qrsc flux div (step {s})", out["qrsc"][r_out],
                  div_sw_c, rtol=1e-6, atol=1e-6, what="qrsc vs clear-sky div")
        qrl_max = max(qrl_max, chk.close(f"out qrl flux div (step {s})",
                                         out["qrl"][r_out], div_lw,
                                         rtol=1e-6, atol=1e-6,
                                         what="qrl vs g*(fnl(k+1)-fnl(k))/dp"))
        chk.close(f"out qrlc flux div (step {s})", out["qrlc"][r_out],
                  div_lw_c, rtol=1e-6, atol=1e-6, what="qrlc vs clear-sky div")
        # column integral: sum qrs*dp/g = fsnt - fsns
        rs = steps_sc[s]
        tot_sw = float(np.sum(qo * dplev / G))
        net_sw = float(sc["fsnt"][rs][0] - sc["fsns"][rs][0])
        chk.close(f"out qrs column integral (step {s})", tot_sw, net_sw,
                  rtol=1e-6, atol=1e-6,
                  what=f"sum(qrs dp/g)={tot_sw:.6f} vs fsnt-fsns={net_sw:.6f}")
        tot_lw = float(np.sum(out["qrl"][r_out] * dplev / G))
        net_lw = float(sc["flns"][rs][0] - sc["flnt"][rs][0])
        chk.close(f"out qrl column integral (step {s})", tot_lw, net_lw,
                  rtol=1e-6, atol=1e-6,
                  what=f"sum(qrl dp/g)={tot_lw:.6f} vs flns-flnt={net_lw:.6f}")
    chk.info_note(f"max |qrs - flux-div| = {qrs_max:.3e} W/kg, "
                  f"max |qrl - flux-div| = {qrl_max:.3e} W/kg")

    # ---- cam_rad_cldopt_sw.csv ---------------------------------------------
    _, cs = D["cam_rad_cldopt_sw.csv"]
    steps_cs, _ = by_step(cs)
    rows = step_rows(steps_cs, RAD_STEPS)
    tau, tauw, tauwg, tauwf = (cs["tau"][rows], cs["tau_w"][rows],
                               cs["tau_w_g"][rows], cs["tau_w_f"][rows])
    ssa, asm, fsfc = cs["ssa"][rows], cs["asm"][rows], cs["fsfc"][rows]
    chk.check("cldopt_sw band range",
              cs["band"][rows].min() == 1 and cs["band"][rows].max() == 14,
              f"band {cs['band'][rows].min()}..{cs['band'][rows].max()}")
    tz = tau > 0
    chk.close("cldopt_sw ssa", ssa[tz], tauw[tz] / tau[tz],
              rtol=1e-12, atol=1e-15, what="ssa vs tau_w/tau (tau>0)")
    chk.close("cldopt_sw ssa (tau=0)", ssa[~tz], np.ones(int((~tz).sum())),
              rtol=0, atol=1e-15, what="ssa=1 when tau=0")
    wz = tauw > 0
    chk.close("cldopt_sw asm", asm[wz], tauwg[wz] / tauw[wz],
              rtol=1e-12, atol=1e-15, what="asm vs tau_w_g/tau_w")
    chk.close("cldopt_sw asm (tau_w=0)", asm[~wz], np.zeros(int((~wz).sum())),
              rtol=0, atol=1e-15, what="asm=0 when tau_w=0")
    chk.close("cldopt_sw fsfc", fsfc[wz], tauwf[wz] / tauw[wz],
              rtol=1e-12, atol=1e-15, what="fsfc vs tau_w_f/tau_w")
    chk.close("cldopt_sw fsfc (tau_w=0)", fsfc[~wz], np.zeros(int((~wz).sum())),
              rtol=0, atol=1e-15, what="fsfc=0 when tau_w=0")
    chk.check("cldopt_sw bounds", np.all(tau >= 0) and np.all(tauw >= 0) and
              np.all(tauwg >= 0) and np.all(tauwf >= 0) and
              np.all(tauw <= tau + 1e-15) and
              np.all(tauwg <= tauw + 1e-15) and np.all(tauwf <= tauw + 1e-15),
              "tau/tau_w/tau_w_g/tau_w_f bounds violated")
    cz = cs["cld"][rows] == 0
    chk.check("cldopt_sw tau=0 when cld=0",
              np.all(tau[cz] == 0),
              f"{int((tau[cz] != 0).sum())} nonzero tau with cld=0")

    # ---- cam_rad_cldopt_lw.csv ---------------------------------------------
    _, cl = D["cam_rad_cldopt_lw.csv"]
    steps_cl, _ = by_step(cl)
    rows = step_rows(steps_cl, RAD_STEPS)
    taulw = cl["tau_lw"][rows]
    chk.check("cldopt_lw band range",
              cl["band"][rows].min() == 1 and cl["band"][rows].max() == 16,
              f"band {cl['band'][rows].min()}..{cl['band'][rows].max()}")
    chk.check("cldopt_lw tau_lw >= 0", np.all(taulw >= 0),
              f"min {taulw.min():.3e}")
    cz = cl["cld"][rows] == 0
    chk.check("cldopt_lw tau_lw=0 when cld=0", np.all(taulw[cz] == 0),
              f"{int((taulw[cz] != 0).sum())} nonzero tau_lw with cld=0")

    # ---- cam_rad_solver_sw.csv ---------------------------------------------
    _, ss = D["cam_rad_solver_sw.csv"]
    steps_ss, _ = by_step(ss)
    rows = step_rows(steps_ss, RAD_STEPS)
    gpt = ss["gpt"][rows].astype(int)
    csto = ss["cld_sto"][rows]
    chk.check("solver_sw gpt range", gpt.min() == 1 and gpt.max() == 112,
              f"gpt {gpt.min()}..{gpt.max()}")
    chk.check("solver_sw cld_sto in {0,1}",
              np.all((csto == 0) | (csto == 1)),
              f"values {np.unique(csto)}")
    # join band values from cldopt_sw by (step, band, lev); rebuild lookups
    # explicitly instead of assuming row order across files:
    tau_band = np.zeros_like(csto)
    ssa_band = np.zeros_like(csto)
    asm_band = np.zeros_like(csto)
    fsf_band = np.zeros_like(csto)
    cld_val = np.zeros_like(csto)
    for s in RAD_STEPS:
        rcs = steps_cs[s]
        rss = steps_ss[s]
        key_cs = (cs["band"][rcs].astype(int), cs["lev"][rcs].astype(int))
        d_cs = {kc: i for i, kc in enumerate(zip(*key_cs))}
        gpts = ss["gpt"][rss].astype(int)
        bands = SW_GPT_BAND[gpts - 1]
        levs = ss["lev"][rss].astype(int)
        j = np.array([d_cs[(b, l)] for b, l in zip(bands, levs)])
        tau_band[rss] = cs["tau"][rcs[j]]
        ssa_band[rss] = cs["ssa"][rcs[j]]
        asm_band[rss] = cs["asm"][rcs[j]]
        fsf_band[rss] = cs["fsfc"][rcs[j]]
        cld_val[rss] = cs["cld"][rcs[j]]
    cloudy = csto == 1
    chk.close("solver_sw tauc_sto", ss["tauc_sto"][rows][cloudy],
              tau_band[cloudy], rtol=1e-12, atol=1e-15,
              what="tauc_sto vs tau(band) when cld_sto=1")
    chk.close("solver_sw tauc_sto clear", ss["tauc_sto"][rows][~cloudy],
              np.zeros(int((~cloudy).sum())), rtol=0, atol=1e-15,
              what="tauc_sto=0 when cld_sto=0")
    chk.close("solver_sw ssac_sto", ss["ssac_sto"][rows][cloudy],
              ssa_band[cloudy], rtol=1e-12, atol=1e-15,
              what="ssac_sto vs ssa(band) when cld_sto=1")
    # clear-gpt ssac_sto=1 only applies on day steps: on night steps the
    # dump records the all-zero SW boundary (flush zero-fill)
    day_ss = np.concatenate([steps_ss[s] for s in day_steps])
    clear_day = ~cloudy[np.isin(rows, day_ss)]
    chk.close("solver_sw ssac_sto clear", ss["ssac_sto"][rows][np.isin(rows, day_ss)][clear_day],
              np.ones(int(clear_day.sum())), rtol=0, atol=1e-15,
              what="ssac_sto=1 when cld_sto=0 (day steps)")
    night_ss = np.concatenate([steps_ss[s] for s in night_steps])
    chk.close("solver_sw ssac_sto night zero", ss["ssac_sto"][night_ss],
              np.zeros(len(night_ss)), rtol=0, atol=1e-15,
              what="ssac_sto=0 on night steps (all-zero SW boundary)")
    chk.close("solver_sw asmc_sto", ss["asmc_sto"][rows][cloudy],
              asm_band[cloudy], rtol=1e-12, atol=1e-15,
              what="asmc_sto vs asm(band) when cld_sto=1")
    chk.close("solver_sw asmc_sto clear", ss["asmc_sto"][rows][~cloudy],
              np.zeros(int((~cloudy).sum())), rtol=0, atol=1e-15,
              what="asmc_sto=0 when cld_sto=0")
    chk.close("solver_sw fsfc_sto", ss["fsfc_sto"][rows][cloudy],
              fsf_band[cloudy], rtol=1e-12, atol=1e-15,
              what="fsfc_sto vs fsfc(band) when cld_sto=1")
    chk.close("solver_sw fsfc_sto clear", ss["fsfc_sto"][rows][~cloudy],
              np.zeros(int((~cloudy).sum())), rtol=0, atol=1e-15,
              what="fsfc_sto=0 when cld_sto=0")
    # cld_sto=1 only where cld>0
    chk.check("solver_sw cld_sto only where cld>0",
              np.all(cld_val[cloudy] > 0),
              f"{int((cld_val[cloudy] <= 0).sum())} cloudy gpts with cld=0")
    # Monte-Carlo consistency: mean(cld_sto) ~ cld per (step, lev)
    for s in RAD_STEPS:
        rss = steps_ss[s]
        rcs = steps_cs[s]
        for l in range(1, PVER + 1):
            m = np.mean(ss["cld_sto"][rss][ss["lev"][rss] == l])
            c = cs["cld"][rcs][cs["lev"][rcs] == l][0]
            if abs(m - c) > 0.25:
                chk.info_note(f"step {s} lev {l}: mean cld_sto {m:.3f} vs "
                              f"cld {c:.3f} (McICA sampling spread)")

    # ---- cam_rad_solver_lw.csv ---------------------------------------------
    _, sl = D["cam_rad_solver_lw.csv"]
    steps_sl, _ = by_step(sl)
    rows = step_rows(steps_sl, RAD_STEPS)
    gpt = sl["gpt"][rows].astype(int)
    csto = sl["cld_sto"][rows]
    chk.check("solver_lw gpt range", gpt.min() == 1 and gpt.max() == 140,
              f"gpt {gpt.min()}..{gpt.max()}")
    chk.check("solver_lw cld_sto in {0,1}",
              np.all((csto == 0) | (csto == 1)),
              f"values {np.unique(csto)}")
    taulw_band = np.zeros_like(csto)
    cld_val = np.zeros_like(csto)
    for s in RAD_STEPS:
        rcl = steps_cl[s]
        rsl = steps_sl[s]
        key_cl = (cl["band"][rcl].astype(int), cl["lev"][rcl].astype(int))
        d_cl = {kc: i for i, kc in enumerate(zip(*key_cl))}
        gpts = sl["gpt"][rsl].astype(int)
        bands = LW_GPT_BAND[gpts - 1]
        levs = sl["lev"][rsl].astype(int)
        j = np.array([d_cl[(b, l)] for b, l in zip(bands, levs)])
        taulw_band[rsl] = cl["tau_lw"][rcl[j]]
        cld_val[rsl] = cl["cld"][rcl[j]]
    cloudy = csto == 1
    chk.close("solver_lw tauc_sto", sl["tauc_sto"][rows][cloudy],
              taulw_band[cloudy], rtol=1e-12, atol=1e-15,
              what="tauc_sto vs tau_lw(band) when cld_sto=1")
    chk.close("solver_lw tauc_sto clear", sl["tauc_sto"][rows][~cloudy],
              np.zeros(int((~cloudy).sum())), rtol=0, atol=1e-15,
              what="tauc_sto=0 when cld_sto=0")
    chk.check("solver_lw cld_sto only where cld>0",
              np.all(cld_val[cloudy] > 0),
              f"{int((cld_val[cloudy] <= 0).sum())} cloudy gpts with cld=0")

    # ---- summary -----------------------------------------------------------
    print(f"checked {args.dumps_dir}")
    if chk.failures:
        print(f"FAIL: {len(chk.failures)} check(s) failed:")
        for f in chk.failures[:50]:
            print("  -", f)
        if len(chk.failures) > 50:
            print(f"  ... and {len(chk.failures) - 50} more")
        sys.exit(1)
    print("PASS: all hard checks passed")
    if chk.info_notes:
        print("info notes:")
        for n in chk.info_notes:
            print("  -", n)
    sys.exit(0)


if __name__ == "__main__":
    main()
