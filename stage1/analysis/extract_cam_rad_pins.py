#!/usr/bin/env python3
"""Extract stage-2 test tuples for stage2/cam_radiation.esm from the CAM
radiation (RRTMG) stage-1 dumps (run 12255018).

Reads the 8 cam_rad_*.csv files (cam_dump_mod.F90 cam_dump_rad_flush,
radiation.F90 radiation_tend, climate icall=0 call) and emits, as JSON on
stdout, the test tuples for the six models of cam_radiation.esm:

  CamRadState       rrtmg_state laws (rrtmg_state.F90:60-247): pmid branch
                    (k=1: 0.5*pint_hi(1), k>=2: 0.5*(pint_lo+pint_hi)),
                    the shared edges pint_hi(k)=pint_lo(k+1) and
                    T_ifc_hi(k)=T_ifc_lo(k+1) (k=1..32, the surface values
                    pint_hi(33)/T_ifc_hi(33) are inputs), h2ovmr =
                    (q/(1-q))*amdw, cfc22vmr = ccl4vmr = 0 (F90:222-223).
                    Dump columns: cam_rad_input.csv pmid, pint_lo, pint_hi,
                    T_ifc_lo, T_ifc_hi, q_H2O, h2ovmr, cfc22vmr, ccl4vmr.
  CamRadSolin       solin = sum(sfac_b*solar_band_irrad)*eccf*coszrs on day
                    steps, 0 at night (radsw.F90:395-397; solar_band_irrad
                    the radconstants.F90:34-40 constants). Dump column:
                    cam_rad_scal.csv solin (+ sfac_b1..b14, eccf, coszrs).
  CamRadSolverInputs  ssa = tau_w/tau (tau>0) else 1; asm = tau_w_g/tau_w,
                    fsfc = tau_w_f/tau_w (tau_w>0) else 0 (radsw.F90:437-475,
                    the max(x,1e-80) guards are vacuous for the dump's
                    nonzero ratios and the else branches pick the exact
                    zeros); cldfprime = max(cld, cldfsnow)
                    (radiation.F90:1012-1016). Dump columns:
                    cam_rad_cldopt_sw.csv tau, tau_w, tau_w_g, tau_w_f, ssa,
                    asm, fsfc; cam_rad_input.csv cld, cldfsnow (rows
                    k=2..33 = CAM lev 1..32).
  CamRadMcicaSw     tauc_sto = cld_sto*tau(band(gpt)), ssac_sto =
                    cld_sto*ssa(band(gpt)) + (1-cld_sto), asmc_sto =
                    cld_sto*asm(band(gpt)), fsfc_sto = cld_sto*fsfc(band(gpt))
                    (mcica_subcol_gen_sw.f90:439-444); band(gpt) the SW
                    k-distribution map (rrtmg_sw_init.f90:266-280, ngptsw=112,
                    gpt counts 6,12,8,8,10,10,2,10,8,6,6,8,6,12). Dump
                    columns: cam_rad_solver_sw.csv cld_sto, tauc_sto,
                    ssac_sto, asmc_sto, fsfc_sto joined to the cldopt bands.
  CamRadMcicaLw     tauc_sto = cld_sto*tau_lw(band(gpt))
                    (mcica_subcol_gen_lw.f90:408-416); band(gpt) the LW
                    k-distribution map (rrtmg_lw_init.f90:308, ngptlw=140,
                    gpt counts 10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2).
                    Dump columns: cam_rad_solver_lw.csv joined to
                    cam_rad_cldopt_lw.csv.
  CamRadHeating     qrs = -g*(fns(k+1)-fns(k))/dp and qrl =
                    +g*(fnl(k+1)-fnl(k))/dp as flux divergences on the esd
                    column (esd lev 1 = surface = CAM lev 32; the law reads
                    qrs = g*D(fns), qrl = -g*D(fnl) on the ze = -pint
                    coordinate), the scalar history fluxes fsnt = fns(1) etc.
                    as index picks and the column integrals sum(qrs dp/g) =
                    fsnt - fsns, sum(qrl dp/g) = flns - flnt (radsw.F90:565-575,
                    radlw.F90:237-245). Dump columns: cam_rad_out_ifc.csv
                    fds, fus, fdl, ful (esd-reversed), cam_rad_input.csv
                    pint_hi (ze = -pint_hi reversed), cam_rad_out.csv qrs/qrl,
                    cam_rad_scal.csv fsnt..flwds.

Every emitted assertion is copy-checked against the test's own input
overrides (|expected - input| must exceed rel*|expected| at the assertion
coordinate, the tolerance-aware form of the house rule): an assertion that
cannot fail is dropped with a warning. Steps pinned: day 3 (clear, solin
172.5), 25 (clear, surface absorption fsns = 104.5), 37 (thin cloud
cld ~ 0.12), 47 (dense cloud cld ~ 1.0) and night 21 (clear, LW only) for
the LW heating; the night SW boundary is all-zero (radsw.F90 early return +
flush zero-fill) and is never pinned.

Usage: python3 extract_cam_rad_pins.py [--dumps-dir DIR] [--render] > out.json
Without --render the output is the tuple summary (inputs, assertions) per
model; with --render the output is the ready-to-paste `tests` array entries
for stage2/cam_radiation.esm (the tests block skeleton and every
description are hand-authored; only the numeric tuples are extracted here,
as the house rules allow). The rendered numbers are pasted verbatim into
the hand-written tests blocks; this script never writes .esm files.
"""
import argparse
import json
import sys

import numpy as np

G = 9.80616          # gravit (physconst.F90)
AMDW = 1.607793      # amdw (rrtmg_state.F90:47)
SOLAR_REF = np.array([12.11, 20.3600000000001, 23.73, 22.43, 55.63, 102.93,
                      24.29, 345.74, 218.19, 347.20, 129.49, 50.15, 3.08,
                      12.89])  # W/m2 (radconstants.F90:34-40, sum 1368.22)
SW_GPT_COUNTS = [6, 12, 8, 8, 10, 10, 2, 10, 8, 6, 6, 8, 6, 12]   # ngptsw=112
LW_GPT_COUNTS = [10, 12, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2]  # 140
SW_GPT_BAND = np.repeat(np.arange(1, 15), SW_GPT_COUNTS)   # band of gpt 1..112
LW_GPT_BAND = np.repeat(np.arange(1, 17), LW_GPT_COUNTS)   # band of gpt 1..140

PVER = 32            # CAM levels
RAD_STEPS = list(range(3, 49, 2))  # iradsw=iradlw=2 radiation steps
DAY_STEPS = [3, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47]
NIGHT_STEPS = [s for s in RAD_STEPS if s not in DAY_STEPS]

REL = 1e-12          # assertion tolerance used by the tests


def read_csv(path):
    cols = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# description"):
                continue
            if line.startswith("# columns"):
                cols = line[len("# columns ="):].strip().split(",")
                continue
            rows.append([x.strip() for x in line.split(",")])
    if cols is None:
        raise ValueError(f"{path}: no '# columns' header")
    return {c: np.array([float(r[i]) for r in rows]) for i, c in enumerate(cols)}


def step_rows(data, s):
    return np.where(data["step"] == s)[0]


def esd(v):
    """Reverse a CAM top-down column (lev 1..32 or ifc 1..33) to esd order
    (lev 1 = surface)."""
    return v[::-1]


class Tuples:
    def __init__(self):
        self.models = {}
        self.warnings = []

    def add_test(self, model, test_id, description, inputs, assertions):
        """inputs: {name: (value, shape_or_None)}; assertions: list of
        (variable, coords_or_None, expected). Copy-checks each assertion."""
        keep = []
        for variable, coords, expected in assertions:
            bad = False
            for name, (val, shape) in inputs.items():
                if shape is None:
                    if coords is None and abs(val - expected) <= REL * abs(expected):
                        bad = True
                        self.warnings.append(
                            f"{model}/{test_id}: {variable}={expected:.17e} is a copy "
                            f"of scalar input {name}={val:.17e}")
                elif coords is not None and set(shape) == set(coords.keys()):
                    idx = tuple(coords[ax] - 1 for ax in shape)
                    if abs(val[idx] - expected) <= REL * abs(expected):
                        bad = True
                        self.warnings.append(
                            f"{model}/{test_id}: {variable}[{coords}]={expected:.17e} "
                            f"is a copy of input {name}[{coords}]={val[idx]:.17e}")
            if not bad:
                keep.append((variable, coords, expected))
        if keep:
            self.models.setdefault(model, {})[test_id] = {
                "description": description,
                "inputs": {n: (v.tolist() if s is not None else v)
                           for n, (v, s) in inputs.items()},
                "assertions": keep,
            }
        else:
            self.warnings.append(f"{model}/{test_id}: ALL assertions dropped")


def render(tuples):
    """Render the tuple summary as the ready-to-paste tests entries:
    [{"id", "description", "time_span", "initial_conditions",
      "parameter_overrides", "assertions"}], one per model test."""
    out = []
    for model, tests in tuples.items():
        for tid, t in tests.items():
            out.append({
                "id": tid,
                "description": t["description"],
                "time_span": {"start": 0.0, "end": 0.0},
                "initial_conditions": {},
                "parameter_overrides": t["inputs"],
                "assertions": [
                    {"variable": v, "time": 0.0,
                     **({} if c is None else {"coords": c}),
                     "expected": e, "tolerance": {"rel": REL}}
                    for v, c, e in t["assertions"]],
            })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dumps-dir", default="stage1/dumps/cam")
    ap.add_argument("--render", action="store_true",
                    help="render the ready-to-paste tests entries")
    args = ap.parse_args()
    d = args.dumps_dir
    inp = read_csv(f"{d}/cam_rad_input.csv")
    sc = read_csv(f"{d}/cam_rad_scal.csv")
    out = read_csv(f"{d}/cam_rad_out.csv")
    oi = read_csv(f"{d}/cam_rad_out_ifc.csv")
    cs = read_csv(f"{d}/cam_rad_cldopt_sw.csv")
    cl = read_csv(f"{d}/cam_rad_cldopt_lw.csv")
    ss = read_csv(f"{d}/cam_rad_solver_sw.csv")
    sl = read_csv(f"{d}/cam_rad_solver_lw.csv")

    t = Tuples()

    # ---------------- CamRadState ----------------
    for s in (3, 41, 47):
        r = step_rows(inp, s)
        pint_lo = inp["pint_lo"][r]           # k = 1..33
        pint_hi = inp["pint_hi"][r]
        tifc_lo = inp["T_ifc_lo"][r]
        tifc_hi = inp["T_ifc_hi"][r]
        q = inp["q_H2O"][r]
        inputs = {
            "pint_lo": (pint_lo, ("k",)),
            "pint_sfc": (pint_hi[32], None),
            "T_ifc_lo": (tifc_lo, ("k",)),
            "T_ifc_sfc": (tifc_hi[32], None),
            "q_H2O": (q, ("k",)),
        }
        asserts = []
        for k in (1, 2, 17, 33):
            asserts.append(("pmid", {"k": k}, inp["pmid"][r][k - 1]))
        for k in (1, 2, 17, 32):
            asserts.append(("pint_hi", {"k": k}, pint_hi[k - 1]))
        for k in (17, 32):
            asserts.append(("T_ifc_hi", {"k": k}, tifc_hi[k - 1]))
        for k in (1, 17, 33):
            asserts.append(("h2ovmr", {"k": k}, inp["h2ovmr"][r][k - 1]))
        asserts.append(("cfc22vmr", {"k": 17}, 0.0))
        asserts.append(("ccl4vmr", {"k": 17}, 0.0))
        t.add_test(
            "CamRadState", f"step{s}",
            f"cam_rad_input.csv rows (step {s}, k = 1..33, RRTMG layer "
            f"top-down): the rrtmg_state_create laws of rrtmg_state.F90:60-247 "
            f"pinned at the extra layer (k=1), top CAM layer (k=2), mid column "
            f"(k=17) and surface layer (k=33): pmid(1) = 0.5*pint_hi(1), "
            f"pmid(k) = 0.5*(pint_lo+pint_hi) (F90:126,140), the shared edges "
            f"pint_hi(k) = pint_lo(k+1) and T_ifc_hi(k) = T_ifc_lo(k+1) "
            f"(k <= 32; the surface edges pint_hi(33)/T_ifc_hi(33) are "
            f"inputs), h2ovmr = (q/(1-q))*amdw (F90:47,198) and the "
            f"cfc22vmr = ccl4vmr = 0 configuration zeroing (F90:222-223).",
            inputs, asserts)

    # ---------------- CamRadSolin ----------------
    for s in (3, 25, 37, 47):
        r = step_rows(sc, s)[0]
        sfac = np.array([sc[f"sfac_b{i}"][r] for i in range(1, 15)])
        inputs = {"sfac": (sfac, ("band",)), "eccf": (sc["eccf"][r], None),
                  "coszrs": (sc["coszrs"][r], None)}
        t.add_test(
            "CamRadSolin", f"step{s}",
            f"cam_rad_scal.csv row (step {s}, day, coszrs = "
            f"{sc['coszrs'][r]:.4f}): solin = sum_b sfac_b*solar_band_irrad_b "
            f"* eccf * coszrs (radsw.F90:395-397) with the solar band "
            f"irradiances of radconstants.F90:34-40 (sum = 1368.22 W/m2); "
            f"the band sum is a faq sum over the band axis.",
            inputs, [("solin", None, sc["solin"][r])])

    # ---------------- CamRadSolverInputs ----------------
    for s in (3, 33, 47):
        r = step_rows(cs, s)
        band = cs["band"][r].astype(int)
        lev = cs["lev"][r].astype(int)
        key = {(b, l): i for i, (b, l) in enumerate(zip(band, lev))}
        tau = np.zeros((14, 32)); tau_w = np.zeros((14, 32))
        tau_wg = np.zeros((14, 32)); tau_wf = np.zeros((14, 32))
        ssa = np.zeros((14, 32)); asm = np.zeros((14, 32)); fsfc = np.zeros((14, 32))
        for b in range(1, 15):
            for l in range(1, 33):
                i = key[(b, l)]
                tau[b - 1, l - 1] = cs["tau"][r][i]
                tau_w[b - 1, l - 1] = cs["tau_w"][r][i]
                tau_wg[b - 1, l - 1] = cs["tau_w_g"][r][i]
                tau_wf[b - 1, l - 1] = cs["tau_w_f"][r][i]
                ssa[b - 1, l - 1] = cs["ssa"][r][i]
                asm[b - 1, l - 1] = cs["asm"][r][i]
                fsfc[b - 1, l - 1] = cs["fsfc"][r][i]
        # cld / cldfsnow per CAM lev from cam_rad_input rows k = 2..33
        ri = step_rows(inp, s)
        cld = inp["cld"][ri][1:]
        cldfsnow = inp["cldfsnow"][ri][1:]
        inputs = {"cld": (cld, ("lev",)), "cldfsnow": (cldfsnow, ("lev",)),
                  "tau": (tau, ("band", "lev")), "tau_w": (tau_w, ("band", "lev")),
                  "tau_w_g": (tau_wg, ("band", "lev")),
                  "tau_w_f": (tau_wf, ("band", "lev"))}
        asserts = []
        if s == 3:
            # clear step: the tau = 0 else branch (ssa = 1); asm/fsfc = 0
            # there are copies of the zero tau_w_g/tau_w_f inputs and are not
            # emitted (the law's else branch is pinned only by ssa)
            for b in (1, 8, 14):
                asserts.append(("ssa", {"band": b, "lev": 1}, 1.0))
            asserts.append(("ssa", {"band": 8, "lev": 17}, 1.0))
        else:
            levs = {33: 15, 47: 20}[s]   # the cloudy levels of the step
            for (b, l) in ((1, levs), (8, levs), (14, levs)):
                asserts.append(("ssa", {"band": b, "lev": l}, ssa[b - 1, l - 1]))
                asserts.append(("asm", {"band": b, "lev": l}, asm[b - 1, l - 1]))
                asserts.append(("fsfc", {"band": b, "lev": l}, fsfc[b - 1, l - 1]))
        t.add_test(
            "CamRadSolverInputs", f"step{s}",
            f"cam_rad_cldopt_sw.csv rows (step {s}, band x lev) and "
            f"cam_rad_input.csv cld/cldfsnow (rows k = 2..33 = CAM lev "
            f"1..32): the solver-input derivations of radsw.F90:437-475 - "
            f"ssa = tau_w/tau when tau > 0 else 1, asm = tau_w_g/tau_w and "
            f"fsfc = tau_w_f/tau_w when tau_w > 0 else 0 (the Fortran "
            f"max(x,1e-80) guards are vacuous for the dump's nonzero ratios) "
            f"- and cldfprime = max(cld, cldfsnow) (radiation.F90:1012-1016; "
            f"cldfsnow = cld bitwise in the window, so cldfprime is a copy of "
            f"an input and is not pinned).",
            inputs, asserts)

    # ---------------- CamRadMcicaSw ----------------
    for s, lvl in ((37, 15), (47, 20)):
        rss = step_rows(ss, s)
        gpt = ss["gpt"][rss].astype(int)
        lev = ss["lev"][rss].astype(int)
        key = {(g, l): i for i, (g, l) in enumerate(zip(gpt, lev))}
        cld_sto = np.zeros((112, 32))
        for g in range(1, 113):
            for l in range(1, 33):
                cld_sto[g - 1, l - 1] = ss["cld_sto"][rss][key[(g, l)]]
        rcs = step_rows(cs, s)
        ckey = {(int(cs["band"][rcs][i]), int(cs["lev"][rcs][i])): i
                for i in range(len(rcs))}
        tau = np.zeros((14, 32)); ssa = np.zeros((14, 32))
        asm = np.zeros((14, 32)); fsfc = np.zeros((14, 32))
        for b in range(1, 15):
            for l in range(1, 33):
                i = ckey[(b, l)]
                tau[b - 1, l - 1] = cs["tau"][rcs][i]
                ssa[b - 1, l - 1] = cs["ssa"][rcs][i]
                asm[b - 1, l - 1] = cs["asm"][rcs][i]
                fsfc[b - 1, l - 1] = cs["fsfc"][rcs][i]
        inputs = {"cld_sto": (cld_sto, ("gpt", "lev")),
                  "tau": (tau, ("band", "lev")), "ssa": (ssa, ("band", "lev")),
                  "asm": (asm, ("band", "lev")), "fsfc": (fsfc, ("band", "lev")),
                  "band_of_gpt": (SW_GPT_BAND, ("gpt",))}
        asserts = []
        # first cloudy subcolumn of the chosen bands at the step's cloudy
        # level (gpt -> band via the SW k-distribution map)
        bands_want = {37: [2, 6, 14], 47: [2, 3, 6, 14]}[s]
        for b in bands_want:
            g = int(np.where((SW_GPT_BAND == b) & (cld_sto[:, lvl - 1] == 1))[0][0]) + 1
            i = key[(g, lvl)]
            asserts.append(("tauc_sto", {"gpt": g, "lev": lvl},
                            ss["tauc_sto"][rss][i]))
            asserts.append(("ssac_sto", {"gpt": g, "lev": lvl},
                            ss["ssac_sto"][rss][i]))
            asserts.append(("asmc_sto", {"gpt": g, "lev": lvl},
                            ss["asmc_sto"][rss][i]))
            asserts.append(("fsfc_sto", {"gpt": g, "lev": lvl},
                            ss["fsfc_sto"][rss][i]))
        # the cld_sto = 0 fill branch (ssac_sto = 1) at the cloudy level
        g0 = int(np.where(cld_sto[:, lvl - 1] == 0)[0][0]) + 1
        asserts.append(("ssac_sto", {"gpt": g0, "lev": lvl}, 1.0))
        t.add_test(
            "CamRadMcicaSw", f"step{s}",
            f"cam_rad_solver_sw.csv rows (step {s}, day, gpt x lev) joined to "
            f"the cam_rad_cldopt_sw.csv band optics: the McICA subcolumn "
            f"masking of mcica_subcol_gen_sw.f90:439-444 - tauc_sto = "
            f"cld_sto*tau(band(gpt)), ssac_sto = cld_sto*ssa(band(gpt)) + "
            f"(1-cld_sto), asmc_sto = cld_sto*asm(band(gpt)), fsfc_sto = "
            f"cld_sto*fsfc(band(gpt)) - with band(gpt) the SW k-distribution "
            f"map (rrtmg_sw_init.f90:266-280, counts 6,12,8,8,10,10,2,10,"
            f"8,6,6,8,6,12). The stochastic cld_sto itself is an input: it is "
            f"generated by the deterministic RNG of mcica_subcol_gen_sw.f90 "
            f"(permuteseed = 1, single column), which is out of scope.",
            inputs, asserts)

    # ---------------- CamRadMcicaLw ----------------
    for s, lvl in ((41, 20), (47, 20)):
        rsl = step_rows(sl, s)
        gpt = sl["gpt"][rsl].astype(int)
        lev = sl["lev"][rsl].astype(int)
        key = {(g, l): i for i, (g, l) in enumerate(zip(gpt, lev))}
        cld_sto = np.zeros((140, 32))
        for g in range(1, 141):
            for l in range(1, 33):
                cld_sto[g - 1, l - 1] = sl["cld_sto"][rsl][key[(g, l)]]
        rcl = step_rows(cl, s)
        ckey = {(int(cl["band"][rcl][i]), int(cl["lev"][rcl][i])): i
                for i in range(len(rcl))}
        taulw = np.zeros((16, 32))
        for b in range(1, 17):
            for l in range(1, 33):
                taulw[b - 1, l - 1] = cl["tau_lw"][rcl][ckey[(b, l)]]
        inputs = {"cld_sto": (cld_sto, ("gpt_lw", "lev")),
                  "tau_lw": (taulw, ("band_lw", "lev")),
                  "band_of_gpt_lw": (LW_GPT_BAND, ("gpt_lw",))}
        asserts = []
        # first cloudy subcolumn of the chosen bands at the cloudy level
        bands_want = {41: [4, 8, 14], 47: [4, 8, 14]}[s]
        for b in bands_want:
            g = int(np.where((LW_GPT_BAND == b) & (cld_sto[:, lvl - 1] == 1))[0][0]) + 1
            i = key[(g, lvl)]
            asserts.append(("tauc_sto_lw", {"gpt_lw": g, "lev": lvl},
                            sl["tauc_sto"][rsl][i]))
        t.add_test(
            "CamRadMcicaLw", f"step{s}",
            f"cam_rad_solver_lw.csv rows (step {s}, gpt x lev) joined to the "
            f"cam_rad_cldopt_lw.csv band optics: the McICA longwave subcolumn "
            f"masking of mcica_subcol_gen_lw.f90:408-416 - tauc_sto = "
            f"cld_sto*tau_lw(band(gpt)) - with band(gpt) the LW "
            f"k-distribution map (rrtmg_lw_init.f90:308, counts 10,12,16,14,"
            f"16,8,12,8,12,6,8,8,4,2,2,2). The stochastic cld_sto is an input "
            f"(deterministic RNG, permuteseed = 150).",
            inputs, asserts)

    # ---------------- CamRadHeating ----------------
    for s in (3, 25, 37, 47, 21):
        r = step_rows(oi, s)
        fds = esd(oi["fds"][r])     # esd ifc 1 = surface = CAM ifc 33
        fus = esd(oi["fus"][r])
        fdl = esd(oi["fdl"][r])
        ful = esd(oi["ful"][r])
        ri = step_rows(inp, s)
        ze = -esd(inp["pint_hi"][ri])   # esd coordinate = -pint, increasing up
        rs = step_rows(sc, s)[0]
        rq = step_rows(out, s)
        lev = out["lev"][rq].astype(int)
        qkey = {l: i for i, l in enumerate(lev)}
        inputs = {"fds": (fds, ("lev_nodes",)), "fus": (fus, ("lev_nodes",)),
                  "fdl": (fdl, ("lev_nodes",)), "ful": (ful, ("lev_nodes",)),
                  "ze": (ze, ("lev_nodes",))}
        asserts = []
        day = s in DAY_STEPS
        for l in (1, 17, 32):        # esd lev = CAM lev 33-l
            i = qkey[33 - l]
            if day:
                asserts.append(("qrs", {"lev": l}, out["qrs"][rq][i]))
            asserts.append(("qrl", {"lev": l}, out["qrl"][rq][i]))
        for i in (1, 17, 33):        # esd ifc = CAM ifc 34-i
            fnl_i = oi["fnl"][r][33 - i]
            asserts.append(("fnl", {"lev_nodes": i}, fnl_i))
            if day:
                fns_i = oi["fns"][r][33 - i]
                # skip the surface point where fns = 0 exactly (albedo = 1
                # on steps 3/41/47, an identically-zero pin)
                if i == 1 and fns_i == 0.0:
                    continue
                asserts.append(("fns", {"lev_nodes": i}, fns_i))
        if day:
            asserts.append(("fsnt", None, sc["fsnt"][rs]))
            asserts.append(("fsds", None, sc["fsds"][rs]))
            asserts.append(("sw_net", None, sc["fsnt"][rs] - sc["fsns"][rs]))
            if abs(sc["fsns"][rs]) > 1e-6:
                asserts.append(("fsns", None, sc["fsns"][rs]))
        asserts.append(("flnt", None, sc["flnt"][rs]))
        asserts.append(("flns", None, sc["flns"][rs]))
        asserts.append(("flut", None, sc["flut"][rs]))
        asserts.append(("flwds", None, sc["flwds"][rs]))
        asserts.append(("lw_net", None, sc["flns"][rs] - sc["flnt"][rs]))
        t.add_test(
            "CamRadHeating", f"step{s}",
            f"cam_rad_out_ifc.csv rows (step {s}, ifc 1..33, esd-reversed), "
            f"cam_rad_input.csv pint_hi (ze = -pint_hi reversed) and "
            f"cam_rad_out.csv / cam_rad_scal.csv: the RRTMG heating as the "
            f"layer flux divergence on the esd column - qrs = -g*(fns(k+1)-"
            f"fns(k))/dp reads qrs = g*D(fns) and qrl = +g*(fnl(k+1)-"
            f"fnl(k))/dp reads qrl = -g*D(fnl) on the ze coordinate "
            f"(radsw.F90:573-574, radlw.F90:244-245, verified to 3e-16 W/kg "
            f"by the stage-1 validator) - fns = fds - fus and fnl = ful - fdl "
            f"(radsw.F90:565-567, radlw.F90:237-238), the history scalars as "
            f"index picks (fsnt = fns(1) = fns[lev_nodes 33] etc., "
            f"radsw.F90:552-562, radlw.F90:214-228) and the column integrals "
            f"sum(qrs dp/g) = fsnt - fsns, sum(qrl dp/g) = flns - flnt.",
            inputs, asserts)

    for w in t.warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    if args.render:
        print(json.dumps(render(t.models), indent=1))
    else:
        print(json.dumps({"models": t.models}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
