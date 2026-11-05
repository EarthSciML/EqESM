#!/usr/bin/env python3
"""Extract the W3 soiltemp test tuples (column 1, step 1) from the stage-1 dump.

The tests block of stage2/clm_soiltemp.esm pins the SoilTemperature PDE model
against dumped CLM state. Its numeric tuples (parameter_overrides and assertion
expected values) are extracted from the dump by this script -- the sanctioned
"extract numeric test tuples from C++ dumps into an existing hand-authored
tests block" workflow (AGENTS.md). The hand-authored structure lives in the
.esm file; this script only produces the numbers.

The derived values reproduce, at full double precision, what the .esm model
computes (its equations are the same math, see README-clm-soiltemp.md):

  * ze:      interface coordinates, ze[n] = -zi_bot(26-n), ze[35] = -zi_top(-8)
             (the ESD lev axis runs bottom-up; lev n = CLM j = 26-n)
  * K:       interface conductivities, K[n] = dumped tk(26-n); K[1] = 0 (bottom
             boundary, tk(nlevgrnd) = 0), K[35] = 0 dummy (no interface above
             the top snow layer)
  * cv, fn:  dumped (fn is the ENTRY-state interface flux -G(t_entry) written
             at the exit stage, verified numerically to 1e-13)
  * dTdt:    the instantaneous PDE tendency (dz*D(K*D(T,lev),lev) + sabg)/cv
             evaluated at the entry and at the exit state -- derived values,
             computed here from the dump; K is frozen from t_entry as in CLM
  * cn_balance: the derived Crank-Nicolson identity
             0.5*(F(t_entry) + F(t_exit)) - cv*(t_exit - t_entry)/dt with
             F(u) = dz*D(K*D(u,lev),lev) + sabg and the capr top-row measure
             M = 0.5*(ze[35]-zc[34] + capr*(ze[35]-zc[33])); the interior and
             bottom rows evaluate to the CLM solver's round-off (max |r| =
             6.3e-10 W/m2 for lev 1..33), the top row deviates by the
             fully-implicit dhsdT term of SetMatrix_Snow and is not pinned

Usage:

    python3 stage1/analysis/extract_soiltemp_testdata.py            # print JSON
    python3 stage1/analysis/extract_soiltemp_testdata.py --check    # verify the
                             # tuples embedded in stage2/clm_soiltemp.esm match
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dumpcsv import load, select

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
# The dumps are gitignored and live in the main checkout (stage1/dumps/); the
# worktree copy of the same checkout also works when present.
_DUMP_CANDIDATES = [
    REPO / "stage1" / "dumps" / "clm" / "clm_soiltemp.csv",
    REPO.parent / "EqESM" / "stage1" / "dumps" / "clm" / "clm_soiltemp.csv",
]
DUMP = next((p for p in _DUMP_CANDIDATES if p.exists()), _DUMP_CANDIDATES[0])
MODEL = REPO / "stage2" / "clm_soiltemp.esm"
NLEV = 34
CAPR = 0.34
DT = 1800.0  # CLM land-model time step of the traced run (dtime)


def rev(x):
    """CLM j (-8..25, top down) -> ESD lev (1..34, bottom up): lev n = 26-j."""
    return [x[33 - k] for k in range(34)]


def load_dump():
    meta, rows = load(str(DUMP))
    ent = select(rows, step=1, c=1, stage="entry")
    exc = select(rows, step=1, c=1, stage="exit")
    assert len(ent) == NLEV and len(exc) == NLEV, "dump must have 34 rows per stage"
    t_ent = rev([r["t_soisno"] for r in ent])
    t_exc = rev([r["t_soisno"] for r in exc])
    dz = rev([r["dz"] for r in exc])
    cv = rev([r["cv"] for r in exc])
    tk = rev([r["tk"] for r in exc])
    fn = rev([r["fn"] for r in exc])
    h2osoi_liq = rev([r["h2osoi_liq"] for r in ent])
    h2osoi_ice = rev([r["h2osoi_ice"] for r in ent])
    watsat = rev([r["watsat"] if r["watsat"] == r["watsat"] else 0.0 for r in ent])
    tkmg = rev([r["tkmg"] if r["tkmg"] == r["tkmg"] else 0.0 for r in ent])
    tkdry = rev([r["tkdry"] if r["tkdry"] == r["tkdry"] else 0.0 for r in ent])
    csol = rev([r["csol"] if r["csol"] == r["csol"] else 0.0 for r in ent])
    sabg = rev([r["sabg_lyr_col"] if r["sabg_lyr_col"] == r["sabg_lyr_col"] else 0.0 for r in exc])
    sabg[33] = 0.0  # top snow layer (j = -8): its 57.79 W/m2 absorption is inside hs_top_snow
    hs = exc[0]["hs_top_snow"]
    dhsdT = exc[0]["dhsdT"]
    eflx_bot = exc[0]["eflx_bot"]
    frac_sno = exc[0]["frac_sno_fluxes"]
    zi_bot = rev([r["zi_bot"] for r in exc])
    ze = [-zi_bot[n] for n in range(NLEV)] + [-ent[0]["zi_top"]]
    zc = [(ze[n] + ze[n + 1]) / 2 for n in range(NLEV)]
    M = 0.5 * (ze[34] - zc[33] + CAPR * (ze[34] - zc[32]))
    K = [tk[n] for n in range(NLEV)] + [0.0]
    return dict(t_ent=t_ent, t_exc=t_exc, dz=dz, cv=cv, K=K, fn=fn,
                h2osoi_liq=h2osoi_liq, h2osoi_ice=h2osoi_ice, watsat=watsat,
                tkmg=tkmg, tkdry=tkdry, csol=csol, sabg=sabg, hs=hs, dhsdT=dhsdT,
                eflx_bot=eflx_bot, frac_sno=frac_sno, ze=ze, zc=zc, M=M)


def compute(d):
    t_ent, t_exc, dz, cv, K = d["t_ent"], d["t_exc"], d["dz"], d["cv"], d["K"]
    zc, ze, M = d["zc"], d["ze"], d["M"]
    hs, dhsdT, eflx_bot, sabg = d["hs"], d["dhsdT"], d["eflx_bot"], d["sabg"]

    def G(u, n):
        """Interface flux at node n (between lev n-1 and n), +up."""
        return K[n - 1] * (u[n - 1] - u[n - 2]) / (zc[n - 1] - zc[n - 2])

    def kdudz_top(u):
        return hs + dhsdT * (u[33] - t_ent[33])

    def L(u):
        """Lowered divergence D(K*D(u,lev),lev): capr top-row measure."""
        out = [0.0] * NLEV
        for n in range(NLEV):
            n1 = n + 1
            if n1 == 1:
                out[n] = (G(u, 2) - eflx_bot) / dz[0]
            elif n1 == NLEV:
                out[n] = (kdudz_top(u) - G(u, NLEV)) / M
            else:
                out[n] = (G(u, n1 + 1) - G(u, n1)) / dz[n]
        return out

    def F(u):
        """Per-layer heat input (W/m2): dz*L + sabg."""
        return [dz[n] * L(u)[n] + sabg[n] for n in range(NLEV)]

    dTdt_ent = [F(t_ent)[n] / cv[n] for n in range(NLEV)]
    dTdt_exc = [F(t_exc)[n] / cv[n] for n in range(NLEV)]
    cn_res = [0.5 * (F(t_ent)[n] + F(t_exc)[n]) - cv[n] * (t_exc[n] - t_ent[n]) / DT
              for n in range(NLEV)]
    return dict(dTdt_ent=dTdt_ent, dTdt_exc=dTdt_exc, cn_res=cn_res,
                kdudz_top_exit=kdudz_top(t_exc))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="verify the tuples embedded in stage2/clm_soiltemp.esm")
    ap.add_argument("--dump", default=str(DUMP), help="path to clm_soiltemp.csv")
    a = ap.parse_args(argv)

    d = load_dump()
    c = compute(d)

    if not a.check:
        out = {
            "ze": d["ze"], "dz": d["dz"], "t_entry": d["t_ent"], "t_exit": d["t_exc"],
            "h2osoi_liq": d["h2osoi_liq"], "h2osoi_ice": d["h2osoi_ice"],
            "watsat": d["watsat"], "tkmg": d["tkmg"], "tkdry": d["tkdry"],
            "csol": d["csol"], "sabg": d["sabg"], "hs": d["hs"], "dhsdT": d["dhsdT"],
            "eflx_bot": d["eflx_bot"], "frac_sno_fluxes": d["frac_sno"],
            "K": d["K"], "cv": d["cv"], "fn": d["fn"],
            "dTdt_entry": c["dTdt_ent"], "dTdt_exit": c["dTdt_exc"],
            "cn_residual": c["cn_res"],
        }
        print(json.dumps(out, indent=1, allow_nan=False))
        print(f"# M = {d['M']!r}  A = M/dz[34] = {d['M']/d['dz'][33]!r}", file=sys.stderr)
        print(f"# cn residual: max|r| rows 1..33 = "
              f"{max(abs(r) for r in c['cn_res'][:33])!r} (top row = {c['cn_res'][33]!r})",
              file=sys.stderr)
        return 0

    # --check: compare against the tuples embedded in the .esm tests block
    doc = json.loads(MODEL.read_text())
    tests = doc["models"]["SoilTemperature"]["tests"]
    by_id = {t["id"]: t for t in tests}
    t1, t2, t3 = by_id["clm_soiltemp_entry"], by_id["clm_soiltemp_exit"], by_id["clm_soiltemp_cn_balance"]

    def coords_vals(assertions):
        out = {}
        for a in assertions:
            key = (a["variable"], tuple(sorted(a.get("coords", {}).items())))
            out[key] = a["expected"]
        return out

    a1 = coords_vals(t1["assertions"])
    a2 = coords_vals(t2["assertions"])
    a3 = coords_vals(t3["assertions"])
    ok = True

    def cmp(name, got, expected):
        nonlocal ok
        if got != expected:
            ok = False
            print(f"  MISMATCH {name}: esm {expected!r} vs dump-derived {got!r}")

    ov = t1["parameter_overrides"]
    for key, got in [("ze", d["ze"]), ("dz", d["dz"]), ("t_entry", d["t_ent"]),
                     ("t_exit", d["t_exc"]), ("h2osoi_liq", d["h2osoi_liq"]),
                     ("h2osoi_ice", d["h2osoi_ice"]), ("watsat", d["watsat"]),
                     ("tkmg", d["tkmg"]), ("tkdry", d["tkdry"]), ("csol", d["csol"]),
                     ("sabg", d["sabg"])]:
        cmp(f"override {key}", got, ov[key])
    for key in ("hs", "dhsdT", "eflx_bot"):
        cmp(f"override {key}", d[key], ov[key])
    cmp("override frac_sno_fluxes", d["frac_sno"], ov["frac_sno_fluxes"])
    assert t1["initial_conditions"]["T"] == d["t_ent"]
    assert t2["initial_conditions"]["T"] == d["t_exc"]
    for n in range(1, 36):
        cmp(f"K[lev_nodes {n}]", d["K"][n - 1], a1[("K", (("lev_nodes", n),))])
    for n in range(1, 35):
        cmp(f"cv[lev {n}]", d["cv"][n - 1], a1[("cv", (("lev", n),))])
        cmp(f"fn_esd[lev {n}]", d["fn"][n - 1], a1[("fn_esd", (("lev", n),))])
        cmp(f"dTdt entry[lev {n}]", c["dTdt_ent"][n - 1], a1[("dTdt", (("lev", n),))])
        cmp(f"dTdt exit[lev {n}]", c["dTdt_exc"][n - 1], a2[("dTdt", (("lev", n),))])
    cmp("kdudz_top exit", c["kdudz_top_exit"], a2[("kdudz_top", ())])
    for n in range(1, 34):
        cmp(f"cn_balance[lev {n}] expected", 0.0, a3[("cn_balance", (("lev", n),))])
    print("OK: every tuple embedded in stage2/clm_soiltemp.esm matches the dump"
          if ok else "MISMATCHES FOUND")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
