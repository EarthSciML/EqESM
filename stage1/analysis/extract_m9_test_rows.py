#!/usr/bin/env python3
"""Extract and verify the numeric test tuples of the M9 datm->CLM subassembly
(stage2/clm_atm_forcing.esm and stage2/m9_datm_clm.esm) from the stage-1
dumps, and compute the two values that are NOT dump columns.

The overlapping (cell, step) rows of the I2000Clm60Sp reference run:

  datm stream  rows (step t,   cell n=1..5) of cdeps/datm_gswp3_streams.csv
  datm export  rows (step t,   cell n=1..5) of cdeps/datm_gswp3_exports.csv
  CLM forcing  rows (step t+1, col  c=1..5) of clm/clm_bareground.csv,
               clm/clm_canopyflux.csv (stage=ent rows), clm/clm_surfrad.csv

with the one-coupling-interval offset that is the CESM data-atmosphere
sequence itself: within each coupling interval the data atmosphere runs
AFTER the land, so the land is forced by the atmosphere of the previous
interval.  (Verified: every forc_* dump value at (t+1, c) equals its
one-step-earlier datm export at (t, n) bit-exact; the five traced cells are
value-identical, so the tests pin cell 1 and cover all of them.)

Laws verified here for all 12 overlapping steps (all to round-off):

  wind split (datm_datamode_clmncep_mod.F90:402-403): Sa_u = Sa_v = wind/sqrt(2)
  import copies (lnd_import_export.F90:542-577): forc_t = Sa_tbot,
    forc_pbot = Sa_pbot, forc_u = Sa_u, forc_v = Sa_v, forc_q = Sa_shum,
    forc_solad_vis = Faxa_swvdr, forc_solad_nir = Faxa_swndr,
    forc_solai_vis = Faxa_swvdf, forc_solai_nir = Faxa_swndf
  rain/snow totals (lnd_import_export_utils.F90:85-86): forc_rain = rainc+rainl,
    forc_snow = snowc+snowl (both Faxa_rainc/Faxa_rainl are 0 on every traced
    step, so forc_rain = 0 and forc_snow = Faxa_snowc + Faxa_snowl)
  wind speed (lnd_import_export_utils.F90:77): forc_wind = sqrt(forc_u^2+forc_v^2)
    -- DERIVED (not a dump column)
  total solar (lnd_import_export_utils.F90:79-83): forc_solar = solad_vis +
    solad_nir + solai_vis + solai_nir -- DERIVED (not a dump column)
  density (lnd_import_export_utils.F90:72-73 / cdeps density template):
    forc_rho = (pbot - 0.378*vp)/(rair*t), vp = q*pbot/(0.622 + 0.378*q)

The .esm tests' parameter_overrides and assertions are then checked against
these dump-derived values at rel 1e-12, so every pinned number in the two
stage-2 files is traceable to a dump row (or to one of the two derived
values above, which are marked "derived (not a dump column)" in the test
descriptions and cite this script).

Usage:
  python3 stage1/analysis/extract_m9_test_rows.py [dumps_dir] [esm_dir]
where dumps_dir defaults to <main checkout>/stage1/dumps (dumps are
gitignored and live in the main checkout, not in the worktree; the
environment variable EQESM_MAIN overrides the main-checkout location) and
esm_dir defaults to <this script's repo>/stage2.
"""
import csv, json, math, os, sys

EQESM_MAIN = os.environ.get("EQESM_MAIN", "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM")
DUMPS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(EQESM_MAIN, "stage1/dumps")
ESM_DIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "stage2")

TKFRZ = 273.15
RDAIR = 287.0423113650487
EPS = 0.622
EPS1 = 1.0 - EPS
STEPS = range(1, 13)          # CLM dump steps s = 1..12 (South Pole cell)
CELL = 1                      # traced cell: the five traced cells are value-identical


def load_dump(path):
    cols, data = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("# columns"):
                cols = [c.strip() for c in line.split("=", 1)[1].strip().split(",")]
            elif not line.startswith("#"):
                data.append(line)
    return list(csv.DictReader(data, fieldnames=cols))


def get(rows, **keys):
    for r in rows:
        if all(r[k].strip() == str(v) for k, v in keys.items()):
            return r
    raise KeyError(keys)


def f(x):
    return float(x.strip())


def verify_chain():
    """Verify the Fortran chain on all 12 overlapping steps; return the
    (step -> {override values, expectation values, derived values}) table."""
    exp_rows = load_dump(os.path.join(DUMPS, "cdeps/datm_gswp3_exports.csv"))
    strm_rows = load_dump(os.path.join(DUMPS, "cdeps/datm_gswp3_streams.csv"))
    bg = load_dump(os.path.join(DUMPS, "clm/clm_bareground.csv"))
    cf = load_dump(os.path.join(DUMPS, "clm/clm_canopyflux.csv"))
    sr = load_dump(os.path.join(DUMPS, "clm/clm_surfrad.csv"))

    table = {}
    worst = 0.0
    for s in STEPS:
        t = s - 1  # datm step
        ex = get(exp_rows, step=t, n=CELL)
        st = get(strm_rows, step=t, n=CELL)
        b = get(bg, step=s, c=CELL)
        c = get(cf, step=s, stage="ent", c=CELL)
        r = get(sr, step=s, c=CELL)

        # datm datamode wind split: Sa_u = Sa_v = wind/sqrt(2)
        assert abs(f(ex["Sa_u"]) - f(st["strm_wind"]) / math.sqrt(2)) <= 1e-12 * f(st["strm_wind"]), s
        assert abs(f(ex["Sa_v"]) - f(ex["Sa_u"])) <= 1e-12 * f(ex["Sa_u"]), s

        # import copies: datm export (t) == CLM forcing input (t+1), bit-exact
        for excol, bgcol in [("Sa_tbot", "forc_t"), ("Sa_pbot", "forc_pbot"),
                             ("Sa_u", "forc_u"), ("Sa_v", "forc_v"),
                             ("Sa_shum", "forc_q")]:
            assert f(ex[excol]) == f(b[bgcol]), (s, excol, ex[excol], b[bgcol])
        for excol, srcol in [("Faxa_swvdr", "forc_solad_vis"), ("Faxa_swndr", "forc_solad_nir"),
                             ("Faxa_swvdf", "forc_solai_vis"), ("Faxa_swndf", "forc_solai_nir")]:
            assert f(ex[excol]) == f(r[srcol]), (s, excol, ex[excol], r[srcol])

        # rain/snow totals; rain is 0 on every traced step
        rain = f(ex["Faxa_rainc"]) + f(ex["Faxa_rainl"])
        snow = f(ex["Faxa_snowc"]) + f(ex["Faxa_snowl"])
        assert rain == f(c["forc_rain"]) and rain == 0.0, s
        assert snow == f(c["forc_snow"]), (s, snow, c["forc_snow"])

        # derived (not dump columns)
        wind = math.sqrt(f(ex["Sa_u"]) ** 2 + f(ex["Sa_v"]) ** 2)
        solar = (f(r["forc_solad_vis"]) + f(r["forc_solad_nir"])
                 + f(r["forc_solai_vis"]) + f(r["forc_solai_nir"]))
        assert abs(wind - f(st["strm_wind"])) <= 1e-12 * f(st["strm_wind"]), s

        # density law: forc_rho = (pbot - 0.378*vp)/(rair*t), vp = q*pbot/(0.622+0.378*q)
        q, pbot, tbot = f(ex["Sa_shum"]), f(ex["Sa_pbot"]), f(ex["Sa_tbot"])
        vp = q * pbot / (EPS + EPS1 * q)
        rho = (pbot - EPS1 * vp) / (RDAIR * tbot)
        err = abs(rho - f(b["forc_rho"])) / f(b["forc_rho"])
        worst = max(worst, err)
        assert err <= 1e-12, (s, err)

        table[s] = {
            "datm_step": t, "clm_step": s, "cell": CELL,
            "overrides": {
                "Sa_z": f(ex["Sa_z"]), "Sa_topo": f(ex["Sa_topo"]),
                "Sa_u": f(ex["Sa_u"]), "Sa_v": f(ex["Sa_v"]),
                "Sa_shum": f(ex["Sa_shum"]), "Sa_ptem": f(ex["Sa_ptem"]),
                "Sa_pbot": f(ex["Sa_pbot"]), "Sa_tbot": f(ex["Sa_tbot"]),
                "Faxa_rainc": f(ex["Faxa_rainc"]), "Faxa_rainl": f(ex["Faxa_rainl"]),
                "Faxa_snowc": f(ex["Faxa_snowc"]), "Faxa_snowl": f(ex["Faxa_snowl"]),
                "Faxa_lwdn": f(ex["Faxa_lwdn"]), "Faxa_swvdr": f(ex["Faxa_swvdr"]),
                "Faxa_swndr": f(ex["Faxa_swndr"]), "Faxa_swvdf": f(ex["Faxa_swvdf"]),
                "Faxa_swndf": f(ex["Faxa_swndf"]),
                "strm_tbot": f(st["strm_tbot"]), "strm_pbot": f(st["strm_pbot"]),
                "strm_wind": f(st["strm_wind"]), "strm_shum": f(st["strm_shum"]),
                "strm_lwdn": f(st["strm_lwdn"]), "strm_swdn": f(st["strm_swdn"]),
                "strm_precn": f(st["strm_precn"]), "strm_topo": f(st["strm_topo"]),
                "strm_z": f(st["strm_z"]), "tbotmax": f(st["tbotmax"]),
            },
            "expectations": {
                "forc_t": f(b["forc_t"]), "forc_q": f(b["forc_q"]),
                "forc_u": f(b["forc_u"]), "forc_v": f(b["forc_v"]),
                "forc_pbot": f(b["forc_pbot"]), "forc_rho": f(b["forc_rho"]),
                "forc_rain": f(c["forc_rain"]), "forc_snow": f(c["forc_snow"]),
                "forc_solad_vis": f(r["forc_solad_vis"]), "forc_solad_nir": f(r["forc_solad_nir"]),
                "forc_solai_vis": f(r["forc_solai_vis"]), "forc_solai_nir": f(r["forc_solai_nir"]),
            },
            "derived": {"forc_wind": wind, "forc_solar": solar},
        }
    print(f"chain laws verified on all {len(STEPS)} overlapping steps; worst density-law rel err = {worst:.3e}")
    return table


def check_esm_tests(table):
    """Check every .esm test override/assertion against the dump-derived table."""
    import re
    problems = 0
    for fname in ("clm_atm_forcing.esm", "m9_datm_clm.esm"):
        path = os.path.join(ESM_DIR, fname)
        if not os.path.exists(path):
            print(f"(skipping {fname}: not found)")
            continue
        doc = json.load(open(path))
        for model in doc["models"].values():
            for test in model.get("tests", []):
                m = re.match(r"step(\d+)_p1", test["id"])
                if not m:
                    continue
                row = table[int(m.group(1))]
                for key, val in test.get("parameter_overrides", {}).items():
                    src = key.split(".")[-1]
                    if src in row["overrides"]:
                        exp = row["overrides"][src]
                    elif key.endswith(".Sa_z") or key.endswith(".strm_z"):
                        exp = row["overrides"]["Sa_z"] if key.endswith(".Sa_z") else row["overrides"]["strm_z"]
                    else:
                        continue  # model flags (needs_downscaling etc.) are run-config, not dump values
                    if abs(val - exp) > 1e-12 * abs(exp):
                        problems += 1
                        print(f"  MISMATCH {fname} {test['id']} override {key}: {val} != dump {exp}")
                for a in test.get("assertions", []):
                    var = a["variable"].split(".")[-1]
                    if var in row["expectations"]:
                        exp = row["expectations"][var]
                    elif var in row["derived"]:
                        exp = row["derived"][var]
                    else:
                        continue
                    if abs(a["expected"] - exp) > 1e-12 * abs(exp):
                        problems += 1
                        print(f"  MISMATCH {fname} {test['id']} assertion {a['variable']}: {a['expected']} != {exp}")
    print(f".esm test tuples checked against the dumps: {'OK' if problems == 0 else f'{problems} PROBLEMS'}")
    return problems


if __name__ == "__main__":
    table = verify_chain()
    print()
    print("derived values (not dump columns), clm step / datm step / cell 1:")
    print(f"  {'step':>4} {'datm':>4} {'forc_wind':>20} {'forc_solar':>20}")
    for s in STEPS:
        d = table[s]["derived"]
        print(f"  {s:>4} {table[s]['datm_step']:>4} {d['forc_wind']:>20.17g} {d['forc_solar']:>20.17g}")
    print()
    sys.exit(1 if check_esm_tests(table) else 0)
