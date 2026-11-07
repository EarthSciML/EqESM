#!/usr/bin/env python3
"""Extract numeric test tuples from the stage-1 MOM6 KPP dumps into JSON
fragments for the hand-authored `tests` blocks of stage2/mom6_kpp.esm.

Emits, for each hand-picked dump row, the `parameter_overrides` + `assertions`
of one `time_span {0,0}` test.  The tests blocks themselves are hand-authored
(the descriptions, the row choice, and the pin selection rules below live in
the .esm file; this script only fills in the numeric tuples at full dump
precision, rel 1e-12).

BLD pins (mom6_kpp_obl.csv rows, pre-smoothing OBLdepth/kOBL), keyed by the
dump's (i, j, tod) columns:
  (7,  82,  3600)  kfirst=3, OBL=6.613560999075406,  kOBL=3.75
  (8,  83,  3600)  kfirst=3, OBL=5.974325344745279,  kOBL=3.25
  (7,  84, 10800)  kfirst=2, OBL=5.0 (floored),      kOBL=2.75
  (11, 84,  3600)  kfirst=1, OBL=5.0 (floored),      kOBL=3.25
  (6,  84, 10800)  kfirst=1, OBL=5.0 (floored),      kOBL=2.75
KPPCoefs pins (mom6_kpp.csv rows, smoothed OBLdepth/kOBL):
  (7, 82, 3600)  OBL=6.1872106097070425  kOBL=3.25  (ktup == kwup-1)
  (8, 82, 3600)  OBL=6.309964392417778   kOBL=3.75  (ktup == kwup)
  (6, 84, 10800) OBL=5.001656526745087   kOBL=2.75  (kwup = 2)

Pin-selection rules (the informative core of each row; the 76-interface
0-below-kwup / copy-of-K_in bulk is deliberately NOT emitted):
- BLD: OBLdepth/kOBL plus kfirst (derived, NOT a dump column: from the
  kpp_ref.py reference reimplementation, which reproduces every dump column
  to <= ~1e-12; the esm test descriptions label it "derived (not a dump
  column)"); BulkRi/deltaU2/N2 at layers 1..kOBL+3; z_inter/z_cell spot
  values at the same 1..kOBL+3 range plus the three deep layers 73..75.
- KPPCoefs: all 7 asserted columns (Kd_heat/Kd_salt/Kv/NLT_heat/Kt_out/
  Ks_out/Kv_out) at interfaces 2..kwup+1 (the KPP-active region plus the
  first restoration interface) and at the deep interface 76.

Usage:
  python3 extract_mom6_kpp_test_rows.py [dumps_dir] > tests_fragments.json
where dumps_dir defaults to <main checkout>/stage1/dumps/mom6 (dumps are
gitignored and live in the main checkout, not in the worktree; the
environment variable EQESM_MAIN overrides the main-checkout location).

The kfirst values need the kpp_ref.py reference (same directory).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dumpcsv import load
from kpp_ref import compute_column

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MAIN = os.environ.get("EQESM_MAIN", ROOT)
NLEV = 75
DEEP_LAYERS = (73, 74, 75)   # the three deepest layers of the 75-layer grid
DEEP_IFACE = 76              # the bottom interface

# hand-picked rows: (i, j, tod) -> test id
BLD_ROWS = [
    ((7, 82, 3600), "cell782_step3600_quad_k3"),
    ((8, 83, 3600), "cell883_step3600_quad_k3"),
    ((7, 84, 10800), "cell784_step10800_two_point_k2"),
    ((11, 84, 3600), "cell1184_step3600_kfirst1"),
    ((6, 84, 10800), "cell684_step10800_kfirst1_shallow"),
]
COEFS_ROWS = [
    ((7, 82, 3600), "cell782_step3600_kOBL325"),
    ((8, 82, 3600), "cell882_step3600_kOBL375"),
    ((6, 84, 10800), "cell684_step10800_kOBL275"),
]


def arr(r, prefix, k0, k1):
    return [r[f"{prefix}_{k:02d}"] for k in range(k0, k1 + 1)]


def bld_test(r, idx, r2):
    """The BLD test for one mom6_kpp_obl.csv row (r is the row dict; r2 the
    matching mom6_kpp.csv row, consumed by the kfirst reference only)."""
    overrides = {
        "h": arr(r, "h", 1, NLEV),
        "dz": arr(r, "dz", 1, NLEV),
        "T": arr(r, "T", 1, NLEV),
        "S": arr(r, "S", 1, NLEV),
        "U_H": arr(r, "U_H", 1, NLEV),
        "V_H": arr(r, "V_H", 1, NLEV),
        "buoyFlux": arr(r, "buoyFlux", 1, NLEV + 1),
        "uStar": r["uStar"],
        "p_surf": r["p_surf"],
    }
    asserts = [
        {"variable": "OBLdepth", "time": 0.0, "expected": r["OBLdepth"]},
        {"variable": "kOBL", "time": 0.0, "expected": r["kOBL"]},
        # kfirst is NOT a dump column: it comes from the kpp_ref.py
        # reimplementation of the OBL_depth_low branch selection.
        {"variable": "kfirst", "time": 0.0, "expected": None},
    ]
    kmax = int(r["kOBL"]) + 3  # layers 1..kOBL+3
    for k in range(1, kmax + 1):
        for col in ("BulkRi", "deltaU2", "N2"):
            asserts.append({"variable": col, "time": 0.0, "expected": r[f"{col}_{k:02d}"], "coords": {"layer": k}})
    for k in list(range(1, kmax + 1)) + list(DEEP_LAYERS):
        asserts.append({"variable": "z_cell", "time": 0.0, "expected": r[f"z_cell_{k:02d}"], "coords": {"layer": k}})
        asserts.append({"variable": "z_inter", "time": 0.0, "expected": r[f"z_inter_{k:02d}"], "coords": {"iface": k}})
    # the kfirst branch selection from the reference (0-based -> 1-based)
    c = compute_column(r, r2)
    asserts[2]["expected"] = c["kfirst"] + 1.0
    return overrides, asserts


def coefs_test(r2, idx):
    """The KPPCoefs test for one mom6_kpp.csv row (r2 is the row dict)."""
    overrides = {
        "OBLdepth_s": r2["OBLdepth"],
        "kOBL_s": r2["kOBL"],
        "uStar_s": r2["uStar"],
        "surfBuoyFlux_s": r2["surfBuoyFlux"],
        "z_inter_s": arr(r2, "z_inter", 1, NLEV + 1),
        "Kt_in": arr(r2, "Kt_in", 1, NLEV + 1),
        "Ks_in": arr(r2, "Ks_in", 1, NLEV + 1),
        "Kv_in": arr(r2, "Kv_in", 1, NLEV + 1),
    }
    cols = (("Kd_heat", "Kd_heat"), ("Kd_salt", "Kd_salt"), ("Kv", "Kv"),
            ("NLT_heat", "NLT_heat"), ("Kt_out", "Kt_out"),
            ("Ks_out", "Ks_out"), ("Kv_out", "Kv_out"))
    kwup = int(r2["kOBL"])  # floor(kOBL_s): the kw loop upper bound
    asserts = []
    for k in list(range(2, kwup + 2)) + [DEEP_IFACE]:
        for var, col in cols:
            asserts.append({"variable": var, "time": 0.0, "expected": r2[f"{col}_{k:02d}"], "coords": {"iface": k}})
    return overrides, asserts


def emit(dumps_dir):
    meta, rows = load(f"{dumps_dir}/mom6_kpp_obl.csv")
    _, rows2 = load(f"{dumps_dir}/mom6_kpp.csv")
    by_key = {}
    for r in rows:
        by_key[(int(r["i"]), int(r["j"]), int(r["tod"]))] = r
    by_key2 = {}
    for r in rows2:
        by_key2[(int(r["i"]), int(r["j"]), int(r["tod"]))] = r

    out = {"bld": {}, "coefs": {}}
    for key, tid in BLD_ROWS:
        ov, asr = bld_test(by_key[key], key, by_key2[key])
        out["bld"][tid] = {"parameter_overrides": ov, "assertions": asr}
    for key, tid in COEFS_ROWS:
        ov, asr = coefs_test(by_key2[key], key)
        out["coefs"][tid] = {"parameter_overrides": ov, "assertions": asr}

    json.dump(out, sys.stdout, indent=1)
    sys.stderr.write(f"wrote {len(out['bld'])} BLD + {len(out['coefs'])} coefs test fragments\n")
    for tid, d in out["bld"].items():
        sys.stderr.write(f"  bld {tid}: {len(d['assertions'])} assertions\n")
    for tid, d in out["coefs"].items():
        sys.stderr.write(f"  coefs {tid}: {len(d['assertions'])} assertions\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dumps_dir", nargs="?", default=f"{MAIN}/stage1/dumps/mom6",
                    help="directory with mom6_kpp_obl.csv and mom6_kpp.csv "
                         "(default: <main checkout>/stage1/dumps/mom6)")
    args = ap.parse_args()
    emit(args.dumps_dir)


if __name__ == "__main__":
    main()
