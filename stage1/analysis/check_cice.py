#!/usr/bin/env python3
"""Re-verify the pre-existing CICE dump laws against the instrument record.

Checks the 13 pre-dEdd CSVs (plus the two new ones for row counts and blank
lines) against the checks recorded in stage1/runs/cice_instrument.json:
  - row counts for every dump file,
  - zero blank lines in every dump file,
  - freezing-temperature law Tf = max(S/(-18.48 + 0.01848*S), -2.0) on all
    384 rows of cice_freezing_temp.csv, identical to run_imports sss/Tf,
  - EVP stress relaxation laws newp/newm/new12 on all 184320 rows,
  - EVP momentum final relation uvel=(cca*cc1+ccb*cc2)/ab2,
    vvel=(cca*cc2-ccb*cc1)/ab2 on all 40320 rows,
  - cice_dedd.csv: 27 identity checks per row (delegated to
    check_cice_dedd.py) plus the dedd structural row counts.

All identities recompute the same float64 operations in the same order as
the Fortran, so they hold to round-off; tolerance rtol=1e-12, atol=1e-12.
"""
import csv
import re
import subprocess
import sys

# gfortran prints denormal reals without the 'E' (e.g. 5.58...-322); fix that
DENORM = re.compile(r"([0-9])([+-])([0-9]{3})$")

def tofloat(s):
    s = s.strip()
    m = DENORM.search(s)
    if m:
        s = s[: m.start(1) + 1] + "E" + m.group(2) + m.group(3)
    return float(s)

DUMP = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cice"
RTOL, ATOL = 1e-12, 1e-12

# expected data-row counts from the record's row_counts check
EXPECTED = {
    "cice_run_imports.csv": 384,
    "cice_run_exports.csv": 416,   # incl. the step=0 initial-export row
    "cice_freezing_temp.csv": 384,
    "cice_thermo_state.csv": 3840,
    "cice_dyn_state.csv": 3840,
    "cice_transport_state.csv": 3840,
    "cice_transport_tracers.csv": 3840,
    "cice_evp_stress.csv": 184320,
    "cice_evp_momentum.csv": 40320,
    "cice_ice_strength.csv": 384,
    "cice_albedo.csv": 1920,
    "cice_flux.csv": 768,
    "cice_thermo_flux.csv": 384,
    "cice_budget.csv": 1,
    "cice_dedd.csv": 5760,
}

def rows_of(name):
    path = f"{DUMP}/{name}"
    with open(path) as f:
        lines = f.read().splitlines()
    if any(not l.strip() for l in lines[2:]):
        print(f"FAIL {name}: blank data line")
        sys.exit(1)
    cols = lines[1].split(" = ", 1)[1].split(",")
    return list(csv.DictReader(lines[2:], fieldnames=cols))

def close(a, b):
    return abs(a - b) <= RTOL * abs(b) + ATOL

def main():
    # row counts + blank lines
    for name, n in EXPECTED.items():
        rows = rows_of(name)
        if len(rows) != n:
            print(f"FAIL {name}: {len(rows)} rows, expected {n}")
            sys.exit(1)
        print(f"PASS {name}: {len(rows)} rows, no blank lines")

    # freezing temperature law: Tf = max(S/(-18.48 + 0.01848*S), -2.0)
    frz = rows_of("cice_freezing_temp.csv")
    n = 0
    for r in frz:
        s = tofloat(r["sss"])
        t = tofloat(r["Tf"])
        expect = max(s / (-18.48 + 0.01848 * s), -2.0)
        n += 1
        if not close(t, expect):
            print(f"FAIL freezing_temp row {n}: Tf={t} vs {expect} (sss={s})")
            sys.exit(1)
    print(f"PASS freezing_temp_law: {n}/{n} rows "
          f"Tf = max(S/(-18.48+0.01848*S), -2.0)")

    # freezing temp identical to run_imports sss/Tf on the same rows
    imp = rows_of("cice_run_imports.csv")
    key = {(r["step"], r["i"], r["j"]): r for r in imp}
    n = 0
    for r in frz:
        r2 = key[(r["step"], r["i"], r["j"])]
        n += 1
        if not (close(tofloat(r2["sss"]), tofloat(r["sss"]))
                and close(tofloat(r2["Tf"]), tofloat(r["Tf"]))):
            print(f"FAIL freezing_temp==imports row {n}: "
                  f"{r['sss']}/{r['Tf']} vs {r2['sss']}/{r2['Tf']}")
            sys.exit(1)
    print(f"PASS freezing_temp==run_imports sss/Tf: {n}/{n} rows")

    # EVP stress relaxation laws (classic EVP, revp=0)
    evp = rows_of("cice_evp_stress.csv")
    n = 0
    for r in evp:
        arlx1i = tofloat(r["arlx1i"])
        revp = tofloat(r["revp"])
        denom1 = tofloat(r["denom1"])
        f = 1.0 - arlx1i * revp
        for var, law in (("newp", "(oldp*f + arlx1i*(zetax2*divu-rep_prs))*denom1"),
                         ("newm", "(oldm*f + arlx1i*etax2*tension)*denom1"),
                         ("new12", "(old12*f + arlx1i*0.5*etax2*shear)*denom1")):
            n += 1
            expect = {
                "newp": (tofloat(r["oldp"]) * f
                         + arlx1i * (tofloat(r["zetax2"]) * tofloat(r["divu"])
                                     - tofloat(r["rep_prs"]))) * denom1,
                "newm": (tofloat(r["oldm"]) * f
                         + arlx1i * tofloat(r["etax2"]) * tofloat(r["tension"])) * denom1,
                "new12": (tofloat(r["old12"]) * f
                          + arlx1i * 0.5 * tofloat(r["etax2"])
                          * tofloat(r["shear"])) * denom1,
            }[var]
            if not close(tofloat(r[var]), expect):
                print(f"FAIL evp_stress {var} row {n} "
                      f"({r['i']},{r['j']} ksub={r['ksub']} corner={r['corner']})")
                sys.exit(1)
    print(f"PASS evp_stress_law: {n}/{n} checks "
          f"(newp/newm/new12 relaxation laws, revp=0)")

    # EVP momentum final relation
    mom = rows_of("cice_evp_momentum.csv")
    n = 0
    for r in mom:
        cca, ccb = tofloat(r["cca"]), tofloat(r["ccb"])
        ab2 = tofloat(r["ab2"])
        cc1, cc2 = tofloat(r["cc1"]), tofloat(r["cc2"])
        n += 2
        if not close(tofloat(r["uvel"]), (cca * cc1 + ccb * cc2) / ab2):
            print(f"FAIL evp_momentum uvel row {n // 2} ({r['i']},{r['j']} ksub={r['ksub']})")
            sys.exit(1)
        if not close(tofloat(r["vvel"]), (cca * cc2 - ccb * cc1) / ab2):
            print(f"FAIL evp_momentum vvel row {n // 2} ({r['i']},{r['j']} ksub={r['ksub']})")
            sys.exit(1)
    print(f"PASS evp_momentum_law: {n}/{n} checks "
          f"uvel=(cca*cc1+ccb*cc2)/ab2, vvel=(cca*cc2-ccb*cc1)/ab2")

    # dEdd identity laws (delegated)
    subprocess.run([sys.executable,
                    "stage1/analysis/check_cice_dedd.py"], check=True)

if __name__ == "__main__":
    main()
