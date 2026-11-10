#!/usr/bin/env python3
"""Identity checks on cice_dedd.csv (CICE delta-Eddington shortwave dump).

The dump has one row per (cell, category, timestep) triple per spectral band
(ns = 1 visible, 2/3 near-IR); the three rows of a triple carry the running
values of the accumulated outputs, final at ns = 3.  The band loop in
compute_dEdd_3bd (icepack_shortwave.F90:2116-2977) assigns, per band ns,

    wghtns(1) = 1
    wghtns(2) = 0.67 + 0.11*(1 - fnidr)
    wghtns(3) = 1 - wghtns(2)
    swdr/swdf = swvdr/swvdf (ns=1), swidr/swidf (ns=2,3)
    flux_ns(k) = dfdir(k)*swdr + dfdif(k)*swdf      (interface k = 0..13)
    avdr = rupdir(0), avdf = rupdif(0)              (ns = 1 only)
    aidr += rupdir(0)*wghtns(ns), aidf += rupdif(0)*wghtns(ns)   (ns = 2,3)
    fsfc  += (flux(0) - flux(ksrf))   * wghtns(ns)
    fint  += (flux(ksrf) - flux(13))  * wghtns(ns)
    fthru +=  flux(13)                * wghtns(ns)
    fthrul(i) = flux(4+i)                        (ns = 1 only; i = 1..9)
    Sabs(i)  += flux(i) - flux(i+1)              (ns = 1..3, srftyp=1 only)
    Iabs(i)  += flux(km) - flux(kp)              (ns = 1..3)

with ksrf = 1 (snow) / 5 (bare, pond), and the Iabs layer indexing
km/kp = (4+i, 5+i), except for snow-covered ice layer 1 (i=1): km=4, kp=6
(the snow surface scattering layer + ice DL layer 1 are absorbed together).

Each row is checked against 27 identities (5760 rows x 27 = 155520):
  fthrul_1..9  = visible flux at interfaces 5..13 (ns=1 row's dfdir/dfdif),
  avdr/avdf    = ns=1 row's rupdir(0)/rupdif(0),
  aidr/aidf    = ns=2..ns(row) running sum of rupdir(0)/rupdif(0)*wghtns,
  Sabs_1..3, Iabs_1..8, fsfc, fint, fthru = running band sums 1..ns(row)
  of the per-band contributions above (with wghtns applied exactly as the
  Fortran does: (flux(km)*swd + flux(kp)*swd)*wghtns, one rounding each).

Recomputation uses the dumped dfdir/dfdif values (which already carry the
echmod < puny -> 0 clamp) in the same operation order as the Fortran, so
the checks hold to round-off; tolerance rtol=1e-12, atol=1e-12.
"""
import csv
import re
import sys
from collections import defaultdict

# gfortran prints denormal reals without the 'E' (e.g. 5.58...-322); fix that
DENORM = re.compile(r"([0-9])([+-])([0-9]{3})$")

def tofloat(s):
    s = s.strip()
    m = DENORM.search(s)
    if m:
        s = s[: m.start(1) + 1] + "E" + m.group(2) + m.group(3)
    return float(s)

PATH = sys.argv[1] if len(sys.argv) > 1 else \
    "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cice/cice_dedd.csv"

RTOL, ATOL = 1e-12, 1e-12

def close(a, b):
    return abs(a - b) <= RTOL * abs(b) + ATOL

def main():
    with open(PATH) as f:
        f.readline()  # # description
        cols = f.readline().strip().split(" = ", 1)[1].split(",")
        rows = list(csv.DictReader(f, fieldnames=cols))

    n = len(rows)
    if n != 5760:
        print(f"FAIL: expected 5760 rows, got {n}")
        sys.exit(1)

    # group rows into (step, i, j, ncat) triples, one row per band ns
    triples = defaultdict(list)
    for r in rows:
        triples[(r["step"], r["i"], r["j"], r["ncat"])].append(r)

    ncheck = 0
    fails = []
    fthru_eq_fthrul9 = 0

    for key, grp in triples.items():
        if len(grp) != 3:
            fails.append(f"{key}: {len(grp)} rows instead of 3")
            continue
        byns = {int(r["ns"]): r for r in grp}
        if set(byns) != {1, 2, 3}:
            fails.append(f"{key}: bands {sorted(byns)}")
            continue
        r1, r2, r3 = byns[1], byns[2], byns[3]
        srftyp = int(r1["srftyp"])
        if not (r1["srftyp"] == r2["srftyp"] == r3["srftyp"]):
            fails.append(f"{key}: srftyp differs within triple")

        def df(r, k):
            return tofloat(r[f"dfdir_{k}"]), tofloat(r[f"dfdif_{k}"])

        swvdr, swvdf = tofloat(r1["swvdr"]), tofloat(r1["swvdf"])
        swidr, swidf = tofloat(r1["swidr"]), tofloat(r1["swidf"])
        fnidr = tofloat(r1["fnidr"])
        w2 = 0.67 + 0.11 * (1.0 - fnidr)
        w3 = 1.0 - w2
        w = {1: 1.0, 2: w2, 3: w3}

        def flux(r, k, b):
            # flux at interface k in band b, as the Fortran forms it
            ddr, ddf = df(r, k)
            swdr = swvdr if b == 1 else swidr
            swdf = swvdf if b == 1 else swidf
            return (ddr * swdr + ddf * swdf) * w[b]

        def check(cond, label):
            nonlocal ncheck
            ncheck += 1
            if not cond:
                fails.append(f"{key} {label}")

        # ksrf: interface index of surface absorption
        ksrf = 1 if srftyp == 1 else 5

        for r, ns in ((r1, 1), (r2, 2), (r3, 3)):
            # 1) fthrul_1..9 = visible flux at interfaces 5..13
            for i in range(1, 10):
                ddr, ddf = df(r1, i + 4)
                check(close(tofloat(r[f"fthrul_{i}"]), ddr * swvdr + ddf * swvdf),
                      f"fthrul_{i} ns={ns}")
            # 2) avdr/avdf = ns=1 row's rupdir(0)/rupdif(0)
            check(close(tofloat(r["avdr"]), tofloat(r1["rupdir_0"])), f"avdr ns={ns}")
            check(close(tofloat(r["avdf"]), tofloat(r1["rupdif_0"])), f"avdf ns={ns}")
            # 3) aidr/aidf = running sum over bands 2..ns
            aidr = sum(tofloat(byns[b]["rupdir_0"]) * w[b] for b in (2, 3) if b <= ns)
            aidf = sum(tofloat(byns[b]["rupdif_0"]) * w[b] for b in (2, 3) if b <= ns)
            check(close(tofloat(r["aidr"]), aidr), f"aidr ns={ns}")
            check(close(tofloat(r["aidf"]), aidf), f"aidf ns={ns}")
            # 4) Sabs_1..3 running band sums (0 unless snow-covered)
            for i in range(1, 4):
                exp = 0.0 if srftyp != 1 else \
                    sum(flux(byns[b], i, b) - flux(byns[b], i + 1, b)
                        for b in range(1, ns + 1))
                check(close(tofloat(r[f"Sabs_{i}"]), exp), f"Sabs_{i} ns={ns}")
            # 5) Iabs_1..8 running band sums
            for i in range(1, 9):
                if srftyp == 1 and i == 1:
                    km, kp = 4, 6
                else:
                    km, kp = 4 + i, 5 + i
                exp = sum(flux(byns[b], km, b) - flux(byns[b], kp, b)
                          for b in range(1, ns + 1))
                check(close(tofloat(r[f"Iabs_{i}"]), exp), f"Iabs_{i} ns={ns}")
            # 6) fsfc/fint/fthru running band sums
            fsfc = sum(flux(byns[b], 0, b) - flux(byns[b], ksrf, b)
                       for b in range(1, ns + 1))
            fint = sum(flux(byns[b], ksrf, b) - flux(byns[b], 13, b)
                       for b in range(1, ns + 1))
            fthru = sum(flux(byns[b], 13, b) for b in range(1, ns + 1))
            check(close(tofloat(r["fsfc"]), fsfc), f"fsfc ns={ns}")
            check(close(tofloat(r["fint"]), fint), f"fint ns={ns}")
            check(close(tofloat(r["fthru"]), fthru), f"fthru ns={ns}")

        # informational: fthrul_9 (= visible flux at 13) vs total fthru
        if close(tofloat(r3["fthrul_9"]), tofloat(r3["fthru"])):
            fthru_eq_fthrul9 += 1

    print(f"rows: {n}, triples: {len(triples)}")
    if fails:
        print(f"FAIL: {len(fails)} of {ncheck} checks failed, e.g.:")
        for f in fails[:10]:
            print("  ", f)
        sys.exit(1)
    print(f"PASS: {ncheck}/{ncheck} identity checks "
          f"({ncheck // n} per row x {n} rows)")
    print(f"fthrul_9 == fthru on {fthru_eq_fthrul9} of {len(triples)} triples")

if __name__ == "__main__":
    main()
