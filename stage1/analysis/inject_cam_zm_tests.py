#!/usr/bin/env python3
"""Materialize the dump-extracted test tuples into the hand-authored skeletons
stage2/cam_zm_q1q2.esm and stage2/cam_zm_evap.esm (the sanctioned use of
scripts: extracting numeric test tuples from the stage-1 dumps into the
hand-authored tests blocks).  Run from the repo root, after
stage1/analysis/extract_cam_zm_tests.py has written
scratch/q1q2_testdata.json and scratch/evap_testdata.json:
  python3 stage1/analysis/inject_cam_zm_tests.py
The skeleton structure, equations, descriptions and assertion choices are
hand-authored; this script fills in only the numbers: the parameter override
arrays and the expected values, both taken from the dumps (with evpvint and
prec_gen the only derived quantities, see extract_cam_zm_tests.py), and the
estbl table copy (taken verbatim from the committed cam_wv_saturation.esm
table).  It also reports the min-limit binding counts per test column for the
test descriptions: the counts classify which of the four bounds equals the
min-selected evpprec per level (min() returns one of its arguments exactly,
so the classification is exact for the dump-fed chain; it feeds the
description text only, never the assertions).  Idempotent: the tests arrays
and the estbl table entries are always overwritten with the extracted
values, so re-running after a data fix is safe.

The assertion expectations are the DUMP values (never this script's), and
the copy check that the replica reproduces the dump within the assertion
tolerance (rel 1e-12, abs 1e-24 for identically-zero pins) is the hard gate
in extract_cam_zm_tests.py; the tolerance selection here for the zero
prec_evap pins is likewise within the abs tolerance (|prec_evap| < 1e-24),
not an exact == 0 comparison.
"""
import csv
import json

M = '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cam/'
PVER = 32
TOL = {"rel": 1e-12}
TOL0 = {"rel": 1e-12, "abs": 1e-24}
ABS0 = 1e-24

ZERO_FAIL = ("ntsnprd = rprd work2 - evpsnow - snowmlt is a three-term "
             "composite whose zero requires all three branch states: "
             "fsnow_conv = 0 (the dumped snow-production fraction; if the "
             "snow region reached this level, work2 = "
             "max(fsnow_conv, flxsnow/flxprec) > 0), flxsntm = 0 and "
             "snowmlt = 0 (both from t1 > tmelt here - the melting level 22 "
             "is above, so no snow enters or melts in the layer), and the "
             "vanished snow flux flxsnow_hi = 0 below the melt - the pin "
             "fails if any of these misbehaves")

MELT_TAIL = ("and the snow flux vanishes below the melting level "
             "(flxsnow_hi(22) = 0, so no snow reaches the surface; "
             "snow_evap = flxsnow_hi(32)/1000 is consequently 0 and is not "
             "separately pinned - the same zero region downstream of the "
             "flxsnow_hi(22) pin).")


def load(fn):
    with open(fn) as f:
        lines = f.readlines()
    cols = lines[1][len('# columns = '):].strip().split(',')
    return cols, list(csv.reader(lines[2:]))


def binder_counts(step):
    """Min-limit binding counts for one dumped column: which of the four
    bounds equals evpprec at each level (0 = Sundqvist raw, 2 = the deltat
    saturation limit, 3 = the falling-flux limit, 4 = the
    remaining-precipitation limit).  The chain is fed the DUMP flux
    (flxprec_hi) at each level, so the classification describes the true
    column; min() returns one of its arguments exactly, so the counts are
    exact for this chain (two bounds within rounding of each other are
    counted as the one min() selected, which can differ from the
    recomputed-chain histogram of extract_cam_zm_tests.py by a level)."""
    cols, rows = load(M + 'cam_zm.csv')
    scols, srows = load(M + 'cam_zm_scal.csv')
    ci = {c: i for i, c in enumerate(cols)}
    si = {c: i for i, c in enumerate(scols)}
    f = float
    s = [r for r in srows if (int(r[0]), int(r[3])) == (step, 1)][0]
    pr = [r for r in rows if (int(r[0]), int(r[3])) == (step, 1)]
    pr.sort(key=lambda r: int(r[4]))
    estbl = json.load(open('stage2/cam_wv_saturation.esm'))
    ESTBL = estbl['function_tables']['estbl']['data']
    EPS = 0.6219705862045155
    OM = 1 - EPS
    grav = f(s[si['grav']])
    ke = f(s[si['ke']])
    deltat = f(s[si['ztodt']])
    prec_gen = f(s[si['prec_convr']]) * 1000.0
    flxprec = 0.0
    evpvint = 0.0
    counts = {0: 0, 2: 0, 3: 0, 4: 0}
    for k, r in enumerate(pr):
        t = f(r[ci['t1']])
        q = f(r[ci['q1']])
        pm = f(r[ci['pmid']])
        pd = f(r[ci['pdel']])
        cf = f(r[ci['cldfrc']])
        tt = max(min(t, 375.16) - 127.16, 0.0)
        i = int(tt) + 1
        w = tt - int(tt)
        es = (1 - w) * ESTBL[i - 1] + w * ESTBL[i]
        qs = EPS * es / (pm - OM * es)
        raw = ke * (1 - cf) * max(1 - q / qs, 0.0) * (flxprec ** 0.5)
        f2 = max(0.0, (qs - q) / deltat)
        f3 = flxprec * grav / pd
        f4 = (prec_gen - evpvint) * grav / pd
        evpprec = min(min(min(f2, f3), f4), raw)
        bind = 0 if evpprec == raw else 2 if evpprec == f2 else \
            3 if evpprec == f3 else 4
        counts[bind] += 1
        evpvint += evpprec * pd / grav
        flxprec = f(r[ci['flxprec_hi']])
    return counts


def make_q1q2_tests():
    q1q2 = json.load(open('scratch/q1q2_testdata.json'))
    tests = []
    for st in ['18', '40', '48']:
        d = q1q2[st]
        sc = d['scalars']
        a = d['arrays']
        ov = {c: a[c] for c in ['cme', 'mu', 'md', 'du', 'qu', 'su', 'qd',
                                'sd', 'qhat', 'shat', 'qcde', 'pdel']}
        ov.update(jt=sc['jt'], maxg=sc['maxg'], dsubcld=sc['dsubcld'],
                  rl=sc['rl'], cpres=sc['cpres'])
        ex = d['expects']
        pins = [(15, 'dqdt'), (17, 'dqdt'), (25, 'dqdt'), (32, 'dqdt'),
                (15, 'dsdt'), (17, 'dsdt'), (25, 'dsdt'), (32, 'dsdt'),
                (15, 'dlf'), (16, 'dlf')]
        assertions = [
            {"variable": v, "time": 0.0, "expected": ex[v][k - 1],
             "coords": {"lev": k}, "tolerance": TOL} for k, v in pins]
        tests.append({
            "id": f"s{st}_g1",
            "description": (
                f"Step {st}, g = 1 column (cam_zm.csv step {st} rows lev "
                f"1..32 with cam_zm_scal.csv step {st}; jt = {sc['jt']}, "
                f"maxg = {sc['maxg']}, dsubcld = {sc['dsubcld']:.6g} mb): "
                f"the q1q2_pjr tendencies of zm_convr.F90:2584-2622 pinned at "
                f"the convection top lev 15 (= jt, where du > 0 so dl is "
                f"active), the downdraft region lev 17, a deep level lev 25 "
                f"and the subcloud level lev 32 (= maxg = pver, the "
                f"subcloud branch of F90:2611-2619); dl at the two detraining "
                f"levels 15-16. The k > maxg copy-down of F90:2620-2622 "
                f"cannot fire (maxg = pver on all dumped columns)."),
            "time_span": {"start": 0.0, "end": 0.0},
            "parameter_overrides": ov,
            "assertions": assertions})
    return tests


def make_evap_tests():
    evap = json.load(open('scratch/evap_testdata.json'))
    tests = []
    for st in ['18', '22', '40', '48']:
        d = evap[st]
        sc = d['scalars']
        a = d['arrays']
        ov = {c: a[c] for c in ['t1', 'q1', 'pmid', 'pdel', 'cldfrc',
                                'fsnow_conv', 'rprd', 'flxprec_lo',
                                'flxsnow_lo']}
        ov['evpvint'] = d['evpvint']
        # tmelt is deliberately not overridden: the ZMevap model uses the
        # shared tmelt constant fragment of cam_templates.esm (TKFRZ = 273.15
        # K), not a local parameter.
        ov.update(prec_convr=sc['prec_convr'], gravit=sc['gravit'],
                  latvap=sc['latvap'], latice=sc['latice'],
                  cpres=sc['cpres'], ke=sc['ke'], deltat=sc['deltat'])
        ex = {('tend_s_snwprd' if k == 'tend_snwprd' else
               'tend_s_snwevmlt' if k == 'tend_snwevmlt' else k): v
              for k, v in d['expects'].items()}
        cnt = binder_counts(int(st))
        deep = (f"the deep-evaporation region (lev 25, where the net "
                f"production ntprprd goes negative and the snow flux has "
                f"vanished so ntsnprd(25) = 0: {ZERO_FAIL})" if st != '22'
                else f"the deep-evaporation region (levs 24-25, where the "
                f"net production ntprprd goes negative and the snow flux "
                f"has vanished so ntsnprd(24) = 0: {ZERO_FAIL})")
        intro = (f"Step {st}, g = 1 column (cam_zm.csv step {st} rows lev "
                 f"1..32 with cam_zm_scal.csv step {st}; prec_convr = "
                 f"{sc['prec_convr']:.6e} m/s, ke = {sc['ke']:.1e}, deltat = "
                 f"{sc['deltat']:.0f} s): the zm_conv_evap_run chain of "
                 f"zm_conv_evap.F90:118-236 pinned at the first raining "
                 f"levels (levs 16-17), {deep} and the bottom boundary "
                 f"(lev 32, the surface rate prec_evap = flxprec_hi(32)/1000 "
                 f"of F90:234-236). ")
        if st == '18':
            pins = [(17, 'ntprprd'), (25, 'ntprprd'), (31, 'ntprprd'),
                    (16, 'flxprec_hi'), (25, 'flxprec_hi'),
                    (31, 'flxprec_hi', True), (32, 'flxprec_hi', True),
                    (17, 'ntsnprd'), (25, 'ntsnprd'), (17, 'tend_q'),
                    (17, 'tend_s'), (16, 'tend_s_snwprd'),
                    (17, 'tend_s_snwevmlt')]
            note = (f"The precipitation evaporates completely in the column "
                    f"(prec_evap = 0): at lev 31 the unclamped flux "
                    f"flxprec_lo + ntprprd pdel/gravit is a roundoff residual "
                    f"of order -1e-19 kg/m2/s and the max guard of "
                    f"F90:222-224 pins it to 0 (the dumped flxprec_hi(31) = "
                    f"0), so flxprec_hi(32) and prec_evap are 0; lev 17 pins "
                    f"the evaporation of the first raining layer (tend_q = "
                    f"evpprec > 0); lev 25 pins the net production going "
                    f"negative (evaporation of the downdraft-evaporated "
                    f"precipitation exceeds the production) with ntsnprd = 0 "
                    f"below the snow region (fsnow_conv = 0, work2 = 0). "
                    f"The evpprec min chain binds the falling-flux limit "
                    f"F90:173 at {cnt[3]} levels and the "
                    f"remaining-precipitation limit F90:176 at {cnt[4]} "
                    f"levels of this column (never the deltat limit F90:168, "
                    f"which does not bind on any dumped column); the "
                    f"Sundqvist rate F90:164 is never the active bound with "
                    f"flxprec_lo > 0 here (it binds only where the entering "
                    f"flux is zero). ")
        elif st == '22':
            pins = [(17, 'ntprprd'), (24, 'ntprprd'), (31, 'ntprprd'),
                    (16, 'flxprec_hi'), (24, 'flxprec_hi'),
                    (31, 'flxprec_hi', True), (32, 'flxprec_hi', True),
                    (17, 'ntsnprd'), (24, 'ntsnprd'), (24, 'tend_q'),
                    (24, 'tend_s')]
            note = (f"prec_evap = 0 (the column evaporates its precipitation; "
                    f"the roundoff guard of F90:222-224 fires at lev 31, "
                    f"flxprec_hi(31) = 0). lev 24 pins the Sundqvist limit "
                    f"binding NON-degenerately: at levels 23-26 of this "
                    f"column evpprec = evpprec_raw = ke (1 - cldfrc) "
                    f"max(1 - q/qs, 0) sqrt(flxprec_lo) (F90:160-164, 178) "
                    f"with flxprec_lo > 0, so tend_q pins the "
                    f"relative-humidity factor max(1 - q/qs, 0) of F90:158 "
                    f"through the table-interpolated qs. The min chain binds "
                    f"the falling-flux limit at {cnt[3]} levels and the "
                    f"remaining-precipitation limit at {cnt[4]} levels. ")
        else:
            pins = [(17, 'ntprprd'), (25, 'ntprprd'), (32, 'ntprprd'),
                    (17, 'flxprec_hi'), (25, 'flxprec_hi'),
                    (32, 'flxprec_hi'), (17, 'ntsnprd'), (25, 'ntsnprd'),
                    (17, 'tend_q'), (17, 'tend_s'), (17, 'tend_s_snwprd'),
                    (17, 'tend_s_snwevmlt')]
            if st == '40':
                note = (f"prec_evap = {d['prec_evap']:.6e} m/s, the "
                        f"second-window early column (steps 40-48); the "
                        f"evpprec min chain binds the Sundqvist rate "
                        f"evpprec_raw of F90:160-164 at every level with "
                        f"flxprec_lo > 0 (levels 17-32) - the "
                        f"falling-flux and remaining-precipitation limits of "
                        f"F90:173/176 never bind in this column - so tend_q "
                        f"at lev 17 pins the max(1 - q/qs, 0) factor of "
                        f"F90:158 through the table-interpolated qs, and "
                        f"tend_q at lev 30 pins the same factor deeper in "
                        f"the column. ")
                pins = pins + [(30, 'tend_q')]
            else:
                note = (f"prec_evap = {d['prec_evap']:.6e} m/s, the "
                        f"largest of the window; as in step 40 the Sundqvist "
                        f"rate evpprec_raw of F90:160-164 binds at every "
                        f"level with flxprec_lo > 0, so tend_q at levs 30 "
                        f"and 32 pins the max(1 - q/qs, 0) factor of "
                        f"F90:158 through the table-interpolated qs, and "
                        f"tend_s at lev 32 pins the evaporation heating of "
                        f"the bottom layer. ")
                pins = pins + [(30, 'tend_q'), (32, 'tend_q'), (32, 'tend_s')]
        # The melting level: lev 22 of every dumped column has t1 > tmelt
        # with flxsnow_lo > 0, so the melt branch of F90:128-134 fires there
        # (snowmlt > 0, work2 zeroed at F90:212, ntsnprd = -snowmlt) and the
        # snow flux vanishes below (flxsnow_hi(22) = 0).
        pins = pins + [(22, 'ntsnprd'), (22, 'flxsnow_hi', True),
                       (22, 'tend_s_snwevmlt')]
        melt = (f"lev 22 is the melting level of this column (t1(22) = "
                f"{a['t1'][21]:.2f} K > tmelt with flxsnow_lo > 0): the melt "
                f"branch of F90:128-134 fires (flxsntm = 0, snowmlt = "
                f"flxsnow_lo gravit/pdel > 0), work2 is zeroed at F90:212 so "
                f"ntsnprd(22) = -snowmlt (the dumped net snow production goes "
                f"negative), tend_s_snwevmlt(22) = -snowmlt latice is the melt "
                f"cooling, {MELT_TAIL}")
        assertions = []
        for pin in pins:
            k, v = pin[0], pin[1]
            zero = len(pin) > 2
            a_ = {"variable": v, "time": 0.0, "expected": ex[v][k - 1],
                  "coords": {"lev": k}, "tolerance": TOL0 if zero else TOL}
            assertions.append(a_)
        assertions.append({"variable": "prec_evap", "time": 0.0,
                           "expected": d['prec_evap'],
                           "tolerance": TOL0
                           if abs(d['prec_evap']) < ABS0 else TOL})
        tests.append({
            "id": f"s{st}_g1",
            "description": intro + note + melt,
            "time_span": {"start": 0.0, "end": 0.0},
            "parameter_overrides": ov,
            "assertions": assertions})
    return tests


def main():
    wv = json.load(open('stage2/cam_wv_saturation.esm'))
    estbl = wv['function_tables']['estbl']

    q = json.load(open('stage2/cam_zm_q1q2.esm'))
    q['models']['ZMq1q2']['tests'] = make_q1q2_tests()
    json.dump(q, open('stage2/cam_zm_q1q2.esm', 'w'), indent=1)
    json.dump(q, open('scratch/cam_zm_q1q2.esm.materialized', 'w'), indent=1)

    e = json.load(open('stage2/cam_zm_evap.esm'))
    e['function_tables']['estbl']['axes'][0]['values'] = \
        estbl['axes'][0]['values']
    e['function_tables']['estbl']['data'] = estbl['data']
    e['models']['ZMevap']['tests'] = make_evap_tests()
    json.dump(e, open('stage2/cam_zm_evap.esm', 'w'), indent=1)
    json.dump(e, open('scratch/cam_zm_evap.esm.materialized', 'w'), indent=1)

    for st in ['18', '22', '40', '48']:
        print(f"step {st}: binders {binder_counts(int(st))}")
    print("injected")


if __name__ == '__main__':
    main()
