#!/usr/bin/env python3
"""Extract numeric test tuples for cam_zm_q1q2.esm and cam_zm_evap.esm from
the stage-1 ZM dumps (cam_zm.csv / cam_zm_scal.csv), and verify that the
planned esm laws reproduce the dump before the tuples are committed into the
hand-authored tests blocks.  Run from the repo root:
  python3 stage1/analysis/extract_cam_zm_tests.py
writes scratch/q1q2_testdata.json and scratch/evap_testdata.json, which
inject_cam_zm_tests.py then fills into the hand-authored tests skeletons.

For cam_zm_evap.esm the Fortran chain (zm_conv_evap.F90) is reproduced level
by level with the estbl table values taken from the committed
cam_wv_saturation.esm (the same svp_trans-generated table the Fortran
interpolates), including the evpvint running accumulator (derived from the
dumped rprd and ntprprd columns: evpvint[k] = sum_{j<k} (rprd[j] -
ntprprd[j]) pdel[j]/grav), the min() limit chain in the Fortran's nesting,
the snow-fraction work1/work2 rules and the max(,0) roundoff guard.  The
assertion expectations in the tests blocks are the DUMP values (never this
script's), with evpvint/prec_gen the only derived inputs.

COPY CHECK: the replica-versus-dump verification at the end compares with the
assertion tolerance (rel 1e-12, the same TOL the tests use; abs 1e-24 for
identically-zero dump values, the tests' TOL0), never with ==, and asserts
(gates) on it, so a dump or law change that would break the committed pins
fails this script with the exact offending column and level.  The binder
histogram below is description-only (which evplimit term the min chain
selects per level, min() returning one of its arguments exactly, so the
classification is exact for the replica chain); it feeds the test-description
text, not the assertions.

Level convention: CAM levels 1..32 top-down, matching the dump rows; the esm
models declare lev with the same top-down convention (k = 1 model top, k =
32 surface), so dumped arrays enter the tests unchanged.
"""
import csv
import json

M = '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cam/'
PVER = 32
REL = 1e-12       # the tests' assertion tolerance (rel)
ABS0 = 1e-24      # the tests' abs tolerance for identically-zero pins
EPSILO = 0.6219705862045155
OMEPS = 1.0 - EPSILO
TMAX = 375.16
TMIN = 127.16


def load(fn):
    with open(fn) as f:
        lines = f.readlines()
    cols = lines[1][len('# columns = '):].strip().split(',')
    return cols, list(csv.reader(lines[2:]))


def estbl_from_wv_sat():
    d = json.load(open('stage2/cam_wv_saturation.esm'))
    return d['function_tables']['estbl']['data']


def estblf(t, ESTBL):
    t_tmp = max(min(t, TMAX) - TMIN, 0.0)
    i = int(t_tmp) + 1                     # 1-based index
    w = t_tmp - int(t_tmp)
    return (1.0 - w) * ESTBL[i - 1] + w * ESTBL[i]


def svp_to_qsat(es, p):
    return 1.0 if (p - es) <= 0.0 else EPSILO * es / (p - OMEPS * es)


def check(got, want, label):
    """Copy check within the assertion tolerance (never ==): rel 1e-12, abs
    1e-24 for identically-zero dump values (the tests' TOL0 for zero pins)."""
    if isinstance(got, list):
        worst = 0.0
        for g, w in zip(got, want):
            err = (abs(g - w) if w == 0.0
                   else abs(g - w) / abs(w))
            worst = max(worst, err)
        lim = ABS0 if all(w == 0.0 for w in want) else REL
        assert worst <= lim, (
            f"{label}: replica diverges from the dump by rel {worst:.3e} "
            f"> {lim:.0e} (the assertion tolerance) - the committed pins "
            f"would fail")
        return worst
    err = (abs(got - want) if want == 0.0 else abs(got - want) / abs(want))
    lim = ABS0 if want == 0.0 else REL
    assert err <= lim, (
        f"{label}: replica diverges from the dump by rel {err:.3e} "
        f"> {lim:.0e} (the assertion tolerance)")
    return err


def main():
    cols, rows = load(M + 'cam_zm.csv')
    scols, srows = load(M + 'cam_zm_scal.csv')
    ci = {c: i for i, c in enumerate(cols)}
    si = {c: i for i, c in enumerate(scols)}
    ESTBL = estbl_from_wv_sat()

    def f(x):
        return float(x)

    def profile(key):
        pr = [r for r in rows if (int(r[0]), int(r[3])) == key]
        pr.sort(key=lambda r: int(r[4]))
        return pr

    def scal(key):
        for r in srows:
            if (int(r[0]), int(r[3])) == key:
                return r

    steps = [18, 22, 40, 48]
    q1q2_tests = {}
    evap_tests = {}
    for st in steps:
        s = scal((st, 1))
        pr = profile((st, 1))
        assert len(pr) == PVER, (st, len(pr))
        grav = f(s[si['grav']])
        latvap = f(s[si['rl']])
        latice = f(s[si['latice']])
        tmelt = f(s[si['tfreez']])
        cpres = f(s[si['cpres']])
        ke = f(s[si['ke']])
        deltat = f(s[si['ztodt']])
        prec_convr = f(s[si['prec_convr']])
        jt = int(f(s[si['jt']]))
        maxg = int(f(s[si['maxg']]))
        dsubcld = f(s[si['dsubcld']])
        A = {c: [f(r[ci[c]]) for r in pr] for c in [
            't1', 'q1', 'pmid', 'pdel', 'cldfrc', 'fsnow_conv', 'rprd',
            'flxprec_lo', 'flxsnow_lo', 'ntprprd', 'ntsnprd',
            'flxprec_hi', 'flxsnow_hi', 'tend_s', 'tend_q',
            'tend_snwprd', 'tend_snwevmlt', 'cme', 'mu', 'md', 'du',
            'qu', 'su', 'qd', 'sd', 'qhat', 'shat', 'qcde', 'dqdt',
            'dsdt', 'dlf']}

        # ---------- q1q2_pjr chain (zm_convr.F90:2584-2622) -------------
        dqdt = [0.0] * PVER
        dsdt = [0.0] * PVER
        dl = [0.0] * PVER
        for k in range(jt - 1, PVER - 1):          # k = jt..pver-1
            dp = 0.01 * A['pdel'][k]
            emc = -A['cme'][k]
            dsdt[k] = (-latvap / cpres * emc +
                       (A['mu'][k + 1] * (A['su'][k + 1] - A['shat'][k + 1])
                        - A['mu'][k] * (A['su'][k] - A['shat'][k])
                        + A['md'][k + 1] * (A['sd'][k + 1] - A['shat'][k + 1])
                        - A['md'][k] * (A['sd'][k] - A['shat'][k])) / dp)
            dqdt[k] = (emc +
                       (A['mu'][k + 1] * (A['qu'][k + 1] - A['qhat'][k + 1])
                        - A['mu'][k] * (A['qu'][k] - A['qhat'][k])
                        + A['md'][k + 1] * (A['qd'][k + 1] - A['qhat'][k + 1])
                        - A['md'][k] * (A['qd'][k] - A['qhat'][k])) / dp)
            dl[k] = A['du'][k] * A['qcde'][k + 1]
        for k in range(maxg - 1, PVER):            # k = maxg..pver
            if k == maxg - 1:
                dsdt[k] = (1.0 / dsubcld) * (
                    -A['mu'][k] * (A['su'][k] - A['shat'][k])
                    - A['md'][k] * (A['sd'][k] - A['shat'][k]))
                dqdt[k] = (1.0 / dsubcld) * (
                    -A['mu'][k] * (A['qu'][k] - A['qhat'][k])
                    - A['md'][k] * (A['qd'][k] - A['qhat'][k]))
            elif k > maxg - 1:
                dqdt[k] = dqdt[k - 1]
                dsdt[k] = dsdt[k - 1]

        for label, got, want in [('dqdt', dqdt, A['dqdt']),
                                 ('dsdt', dsdt, A['dsdt']),
                                 ('dlf', dl, A['dlf'])]:
            e = check(got, want, f"step {st} q1q2 {label}")
            print(f"step {st} q1q2: max rel err {label} {e:.3e}")
        q1q2_tests[st] = dict(
            arrays={c: A[c] for c in ['cme', 'mu', 'md', 'du', 'qu', 'su',
                                      'qd', 'sd', 'qhat', 'shat', 'qcde',
                                      'pdel']},
            scalars=dict(jt=jt, maxg=maxg, dsubcld=dsubcld, rl=latvap,
                         cpres=cpres),
            expects=dict(dqdt=[dqdt[k] for k in range(PVER)],
                         dsdt=[dsdt[k] for k in range(PVER)],
                         dlf=[dl[k] for k in range(PVER)]))

        # ---------- zm_conv_evap_run chain (F90:118-236) ----------------
        flxprec = 0.0
        flxsnow = 0.0
        evpvint = 0.0
        prec_gen = prec_convr * 1000.0
        out = {}
        binder = []
        for k in range(PVER):
            t = A['t1'][k]
            q = A['q1'][k]
            if t > tmelt:
                flxsntm = 0.0
                snowmlt = flxsnow * grav / A['pdel'][k]
            else:
                flxsntm = flxsnow
                snowmlt = 0.0
            qs = svp_to_qsat(estblf(t, ESTBL), A['pmid'][k])
            evplimit_rh = max(1.0 - q / qs, 0.0)
            evpprec_raw = ke * (1.0 - A['cldfrc'][k]) * evplimit_rh * \
                (flxprec ** 0.5)
            evplimit = max(0.0, (qs - q) / deltat)
            evplimit = min(evplimit, flxprec * grav / A['pdel'][k])
            evplimit = min(evplimit,
                           (prec_gen - evpvint) * grav / A['pdel'][k])
            evpprec = min(evplimit, evpprec_raw)
            if evpprec == evplimit and evplimit == evpprec_raw:
                bind = 0
            elif evpprec == max(0.0, (qs - q) / deltat):
                bind = 2
            elif evpprec == flxprec * grav / A['pdel'][k]:
                bind = 3
            else:
                bind = 4
            binder.append(bind)
            if flxprec > 0.0:
                work1 = min(max(0.0, flxsntm / flxprec), 1.0)
                evpsnow = evpprec * work1
            else:
                evpsnow = 0.0
            out.setdefault('evpvint', []).append(evpvint)
            evpvint += evpprec * A['pdel'][k] / grav
            ntprprd = A['rprd'][k] - evpprec
            if flxprec > 0.0:
                work1b = min(max(0.0, flxsnow / flxprec), 1.0)
            else:
                work1b = 0.0
            work2 = max(A['fsnow_conv'][k], work1b)
            if snowmlt > 0.0:
                work2 = 0.0
            ntsnprd = A['rprd'][k] * work2 - evpsnow - snowmlt
            tend_snwprd = A['rprd'][k] * work2 * latice
            tend_snwevmlt = -(evpsnow + snowmlt) * latice
            flxprec_hi = flxprec + ntprprd * A['pdel'][k] / grav
            flxsnow_hi = flxsnow + ntsnprd * A['pdel'][k] / grav
            flxprec_hi = max(flxprec_hi, 0.0)
            flxsnow_hi = max(flxsnow_hi, 0.0)
            tend_s = -evpprec * latvap + ntsnprd * latice
            tend_q = evpprec
            flxprec = flxprec_hi
            flxsnow = flxsnow_hi
            for key, val in [('ntprprd', ntprprd), ('ntsnprd', ntsnprd),
                             ('flxprec_hi', flxprec_hi),
                             ('flxsnow_hi', flxsnow_hi), ('tend_s', tend_s),
                             ('tend_q', tend_q),
                             ('tend_snwprd', tend_snwprd),
                             ('tend_snwevmlt', tend_snwevmlt)]:
                out.setdefault(key, []).append(val)
        out['prec_evap'] = flxprec / 1000.0
        out['snow_evap'] = flxsnow / 1000.0
        worst = 0.0
        for c in ['ntprprd', 'ntsnprd', 'flxprec_hi', 'flxsnow_hi',
                  'tend_s', 'tend_q', 'tend_snwprd', 'tend_snwevmlt']:
            e = check(out[c], A[c], f"step {st} evap {c}")
            worst = max(worst, e)
        pe = f(s[si['prec_evap']])
        e = check(out['prec_evap'], pe, f"step {st} evap prec_evap")
        worst = max(worst, e)
        print(f"step {st} evap: max rel err {worst:.3e} over the checked "
              f"columns and prec_evap; snow_evap={out['snow_evap']:.3e}")
        print(f"step {st} evap: binder histogram "
              f"{{0: {binder.count(0)}, 2: {binder.count(2)}, "
              f"3: {binder.count(3)}, 4: {binder.count(4)}}}  "
              f"levels: {[k + 1 for k in range(PVER) if binder[k] != 0]}")
        evap_tests[st] = dict(
            arrays={c: A[c] for c in ['t1', 'q1', 'pmid', 'pdel', 'cldfrc',
                                      'fsnow_conv', 'rprd', 'flxprec_lo',
                                      'flxsnow_lo']},
            evpvint=out['evpvint'],
            scalars=dict(prec_convr=prec_convr, gravit=grav, latvap=latvap,
                         latice=latice, cpres=cpres, ke=ke, deltat=deltat,
                         tmelt=tmelt),
            expects={c: out[c] for c in
                     ['ntprprd', 'ntsnprd', 'flxprec_hi', 'flxsnow_hi',
                      'tend_s', 'tend_q', 'tend_snwprd', 'tend_snwevmlt']},
            prec_evap=out['prec_evap'], snow_evap=out['snow_evap'])

    json.dump(q1q2_tests, open('scratch/q1q2_testdata.json', 'w'),
              indent=1, sort_keys=True)
    json.dump(evap_tests, open('scratch/evap_testdata.json', 'w'),
              indent=1, sort_keys=True)
    print('wrote scratch/q1q2_testdata.json and scratch/evap_testdata.json')


if __name__ == '__main__':
    main()
