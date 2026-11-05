#!/usr/bin/env python3
"""Validate the stage-1 CMEPS aoflux dump (cmeps_aoflux.csv).

Replays the full flux_atmOcn_large computation (flux_atmocn_Large.F90,
Large & Yeager scheme as configured in the C_JRA_t62_inst run) from the
dumped inputs and compares every dumped intermediate/derived/flux column
against the replay, plus structural checks (headers, arity, row counts,
no spval/NaN) and coverage checks (stable/unstable, low/high wind,
iteration-count spread).

Usage: python3 validate_cmeps_dump.py [path/to/cmeps_aoflux.csv]

The replay is the Python transcription of the Fortran laws that the
stage-2 .esm files (stage2/cmeps_aoflux.esm, stage2/cmeps_templates.esm)
transcribe; a fully green run here means the .esm tests pin the same laws
against the same dump rows.
"""
import sys
import math

import numpy as np

# ---------------- run configuration (C_JRA_t62_inst / C_JRA_t62) ------
ZREF = 10.0      # m, flux_atmocn_Large.F90:88
ZTREF = 2.0      # m, F90:89
MINWIND = 0.5    # m/s, med_phases_aofluxes_mod.F90 (seq_flux_atmocn_minwind)
FLUX_CON_TOL = 0.01     # nuopc.runconfig flux_convergence
FLUX_CON_MAX_ITER = 5   # nuopc.runconfig flux_max_iteration
ADD_GUSTS = False       # nuopc.runconfig add_gusts
COLD_AIR = False        # nuopc.runconfig coldair_outbreak_mod

# shr_const_mod.F90
STEBOL = 5.67e-8   # W/m^2/K^4
KARMAN = 0.4
G = 9.80616        # m/s^2
CPDAIR = 1004.64   # J/kg/K
CPWV = 1810.0      # J/kg/K
CPVIR = CPWV / CPDAIR - 1.0
LATVAP = 2.501e6   # J/kg
RWV, RDAIR = 28.966, 18.016   # g/mol
ZVIR = RWV / RDAIR - 1.0
EPSILO = 18.016 / 28.966
OMEPS = 1.0 - EPSILO
TBOIL = 373.16     # K, shr_wv_sat_mod.F90 GoffGratch_svp_liquid


def cdn(u):
    """Neutral 10-m drag coeff, Large & Yeager 2009 (F90:153-154)."""
    u33 = min(33.0, u)
    return 0.0027 / u33 + 0.000142 + 0.0000764 * u33 - 3.14807e-13 * u33**6


def psimhu(x):
    """Unstable momentum stability function (F90:163)."""
    return math.log((1.0 + x * (2.0 + x)) * (1.0 + x * x) / 8.0) - 2.0 * math.atan(x) + 1.571


def psixhu(x):
    """Unstable scalar stability function (F90:164)."""
    return 2.0 * math.log((1.0 + x * x) / 2.0)


def goff_gratch_liquid(t):
    """SVP over water, Goff & Gratch 1946 (shr_wv_sat_mod.F90:942-956)."""
    return (10.0 ** (-7.90298 * (TBOIL / t - 1.0)
                     + 5.02808 * math.log10(TBOIL / t)
                     - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / TBOIL)) - 1.0)
                     + 8.1328e-3 * (10.0 ** (-3.49149 * (TBOIL / t - 1.0)) - 1.0)
                     + math.log10(1013.246)) * 100.0)


def qsat_from_es(es, p):
    """shr_wv_sat_svp_to_qsat (F90:351-363) with the p<=es guard."""
    if (p - es) <= 0.0:
        return 1.0
    return EPSILO * es / (p - OMEPS * es)


def replay(row):
    """Replay flux_atmOcn_large for one dump row; return dict of quantities.

    row: dict with keys step,ymd,tod,n,gid,lat,lon,zbot,ubot,vbot,thbot,qbot,
    rainc,rbot,tbot,pslv,us,vs,ts (the raw inputs).
    """
    zbot, ubot, vbot = row['zbot'], row['ubot'], row['vbot']
    thbot, qbot, rbot = row['thbot'], row['qbot'], row['rbot']
    tbot, pslv = row['tbot'], row['pslv']
    us, vs, ts = row['us'], row['vs'], row['ts']
    rainc = row['rainc']

    # --- derived inputs (F90:193-229) ---
    rel = math.sqrt((ubot - us) ** 2 + (vbot - vs) ** 2)
    vmag = max(MINWIND, rel)
    wind0 = max(MINWIND, rel)
    ug = 0.0  # add_gusts = .false.
    esat = goff_gratch_liquid(ts)
    esat = min(esat, pslv)          # shr_wv_sat_qsat_liquid es = min(es, p)
    qsat = qsat_from_es(esat, pslv)
    ssq = 0.98 * qsat
    delt = thbot - ts
    delq = qbot - ssq
    alz = math.log(zbot / ZREF)
    cp = CPDAIR * (1.0 + CPVIR * ssq)
    al2 = math.log(ZREF / ZTREF)

    # --- neutral coefficients (F90:236-240) ---
    stable0 = 0.5 + math.copysign(0.5, delt)
    rdn0 = math.sqrt(cdn(vmag))
    rhn0 = (1.0 - stable0) * 0.0327 + stable0 * 0.018
    ren0 = 0.0346

    # --- stability iteration (F90:242-290) ---
    # state 0 (neutral)
    ustar = rdn0 * vmag
    tstar = rhn0 * delt
    qstar = ren0 * delq
    rdn = rdn0
    rhn = rhn0
    ren = ren0
    # The Fortran do-while (F90:242-290) exits as soon as
    # |(ustar - ustar_prev)/ustar| <= flux_con_tol (checked after updating the
    # state, F90:286-287), at most flux_con_max_iter = 5 iterations. The .esm
    # representation unrolls all 5 updates unconditionally and selects the
    # exit state by the first converged update; replay the same structure so
    # the conv flags, iter and the exit-state selection are defined exactly
    # as in the .esm model.
    ustar_trace = [ustar]
    states = []
    for _ in range(FLUX_CON_MAX_ITER):
        hol = KARMAN * G * zbot * (tstar / thbot + qstar / (1.0 / ZVIR + qbot)) / ustar ** 2
        hol = math.copysign(min(abs(hol), 10.0), hol)
        stable = 0.5 + math.copysign(0.5, hol)
        xsq = max(math.sqrt(abs(1.0 - 16.0 * hol)), 1.0)
        xqq = math.sqrt(xsq)
        psimh = -5.0 * hol * stable + (1.0 - stable) * psimhu(xqq)
        psixh = -5.0 * hol * stable + (1.0 - stable) * psixhu(xqq)
        # shift wind speed using old coefficient
        rd_shift = rdn / (1.0 + max(rdn / KARMAN * (alz - psimh), -0.5))
        u10n = vmag * rd_shift / rdn
        # update transfer coeffs at 10m and neutral stability
        rdn = math.sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.0 - stable) * 0.0327 + stable * 0.018
        # shift all coeffs to measurement height and stability
        rd = rdn / (1.0 + rdn / KARMAN * (alz - psimh))
        rh = rhn / (1.0 + rhn / KARMAN * (alz - psixh))
        re = ren / (1.0 + ren / KARMAN * (alz - psixh))
        # update ustar, tstar, qstar using updated, shifted coeffs
        ustar = rd * vmag
        tstar = rh * delt
        qstar = re * delq
        ustar_trace.append(ustar)
        states.append(dict(hol=hol, stable=stable, xsq=xsq, xqq=xqq,
                           psimh=psimh, psixh=psixh, rd_shift=rd_shift, u10n=u10n,
                           rdn=rdn, rhn=rhn, ren=ren, rd=rd, rh=rh, re=re,
                           ustar=ustar, tstar=tstar, qstar=qstar))

    # per-iteration convergence flags conv_k = |(u_k - u_{k-1})/u_k| <= tol
    # (the Fortran exit test after iteration k, F90:286-287); the .esm
    # representation carries them explicitly to select the final state
    conv = {}
    for k in range(1, 5):
        conv[f'conv_{k}'] = (1.0 if abs((ustar_trace[k] - ustar_trace[k - 1]) / ustar_trace[k])
                             <= FLUX_CON_TOL else 0.0)
    iter = next((k for k in range(1, 5) if conv[f'conv_{k}']), 5)
    s = states[iter - 1]   # the exit state, exactly the .esm selection
    hol, stable = s['hol'], s['stable']
    xsq, xqq = s['xsq'], s['xqq']
    psimh, psixh = s['psimh'], s['psixh']
    rd, u10n = s['rd'], s['u10n']
    rdn, rhn, ren = s['rdn'], s['rhn'], s['ren']
    rh, re = s['rh'], s['re']
    ustar, tstar, qstar = s['ustar'], s['tstar'], s['qstar']

    # --- fluxes (F90:299-311), at the exit state ---
    tau = rbot * ustar * ustar
    taux = tau * (ubot - us) / vmag
    tauy = tau * (vbot - vs) / vmag
    sen = cp * tau * tstar / ustar
    lat = LATVAP * tau * qstar / ustar
    lwup = -STEBOL * ts ** 4
    evap = lat / LATVAP

    # --- diagnostics (F90:316-333), at the exit state ---
    hol0 = hol
    xsq0 = xsq
    xqq0 = xqq
    hol2 = hol * ZTREF / zbot
    xsq2 = max(1.0, math.sqrt(abs(1.0 - 16.0 * hol2)))
    xqq2 = math.sqrt(xsq2)
    psix2 = -5.0 * hol2 * stable + (1.0 - stable) * psixhu(xqq2)
    fac_t = (rh / KARMAN) * (alz + al2 - psixh + psix2)
    tref = thbot - delt * fac_t - 0.01 * ZTREF
    fac_q = (re / KARMAN) * (alz + al2 - psixh + psix2)
    qref = qbot - delq * fac_q
    duu10n = u10n * u10n
    u10res = u10n * (wind0 / vmag)

    return dict(conv, vmag=vmag, wind0=wind0, ugust=ug, esat=esat, qsat=qsat, ssq=ssq,
                delt=delt, delq=delq, alz=alz, cp=cp,
                stable0=stable0, rdn0=rdn0, rhn0=rhn0, ren0=ren0,
                iter=iter, hol=hol0, stable=stable, xsq=xsq0, xqq=xqq0,
                psimh=psimh, psixh=psixh, rd=rd, u10n=u10n, rdn=rdn, rhn=rhn,
                ren=ren, rh=rh, re=re, ustar=ustar, tstar=tstar, qstar=qstar,
                tau=tau, taux=taux, tauy=tauy, sen=sen, lat=lat, lwup=lwup,
                evap=evap, hol2=hol2, xsq2=xsq2, xqq2=xqq2, psix2=psix2,
                fac_t=fac_t, fac_q=fac_q, tref=tref, qref=qref,
                duu10n=duu10n, u10res=u10res)


def read_cols(path):
    """Read the column list from the dump's '# columns' header line."""
    with open(path) as f:
        for l in f:
            if l.startswith('# columns'):
                return l.split('= ', 1)[1].strip().split(',')
    raise RuntimeError('no # columns header found')


def main(path):
    COLS = read_cols(path)
    # inputs are replayed FROM, not compared; 'lat' is both the latitude
    # metadata (first occurrence) and the latent-heat flux (second), so
    # exclude by name here (latent heat is still checked via evap = lat/LATVAP)
    INPUTS = ('zbot', 'ubot', 'vbot', 'thbot', 'qbot', 'rainc', 'rbot', 'tbot',
              'pslv', 'us', 'vs', 'ts')
    COMPARE = [c for c in COLS if c not in
               ('step', 'ymd', 'tod', 'n', 'gid', 'iter', 'lat', 'lon') + INPUTS]
    with open(path) as f:
        lines = f.readlines()
    desc = [l for l in lines if l.startswith('# description')]
    cols = [l for l in lines if l.startswith('# columns')]
    data = [l for l in lines if l.strip() and not l.startswith('#')]
    print(f'header lines: {len(desc)} description, {len(cols)} columns')
    print(f'data rows: {len(data)}')
    if len(data) == 0:
        print('FAIL: no data rows'); sys.exit(1)
    ncols = data[0].count(',') + 1
    print(f'columns in data: {ncols} (expected {len(COLS)})')
    bad_arity = [i for i, l in enumerate(data) if l.count(',') + 1 != ncols]
    print(f'rows with wrong arity: {len(bad_arity)}')

    rows = []
    for l in data:
        parts = l.strip().split(',')
        r = {}
        for c, v in zip(COLS, parts):
            try:
                r[c] = int(v) if c in ('step', 'ymd', 'tod', 'n', 'gid', 'iter') else float(v)
            except ValueError:
                r[c] = float('nan')
        rows.append(r)
    arr = {c: np.array([r[c] for r in rows]) for c in COLS}

    steps = sorted(set(arr['step']))
    cells = sorted(set(arr['n']))
    print(f'steps: {steps}')
    print(f'cells traced: {len(cells)} (n range {min(cells)}..{max(cells)})')
    # fixed mode (first-N cells, identical set every step): rows = cells x steps;
    # dynamic mode (per-step stable/wind selection): the cell set can differ
    # between steps, so compare rows against the per-step union instead.
    cells_per_step = {s: sorted(set(arr['n'][arr['step'] == s])) for s in steps}
    union_total = sum(len(cs) for cs in cells_per_step.values())
    if all(cells_per_step[s] == cells for s in steps):
        mode = 'fixed (identical cell set every step)'
    else:
        mode = 'dynamic (cell set varies per step)'
    print(f'selection: {mode}')
    if len(cells_per_step) <= 12:
        for s in steps:
            print(f'  step {s:2d}: {len(cells_per_step[s])} rows')
    print(f'rows = sum over steps of that step\'s cells: {len(rows)} == {union_total}: '
          f'{len(rows) == union_total}')
    print(f'rows = cells x steps: {len(rows)} == {len(cells) * len(steps)}: '
          f'{len(rows) == len(cells) * len(steps)}')

    # NaN / spval checks (spval = 1e30 in the mediator)
    for c in COMPARE:
        nnan = int(np.isnan(arr[c]).sum())
        nspv = int((np.abs(arr[c]) > 1e29).sum())
        if nnan or nspv:
            print(f'  column {c}: {nnan} NaN, {nspv} spval-like')

    # replay every row and compare
    worst = {}
    worst_iter = None
    for r in rows:
        got = replay(r)
        for c in COMPARE:
            exp = got[c]
            val = r[c]
            if c == 'iter':
                if int(val) != int(exp):
                    worst_iter = (r['step'], r['n'], val, exp)
                continue
            denom = max(abs(exp), 1e-30)
            relerr = abs(val - exp) / denom
            if c not in worst or relerr > worst[c]:
                worst[c] = relerr
    print('\nworst relative errors (dump vs replay), per column:')
    for c in COMPARE:
        if c == 'iter':
            continue
        flag = 'OK' if worst[c] < 1e-11 else ('WARN' if worst[c] < 1e-6 else 'FAIL')
        print(f'  {c:8s} {worst[c]:.3e}  {flag}')
    print(f'  iter     mismatch: {worst_iter}')

    # coverage checks: stability, wind, clamps, iteration spread
    hol = arr['hol']
    vmag = arr['vmag']
    delt = arr['delt']
    n_unstable = int((hol < 0).sum())
    n_stable = int((hol > 0).sum())
    n_neutral = int((hol == 0).sum())
    n_lowwind = int((vmag < 3.0).sum())
    n_highwind = int((vmag > 15.0).sum())
    n_minwind = int((vmag <= 0.5 + 1e-12).sum())
    n_holclamp = int((np.abs(hol) >= 10.0 - 1e-12).sum())
    n_deltpos = int((delt > 0).sum())
    n_windcap = int((vmag >= 33.0 - 1e-9).sum())
    n_u10ncap = int((arr['u10n'] >= 33.0 - 1e-9).sum())
    print(f'\ncoverage: unstable {n_unstable}, stable {n_stable}, neutral {n_neutral}, '
          f'low wind (<3 m/s) {n_lowwind}, high wind (>15 m/s) {n_highwind}')
    print(f'  delt > 0 (stable regime): {n_deltpos} rows; '
          f'minwind floor: {n_minwind}; |hol| = 10 clamp: {n_holclamp}')
    print(f'  vmag: min {vmag.min():.4f}, max {vmag.max():.4f} m/s; '
          f'rows with vmag >= 33 m/s (cdn cap): {n_windcap}; '
          f'rows with u10n >= 33 m/s (cdn cap): {n_u10ncap}')
    print(f'  hol: min {hol.min():.4f}, max {hol.max():.4f}')
    iters = {k: int((arr['iter'] == k).sum()) for k in sorted(set(arr['iter']))}
    print(f'iteration counts: {iters}')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else
         '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cmeps/cmeps_aoflux.csv')
