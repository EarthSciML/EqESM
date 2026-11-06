#!/usr/bin/env python3
"""Validate the stage-1 CMEPS ocean-albedo and ocn-export merge dumps
(cmeps_ocnalb.csv, cmeps_merge.csv).

Replays, from the dumped inputs:

* cmeps_ocnalb.csv -- the med_phases_ocnalb_run computation
  (mediator/med_phases_ocnalb_mod.F90, Briegleb et al. 1986 scheme as
  configured in the C_JRA_t62_inst run): the solar declination and
  earth-sun distance factor of shr_orb_decl (share/src/shr_orb_mod.F90,
  Berger 1978; dayspy = 365, vernal equinox at calday 80.5), the cosine
  of the solar zenith angle of shr_orb_cosz, and the four ocean
  albedos, and compares them against the dumped columns.

* cmeps_merge.csv -- the fraction-weighted ocn-export merges of
  med_phases_prep_ocn_accum (mediator/med_phases_prep_ocn_mod.F90:
  med_merge_auto with the merge specs of
  esmFldsExchange_cesm_mod.F90, plus the med_phases_prep_ocn_custom
  shortwave-net computation) replayed from the dumped FBImp/FBMed
  inputs, plus the derived columns the dump writer recomputes with the
  exact custom-computation expressions.

Usage: python3 validate_cmeps_ocnalb_merge_dumps.py [ocnalb.csv] [merge.csv]

A fully green run here means the stage-2 .esm files
(stage2/cmeps_ocnalb.esm, stage2/cmeps_merge.esm,
stage2/cmeps_templates.esm) pin the same laws against the same dump
rows.  The dark-side branch of the albedo scheme (cosz <= 0 -> all
albedos 1.0) and the Taylor 1996 scheme are transcribed in the replay
but unpinned by this dump (all traced cells are day-side, Briegleb
scheme).  The enthalpy fluxes (Foxx_h*, Foxx_hevap etc.) are computed
outside the merges dumped here and are a documented gap.
"""
import sys
import math

import numpy as np

# ---------------- run configuration (C_JRA_t62_inst) -------------------
# Briegleb et al. (1986) ocean albedo scheme (ocean_albedo_scheme = 0),
# no albedo floor, no flux averaging:
USE_BRIEGLEB = True
FLUX_ALBAV = False
# med_phases_prep_ocn_custom shortwave band split (F90:818): the
# Fortran literals c1..c4 are written without an _r8 kind suffix, so
# they are SINGLE-precision literals promoted to r8 (0.285f =
# 0.2849999964237213134765625, 0.215f = 0.2150000035762786865234375);
# the dump pins this (vdr/swnet = 0.28499999642372131 in every row)
C1, C2 = np.float32(0.285).astype(float), np.float32(0.285).astype(float)
C3, C4 = np.float32(0.215).astype(float), np.float32(0.215).astype(float)

# shr_orb_mod.F90
DAYS_PER_YEAR = 365.0     # dayspy
VERNAL_EQUINOX = 80.5     # ve, calday of vernal equinox (Jan 1 = day 1)
PI = math.pi


def shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr):
    """Solar declination [rad] and earth-sun distance factor (Berger 1978,
    shr_orb_mod.F90:736-814)."""
    lambm = lambm0 + (calday - VERNAL_EQUINOX) * 2.0 * PI / DAYS_PER_YEAR
    lmm = lambm - mvelpp
    sinl = math.sin(lmm)
    lamb = (lambm + eccen * (2.0 * sinl + eccen * (1.25 * math.sin(2.0 * lmm)
            + eccen * ((13.0 / 12.0) * math.sin(3.0 * lmm) - 0.25 * sinl))))
    invrho = (1.0 + eccen * math.cos(lamb - mvelpp)) / (1.0 - eccen * eccen)
    delta = math.asin(math.sin(obliqr) * math.sin(lamb))
    eccf = invrho * invrho
    return delta, eccf


def shr_orb_cosz(jday, lat, lon, declin):
    """Cosine of the solar zenith angle (shr_orb_mod.F90:110-159, the
    no-dt_avg branch)."""
    return (math.sin(lat) * math.sin(declin)
            - math.cos(lat) * math.cos(declin)
            * math.cos((jday - math.floor(jday)) * 2.0 * PI + lon))


def briegleb_anidr(cosz):
    """Briegleb et al. (1986) direct visible ocean albedo (F90:476-478)."""
    return (0.026 / (cosz ** 1.7 + 0.065)
            + 0.150 * (cosz - 0.100) * (cosz - 0.500) * (cosz - 1.000))


def taylor_anidr(cosz):
    """Taylor et al. (1996) direct visible ocean albedo (F90:480-482)."""
    return 0.037 / (1.1 * cosz ** 1.4 + 0.15)


def replay_ocnalb(row):
    """Replay the ocnalb phase for one dump row; return the computed
    quantities.  row: dict of the dumped columns."""
    nextsw_cday = row['nextsw_cday']
    eccen, obliqr = row['eccen'], row['obliqr']
    lambm0, mvelpp = row['lambm0'], row['mvelpp']
    albdif, albdir = row['albdif'], row['albdir']
    min_albedo = row['min_albedo']
    use_min_albedo = row['use_min_albedo'] > 0.5
    use_briegleb = row['use_briegleb'] > 0.5
    flux_albav = row['flux_albav'] > 0.5
    rlat, rlon = row['rlat'], row['rlon']

    delta, eccf = shr_orb_decl(nextsw_cday, eccen, mvelpp, lambm0, obliqr)
    if flux_albav:
        # flat mean albedos (F90:453-459)
        return dict(delta=delta, eccf=eccf, cosz=0.0,
                    anidr=albdir, avsdr=albdir, anidf=albdif, avsdf=albdif)
    cosz = shr_orb_cosz(nextsw_cday, rlat, rlon, delta)
    if cosz > 0.0:
        if use_briegleb:
            anidr = briegleb_anidr(cosz)
        else:
            anidr = taylor_anidr(cosz)
        if use_min_albedo:
            anidr = max(anidr, min_albedo)
        avsdr = anidr
        anidf = albdif
        avsdf = albdif
    else:
        # dark side of the earth (F90:489-494)
        anidr = avsdr = anidf = avsdf = 1.0
    return dict(delta=delta, eccf=eccf, cosz=cosz,
                anidr=anidr, avsdr=avsdr, anidf=anidf, avsdf=avsdf)


def replay_merge(row):
    """Replay the fraction-weighted merges and the derived shortwave
    columns for one dump row.  The Fortran accumulates
    dp = dp + src*weight (med_merge_mod.F90:405-412) in the source
    order of esmFldsExchange_cesm_mod.F90, so the replay preserves the
    same term order and association."""
    ifrac, ofrac = row['ifrac'], row['ofrac']
    ifrad, ofrad = row['ifrad'], row['ofrad']

    # med_merge_auto (med_phases_prep_ocn_mod.F90:161-168), merge specs
    # of esmFldsExchange_cesm_mod.F90
    faxa_rain = row['Faxa_rainc'] * ofrac + row['Faxa_rainl'] * ofrac      # sum_with_weights ofrac (L1940-1942)
    faxa_snow = row['Faxa_snowc'] * ofrac + row['Faxa_snowl'] * ofrac      # sum_with_weights ofrac (L1949-1951)
    foxx_lwnet = row['Faox_lwup'] * ofrac + row['Faxa_lwdn'] * ofrac       # merge ofrac x2 (L1759-1762)
    foxx_sen = row['Faox_sen'] * ofrac                                      # merge ofrac (L1966-1967)
    foxx_taux = row['Fioi_taux_ice'] * ifrac + row['Faox_taux'] * ofrac    # merge ifrac + ofrac (L2154-2158)
    foxx_tauy = row['Fioi_tauy_ice'] * ifrac + row['Faox_tauy'] * ofrac    # merge ifrac + ofrac (L2169-2173)
    foxx_evap = row['Faox_evap'] * ofrac                                    # merge ofrac (L1992-1993)
    fioi_melth = row['Fioi_melth_ice'] * ifrac                             # copy_with_weights ifrac (L2204-2206)
    fioi_meltw = row['Fioi_meltw_ice'] * ifrac                             # copy_with_weights ifrac (L2189-2191)
    fioi_salt = row['Fioi_salt_ice'] * ifrac                               # copy_with_weights ifrac (L2218-2220)
    foxx_rofl = row['Forr_rofl']                                            # sum (L2354)
    foxx_rofi = row['Forr_rofi']                                            # sum (L2394)

    # med_phases_prep_ocn_custom (F90:781-825): absorbed shortwave,
    # scaled fractions and the merged net shortwave
    fswabsv = row['Faxa_swvdr'] * (1.0 - row['So_avsdr']) + row['Faxa_swvdf'] * (1.0 - row['So_avsdf'])
    fswabsi = row['Faxa_swndr'] * (1.0 - row['So_anidr']) + row['Faxa_swndf'] * (1.0 - row['So_anidf'])
    frac_sum = ifrac + ofrac
    ifrac_scaled = ifrac / frac_sum if frac_sum != 0.0 else ifrac
    ofrac_scaled = ofrac / frac_sum if frac_sum != 0.0 else ofrac
    frac_sum = ifrad + ofrad
    ifracr_scaled = ifrad / frac_sum if frac_sum != 0.0 else ifrad
    ofracr_scaled = ofrad / frac_sum if frac_sum != 0.0 else ofrad
    swnet = ofracr_scaled * (fswabsv + fswabsi) + ifrac_scaled * row['Fioi_swpen_ice']
    swnet_vdr = C1 * swnet
    swnet_vdf = C2 * swnet
    swnet_idr = C3 * swnet
    swnet_idf = C4 * swnet

    return dict(faxa_rain=faxa_rain, faxa_snow=faxa_snow,
                foxx_lwnet=foxx_lwnet, foxx_sen=foxx_sen,
                foxx_taux=foxx_taux, foxx_tauy=foxx_tauy, foxx_evap=foxx_evap,
                fioi_melth=fioi_melth, fioi_meltw=fioi_meltw, fioi_salt=fioi_salt,
                foxx_rofl=foxx_rofl, foxx_rofi=foxx_rofi,
                fswabsv=fswabsv, fswabsi=fswabsi,
                ifrac_scaled=ifrac_scaled, ofrac_scaled=ofrac_scaled,
                ifracr_scaled=ifracr_scaled, ofracr_scaled=ofracr_scaled,
                swnet=swnet, swnet_vdr=swnet_vdr, swnet_vdf=swnet_vdf,
                swnet_idr=swnet_idr, swnet_idf=swnet_idf)


def read_cols(path):
    """Read the column list from the dump's '# columns' header line."""
    with open(path) as f:
        for l in f:
            if l.startswith('# columns'):
                return l.split('= ', 1)[1].strip().split(',')
    raise RuntimeError(f'no # columns header found in {path}')


def load(path):
    """Structural checks and row dicts for one dump."""
    COLS = read_cols(path)
    with open(path) as f:
        lines = f.readlines()
    desc = [l for l in lines if l.startswith('# description')]
    cols = [l for l in lines if l.startswith('# columns')]
    data = [l for l in lines if l.strip() and not l.startswith('#')]
    print(f'  header lines: {len(desc)} description, {len(cols)} columns')
    print(f'  data rows: {len(data)}')
    if len(data) == 0:
        # an empty dump (headers only, e.g. a phase whose traced rows were
        # never written) is reported by the caller, not a hard error here
        return COLS, [], {}, {}, []
    ncols = data[0].count(',') + 1
    print(f'  columns in data: {ncols} (expected {len(COLS)})')
    bad_arity = [i for i, l in enumerate(data) if l.count(',') + 1 != ncols]
    print(f'  rows with wrong arity: {len(bad_arity)}')
    rows = []
    for l in data:
        parts = l.strip().split(',')
        r = {}
        for c, v in zip(COLS, parts):
            try:
                r[c] = int(v) if c in ('step', 'ymd', 'tod', 'n', 'gid') else float(v)
            except ValueError:
                r[c] = float('nan')
        rows.append(r)
    arr = {c: np.array([r[c] for r in rows]) for c in COLS}
    steps = sorted(set(arr['step']))
    cells = sorted(set(arr['n']))
    cells_per_step = {s: sorted(set(arr['n'][arr['step'] == s])) for s in steps}
    union_total = sum(len(cs) for cs in cells_per_step.values())
    print(f'  steps: {steps}')
    print(f'  cells traced: {len(cells)} (n range {min(cells)}..{max(cells)})')
    print(f'  rows = sum over steps of that step\'s cells: {len(rows)} == {union_total}: '
          f'{len(rows) == union_total}')
    for c in COLS:
        nnan = int(np.isnan(arr[c]).sum())
        nspv = int((np.abs(arr[c]) > 1e29).sum())
        if nnan or nspv:
            print(f'  column {c}: {nnan} NaN, {nspv} spval-like')
    return COLS, rows, arr, cells_per_step, cells


def compare(rows, replay, mapping, label):
    """Replay every row, compare mapped columns; return worst relerr dict."""
    worst = {}
    for r in rows:
        got = replay(r)
        for out, inp in mapping.items():
            exp = got[inp]
            val = r[out]
            denom = max(abs(exp), 1e-30)
            relerr = abs(val - exp) / denom
            if out not in worst or relerr > worst[out]:
                worst[out] = relerr
    print(f'\n  worst relative errors ({label}, dump vs replay), per column:')
    for c in sorted(worst):
        flag = 'OK' if worst[c] < 1e-11 else ('WARN' if worst[c] < 1e-6 else 'FAIL')
        print(f'    {c:12s} {worst[c]:.3e}  {flag}')
    return worst


def main(ocnalb_path, merge_path):
    print('=' * 70)
    print(f'ocnalb dump: {ocnalb_path}')
    print('=' * 70)
    COLS, rows, arr, cps, cells = load(ocnalb_path)
    if len(rows) == 0:
        print('  SKIP: no ocnalb data rows (dump empty; nothing to replay)')
        worst_o = {}
    else:
        # replay: delta, eccf, cosz and the four albedos from the dumped inputs
        ocnalb_map = {'delta': 'delta', 'eccf': 'eccf', 'cosz': 'cosz',
                      'anidr': 'anidr', 'avsdr': 'avsdr', 'anidf': 'anidf',
                      'avsdf': 'avsdf'}
        worst_o = compare(rows, replay_ocnalb, ocnalb_map, 'ocnalb')

    # internal consistency of the dumped rows
    print('\n  consistency checks:')
    if len(rows) == 0:
        print('    (skipped: no rows)')
    else:
        n_day = int((arr['cosz'] > 0).sum())
        n_dark = int((arr['cosz'] <= 0).sum())
        print(f'    day-side rows (cosz > 0): {n_day}; dark-side rows: {n_dark}')
        print(f'    avsdr == anidr everywhere: {bool(np.all(arr["avsdr"] == arr["anidr"]))}')
        print(f'    anidf == avsdf == albdif (day side): '
              f'{bool(np.all(arr["anidf"][arr["cosz"] > 0] == arr["albdif"][arr["cosz"] > 0]))}')
        # the ocnalb phase copies ifrac/ofrac into ifrad/ofrad AFTER the
        # albedo loop (F90:502-521), so at the dump point (inside the loop)
        # ifrad/ofrad hold the PREVIOUS step's ice/open fractions: verify
        # ifrad(step k) == ifrac(step k-1) per cell for k > 1 (step 1 shows
        # the init-time equality ifrad = ifrac)
        lag_ok = True
        prev = {}
        for s in sorted(set(r['step'] for r in rows)):
            cur = {}
            for r in rows:
                if r['step'] != s:
                    continue
                cur[r['n']] = r
            for n, r in cur.items():
                if n in prev and (r['ifrad'] != prev[n]['ifrac'] or r['ofrad'] != prev[n]['ofrac']):
                    lag_ok = False
            prev = cur
        print(f'    ifrad/ofrad lag by one step (ifrad(k) == ifrac(k-1)): {lag_ok}')
        print(f'    ifrac + ofrac == 1 (tol 1e-12): {bool(np.all(np.abs(arr["ifrac"] + arr["ofrac"] - 1.0) < 1e-12))}')
        print(f'    use_briegleb all 1: {bool(np.all(arr["use_briegleb"] == 1.0))}, '
              f'flux_albav all 0: {bool(np.all(arr["flux_albav"] == 0.0))}, '
              f'use_min_albedo all 0: {bool(np.all(arr["use_min_albedo"] == 0.0))}')
        print(f'    cosz range: {arr["cosz"].min():.5f} .. {arr["cosz"].max():.5f}; '
              f'delta range: {arr["delta"].min():.5f} .. {arr["delta"].max():.5f} rad')

    print()
    print('=' * 70)
    print(f'merge dump: {merge_path}')
    print('=' * 70)
    COLS, rows, arr, cps, cells = load(merge_path)
    merge_map = {'Faxa_rain': 'faxa_rain', 'Faxa_snow': 'faxa_snow',
                 'Foxx_lwnet': 'foxx_lwnet', 'Foxx_sen': 'foxx_sen',
                 'Foxx_taux': 'foxx_taux', 'Foxx_tauy': 'foxx_tauy',
                 'Foxx_evap': 'foxx_evap',
                 'Fioi_melth': 'fioi_melth', 'Fioi_meltw': 'fioi_meltw',
                 'Fioi_salt': 'fioi_salt',
                 'Foxx_rofl': 'foxx_rofl', 'Foxx_rofi': 'foxx_rofi',
                 'fswabsv': 'fswabsv', 'fswabsi': 'fswabsi',
                 'ifrac_scaled': 'ifrac_scaled', 'ofrac_scaled': 'ofrac_scaled',
                 'ifracr_scaled': 'ifracr_scaled', 'ofracr_scaled': 'ofracr_scaled',
                 'Foxx_swnet': 'swnet',
                 'Foxx_swnet_vdr': 'swnet_vdr', 'Foxx_swnet_vdf': 'swnet_vdf',
                 'Foxx_swnet_idr': 'swnet_idr', 'Foxx_swnet_idf': 'swnet_idf'}
    worst_m = compare(rows, replay_merge, merge_map, 'merge')

    print('\n  consistency checks:')
    print(f'    ifrac + ofrac == 1 (tol 1e-12): {bool(np.all(np.abs(arr["ifrac"] + arr["ofrac"] - 1.0) < 1e-12))}')
    print(f'    ifrad + ofrad == 1 (tol 1e-12): {bool(np.all(np.abs(arr["ifrad"] + arr["ofrad"] - 1.0) < 1e-12))}')
    print(f'    ifrac_scaled + ofrac_scaled == 1: '
          f'{bool(np.all(np.abs(arr["ifrac_scaled"] + arr["ofrac_scaled"] - 1.0) < 1e-12))}')
    print(f'    ifracr_scaled + ofracr_scaled == 1: '
          f'{bool(np.all(np.abs(arr["ifracr_scaled"] + arr["ofracr_scaled"] - 1.0) < 1e-12))}')
    fracs = np.unique(np.round(arr['ifrac'], 4))
    print(f'    distinct ifrac values: {len(fracs)} (min {fracs.min():.4f}, max {fracs.max():.4f})')
    n_partial = int(((arr['ifrac'] > 0) & (arr['ifrac'] < 1)).sum())
    n_open = int((arr['ifrac'] == 0).sum())
    n_ice = int((arr['ifrac'] == 1).sum())
    print(f'    partial-ice rows: {n_partial}, open-ocean rows: {n_open}, fully-iced rows: {n_ice}')

    ok_o = all(w < 1e-11 for w in worst_o.values())
    ok_m = all(w < 1e-11 for w in worst_m.values())
    print()
    print(f'RESULT: ocnalb {"PASS" if ok_o else "FAIL"}, merge {"PASS" if ok_m else "FAIL"}')
    sys.exit(0 if (ok_o and ok_m) else 1)


if __name__ == '__main__':
    base = '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cmeps/'
    main(sys.argv[1] if len(sys.argv) > 1 else base + 'cmeps_ocnalb.csv',
         sys.argv[2] if len(sys.argv) > 2 else base + 'cmeps_merge.csv')
