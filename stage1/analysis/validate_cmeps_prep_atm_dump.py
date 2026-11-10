#!/usr/bin/env python3
"""Validate the stage-1 CMEPS atmosphere-export merge dump
(cmeps_prep_atm.csv).

Replays, from the dumped inputs, the fraction-weighted atmosphere
exports of med_phases_prep_atm (mediator/med_phases_prep_atm_mod.F90,
as configured in the F2000climo_f19_inst run: coupling_mode = cesm,
aoflux_grid = ogrid, component_computes_enthalpy_flux = none,
mrg_fracname_lnd2atm_state = 'lfrac', mrg_fracname_lnd2atm_flux =
'lfrac', samegrid_atmlnd = .true.):

* the med_merge_auto accumulation dp = dp + src*weight of
  mediator/med_merge_mod.F90:405-412 applied with the merge lists of
  mediator/esmFldsExchange_cesm_mod.F90, one replay per dumped export
  field: Sx_* = Sl_* * lfrac + Si_* * ifrac + So_* * ofrac (source
  order lnd, ice, ocn, weights lfrac/ifrac/ofrac) for the albedo
  group (Sx_avsdr/Sx_avsdf/Sx_anidr/Sx_anidf from the Sl_* land
  fields, the Si_* ice fields and the So_* ocean albedos of
  FBMed_ocnalb_a; esmFldsExchange L875-887, 902-914, 929-941,
  956-968), the reference group (Sx_tref/Sx_u10/Sx_qref with the So_*
  of FBMed_aoflux_a; L988-1002, 1014-1029, 1045-1063), the flux group
  (Faxx_taux/Faxx_tauy [N/m2], Faxx_lat/Faxx_sen/Faxx_lwup [W/m2],
  Faxx_evap [kg/m2/s] with the Fall_* land, Faii_* ice and Faox_*
  aoflux sources; L1085-1099, 1112-1126, 1139-1153, 1166-1180,
  1196-1211, 1236-1250) and the surface-temperature group (Sx_t with
  the So_t ocn source mapped to the atm grid; L1267-1279), plus the
  copy export So_t = So_t source (mrg_type = 'copy', L1282).

The enthalpy branches of the phase (component_computes_enthalpy_flux
= 'med'/'atm') and the Faoo_f* ocean-gas exports (scaled by ifrac +
ofrac in the custom section, F90:220-248) are not exercised by the
traced run (CAM's F2000 import list carries neither) and are a
documented gap.  The dump writes the final export values AFTER the
custom calculations (the fraction copies So_ofrac/Si_ifrac/Sl_lfrac
of F90:169-218 write the same numbers the dump reads from FBFrac), so
the replayed identities close exactly.

Usage: python3 validate_cmeps_prep_atm_dump.py [prep_atm.csv]

A fully green run here means the stage-2 .esm file
(stage2/cmeps_prep_atm.esm, which applies the fraction_weighted_sum
template of cmeps_templates.esm per field) pins the same laws against
the same dump rows.
"""
import sys

import numpy as np

# the fraction-weighted exports of the cesm coupling mode, as dumped:
# for each field the land/ice/ocn source columns and the merged export
# column.  The replay preserves the F90 accumulation order
# (dp = dp + src*weight over the esmFldsExchange merge lists, source
# order lnd, ice, ocn) so the comparison closes to round-off.
MERGE_GROUPS = {
    'avsdr': ('Sl_avsdr', 'Si_avsdr', 'So_avsdr', 'Sx_avsdr'),
    'avsdf': ('Sl_avsdf', 'Si_avsdf', 'So_avsdf', 'Sx_avsdf'),
    'anidr': ('Sl_anidr', 'Si_anidr', 'So_anidr', 'Sx_anidr'),
    'anidf': ('Sl_anidf', 'Si_anidf', 'So_anidf', 'Sx_anidf'),
    'tref': ('Sl_tref', 'Si_tref', 'So_tref', 'Sx_tref'),
    'u10': ('Sl_u10', 'Si_u10', 'So_u10', 'Sx_u10'),
    'qref': ('Sl_qref', 'Si_qref', 'So_qref', 'Sx_qref'),
    'taux': ('Fall_taux', 'Faii_taux', 'Faox_taux', 'Faxx_taux'),
    'tauy': ('Fall_tauy', 'Faii_tauy', 'Faox_tauy', 'Faxx_tauy'),
    'lat': ('Fall_lat', 'Faii_lat', 'Faox_lat', 'Faxx_lat'),
    'sen': ('Fall_sen', 'Faii_sen', 'Faox_sen', 'Faxx_sen'),
    'evap': ('Fall_evap', 'Faii_evap', 'Faox_evap', 'Faxx_evap'),
    'lwup': ('Fall_lwup', 'Faii_lwup', 'Faox_lwup', 'Faxx_lwup'),
    't': ('Sl_t', 'Si_t', 'So_t', 'Sx_t'),
}
# fraction columns of the dump, in the weight order of the merges
FRAC_COLS = ('lfrac', 'ifrac', 'ofrac')


def read_cols(path):
    """Read the column list from the dump's '# columns' header line,
    renaming the duplicated 'So_t' (the source and the copy export) to
    'So_t_exp'."""
    with open(path) as f:
        for l in f:
            if l.startswith('# columns'):
                cols = l.split('= ', 1)[1].strip().split(',')
                seen = set()
                out = []
                for c in cols:
                    if c in seen:
                        out.append(c + '_exp')
                    else:
                        seen.add(c)
                        out.append(c)
                return out
    raise RuntimeError(f'no # columns header found in {path}')


def load(path):
    """Structural checks and row dicts for the dump."""
    COLS = read_cols(path)
    with open(path) as f:
        lines = f.readlines()
    desc = [l for l in lines if l.startswith('# description')]
    cols = [l for l in lines if l.startswith('# columns')]
    data = [l for l in lines if l.strip() and not l.startswith('#')]
    print(f'  header lines: {len(desc)} description, {len(cols)} columns')
    print(f'  data rows: {len(data)}')
    if len(data) == 0:
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


def replay_merge(row):
    """Replay the fraction-weighted merges of one dump row (the F90
    accumulation dp = dp + src*weight in the source order of
    esmFldsExchange_cesm_mod.F90: lnd, ice, ocn)."""
    out = {}
    for suffix, (sl, si, so, sx) in MERGE_GROUPS.items():
        out[sx] = (row[sl] * row['lfrac'] + row[si] * row['ifrac']) + row[so] * row['ofrac']
    return out


def compare(rows, mapping, label):
    """Replay every row, compare mapped columns; return worst relerr dict."""
    worst = {}
    for r in rows:
        got = replay_merge(r)
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


def main(path):
    print('=' * 70)
    print(f'prep_atm dump: {path}')
    print('=' * 70)
    COLS, rows, arr, cps, cells = load(path)
    if len(rows) == 0:
        print('  SKIP: no data rows (dump empty; nothing to replay)')
        return 1

    mapping = {sx: sx for (_, _, _, sx) in MERGE_GROUPS.values()}
    worst = compare(rows, mapping, 'prep_atm')

    print('\n  consistency checks:')
    print(f'    lfrac + ifrac + ofrac == 1 (tol 1e-12): '
          f'{bool(np.all(np.abs(arr["lfrac"] + arr["ifrac"] + arr["ofrac"] - 1.0) < 1e-12))}')
    fracs = np.unique(np.round(arr['ifrac'], 4))
    print(f'    distinct ifrac values: {len(fracs)} (min {fracs.min():.4f}, max {fracs.max():.4f})')
    n_partial_ice = int(((arr['ifrac'] > 0) & (arr['ifrac'] < 1)).sum())
    n_partial_lnd = int(((arr['lfrac'] > 0) & (arr['lfrac'] < 1)).sum())
    n_open = int((arr['ifrac'] == 0).sum())
    n_ice = int((arr['ifrac'] == 1).sum())
    print(f'    partial-ice rows: {n_partial_ice}, open-ocean rows: {n_open}, '
          f'fully-iced rows: {n_ice}, partial-land rows: {n_partial_lnd}')
    print(f'    every step traces the base dump_first_cells window '
          f'(n = 1..12 present in every step; {min(len(c) for c in cps.values())} '
          f'cells/step): '
          f'{bool(all(set(range(1, 13)).issubset(set(c)) for c in cps.values()))}')
    # the So_t export is a copy of the So_t source (esmFldsExchange L1282)
    print(f'    So_t export == So_t source (copy) everywhere: '
          f'{bool(np.all(arr["So_t_exp"] == arr["So_t"]))}')
    n_partial_sel = int(((arr['ifrac'] > 0) & (arr['ifrac'] < 1)).sum()) + \
        int(((arr['lfrac'] > 0) & (arr['lfrac'] < 1)).sum())
    per_step_trig = all(
        ((a['ifrac'] > 0) & (a['ifrac'] < 1)).any() or
        ((a['lfrac'] > 0) & (a['lfrac'] < 1)).any()
        for s in cps
        for a in [{c: arr[c][arr['step'] == s] for c in arr}])
    print(f'    at least one partial-ice or partial-land row per step (selection trigger): '
          f'{bool(per_step_trig)} ({n_partial_sel} such rows total)')
    # lat/lon sanity on the atm grid (f19: -90..90, -180..180)
    print(f'    lat range {arr["lat"].min():.4f}..{arr["lat"].max():.4f}, '
          f'lon range {arr["lon"].min():.4f}..{arr["lon"].max():.4f}')
    # the same traced cells are traced every step (the per-step extras
    # are picked from the same PET-0 local set, so the (n,gid) pairing
    # is stable across steps)
    meta = {(int(n), int(g)) for n, g in zip(arr['n'], arr['gid'])}
    print(f'    distinct (n, gid) pairs across all steps: {len(meta)} '
          f'(cells {len(cells)}; a stable pairing means {len(meta) == len(cells)})')

    ok = all(w < 1e-11 for w in worst.values())
    print()
    print(f'RESULT: prep_atm {"PASS" if ok else "FAIL"}')
    return 0 if ok else 1


if __name__ == '__main__':
    base = '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cmeps/'
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else base + 'cmeps_prep_atm.csv'))
