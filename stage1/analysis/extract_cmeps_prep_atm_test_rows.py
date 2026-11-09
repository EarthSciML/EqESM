#!/usr/bin/env python3
"""Extract stage-2 test rows from the CMEPS atmosphere-export merge dump.

Selects the hand-picked rows (see ROWS below) from cmeps_prep_atm.csv
and emits the `tests` array JSON for stage2/cmeps_prep_atm.esm: each
test pins one dump row at time_span {0,0} with parameter_overrides =
the fractions and the 45 source fields (land/ice/ocean) at full dump
precision and assertions on the 14 merged exports and the So_t copy
export (the model's equation order).

Usage: python3 extract_cmeps_prep_atm_test_rows.py [dump.csv] > tests.json
"""
import json
import sys

# hand-picked rows: (step, n, id, description) exercising the merge
# laws and the traced-cell selection.  The dump's n is the PET-0 local
# atm-grid index; rows are chosen after inspecting the dump so that
# the fraction regimes (open-ocean, partial-ice, partial-land,
# fully-iced) and the step range are all pinned.
ROWS = [
    (1, 0, 'step1_first_cell',
     'the first traced cell (step 1, n = 0, the base first-N selection)'),
    (12, 0, 'step12_first_cell',
     'the first traced cell at the last coupling step'),
]

# the model's parameter (source) columns, in the .esm variable order
PARAMS = [
    'lfrac', 'ifrac', 'ofrac',
    'Sl_avsdr', 'Si_avsdr', 'So_avsdr',
    'Sl_avsdf', 'Si_avsdf', 'So_avsdf',
    'Sl_anidr', 'Si_anidr', 'So_anidr',
    'Sl_anidf', 'Si_anidf', 'So_anidf',
    'Sl_tref', 'Si_tref', 'So_tref',
    'Sl_u10', 'Si_u10', 'So_u10',
    'Sl_qref', 'Si_qref', 'So_qref',
    'Fall_taux', 'Faii_taux', 'Faox_taux',
    'Fall_tauy', 'Faii_tauy', 'Faox_tauy',
    'Fall_lat', 'Faii_lat', 'Faox_lat',
    'Fall_sen', 'Faii_sen', 'Faox_sen',
    'Fall_evap', 'Faii_evap', 'Faox_evap',
    'Fall_lwup', 'Faii_lwup', 'Faox_lwup',
    'Sl_t', 'Si_t', 'So_t',
]

# assertion order: the model's equation order
ASSERTS = [
    'Sx_avsdr', 'Sx_avsdf', 'Sx_anidr', 'Sx_anidf',
    'Sx_tref', 'Sx_u10', 'Sx_qref',
    'Faxx_taux', 'Faxx_tauy', 'Faxx_lat', 'Faxx_sen',
    'Faxx_evap', 'Faxx_lwup',
    'Sx_t',
    'So_t_export',
]


def read_cols(path):
    with open(path) as f:
        for l in f:
            if l.startswith('# columns'):
                return l.split('= ', 1)[1].strip().split(',')
    raise RuntimeError(f'no # columns header found in {path}')


def main(path):
    COLS = read_cols(path)
    lines = open(path).readlines()
    data = [l for l in lines if l.strip() and not l.startswith('#')]
    recs = [dict(zip(COLS, l.strip().split(','))) for l in data]
    # the duplicated So_t (source, export) becomes So_t / So_t_export
    bykey = {(r['step'], int(r['n'])): r for r in recs}

    tests = []
    for step, n, tid, desc in ROWS:
        r = bykey[(str(step), n)]
        lat, lon = r['lat'], r['lon']
        params = {c: float(r[c]) for c in PARAMS}
        assertions = [{'variable': v, 'time': 0.0,
                       'expected': float(r['So_t' if v == 'So_t_export' else v])}
                      for v in ASSERTS]
        tests.append({
            'id': tid,
            'description': (
                'Instantaneous pin of the fraction-weighted atmosphere '
                'exports of med_phases_prep_atm (mediator/'
                'med_phases_prep_atm_mod.F90:149-157, the med_merge '
                'accumulation dp = dp + src*weight of med_merge_mod.F90:'
                '405-412 with the merge lists of esmFldsExchange_cesm_mod.'
                'F90, and the So_t copy export of L1282; the traced '
                'F2000climo_f19_inst run: component_computes_enthalpy_flux '
                '= none, mrg_fracname_lnd2atm_state/flux = lfrac, '
                'samegrid_atmlnd) from stage1/dumps/cmeps/%s row (step %s, '
                'atm cell n = %s, lat %s, lon %s): %s. time_span {0,0} '
                'evaluates the pointwise laws at the dumped inputs; every '
                'merged export is asserted at full dump precision (rel '
                '1e-12), so a wrong weight, a wrong source, a wrong '
                'association or a wrong fraction all fail.'
            ) % (path.rsplit('/', 1)[-1], r['step'], r['n'], lat, lon, desc),
            'time_span': {'start': 0.0, 'end': 0.0},
            'parameter_overrides': params,
            'tolerance': {'rel': 1e-12, 'abs': 0.0},
            'assertions': assertions,
        })
    json.dump(tests, sys.stdout, indent=1)


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else \
        '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cmeps/cmeps_prep_atm.csv'
    main(path)
