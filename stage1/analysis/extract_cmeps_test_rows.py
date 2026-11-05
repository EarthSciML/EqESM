#!/usr/bin/env python3
"""Extract stage-2 test rows from the CMEPS aoflux dump.

Selects the hand-picked rows (see ROWS / ROWS_WIDE below) from
cmeps_aoflux.csv (fixed first-N selection) or cmeps_aoflux_wide.csv
(per-step stable/wind dynamic selection) and emits the `tests` array
JSON for stage2/cmeps_aoflux.esm: each test pins one dump row at
time_span {0,0} with parameter_overrides = the 12 raw inputs (zbot,
ubot, vbot, thbot, qbot, rainc, rbot, tbot, pslv, us, vs, ts) at full
dump precision and assertions on all 47 dumped columns plus the
convergence flags conv_1..4 and iter (computed by the validate_cmeps_dump
replay, which reproduces the dump to <= 4.4e-15 relative on every column).

Usage: python3 extract_cmeps_test_rows.py [dump.csv] > tests.json
       python3 extract_cmeps_test_rows.py --wide > tests.json
"""
import json
import sys

from validate_cmeps_dump import read_cols, replay

# hand-picked rows: (step, gid) with the law each exercises
ROWS = [
    (1, 41583, 'iter2_typical_unstable',
     'the typical traced state (step 1, gid 41583): convergence at iteration 2 '
     '(145/288 rows), vmag 2.72 m/s, hol -1.84, no clamp active'),
    (1, 42662, 'iter3_cold_cell',
     'the coldest traced cell (step 1, gid 42662, ~69 S): convergence at '
     'iteration 3, vmag 1.33 m/s, hol -4.12'),
    (3, 42662, 'iter4_cold_cell',
     'cell 42662 at step 3: convergence at iteration 4, vmag 0.80 m/s, '
     'hol -7.20 (approaching the |hol| <= 10 clamp)'),
    (8, 42662, 'iter5_hol_clamp',
     'cell 42662 at step 8: all 5 iterations run (no convergence by the '
     'flux_convergence = 0.01 test) with the stability parameter pinned at '
     'the |hol| <= 10 clamp (F90:254); vmag 0.52 m/s near the minwind floor'),
    (22, 42662, 'iter5_minwind_clamp',
     'cell 42662 at step 22: all 5 iterations run with the wind pinned at '
     'the minwind = 0.5 m/s floor (F90:193-196, seq_flux_atmocn_minwind) '
     'while hol = -9.76 stays just off the clamp'),
    (24, 41583, 'iter2_late_step',
     'cell 41583 at the last traced step 24: convergence at iteration 2 '
     'again, vmag 2.22 m/s, hol -2.58 (the same cell as iter2_typical_'
     'unstable, 23 coupling steps later)'),
]

# wide dump rows: hand-picked (step, gid) pairs from the dynamic
# (per-step stable/wind ranked) selection of cmeps_aoflux_wide.csv.
# The stable branch and the cdn 33 m/s cap are NOT in this window (no
# PET-0 cell has delt > 0 and the domain wind maximum is 4.8 m/s; the
# selection ranks all 12 masked PET-0 cells every step, so this is
# exhaustive, not a sampling artifact) -- the pins below are the
# iteration-count variants the wide rows produce (2/3/4) plus the
# strongest-wind row the run offers.
ROWS_WIDE = [
    (4, 42667, 'wide_maxwind_cell',
     'the strongest-wind row of the wide dump (step 4, gid 42667, n 7117, '
     '~69.0 S): the wind-ranked selection\'s domain maximum vmag 4.81 m/s, '
     'the strongest wind this run\'s PET-0 slice offers (the cdn 33 m/s cap '
     'of F90:153-154 is not reached anywhere in the run), converged at '
     'iteration 2 with hol -1.16'),
    (13, 42126, 'wide_iter3_coldest_cell',
     'wide dump row (step 13, gid 42126, n 7026, ~69.3 S): convergence at '
     'iteration 3 with delt -4.54 K, vmag 1.75 m/s, hol -4.21, at one of the '
     'domain\'s coldest cells (its delt reaches -5.86 K over the 24 steps)'),
    (19, 41583, 'wide_iter4_holclamp_adj',
     'wide dump row (step 19, gid 41583, n 6933, ~69.5 S): convergence at '
     'iteration 4 with hol -9.901, a tenth off the |hol| <= 10 clamp '
     '(F90:254), at the wind-ranked selection\'s weakest-wind step '
     '(vmag 0.64 m/s; the same cell as iter2_typical_unstable at step 1 of '
     'the fixed dump)'),
]


def main(path, rows):
    COLS = read_cols(path)
    lines = open(path).readlines()
    data = [l for l in lines if l.strip() and not l.startswith('#')]
    INUM = ('step', 'ymd', 'tod', 'n', 'gid', 'iter')
    recs = [dict(zip(COLS, l.strip().split(','))) for l in data]
    bykey = {(r['step'], r['gid'].strip()): r for r in recs}
    # 'lat' names both the latitude metadata (first occurrence) and the
    # latent-heat flux (second); the dict above keeps only the last, so the
    # true cell coordinates for the description come from the raw fields
    lat_idx = COLS.index('lat')
    lon_idx = COLS.index('lon')
    meta = {(r['step'], r['gid'].strip()): (r['n'].strip(),
            l.strip().split(',')[lat_idx], l.strip().split(',')[lon_idx])
            for r, l in zip(recs, data)}

    # assertion order: the model's equation order
    ASSERT_VARS = [
        # derived inputs
        'vmag', 'wind0', 'ugust', 'esat', 'qsat', 'ssq', 'delt', 'delq',
        'alz', 'cp',
        # neutral coefficients
        'stable0', 'rdn0', 'rhn0', 'ren0',
        # exit state of the stability iteration
        'hol', 'stable', 'xsq', 'xqq', 'psimh', 'psixh', 'rd', 'u10n',
        'rdn', 'rhn', 'ren', 'rh', 're', 'ustar', 'tstar', 'qstar',
        # convergence representation
        'conv_1', 'conv_2', 'conv_3', 'conv_4', 'iter',
        # fluxes
        'tau', 'taux', 'tauy', 'sen', 'lat', 'lwup', 'evap',
        # diagnostics
        'hol2', 'xsq2', 'xqq2', 'psix2', 'fac_t', 'fac_q', 'tref', 'qref',
        'duu10n', 'u10res',
    ]
    INPUTS = ['zbot', 'ubot', 'vbot', 'thbot', 'qbot', 'rainc', 'rbot',
              'tbot', 'pslv', 'us', 'vs', 'ts']

    tests = []
    for step, gid, tid, desc in rows:
        r = bykey[(str(step), str(gid))]
        got = replay({c: float(r[c]) for c in INPUTS})
        params = {c: float(r[c]) for c in INPUTS}
        assertions = []
        for v in ASSERT_VARS:
            exp = float(r[v]) if v in r else got[v]
            assertions.append({'variable': v, 'time': 0.0, 'expected': exp})
        tests.append({
            'id': tid,
            'description': (
                'Instantaneous pin of the full flux_atmOcn_large chain '
                '(flux_atmocn_Large.F90, Large & Yeager scheme as configured '
                'in the C_JRA_t62_inst run: add_gusts = .false., '
                'aofluxes_use_shr_wv_sat = .true., minwind = 0.5 m/s, '
                'flux_convergence = 0.01, flux_max_iteration = 5) from '
                'stage1/dumps/cmeps/%s row (step %s, gid %s, '
                'PET-0 local ogrid index n = %s, lat %s, lon %s): %s. '
                'time_span {0,0} evaluates the pointwise laws at the dumped '
                'inputs; every dumped column is asserted at full dump '
                'precision, and the convergence flags conv_k = |(u_k - '
                'u_{k-1})/u_k| <= 0.01 (k = 1..4) and the iteration count '
                'iter = first converged k, else 5, are asserted from the '
                'replay (which reproduces the dump rows to <= 4.4e-15 '
                'relative; see stage1/analysis/validate_cmeps_dump.py).'
            ) % (path.rsplit('/', 1)[-1], r['step'], r['gid'].strip(),
                 meta[(r['step'], r['gid'].strip())][0],
                 meta[(r['step'], r['gid'].strip())][1],
                 meta[(r['step'], r['gid'].strip())][2], desc),
            'time_span': {'start': 0.0, 'end': 0.0},
            'parameter_overrides': params,
            'tolerance': {'rel': 1e-12, 'abs': 0.0},
            'assertions': assertions,
        })
    json.dump(tests, sys.stdout, indent=1)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--wide':
        path = '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cmeps/cmeps_aoflux_wide.csv'
        rows = ROWS_WIDE
    else:
        path = sys.argv[1] if len(sys.argv) > 1 else \
            '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cmeps/cmeps_aoflux.csv'
        rows = ROWS
    main(path, rows)
