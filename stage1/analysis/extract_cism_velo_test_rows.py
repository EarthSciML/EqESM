#!/usr/bin/env python3
"""Extract numeric test tuples from the stage-1 CISM DIVA velocity dumps into
JSON fragments for the hand-authored `tests` blocks of stage2/cism_velo.esm.

The dumps (stage1/dumps/cism/cism_velo.csv, cism_velo_nbr.csv,
cism_velo_stag.csv, cism_basal_traction.csv, cism_velo_marine_exit.csv,
cism_velo_marine_pass.csv; gitignored, living in the main checkout) hold 4
steps x 32 traced cells x (entry, exit) rows, the per-iteration basal-traction
rows of the same 32 vertices, and the 4-step x 8-cell marine-margin boundary.
This script first RE-VERIFIES every law the stage-2 models encode against all
128 solve-exit rows (failing loudly on any mismatch), then emits, for each
hand-picked row, the `parameter_overrides` + `assertions` of one
`time_span {0,0}` test at full dump precision.

Laws verified here (all on the exit rows, aligned by (cism_step, cell)):

  DivaGeometry:
    z_mean = usrf - thck/2 at each of the 4 stencil cells
    dusrf_dx = clamp(0.5*(edge_s + edge_n), +/-max_slope), the HO_GRADIENT_MARGIN_HYBRID
      edge gradients (thklim = 1 m, thck_gradient_ramp = 50 m; land-terminating
      margins damped by the ramp factor; ocean/nunatak edges zeroed)
    dusrf_dy likewise from the west/east edges
    dz_dx, dz_dy = PLAIN centered vertex gradients of z_mean (no hybrid, no clamp:
      the glissade_slope_angle_staggered construction, slope_mask = 1 everywhere)
    theta_x = atan(|dz_dx|), F_x = cos^2(theta_x) + 4 sin^2(theta_x) (and y)
    stagthck = 4-cell mean of thck
    taudx = -rhoi*g*stagthck*dusrf_dx, taudy likewise
    effecpress_stag = rhoi*g*stagthck (bwat = 0 everywhere)

  DivaIntegrals:
    efvs = arithmetic mean of the 4 quad-point viscosities (the Fortran
      efvs(:,i,j) = sum_p efvs_qp(:,p)/nQuadPoints_2d)
    omega_k[kn] = mean_p sum_{k'>=kn} sig_mid(k')/efvs_qp(k',p) * dsigma(k') * thck
      (the Goldberg 2011 single integral, F90:5627-5653)
    omega = sum_k 0.5*(omega_k[k] + omega_k[k+1]) * dsigma(k)  (trapezoid, no thck)
    stag_omega = the ice-masked 4-cell mean of omega (cism_stagger_masked_mean,
      glissade_stagger with stagger_margin = 1, glissade_grid_operators.F90:
      135-148; ice_plus_land_mask of glissade_velo_higher.F90:1532-1537) -- on
      every traced row the 4 masks are 1, so the masked mean equals the plain
      4-cell mean 0.25*(omega + omega_se + omega_nw + omega_ne) exactly
    stag_omega_k[k] = the same masked 4-cell mean applied per velocity level to
      the omega_k of the 4 stencil cells (the per-level stagger of
      glissade_velo_higher.F90:2561-2564; the dump's omega_k_sw column is
      cross-checked against the traced cell's own omega_k)
    inactive cells (ice_mask = 0): efvs = omega_k = omega = 0 (the whole
      block lives inside the active_cell guard, F90:5091, and the dumped
      efvs_q of a cell that lost its ice mid-run hold the stale previous-step
      values while efvs/omega/omega_k read 0)

  DivaVelocity:
    beta_internal = f_ground * beta -- the pseudo-plastic law of the mounted
      BasalTraction, evaluated at the (u, v, bed, N) of the LAST ACCEPTED
      calcbeta pass: the dumped beta_internal is the accepted pass's
      f_ground*beta (rel 0.0 on all 128 rows; the pass is the single
      cism_basal_traction.csv row whose beta equals beta_internal, velo_iter
      = 38 of 39 on every row, and its bed/N equal the exit row's topg - eus /
      effecpress_stag exactly).  The pass velocity iterate differs from the
      exit basal velocity by the final Picard correction, so feeding the exit
      velocity instead would leave a ~6e-7 residual; feeding the PASS velocity
      makes the whole chain exact (rel 0.0 on the law, <= 2.9e-16 downstream)
      -- the reason the stage-2 tests feed the pass row and assert at 1e-12.
    beta_eff_x = F_x*beta_internal/(F_x + beta_internal*stag_omega*cos(theta_x))
      (and y) -- EXACT given the dumped beta_internal
    uvel_basal = uvel_2d*beta_eff_x/beta_internal (Goldberg eq. 34 recovery)
    btractx = beta_internal*uvel_basal = beta_eff_x*uvel_2d (F90:3874, 5780)
    uvel[k] = uvel_basal + btractx*stag_omega_k[k]*cos(theta_x)/F_x, k = 1..10,
      uvel(11) = uvel_basal (stag_omega_k(11) = 0) -- and the same for vvel
      with the y-direction quantities (Goldberg eq. 32 + slope correction)

  Marine-margin boundary (cism_velo_marine_exit.csv / _pass.csv; the first 8
    cells with initial ice over a marine bed, topg < eus): the same coupled
    chain runs on the marine rows, with the same matching-pass beta link.  The
    reference run has NO partially grounded vertices (f_ground = 1 on all 32
    exit rows and all 32 marine exit rows), so the f_ground < 1 branch of
    beta_internal = f_ground*beta is not exercised by any dump row; the marine
    tests pin the fully grounded marine-bed configuration instead.

Usage:
  python3 stage1/analysis/extract_cism_velo_test_rows.py [dumps_dir] > tests_fragments.json
where dumps_dir defaults to <main checkout>/stage1/dumps/cism (dumps are
gitignored and live in the main checkout, not in the worktree; the
environment variable EQESM_MAIN overrides the main-checkout location).
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dumpcsv import load

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MAIN = os.environ.get("EQESM_MAIN", ROOT)
SIBLING_MAIN = os.path.normpath(os.path.join(ROOT, "..", "EqESM"))  # the main checkout

RHOI, G, DX = 917.0, 9.80616, 20000.0
THKLIM, RAMP, MAXSLOPE = 1.0, 50.0, 0.5
NL, NZ = 10, 11          # nlayers, velocity levels

# Hand-picked test rows, keyed by (cism_step, i, j) per model.
GEOM_ROWS = [(1, 25, 11), (1, 27, 13), (3, 24, 14), (3, 28, 13)]
INTEGRAL_ROWS = [(1, 25, 11), (1, 26, 12), (1, 27, 13), (3, 24, 14)]
# ice_mask = 0 with NONZERO dumped efvs_q: (3,17,14) lost its ice between
# steps 2 and 3 (thck 97.3 -> 0.986 m < thklim), so the dumped efvs_q columns
# hold the stale previous-step values while the dumped efvs/omega/omega_k are
# 0 (the active_cell block, F90:5091, is skipped) - the guard must zero the
# integrals and the cell-mean viscosity against nonzero inputs.
INTEGRAL_INACTIVE_ROWS = [(3, 17, 14)]
VELO_ROWS = [(1, 25, 11), (1, 22, 13), (1, 27, 13), (3, 24, 14)]
# Marine-margin cells (first dump_marine_cells = 8 cells with initial ice over
# a marine bed, topg < eus): a thin-bed cell and a deep-bed cell.
MARINE_ROWS = [(1, 27, 49), (3, 39, 51)]


def hybrid_edge_gradient(usrf_w, usrf_e, thck_w, thck_e, ice_w, ice_e, land_w, land_e):
    """ds/dx on one east edge, HO_GRADIENT_MARGIN_HYBRID (grid_operators F90:999-1079)."""
    if usrf_e >= usrf_w:
        iu, il, sign = 'e', 'w', 1.0
    else:
        iu, il, sign = 'w', 'e', -1.0
    u = dict(e=usrf_e, w=usrf_w)[iu]
    l = dict(e=usrf_e, w=usrf_w)[il]
    ice_u = dict(e=ice_e, w=ice_w)[iu]
    ice_l = dict(e=ice_e, w=ice_w)[il]
    land_u = dict(e=land_e, w=land_w)[iu]
    land_l = dict(e=land_e, w=land_w)[il]
    thck_u = dict(e=thck_e, w=thck_w)[iu]
    if land_u == 1 and RAMP > 0.0:
        edge_factor = max(0.0, min(1.0, (thck_u - THKLIM) / RAMP))
    else:
        edge_factor = 1.0
    if ice_u == 1 and (ice_l == 1 or land_l == 1):
        return edge_factor * sign * (u - l) / DX
    return 0.0


def clamp_slope(x):
    return max(-MAXSLOPE, min(MAXSLOPE, x))


def basal_beta(u, v, bed, N_stag, q=0.5, u0=100.0, eps=0.01,
               phimin=5.0, phimax=40.0, bedmin=-300.0, bedmax=700.0):
    """The PSEUDO_PLASTIC_OLD law of glissade_calcbeta (F90:250-292)."""
    if bed <= bedmin:
        phi = phimin
    elif bed >= bedmax:
        phi = phimax
    else:
        phi = phimin + (bed - bedmin) / (bedmax - bedmin) * (phimax - phimin)
    tanphi = math.tan(phi * math.pi / 180.0)
    tau_c = tanphi * N_stag
    speed = math.sqrt(u * u + v * v + eps * eps)
    return tau_c / (u0 ** q * speed ** (1.0 - q))


def rel(a, b):
    return abs(a - b) / abs(b) if b != 0.0 else (0.0 if a == 0.0 else float('inf'))


def masked_mean(vals, masks):
    """The cism_stagger_masked_mean law (glissade_grid_operators.F90:135-148):
    sum(mask*val)/max(sum(mask), 1), with the vertex stencil cells ordered
    (sw, se, nw, ne)."""
    svar = sum(m * v for m, v in zip(masks, vals))
    smask = sum(masks)
    return svar / max(smask, 1.0)


def main():
    ap = argparse.ArgumentParser()
    default_dir = os.path.join(MAIN, "stage1", "dumps", "cism")
    if not os.path.isdir(default_dir) and os.path.isdir(os.path.join(SIBLING_MAIN, "stage1", "dumps", "cism")):
        default_dir = os.path.join(SIBLING_MAIN, "stage1", "dumps", "cism")
    ap.add_argument("dumps_dir", nargs="?", default=default_dir)
    a = ap.parse_args()

    _, velo = load(os.path.join(a.dumps_dir, "cism_velo.csv"))
    _, nbr = load(os.path.join(a.dumps_dir, "cism_velo_nbr.csv"))
    _, stag = load(os.path.join(a.dumps_dir, "cism_velo_stag.csv"))
    _, bt = load(os.path.join(a.dumps_dir, "cism_basal_traction.csv"))
    ve = {f"{r['cism_step']}_{r['i']}_{r['j']}": r for r in velo if r['phase'] == 'exit'}
    ne = {f"{r['cism_step']}_{r['i']}_{r['j']}": r for r in nbr if r['phase'] == 'exit'}
    se = {f"{r['cism_step']}_{r['i']}_{r['j']}": r for r in stag if r['phase'] == 'exit'}
    common = sorted(set(ve) & set(ne) & set(se))
    assert len(common) == 128, len(common)
    print(f"# exit rows: velo {len(ve)}, nbr {len(ne)}, stag {len(se)}, common {len(common)}",
          file=sys.stderr)

    # per-vertex basal-traction pass rows (all rows are passes; velo_iter 1 is
    # the effecpress trace and 2..39 the calcbeta iterations)
    bt_by_cell = {}
    for r in bt:
        bt_by_cell.setdefault((r['cism_step'], r['i'], r['j']), []).append(r)

    stats = {}
    worst = {}
    fails = 0

    def check(name, computed, dumped, rtol=1e-9):
        nonlocal fails
        r = rel(computed, dumped)
        stats.setdefault(name, [0, 0.0])
        stats[name][0] += 1
        stats[name][1] = max(stats[name][1], r)
        if r > rtol:
            fails += 1
            if fails <= 6:
                print(f"  FAIL {name}: computed={computed!r} dumped={dumped!r} rel={r:.2e}",
                      file=sys.stderr)

    chain_worst = 0.0
    pass_iter = {}
    for key in common:
        v, n, s = ve[key], ne[key], se[key]
        step, i, j = key.split('_')
        sigma = [n[f'sigma_{k:03d}'] for k in range(1, NZ + 1)]
        dsigma = [sigma[k + 1] - sigma[k] for k in range(NL)]
        sig_mid = [0.5 * (sigma[k] + sigma[k + 1]) for k in range(NL)]
        efvs_q = [[n[f'efvs_q{p}_{k:03d}'] for k in range(1, NL + 1)] for p in range(1, 5)]
        omega_k = [n[f'omega_k_{k:03d}'] for k in range(1, NZ + 1)]
        thck, ice_mask = v['thck'], v['ice_mask']
        active = ice_mask == 1.0
        fg = v['f_ground']

        # --- DivaIntegrals -------------------------------------------------
        # efvs = arithmetic mean of the quad-point viscosities (nbr-dump column);
        # 0 on inactive cells even when the dumped efvs_q are stale nonzero values
        # (the mean is computed inside the active_cell block, F90:5266-5272)
        efvs_c = [0.25 * sum(efvs_q[p][k] for p in range(4)) for k in range(NL)]
        check('efvs_mean', efvs_c[0] if active else 0.0, n['efvs_001'], 1e-12)
        check('efvs_mean', efvs_c[5] if active else 0.0, n['efvs_006'], 1e-12)
        # omega_k single integral
        omega_k_c = []
        for kk in range(NZ):
            vals = []
            for p in range(4):
                ssum = 0.0
                for k in range(kk, NL):
                    if active:
                        ssum += sig_mid[k] / efvs_q[p][k] * dsigma[k] * thck
                vals.append(ssum)
            omega_k_c.append(0.25 * sum(vals))
        if active:
            for kk in (0, 5, 10):
                check('omega_k_int', omega_k_c[kk], omega_k[kk], 1e-10)
        else:
            assert all(x == 0.0 for x in omega_k), (key, 'inactive omega_k nonzero')
        # omega trapezoid
        omega_c = sum(0.5 * (omega_k[k] + omega_k[k + 1]) * dsigma[k] for k in range(NL))
        check('omega_trap', omega_c, n['omega'], 1e-12)
        # the stag dump's SW-cell omega_k must be the traced cell's own omega_k
        ok_sw = [s[f'omega_k_sw_{k:03d}'] for k in range(1, NZ + 1)]
        for kk in (1, 6, 11):
            check('stag_omega_k_sw', ok_sw[kk - 1], n[f'omega_k_{kk:03d}'], 1e-12)
        # masked 4-cell staggers (cism_stagger_masked_mean)
        masks = [s[f'ice_plus_land_mask_{t}'] for t in ('sw', 'se', 'nw', 'ne')]
        ok_vals = {t: [s[f'omega_k_{t}_{k:03d}'] for k in range(1, NZ + 1)]
                   for t in ('sw', 'se', 'nw', 'ne')}
        check('stag_omega_masked',
              masked_mean([n['omega'], n['omega_se'], n['omega_nw'], n['omega_ne']], masks),
              n['stag_omega'], 1e-12)
        for kk in (1, 6, 11):
            check('stag_omega_k_masked',
                  masked_mean([ok_vals[t][kk - 1] for t in ('sw', 'se', 'nw', 'ne')], masks),
                  n[f'stag_omega_k_{kk:03d}'], 1e-12)

        # --- DivaGeometry --------------------------------------------------
        c = {}
        for tag in ('sw', 'se', 'nw', 'ne'):
            c[tag] = dict(usrf=n[f'usrf_{tag}'], thck=n[f'thck_{tag}'],
                          ice=n[f'ice_mask_{tag}'], land=n[f'land_mask_{tag}'])
        zmean = {t: c[t]['usrf'] - 0.5 * c[t]['thck'] for t in c}
        edge_s = hybrid_edge_gradient(c['sw']['usrf'], c['se']['usrf'], c['sw']['thck'],
                                      c['se']['thck'], c['sw']['ice'], c['se']['ice'],
                                      c['sw']['land'], c['se']['land'])
        edge_n = hybrid_edge_gradient(c['nw']['usrf'], c['ne']['usrf'], c['nw']['thck'],
                                      c['ne']['thck'], c['nw']['ice'], c['ne']['ice'],
                                      c['nw']['land'], c['ne']['land'])
        edge_w = hybrid_edge_gradient(c['sw']['usrf'], c['nw']['usrf'], c['sw']['thck'],
                                      c['nw']['thck'], c['sw']['ice'], c['nw']['ice'],
                                      c['sw']['land'], c['nw']['land'])
        edge_e = hybrid_edge_gradient(c['se']['usrf'], c['ne']['usrf'], c['se']['thck'],
                                      c['ne']['thck'], c['se']['ice'], c['ne']['ice'],
                                      c['se']['land'], c['ne']['land'])
        ddx_c = clamp_slope(0.5 * (edge_s + edge_n))
        ddy_c = clamp_slope(0.5 * (edge_w + edge_e))
        check('dusrf_dx', ddx_c, n['dusrf_dx'], 1e-9)
        check('dusrf_dy', ddy_c, n['dusrf_dy'], 1e-9)
        stagthck_c = 0.25 * sum(c[t]['thck'] for t in c)
        check('stagthck', stagthck_c, v['stagthck'], 1e-12)
        check('taudx', -RHOI * G * stagthck_c * ddx_c, v['taudx'], 1e-9)
        check('taudy', -RHOI * G * stagthck_c * ddy_c, v['taudy'], 1e-9)
        check('effecpress_stag', RHOI * G * stagthck_c, v['effecpress_stag'], 1e-12)
        dz_dx_c = 0.5 * ((zmean['ne'] - zmean['nw']) / DX + (zmean['se'] - zmean['sw']) / DX)
        dz_dy_c = 0.5 * ((zmean['ne'] - zmean['se']) / DX + (zmean['nw'] - zmean['sw']) / DX)
        check('theta_x', math.atan(abs(dz_dx_c)), n['stag_theta_slope_x'], 1e-10)
        check('theta_y', math.atan(abs(dz_dy_c)), n['stag_theta_slope_y'], 1e-10)
        for tag, th in (('x', n['stag_theta_slope_x']), ('y', n['stag_theta_slope_y'])):
            F_c = math.cos(th) ** 2 + 4.0 * math.sin(th) ** 2
            check(f'F_{tag}', F_c, n[f'stag_diva_slope_factor_{tag}'], 1e-12)

        # --- DivaVelocity --------------------------------------------------
        beta_int = n['beta_internal']
        for tag in ('x', 'y'):
            th = n[f'stag_theta_slope_{tag}']
            Fv = n[f'stag_diva_slope_factor_{tag}']
            be_c = Fv * beta_int / (Fv + beta_int * n['stag_omega'] * math.cos(th))
            check(f'beta_eff_{tag}', be_c, n[f'beta_eff_{tag}'], 1e-12)
        u2d, v2d = v['uvel_2d'], v['vvel_2d']
        u_b, v_b = v['uvel_011'], v['vvel_011']
        check('uvel_basal_rec', u2d * n['beta_eff_x'] / beta_int, u_b, 1e-12)
        check('vvel_basal_rec', v2d * n['beta_eff_y'] / beta_int, v_b, 1e-12)
        check('btractx1', beta_int * u_b, v['btractx'], 1e-12)
        check('btracty1', beta_int * v_b, v['btracty'], 1e-12)
        check('btractx2', n['beta_eff_x'] * u2d, v['btractx'], 1e-12)
        check('btracty2', n['beta_eff_y'] * v2d, v['btracty'], 1e-12)
        for tag, bb, uu in (('x', u_b, 'uvel'), ('y', v_b, 'vvel')):
            th = n[f'stag_theta_slope_{tag}']
            Fv = n[f'stag_diva_slope_factor_{tag}']
            corr = math.cos(th) / Fv
            bt_ = v[f'btract{"x" if tag == "x" else "y"}']
            for kk in range(NZ - 1):
                check(f'{uu}_profile', bb + bt_ * n[f'stag_omega_k_{kk + 1:03d}'] * corr,
                      v[f'{uu}_{kk + 1:03d}'], 1e-12)

        # beta chain at the ACCEPTED calcbeta pass (the exact link): the pass
        # row whose beta equals the dumped beta_internal (f_ground = 1)
        rows = bt_by_cell[(v['cism_step'], v['i'], v['j'])]
        hits = [r for r in rows if rel(r['beta'], beta_int) < 1e-14]
        assert len(hits) == 1, (key, len(hits))
        p = hits[0]
        pass_iter[p['velo_iter']] = pass_iter.get(p['velo_iter'], 0) + 1
        beta_law = fg * basal_beta(p['u'], p['v'], p['bed'], p['N'])
        chain_worst = max(chain_worst, rel(beta_law, beta_int))
        check('beta_at_pass', beta_law, beta_int, 1e-12)
        check('pass_bed', p['bed'], v['topg'] - v['eus'], 1e-12)
        check('pass_N', p['N'], v['effecpress_stag'], 1e-12)

    print(f"# accepted-pass velo_iter distribution: {dict(pass_iter)}", file=sys.stderr)
    print(f"# chain worst rel (law at the accepted pass vs dumped beta_internal): "
          f"{chain_worst:.3e}", file=sys.stderr)
    for name in sorted(stats):
        print(f"# {name:16s} {stats[name][0]:4d} checks, worst rel {stats[name][1]:.3e}",
              file=sys.stderr)
    assert fails == 0, f"{fails} law mismatches"
    print("# ALL LAWS VERIFIED on 128 exit rows", file=sys.stderr)

    # ---- marine-margin boundary -------------------------------------------
    me = {}
    mp = {}
    try:
        _, mex = load(os.path.join(a.dumps_dir, "cism_velo_marine_exit.csv"))
        _, mpas = load(os.path.join(a.dumps_dir, "cism_velo_marine_pass.csv"))
    except FileNotFoundError:
        mex, mpas = [], []
    for r in mex:
        me[f"{r['cism_step']}_{r['i']}_{r['j']}"] = r
    for r in mpas:
        mp.setdefault((r['cism_step'], r['i'], r['j']), []).append(r)
    if me:
        print(f"# marine exit rows {len(me)}, pass rows {len(mpas)}", file=sys.stderr)
        mfg = sorted(set(r['f_ground'] for r in me.values()))
        print(f"# marine exit f_ground values: {mfg}", file=sys.stderr)
        for key, m in me.items():
            rows = mp[(m['cism_step'], m['i'], m['j'])]
            hits = [r for r in rows if rel(r['f_ground'] * r['beta'], m['beta_internal']) < 1e-14]
            assert len(hits) == 1, (key, len(hits))
            p = hits[0]
            beta_law = m['f_ground'] * basal_beta(p['u'], p['v'], p['bed'], p['N'])
            assert rel(beta_law, m['beta_internal']) < 1e-12, (key, rel(beta_law, m['beta_internal']))
    else:
        print("# marine dumps absent; marine fragments skipped", file=sys.stderr)

    # ---- emit test fragments ----------------------------------------------
    out = {"geom": [], "integrals": [], "integrals_inactive": [], "velo": [], "marine": []}

    def sigma_of(n):
        return [n[f'sigma_{k:03d}'] for k in range(1, NZ + 1)]

    def efvs_of(n):
        return {p: [n[f'efvs_q{p}_{k:03d}'] for k in range(1, NL + 1)] for p in range(1, 5)}

    def integ_overrides(v, n, s):
        ov = {"sigma": sigma_of(n)}
        for p in range(1, 5):
            ov[f"efvs_q{p}"] = efvs_of(n)[p]
        ov.update({
            "thck": v["thck"], "ice_mask": v["ice_mask"],
            "omega_se": n["omega_se"], "omega_nw": n["omega_nw"], "omega_ne": n["omega_ne"],
            "ipl_mask_sw": s["ice_plus_land_mask_sw"], "ipl_mask_se": s["ice_plus_land_mask_se"],
            "ipl_mask_nw": s["ice_plus_land_mask_nw"], "ipl_mask_ne": s["ice_plus_land_mask_ne"],
            "omega_k_se": [s[f'omega_k_se_{k:03d}'] for k in range(1, NZ + 1)],
            "omega_k_nw": [s[f'omega_k_nw_{k:03d}'] for k in range(1, NZ + 1)],
            "omega_k_ne": [s[f'omega_k_ne_{k:03d}'] for k in range(1, NZ + 1)],
        })
        return ov

    def geom_overrides(n):
        ov = {}
        for tag in ('sw', 'se', 'nw', 'ne'):
            for q in ('usrf', 'thck', 'ice_mask', 'land_mask'):
                ov[f'{q}_{tag}'] = n[f'{q}_{tag}']
        return ov

    # DivaGeometry fragments
    for step, i, j in GEOM_ROWS:
        key = f"{step}_{i}_{j}"
        v, n = ve[key], ne[key]
        src = {"dusrf_dx": n["dusrf_dx"], "dusrf_dy": n["dusrf_dy"], "stagthck": v["stagthck"],
               "taudx": v["taudx"], "taudy": v["taudy"],
               "theta_x": n["stag_theta_slope_x"], "theta_y": n["stag_theta_slope_y"],
               "F_x": n["stag_diva_slope_factor_x"], "F_y": n["stag_diva_slope_factor_y"]}
        # The hybrid edge-gradient values (dusrf_dx/dy and the taudx/taudy that
        # scale them) are differences of ~1800 m surface elevations over 20 km,
        # so the es26.17 dump round-trip puts them at up to 6.6e-13 relative
        # (the worst over all 128 exit rows); theta/F/stagthck stay at 1e-12.
        tol = {"dusrf_dx": 1e-11, "dusrf_dy": 1e-11, "taudx": 1e-11, "taudy": 1e-11}
        out["geom"].append({
            "id": f"geom_step{step}_cell{i}{j}",
            "parameter_overrides": geom_overrides(n),
            "assertions": [{"variable": nm, "time": 0.0, "expected": src[nm],
                            "tolerance": {"rel": tol.get(nm, 1e-12)}} for nm in src],
        })

    # DivaIntegrals fragments
    for step, i, j in INTEGRAL_ROWS:
        key = f"{step}_{i}_{j}"
        v, n, s = ve[key], ne[key], se[key]
        assertions = [
            {"variable": "efvs", "time": 0.0, "expected": n[f'efvs_{k:03d}'],
             "coords": {"lev": k}, "tolerance": {"rel": 1e-12}}
            for k in (1, 6)
        ]
        assertions += [
            {"variable": "omega_k", "time": 0.0, "expected": n[f'omega_k_{k:03d}'],
             "coords": {"lev_nodes": k}, "tolerance": {"rel": 1e-12}}
            for k in (1, 6, 11)
        ]
        assertions += [
            {"variable": "omega", "time": 0.0, "expected": n["omega"],
             "tolerance": {"rel": 1e-12}},
            {"variable": "stag_omega", "time": 0.0, "expected": n["stag_omega"],
             "tolerance": {"rel": 1e-12}},
        ]
        assertions += [
            {"variable": "stag_omega_k", "time": 0.0, "expected": n[f'stag_omega_k_{k:03d}'],
             "coords": {"lev_nodes": k}, "tolerance": {"rel": 1e-12}}
            for k in (1, 6, 11)
        ]
        out["integrals"].append({
            "id": f"int_step{step}_cell{i}{j}",
            "parameter_overrides": integ_overrides(v, n, s),
            "assertions": assertions,
        })
    for step, i, j in INTEGRAL_INACTIVE_ROWS:
        key = f"{step}_{i}_{j}"
        v, n, s = ve[key], ne[key], se[key]
        masks = [s[f'ice_plus_land_mask_{t}'] for t in ('sw', 'se', 'nw', 'ne')]
        ok_vals = {t: [s[f'omega_k_{t}_{k:03d}'] for k in range(1, NZ + 1)]
                   for t in ('sw', 'se', 'nw', 'ne')}
        out["integrals_inactive"].append({
            "id": f"int_step{step}_cell{i}{j}_inactive",
            "parameter_overrides": integ_overrides(v, n, s),
            "assertions": [
                {"variable": "efvs", "time": 0.0, "expected": 0.0,
                 "coords": {"lev": 1}, "tolerance": {"rel": 1e-12}},
                {"variable": "omega_k", "time": 0.0, "expected": 0.0,
                 "coords": {"lev_nodes": 1}, "tolerance": {"rel": 1e-12}},
                {"variable": "omega", "time": 0.0, "expected": 0.0,
                 "tolerance": {"rel": 1e-12}},
                {"variable": "stag_omega", "time": 0.0, "expected": n["stag_omega"],
                 "tolerance": {"rel": 1e-12}},
                {"variable": "stag_omega_k", "time": 0.0,
                 "expected": masked_mean([ok_vals[t][0] for t in ('sw', 'se', 'nw', 'ne')], masks),
                 "coords": {"lev_nodes": 1}, "tolerance": {"rel": 1e-12}},
            ],
        })

    # DivaVelocity fragments
    for step, i, j in VELO_ROWS:
        key = f"{step}_{i}_{j}"
        v, n, s = ve[key], ne[key], se[key]
        p = [r for r in bt_by_cell[(v['cism_step'], v['i'], v['j'])]
             if rel(r['beta'], n['beta_internal']) < 1e-14][0]
        ov = {"BasalTraction.u": p["u"], "BasalTraction.v": p["v"],
              "BasalTraction.bed": p["bed"], "BasalTraction.N_stag": p["N"],
              "uvel_2d": v["uvel_2d"], "vvel_2d": v["vvel_2d"],
              "uvel_basal": v["uvel_011"], "vvel_basal": v["vvel_011"],
              "f_ground": v["f_ground"]}
        ov.update({f"DivaGeometry.{k}": val for k, val in geom_overrides(n).items()})
        ov.update({f"DivaIntegrals.{k}": val for k, val in integ_overrides(v, n, s).items()})
        src = {"beta_internal": n["beta_internal"], "beta_eff_x": n["beta_eff_x"],
               "beta_eff_y": n["beta_eff_y"], "btractx": v["btractx"], "btracty": v["btracty"],
               "uvel_basal_recovered": v["uvel_011"], "vvel_basal_recovered": v["vvel_011"]}
        assertions = [{"variable": nm, "time": 0.0, "expected": src[nm],
                       "tolerance": {"rel": 1e-12}} for nm in src]
        prof_levels = range(1, NZ + 1) if (step, i, j) == (1, 22, 13) else (1, 6, 11)
        for kk in prof_levels:
            assertions.append({"variable": "uvel", "time": 0.0, "expected": v[f'uvel_{kk:03d}'],
                               "coords": {"lev_nodes": kk}, "tolerance": {"rel": 1e-12}})
            assertions.append({"variable": "vvel", "time": 0.0, "expected": v[f'vvel_{kk:03d}'],
                               "coords": {"lev_nodes": kk}, "tolerance": {"rel": 1e-12}})
        out["velo"].append({"id": f"velo_step{step}_cell{i}{j}",
                            "parameter_overrides": ov, "assertions": assertions})

    # marine fragments (same coupled-chain shape as the velo fragments)
    for step, i, j in MARINE_ROWS:
        key = f"{step}_{i}_{j}"
        if key not in me:
            print(f"# marine row {key} absent; skipped", file=sys.stderr)
            continue
        m = me[key]
        p = [r for r in mp[(m['cism_step'], m['i'], m['j'])]
             if rel(r['f_ground'] * r['beta'], m['beta_internal']) < 1e-14][0]
        n = m  # the marine exit row carries the nbr + stag columns
        s = m
        v = m
        ov = {"BasalTraction.u": p["u"], "BasalTraction.v": p["v"],
              "BasalTraction.bed": p["bed"], "BasalTraction.N_stag": p["N"],
              "uvel_2d": v["uvel_2d"], "vvel_2d": v["vvel_2d"],
              "uvel_basal": v["uvel_011"], "vvel_basal": v["vvel_011"],
              "f_ground": v["f_ground"]}
        ov.update({f"DivaGeometry.{k}": val for k, val in geom_overrides(n).items()})
        ov.update({f"DivaIntegrals.{k}": val for k, val in integ_overrides(v, n, s).items()})
        src = {"beta_internal": n["beta_internal"], "beta_eff_x": n["beta_eff_x"],
               "beta_eff_y": n["beta_eff_y"], "btractx": v["btractx"], "btracty": v["btracty"],
               "uvel_basal_recovered": v["uvel_011"], "vvel_basal_recovered": v["vvel_011"]}
        assertions = [{"variable": nm, "time": 0.0, "expected": src[nm],
                       "tolerance": {"rel": 1e-12}} for nm in src]
        for kk in (1, 6, 11):
            assertions.append({"variable": "uvel", "time": 0.0, "expected": v[f'uvel_{kk:03d}'],
                               "coords": {"lev_nodes": kk}, "tolerance": {"rel": 1e-12}})
            assertions.append({"variable": "vvel", "time": 0.0, "expected": v[f'vvel_{kk:03d}'],
                               "coords": {"lev_nodes": kk}, "tolerance": {"rel": 1e-12}})
        out["marine"].append({"id": f"marine_step{step}_cell{i}{j}",
                              "parameter_overrides": ov, "assertions": assertions})

    json.dump(out, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
