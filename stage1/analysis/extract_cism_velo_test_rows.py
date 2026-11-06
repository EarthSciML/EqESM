#!/usr/bin/env python3
"""Extract numeric test tuples from the stage-1 CISM DIVA velocity dumps into
JSON fragments for the hand-authored `tests` blocks of stage2/cism_velo.esm.

The dumps (stage1/dumps/cism/cism_velo.csv, cism_velo_nbr.csv; gitignored,
living in the main checkout) hold 4 steps x 32 traced cells x (entry, exit)
rows.  This script first RE-VERIFIES every law the stage-2 models encode
against all 128 solve-exit rows (failing loudly on any mismatch), then emits,
for each hand-picked row, the `parameter_overrides` + `assertions` of one
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
    stag_omega = 4-cell mean of omega (stagger_margin = 1)
    inactive cells (ice_mask = 0): efvs = omega_k = omega = 0 (the whole
      block lives inside the active_cell guard, F90:5091, and the dumped
      efvs_q of a cell that lost its ice mid-run hold the stale previous-step
      values while efvs/omega/omega_k read 0)

  DivaVelocity:
    beta_internal = f_ground * beta(u_basal, v_basal) -- the pseudo-plastic law
      of the mounted BasalTraction, evaluated at the EXIT basal velocity.  The
      dumped beta_internal is the last ACCEPTED calcbeta pass of the Picard
      loop (its velocity iterate differs from the exit basal velocity by the
      final Picard correction), so this chain is exact only to ~6e-7 relative
      -- the script reports the worst chain error over all rows, which is the
      stated numerical reason for the rel 1e-5 test tolerance.
    beta_eff_x = F_x*beta_internal/(F_x + beta_internal*stag_omega*cos(theta_x))
      (and y) -- EXACT given the dumped beta_internal
    uvel_basal = uvel_2d*beta_eff_x/beta_internal (Goldberg eq. 34 recovery)
    btractx = beta_internal*uvel_basal = beta_eff_x*uvel_2d (F90:3874, 5780)
    uvel[k] = uvel_basal + btractx*stag_omega_k[k]*cos(theta_x)/F_x, k = 1..10,
      uvel(11) = uvel_basal (stag_omega_k(11) = 0) -- and the same for vvel
      with the y-direction quantities (Goldberg eq. 32 + slope correction)

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


def main():
    ap = argparse.ArgumentParser()
    default_dir = os.path.join(MAIN, "stage1", "dumps", "cism")
    if not os.path.isdir(default_dir) and os.path.isdir(os.path.join(SIBLING_MAIN, "stage1", "dumps", "cism")):
        default_dir = os.path.join(SIBLING_MAIN, "stage1", "dumps", "cism")
    ap.add_argument("dumps_dir", nargs="?", default=default_dir)
    a = ap.parse_args()

    _, velo = load(os.path.join(a.dumps_dir, "cism_velo.csv"))
    _, nbr = load(os.path.join(a.dumps_dir, "cism_velo_nbr.csv"))
    ve = {f"{r['cism_step']}_{r['i']}_{r['j']}": r for r in velo if r['phase'] == 'exit'}
    ne = {f"{r['cism_step']}_{r['i']}_{r['j']}": r for r in nbr if r['phase'] == 'exit'}
    common = sorted(set(ve) & set(ne))
    assert len(common) == 128, len(common)
    print(f"# exit rows: velo {len(ve)}, nbr {len(ne)}, common {len(common)}", file=sys.stderr)

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
    for key in common:
        v, n = ve[key], ne[key]
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
                s = 0.0
                for k in range(kk, NL):
                    if active:
                        s += sig_mid[k] / efvs_q[p][k] * dsigma[k] * thck
                vals.append(s)
            omega_k_c.append(0.25 * sum(vals))
        if active:
            for kk in (0, 5, 10):
                check('omega_k_int', omega_k_c[kk], omega_k[kk], 1e-10)
        else:
            assert all(x == 0.0 for x in omega_k), (key, 'inactive omega_k nonzero')
            assert v['omega'] == 0.0 if 'omega' in v else True
        # omega trapezoid
        omega_c = sum(0.5 * (omega_k[k] + omega_k[k + 1]) * dsigma[k] for k in range(NL))
        check('omega_trap', omega_c, n['omega'], 1e-12)
        # stag_omega 4-cell mean
        check('stag_omega_avg',
              0.25 * (n['omega_sw'] + n['omega_se'] + n['omega_nw'] + n['omega_ne']),
              n['stag_omega'], 1e-12)

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
            check(f'beta_eff_{tag}', be_c, n[f'beta_eff_{tag}'], 1e-10)
        u2d, v2d = v['uvel_2d'], v['vvel_2d']
        u_b, v_b = v['uvel_011'], v['vvel_011']
        check('uvel_basal_rec', u2d * n['beta_eff_x'] / beta_int, u_b, 1e-10)
        check('vvel_basal_rec', v2d * n['beta_eff_y'] / beta_int, v_b, 1e-10)
        check('btractx1', beta_int * u_b, v['btractx'], 1e-10)
        check('btracty1', beta_int * v_b, v['btracty'], 1e-10)
        check('btractx2', n['beta_eff_x'] * u2d, v['btractx'], 1e-10)
        check('btracty2', n['beta_eff_y'] * v2d, v['btracty'], 1e-10)
        for tag, bb, uu in (('x', u_b, 'uvel'), ('y', v_b, 'vvel')):
            th = n[f'stag_theta_slope_{tag}']
            Fv = n[f'stag_diva_slope_factor_{tag}']
            corr = math.cos(th) / Fv
            bt = v[f'btract{"x" if tag == "x" else "y"}']
            for kk in range(NZ - 1):
                check(f'{uu}_profile', bb + bt * n[f'stag_omega_k_{kk + 1:03d}'] * corr,
                      v[f'{uu}_{kk + 1:03d}'], 1e-9)
        # beta chain at the exit basal velocity (the approximate link)
        beta_law = fg * basal_beta(u_b, v_b, v['topg'] - v['eus'], v['effecpress_stag'])
        chain_worst = max(chain_worst, rel(beta_law, beta_int))

    print(f"# chain worst rel (law at exit basal velocity vs dumped beta_internal): "
          f"{chain_worst:.3e}", file=sys.stderr)
    for name in sorted(stats):
        print(f"# {name:16s} {stats[name][0]:4d} checks, worst rel {stats[name][1]:.3e}",
              file=sys.stderr)
    assert fails == 0, f"{fails} law mismatches"
    print("# ALL LAWS VERIFIED on 128 exit rows", file=sys.stderr)

    # ---- emit test fragments ----------------------------------------------
    def overrides(row, names):
        return {k: row[k] for k in names}

    # DivaGeometry fragments
    out = {"geom": [], "integrals": [], "integrals_inactive": [], "velo": []}
    for step, i, j in GEOM_ROWS:
        key = f"{step}_{i}_{j}"
        v, n = ve[key], ne[key]
        ov = {}
        for tag in ('sw', 'se', 'nw', 'ne'):
            for q in ('usrf', 'thck', 'ice_mask', 'land_mask'):
                ov[f'{q}_{tag}'] = n[f'{q}_{tag}']
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
            "parameter_overrides": ov,
            "assertions": [{"variable": nm, "time": 0.0, "expected": src[nm],
                            "tolerance": {"rel": tol.get(nm, 1e-12)}} for nm in src],
        })

    # DivaIntegrals fragments
    for step, i, j in INTEGRAL_ROWS:
        key = f"{step}_{i}_{j}"
        v, n = ve[key], ne[key]
        ov = {"sigma": [n[f'sigma_{k:03d}'] for k in range(1, NZ + 1)],
              "efvs_q1": [n[f'efvs_q1_{k:03d}'] for k in range(1, NL + 1)],
              "efvs_q2": [n[f'efvs_q2_{k:03d}'] for k in range(1, NL + 1)],
              "efvs_q3": [n[f'efvs_q3_{k:03d}'] for k in range(1, NL + 1)],
              "efvs_q4": [n[f'efvs_q4_{k:03d}'] for k in range(1, NL + 1)],
              "thck": v["thck"], "ice_mask": v["ice_mask"],
              "omega_se": n["omega_se"], "omega_nw": n["omega_nw"], "omega_ne": n["omega_ne"]}
        # efvs lives in the nbr dump (cell-centered mean over quad points)
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
        out["integrals"].append({
            "id": f"int_step{step}_cell{i}{j}",
            "parameter_overrides": ov,
            "assertions": assertions,
        })
    for step, i, j in INTEGRAL_INACTIVE_ROWS:
        key = f"{step}_{i}_{j}"
        v, n = ve[key], ne[key]
        ov = {"sigma": [n[f'sigma_{k:03d}'] for k in range(1, NZ + 1)],
              "efvs_q1": [n[f'efvs_q1_{k:03d}'] for k in range(1, NL + 1)],
              "efvs_q2": [n[f'efvs_q2_{k:03d}'] for k in range(1, NL + 1)],
              "efvs_q3": [n[f'efvs_q3_{k:03d}'] for k in range(1, NL + 1)],
              "efvs_q4": [n[f'efvs_q4_{k:03d}'] for k in range(1, NL + 1)],
              "thck": v["thck"], "ice_mask": v["ice_mask"],
              "omega_se": n["omega_se"], "omega_nw": n["omega_nw"], "omega_ne": n["omega_ne"]}
        out["integrals_inactive"].append({
            "id": f"int_step{step}_cell{i}{j}_inactive",
            "parameter_overrides": ov,
            "assertions": [
                {"variable": "efvs", "time": 0.0, "expected": 0.0,
                 "coords": {"lev": 1}, "tolerance": {"rel": 1e-12}},
                {"variable": "omega_k", "time": 0.0, "expected": 0.0,
                 "coords": {"lev_nodes": 1}, "tolerance": {"rel": 1e-12}},
                {"variable": "omega", "time": 0.0, "expected": 0.0,
                 "tolerance": {"rel": 1e-12}},
                {"variable": "stag_omega", "time": 0.0, "expected": n["stag_omega"],
                 "tolerance": {"rel": 1e-12}},
            ],
        })

    # DivaVelocity fragments
    for step, i, j in VELO_ROWS:
        key = f"{step}_{i}_{j}"
        v, n = ve[key], ne[key]
        ov = {
            "BasalTraction.u": v["uvel_011"], "BasalTraction.v": v["vvel_011"],
            "BasalTraction.bed": v["topg"] - v["eus"], "BasalTraction.N_stag": v["effecpress_stag"],
            "uvel_2d": v["uvel_2d"], "vvel_2d": v["vvel_2d"],
            "uvel_basal": v["uvel_011"], "vvel_basal": v["vvel_011"],
            "f_ground": v["f_ground"],
            "omega": n["stag_omega"],
            "stag_omega_k": [n[f'stag_omega_k_{k:03d}'] for k in range(1, NZ + 1)],
            "theta_x": n["stag_theta_slope_x"], "theta_y": n["stag_theta_slope_y"],
            "F_x": n["stag_diva_slope_factor_x"], "F_y": n["stag_diva_slope_factor_y"],
        }
        src = {"beta_internal": n["beta_internal"], "beta_eff_x": n["beta_eff_x"],
               "beta_eff_y": n["beta_eff_y"], "btractx": v["btractx"], "btracty": v["btracty"],
               "uvel_basal_recovered": v["uvel_011"], "vvel_basal_recovered": v["vvel_011"]}
        assertions = [{"variable": nm, "time": 0.0, "expected": src[nm],
                       "tolerance": {"rel": 1e-5}} for nm in src]
        prof_levels = range(1, NZ + 1) if (step, i, j) == (1, 22, 13) else (1, 6, 11)
        for kk in prof_levels:
            assertions.append({"variable": "uvel", "time": 0.0, "expected": v[f'uvel_{kk:03d}'],
                               "coords": {"lev_nodes": kk}, "tolerance": {"rel": 1e-5}})
            assertions.append({"variable": "vvel", "time": 0.0, "expected": v[f'vvel_{kk:03d}'],
                               "coords": {"lev_nodes": kk}, "tolerance": {"rel": 1e-5}})
        out["velo"].append({"id": f"velo_step{step}_cell{i}{j}",
                            "parameter_overrides": ov, "assertions": assertions})

    json.dump(out, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
