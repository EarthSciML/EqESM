#!/usr/bin/env python3
"""Validate the MOM6 stage-1 dynamics boundary dumps (mom6_continuity.csv,
mom6_coradv.csv, mom6_pgf.csv) from the instrumented C_JRA_t62_mom6_dyn_inst
run (see stage1/runs/MOM6_LANE.md).

Usage: python3 validate_mom6_dyn.py DUMP_DIR

The dumps trace the first DUMP_DYN_STEPS=12 dynamics steps of the run
(DT=1800 s, FIRST_DIRECTION=0 => x_first; 12 ocean cells of PET 0 x 75
layers; continuity calls 1 (the BT_cont pure-PPM call) and 3 (the final
call); CorAdCalc corrector call 2; PressureForce_FV_Bouss call 1).

Land masks: land cells have Area_h == 0 (mom6_coradv.csv, exact) or, for
cells outside the coradv stencil, a constant thickness of 0.001 m in every
layer at every step (the model's land thickness fill).  The traced cells
are Antarctic-coastal (gi 2-7, gj 78-80); the land cells adjacent to the
edge_ok=1 cells (3,78) and (4,78) are (1,78),(2,78),(5,78),(6,78) and
(3,76),(3,77),(4,76),(4,77).

Laws pinned (all in physical units; conversion factors all 1.0 here):

1. mom6_continuity.csv (69 cols, 21600 rows = 12x75x12x2):
   - zonal convergence (all rows): h_inter = max(hin - dt*IareaT*
     (uh_E - uh_W), 0)  (continuity_zonal_convergence with hmin = 0).
   - meridional convergence (all rows): h_out = max(h_inter - dt*IareaT*
     (vh_N - vh_S), h_min) with h_min = GV%Angstrom_H (the dumped h_min).
     For call-3 rows the dumped hin/h_inter columns are aliased (the call
     passes the same array as hin and h, so the dump reads the final
     thickness); the call-3 input thickness is bit-identical to the call-1
     hin of the same (cell, layer, step) - h is not updated between the
     two calls - so the call-3 identities use the call-1 hin as the base.
   - zonal PPM reconstruction of the traced cell (call-1 rows with
     edge_ok=1): slp = 0.5*(h_{m+1} - h_{m-1}) limited by
     sign*min(|slp|, 2*min(dMx,dMn)), with slp = 0 when any of the 3
     stencil cells is land, and land neighbors replaced by the cell
     value; then h_W = 0.5*(h_im1+h) + (1/6)*(slp_{m-1}-slp_m),
     h_E = 0.5*(h_ip1+h) + (1/6)*(slp_m-slp_{m+1}); then PPM_limit_pos
     with h_min_recon = 2*Angstrom_H.  Same for cells (i+1,j), (i,j-1),
     (i,j+1) using the dumped hin_im4..ip3 / h_inter_jm4..jp3 stencils.
   - zonal flux law (call-1 rows, edge_ok=1): for u_E > 0,
     uh_E = dy_Cu_E*u_E*(h_E_c + CFL*(0.5*(h_W_c-h_E_c) + curv3*(CFL-1.5)))
     with CFL = u_E*dt*IdxT, curv3 = h_W_c+h_E_c-2*hin; the mirror with
     cell (i+1,j) for u_E < 0 (h_W_e/h_E_e, IdxT_e, hin_ip1); the west
     face with cell (i-1,j) edges reconstructed from the hin stencil for
     u_W > 0 (IdxT_w, hin_im1) and the dumped edges for u_W < 0; the
     meridional faces with h_S/h_N edges of the traced cell and (i,j+1)
     (h_S_c/h_N_c, h_S_n/h_N_n, IdyT, IdyT_n, h_inter, h_inter_jp1) and
     the (i,j-1) edges reconstructed from h_inter_jm4..jp1 for v_S > 0
     (IdyT_s, h_inter_jm1).  uh = 0 when the face velocity is exactly 0.
2. mom6_coradv.csv (137 cols, 10800 rows = 12x75x12; CorAdCalc call 2):
   - masks from Area_h == 0 (exact land indicator).
   - dvdx/dudy at the 4 corners (from the 6 dumped v/dyCv and u/dxCu
     faces), rel_vort = mask2dBu*(dvdx-dudy)*IareaBu (mask2dBu = product
     of the 4 surrounding cell masks), abs_vort = f + rel_vort,
     q = abs_vort*Ih_q.
   - Ih_q law at the NE corner (the only corner with all 4 hArea terms
     dumped): Ihq_ne = Aq_ne/((hAu_e+hAu_en)+(hAv_n+hAv_ne)+vol_neglect);
     the 8 hArea_u/hArea_v laws
     (hArea_u(I,j) = 0.5*(Ah(i,j)*h(i,j) + Ah(i+1,j)*h(i+1,j)) and the
     v mirror).
   - KE_c = 0.25*IareaT*(areaCu_w*u_w^2 + areaCu_e*u_e^2 + areaCv_s*v_s^2
     + areaCv_n*v_n^2)  (KE_ARAKAWA, ke_scheme 10).
   - CAu (SADOURNY75_ENERGY, coriolis_scheme 1):
     CAu_e_raw = 0.25*((q_ne*(vh_ne+vh_n)) + (q_se*(vh_s+vh_se)))*IdxCu_e
     clamped by bound_coriolis to [min(fv1..fv4), max(fv1..fv4)] with
     fv1..fv4 = av*v at the 4 surrounding q-points, then
     CAu_e = clamp - KEx_e; the CAu_w mirror.  CAv:
     CAv_n_raw = -0.25*((q_nw*(uh_w+uh_wn)) + (q_ne*(uh_e+uh_en)))*IdyCv_n
     clamped with fu1..fu4 = -av*u, then CAv_n = clamp - KEy_n; the
     CAv_s mirror.
   - grid completeness: areaT*IareaT = 1 and IdxCu*dxCu = IdyCv*dyCv = 1
     at the 4 faces (IareaBu*dxCu*dyCv != 1 on the curvilinear tripole
     grid, so the q-cell metric pairs are not checked).
3. mom6_pgf.csv (46 cols, 10800 rows = 12x75x12; PressureForce_FV_Bouss
   call 1):
   - e_bot = e_top - h*H_to_Z at the c/e/n cells (all rows);
   - e_top(k=1) == e_sfc; cross-layer e continuity e_top(k+1) == e_bot(k);
   - pa_top(k=1) = GxRho_ref*(e_sfc - Z_ref) + p_atm_c;
   - cross-layer pa continuity pa_top(k+1) == pa_bot(k);
   - sanity: H_to_Z = 1, h_neglect = 1e-30 m (H_subroundoff =
     1e-20*Angstrom_H), GxRho_ref = 9.8*rho0 (g_Earth = 9.8 in this
     run), I_Rho0 = 1/rho0, Z_ref = 0.

Tolerances: relative 1e-12 for the computed laws; exact for the
cross-layer/identity rows.
"""

import argparse
import math
import sys

RTOL = 1e-12
KE = 75                # layers of the HYCOM1 75-layer grid
NCELLS = 12            # PET 0 ocean cells traced (sanity only)
NSTEPS = 12            # dynamics steps traced (DUMP_DYN_STEPS)
RHO0 = 1035.0          # RHO_0 default
G_EARTH = 9.80665      # g_Earth default (m s-2)
ANGSTROM_H = 1.0e-10   # GV%Angstrom_H floor [m]


def read_csv(path):
    """Return (description, columns, rows) with the header lines stripped."""
    desc = None
    cols = None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("# description = "):
                desc = line[len("# description = "):]
            elif line.startswith("# columns = "):
                cols = line[len("# columns = "):].rstrip(",").split(",")
            elif line.strip() == "":
                raise ValueError(f"blank line in {path}")
            else:
                rows.append([float(x) for x in line.rstrip(",").split(",")])
    if desc is None or cols is None:
        raise ValueError(f"missing header in {path}")
    for n, row in enumerate(rows):
        if len(row) != len(cols):
            raise ValueError(f"{path} row {n}: {len(row)} fields != {len(cols)} columns")
    return desc, cols, rows


def col_ix(cols, name):
    return cols.index(name)


def rel_err(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def check(cond, msg, failures):
    if not cond:
        failures.append(msg)


# --------------------------------------------------------------------------
# PPM reconstruction (MOM_continuity_PPM.F90 PPM_reconstruction_x/y and
# PPM_limit_pos, with the monotonic=False / simple_2nd=False branches used
# by the reference run).
# --------------------------------------------------------------------------

def ppm_slp(h_at, mask_at, o):
    """Slope at offset o; needs h and mask at o-1, o, o+1."""
    h0 = h_at(o)
    hm1, hp1 = h_at(o - 1), h_at(o + 1)
    if not (mask_at(o - 1) and mask_at(o) and mask_at(o + 1)):
        return 0.0
    s = 0.5 * (hp1 - hm1)
    dMx = max(hp1, hm1, h0) - h0
    dMn = h0 - min(hp1, hm1, h0)
    if s == 0.0:
        return 0.0
    return math.copysign(1.0, s) * min(abs(s), 2.0 * min(dMx, dMn))


def ppm_edges(h_at, mask_at, h_min_recon):
    """PPM reconstruction + PPM_limit_pos for the cell at offset 0.
    Returns the limited (h_L, h_R) edge thicknesses."""
    h0 = h_at(0)
    slp = {o: ppm_slp(h_at, mask_at, o) for o in (-1, 0, 1)}
    him1 = mask_at(-1) * h_at(-1) + (1.0 - mask_at(-1)) * h0
    hip1 = mask_at(1) * h_at(1) + (1.0 - mask_at(1)) * h0
    hL = 0.5 * (him1 + h0) + (1.0 / 6.0) * (slp[-1] - slp[0])
    hR = 0.5 * (hip1 + h0) + (1.0 / 6.0) * (slp[0] - slp[1])
    # PPM_limit_pos
    curv = 3.0 * ((hL + hR) - 2.0 * h0)
    if curv > 0.0:
        dh = hR - hL
        if abs(dh) < curv:
            if h0 <= h_min_recon:
                hL = h0
                hR = h0
            elif 12.0 * curv * (h0 - h_min_recon) < (curv ** 2 + 3.0 * dh ** 2):
                scale = 12.0 * curv * (h0 - h_min_recon) / (curv ** 2 + 3.0 * dh ** 2)
                hL = h0 + scale * (hL - h0)
                hR = h0 + scale * (hR - h0)
    return hL, hR


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    args = ap.parse_args()
    dd = args.dump_dir
    failures = []

    d_cn, c_cn, r_cn = read_csv(f"{dd}/mom6_continuity.csv")
    d_cr, c_cr, r_cr = read_csv(f"{dd}/mom6_coradv.csv")
    d_pg, c_pg, r_pg = read_csv(f"{dd}/mom6_pgf.csv")
    ix_cr = {n: col_ix(c_cr, n) for n in c_cr}

    # ------------------------------------------------------------------ cont
    print("=== mom6_continuity.csv")
    ix = {n: col_ix(c_cn, n) for n in c_cn}
    cells = sorted({(int(r[ix["gi"]]), int(r[ix["gj"]])) for r in r_cn})
    steps = sorted({int(r[ix["dyn_step"]]) for r in r_cn})
    calls = sorted({int(r[ix["call"]]) for r in r_cn})
    check(len(r_cn) == NCELLS * KE * NSTEPS * 2,
          f"cont rows {len(r_cn)} != {NCELLS * KE * NSTEPS * 2}", failures)
    check(len(cells) == NCELLS, f"cont cells {len(cells)} != {NCELLS}", failures)
    check(steps == list(range(1, NSTEPS + 1)),
          f"cont dyn_steps {steps}", failures)
    check(calls == [1, 3], f"cont calls {calls} != [1, 3]", failures)
    check(all(int(r[ix["dt"]]) == 1800 for r in r_cn), "cont dt != 1800", failures)
    check(all(int(r[ix["x_first"]]) == 1 for r in r_cn),
          "cont x_first != 1 (FIRST_DIRECTION)", failures)
    check(all(r[ix["h_min"]] == ANGSTROM_H for r in r_cn),
          "cont h_min != GV%Angstrom_H", failures)
    check(all(r[ix["h_min_recon"]] == 2.0 * ANGSTROM_H for r in r_cn),
          "cont h_min_recon != 2*Angstrom_H", failures)
    # edge_ok: only the cells >= 2 from the PE compute-domain edges qualify
    # for the reconstruction/flux pins; the local domain of PET 0 ends at
    # local j=84, so only (3,78) and (4,78) qualify.
    eo = sorted({(int(r[ix["gi"]]), int(r[ix["gj"]]), int(r[ix["edge_ok"]]))
                 for r in r_cn})
    print(f"rows {len(r_cn)}; cells {len(cells)}; steps {steps[0]}..{steps[-1]}; "
          f"calls {calls}; edge_ok cells: " +
          ", ".join(f"({i},{j})" for i, j, ok in eo if ok))

    # masks: primary = Area_h == 0 from the coradv dump; extension = cells
    # outside the coradv stencil with a constant thickness (0.001 m) in
    # every row (the model's land thickness fill).
    land = set()
    area_map = {}
    for r in r_cr:
        gi, gj = int(r[ix_cr["gi"]]), int(r[ix_cr["gj"]])
        for name, di, dj in (("Ah_sw", -1, -1), ("Ah_s", 0, -1), ("Ah_se", 1, -1),
                             ("Ah_w", -1, 0), ("Ah_c", 0, 0), ("Ah_e", 1, 0),
                             ("Ah_nw", -1, 1), ("Ah_n", 0, 1), ("Ah_ne", 1, 1)):
            key = (gi + di, gj + dj)
            area_map[key] = max(area_map.get(key, 0.0), r[ix_cr[name]])
    for k, v in area_map.items():
        if v == 0.0:
            land.add(k)
    hmm_map = {}
    stencil_cols = ["hin"] + [f"hin_im{d}" for d in (4, 3, 2, 1)] + \
                   [f"hin_ip{d}" for d in (1, 2, 3)] + \
                   [f"hin_jm{d}" for d in (4, 3, 2, 1)] + \
                   [f"hin_jp{d}" for d in (1, 2, 3)]
    for r in r_cn:
        gi, gj = int(r[ix["gi"]]), int(r[ix["gj"]])
        for name in stencil_cols:
            di, dj = 0, 0
            if name.startswith("hin_im"): di = -int(name[6:])
            elif name.startswith("hin_ip"): di = int(name[6:])
            elif name.startswith("hin_jm"): dj = -int(name[6:])
            elif name.startswith("hin_jp"): dj = int(name[6:])
            key = (gi + di, gj + dj)
            lo, hi = hmm_map.get(key, (1e300, -1e300))
            hmm_map[key] = (min(lo, r[ix[name]]), max(hi, r[ix[name]]))
    for k, (lo, hi) in hmm_map.items():
        if k not in area_map and lo == hi:
            land.add(k)   # constant thickness => land fill value
    mask2dT = {}
    for k in set(area_map) | set(hmm_map):
        mask2dT[k] = 0.0 if k in land else 1.0
    print("land cells: " + ", ".join(f"({i},{j})" for i, j in sorted(land)))
    check(all(mask2dT.get((gi, gj), 1.0) == 1.0 for gi, gj in cells),
          "a traced cell inferred land!", failures)

    def m_at(cell):
        return mask2dT.get(cell, 1.0)

    # group rows by (gi, gj, dyn_step, call) -> layer index within the group
    groups = {}
    for r in r_cn:
        key = (int(r[ix["gi"]]), int(r[ix["gj"]]),
               int(r[ix["dyn_step"]]), int(r[ix["call"]]))
        groups.setdefault(key, []).append(r)
    check(len(groups) == NCELLS * NSTEPS * 2,
          f"cont groups {len(groups)} != {NCELLS * NSTEPS * 2}", failures)
    check(all(len(g) == KE for g in groups.values()),
          "cont group size != 75", failures)

    def get(colname, r):
        return r[ix[colname]]

    # ---- convergence identities (all rows)
    n_z = n_m = 0
    hin_c1 = {}
    for key, g in groups.items():
        gi, gj, step, call = key
        for k, r in enumerate(g):
            dt, iareat = get("dt", r), get("IareaT", r)
            uhE, uhW = get("uh_E", r), get("uh_W", r)
            vhN, vhS = get("vh_N", r), get("vh_S", r)
            if call == 1:
                hin = get("hin", r)
                hin_c1[(gi, gj, step, k)] = hin
                h_inter = get("h_inter", r)
            else:
                # call 3: hin/h_inter are aliased artifacts; the true input
                # is bit-identical to the call-1 hin (h is not updated
                # between the two calls of the same step)
                hin = hin_c1[(gi, gj, step, k)]
                h_inter = max(hin - dt * iareat * (uhE - uhW), 0.0)
            hi = max(hin - dt * iareat * (uhE - uhW), 0.0)
            if rel_err(hi, h_inter) > RTOL:
                n_z += 1
            ho = max(h_inter - dt * iareat * (vhN - vhS), get("h_min", r))
            if rel_err(ho, get("h_out", r)) > RTOL:
                n_m += 1
    check(n_z == 0, f"zonal convergence: {n_z} rows differ", failures)
    check(n_m == 0, f"meridional convergence: {n_m} rows differ", failures)
    print(f"convergence identities: zonal {n_z} bad, meridional {n_m} bad "
          f"of {len(groups) * KE} rows")

    # ---- reconstruction + flux-law pins (call-1 rows, edge_ok cells)
    n_rec = n_flux = 0
    n_rec_tot = n_flux_tot = 0
    for key, g in groups.items():
        gi, gj, step, call = key
        if call != 1:
            continue
        for r in g:
            if get("edge_ok", r) != 1.0:
                continue
            hmin_rec = get("h_min_recon", r)
            hin_off = {0: get("hin", r)}
            for d in (1, 2, 3):
                hin_off[d] = get(f"hin_ip{d}", r)
            for d in (-1, -2, -3, -4):
                hin_off[d] = get(f"hin_im{-d}", r)
            hi_off = {0: get("h_inter", r)}
            for d in (1, 2, 3):
                hi_off[d] = get(f"h_inter_jp{d}", r)
            for d in (-1, -2, -3, -4):
                hi_off[d] = get(f"h_inter_jm{-d}", r)

            def h_at(off, table):
                return table[off]

            def m_x(off, gjj):
                return m_at((gi + off, gjj))

            def m_y(off, gii):
                return m_at((gii, gj + off))

            # zonal reconstruction of the traced cell
            n_rec_tot += 2
            hL, hR = ppm_edges(lambda o: h_at(o, hin_off),
                               lambda o: m_x(o, gj), hmin_rec)
            if rel_err(hL, get("h_W_c", r)) > RTOL:
                n_rec += 1
            if rel_err(hR, get("h_E_c", r)) > RTOL:
                n_rec += 1
            # zonal reconstruction of (i+1,j)
            n_rec_tot += 2
            hL, hR = ppm_edges(lambda o: h_at(o + 1, hin_off),
                               lambda o: m_x(o + 1, gj), hmin_rec)
            if rel_err(hL, get("h_W_e", r)) > RTOL:
                n_rec += 1
            if rel_err(hR, get("h_E_e", r)) > RTOL:
                n_rec += 1
            # meridional reconstruction of the traced cell
            n_rec_tot += 2
            hL, hR = ppm_edges(lambda o: h_at(o, hi_off),
                               lambda o: m_y(o, gi), hmin_rec)
            if rel_err(hL, get("h_S_c", r)) > RTOL:
                n_rec += 1
            if rel_err(hR, get("h_N_c", r)) > RTOL:
                n_rec += 1
            # meridional reconstruction of (i,j+1)
            n_rec_tot += 2
            hL, hR = ppm_edges(lambda o: h_at(o + 1, hi_off),
                               lambda o: m_y(o + 1, gi), hmin_rec)
            if rel_err(hL, get("h_S_n", r)) > RTOL:
                n_rec += 1
            if rel_err(hR, get("h_N_n", r)) > RTOL:
                n_rec += 1
            # meridional reconstruction of (i,j-1) (for the vh_S law)
            n_rec_tot += 2
            hL_jm1, hR_jm1 = ppm_edges(lambda o: h_at(o - 1, hi_off),
                                       lambda o: m_y(o - 1, gi), hmin_rec)
            # zonal reconstruction of (i-1,j) (for the uh_W law)
            n_rec_tot += 2
            hL_im1, hR_im1 = ppm_edges(lambda o: h_at(o - 1, hin_off),
                                       lambda o: m_x(o - 1, gj), hmin_rec)

            # ---- flux laws (pure PPM at call 1: no uhbt/vhbt adjustment)
            dt = get("dt", r)
            hin, h_inter = get("hin", r), get("h_inter", r)
            # east face
            u = get("u_E", r)
            n_flux_tot += 1
            if u > 0.0:
                cfl = u * dt * get("IdxT", r)
                curv3 = get("h_W_c", r) + get("h_E_c", r) - 2.0 * hin
                uh = get("dy_Cu_E", r) * u * (get("h_E_c", r) + cfl *
                      (0.5 * (get("h_W_c", r) - get("h_E_c", r)) + curv3 * (cfl - 1.5)))
            elif u < 0.0:
                cfl = -u * dt * get("IdxT_e", r)
                curv3 = get("h_W_e", r) + get("h_E_e", r) - 2.0 * get("hin_ip1", r)
                uh = get("dy_Cu_E", r) * u * (get("h_W_e", r) + cfl *
                      (0.5 * (get("h_E_e", r) - get("h_W_e", r)) + curv3 * (cfl - 1.5)))
            else:
                uh = 0.0
            if rel_err(uh, get("uh_E", r)) > RTOL:
                n_flux += 1
            # west face
            u = get("u_W", r)
            n_flux_tot += 1
            if u > 0.0:
                cfl = u * dt * get("IdxT_w", r)
                curv3 = hL_im1 + hR_im1 - 2.0 * get("hin_im1", r)
                uh = get("dy_Cu_W", r) * u * (hR_im1 + cfl *
                      (0.5 * (hL_im1 - hR_im1) + curv3 * (cfl - 1.5)))
            elif u < 0.0:
                cfl = -u * dt * get("IdxT", r)
                curv3 = get("h_W_c", r) + get("h_E_c", r) - 2.0 * hin
                uh = get("dy_Cu_W", r) * u * (get("h_W_c", r) + cfl *
                      (0.5 * (get("h_E_c", r) - get("h_W_c", r)) + curv3 * (cfl - 1.5)))
            else:
                uh = 0.0
            if rel_err(uh, get("uh_W", r)) > RTOL:
                n_flux += 1
            # north face
            v = get("v_N", r)
            n_flux_tot += 1
            if v > 0.0:
                cfl = v * dt * get("IdyT", r)
                curv3 = get("h_S_c", r) + get("h_N_c", r) - 2.0 * h_inter
                vh = get("dx_Cv_N", r) * v * (get("h_N_c", r) + cfl *
                     (0.5 * (get("h_S_c", r) - get("h_N_c", r)) + curv3 * (cfl - 1.5)))
            elif v < 0.0:
                cfl = -v * dt * get("IdyT_n", r)
                curv3 = get("h_S_n", r) + get("h_N_n", r) - 2.0 * get("h_inter_jp1", r)
                vh = get("dx_Cv_N", r) * v * (get("h_S_n", r) + cfl *
                     (0.5 * (get("h_N_n", r) - get("h_S_n", r)) + curv3 * (cfl - 1.5)))
            else:
                vh = 0.0
            if rel_err(vh, get("vh_N", r)) > RTOL:
                n_flux += 1
            # south face
            v = get("v_S", r)
            n_flux_tot += 1
            if v > 0.0:
                cfl = v * dt * get("IdyT_s", r)
                curv3 = hL_jm1 + hR_jm1 - 2.0 * get("h_inter_jm1", r)
                vh = get("dx_Cv_S", r) * v * (hR_jm1 + cfl *
                     (0.5 * (hL_jm1 - hR_jm1) + curv3 * (cfl - 1.5)))
            elif v < 0.0:
                cfl = -v * dt * get("IdyT", r)
                curv3 = get("h_S_c", r) + get("h_N_c", r) - 2.0 * h_inter
                vh = get("dx_Cv_S", r) * v * (get("h_S_c", r) + cfl *
                     (0.5 * (get("h_N_c", r) - get("h_S_c", r)) + curv3 * (cfl - 1.5)))
            else:
                vh = 0.0
            if rel_err(vh, get("vh_S", r)) > RTOL:
                n_flux += 1
    check(n_rec == 0, f"PPM reconstruction: {n_rec}/{n_rec_tot} edges differ", failures)
    check(n_flux == 0, f"flux laws: {n_flux}/{n_flux_tot} faces differ", failures)
    print(f"reconstruction pins: {n_rec_tot - n_rec}/{n_rec_tot} ok; "
          f"flux-law pins: {n_flux_tot - n_flux}/{n_flux_tot} ok")

    # ----------------------------------------------------------------- coradv
    print("=== mom6_coradv.csv")
    ix = {n: col_ix(c_cr, n) for n in c_cr}
    cells_cr = sorted({(int(r[ix["gi"]]), int(r[ix["gj"]])) for r in r_cr})
    steps_cr = sorted({int(r[ix["dyn_step"]]) for r in r_cr})
    check(len(r_cr) == NCELLS * KE * NSTEPS,
          f"coradv rows {len(r_cr)} != {NCELLS * KE * NSTEPS}", failures)
    check(cells_cr == cells, "coradv cells != cont cells", failures)
    check(steps_cr == steps, "coradv steps != cont steps", failures)
    check(all(int(r[ix["call"]]) == 2 for r in r_cr), "coradv call != 2", failures)
    check(all(int(r[ix["bound_coriolis"]]) == 1 for r in r_cr),
          "coradv bound_coriolis != 1", failures)
    check(all(int(r[ix["coriolis_scheme"]]) == 1 for r in r_cr),
          "coradv coriolis_scheme != 1 (SADOURNY75_ENERGY)", failures)
    check(all(int(r[ix["ke_scheme"]]) == 10 for r in r_cr),
          "coradv ke_scheme != 10 (KE_ARAKAWA)", failures)
    print(f"rows {len(r_cr)}; cells {len(cells_cr)}; steps {steps_cr[0]}..{steps_cr[-1]}")

    def cm(cell):
        return mask2dT.get(cell, 1.0)

    n_rv = n_av = n_q = n_ihq = n_harea = n_ke = n_ca = n_g = n_dv = 0
    n_rv_t = n_av_t = n_q_t = n_ihq_t = n_harea_t = n_ke_t = n_ca_t = n_g_t = n_dv_t = 0
    for r in r_cr:
        # dvdx/dudy at the 4 corners: dvdx(I,J) = v(i+1,J)*dyCv(i+1,J)
        # - v(i,J)*dyCv(i,J), dudy(I,J) = u(I,j+1)*dxCu(I,j+1) - u(I,j)*dxCu(I,j),
        # mapped to the dump's per-cell columns (corner 'se' of cell (i,j) is the
        # q-point (i, j-1), corner 'sw' the q-point (i-1, j-1), etc.).
        dvdx = {
            "sw": r[ix["v_s"]] * r[ix["dyCv_s"]] - r[ix["v_sw"]] * r[ix["dyCv_sw"]],
            "se": r[ix["v_se"]] * r[ix["dyCv_se"]] - r[ix["v_s"]] * r[ix["dyCv_s"]],
            "nw": r[ix["v_n"]] * r[ix["dyCv_n"]] - r[ix["v_nw"]] * r[ix["dyCv_nw"]],
            "ne": r[ix["v_ne"]] * r[ix["dyCv_ne"]] - r[ix["v_n"]] * r[ix["dyCv_n"]],
        }
        dudy = {
            "sw": r[ix["u_w"]] * r[ix["dxCu_w"]] - r[ix["u_ws"]] * r[ix["dxCu_ws"]],
            "se": r[ix["u_e"]] * r[ix["dxCu_e"]] - r[ix["u_es"]] * r[ix["dxCu_es"]],
            "nw": r[ix["u_wn"]] * r[ix["dxCu_wn"]] - r[ix["u_w"]] * r[ix["dxCu_w"]],
            "ne": r[ix["u_en"]] * r[ix["dxCu_en"]] - r[ix["u_e"]] * r[ix["dxCu_e"]],
        }
        gi, gj = int(r[ix["gi"]]), int(r[ix["gj"]])
        # corner masks: mask2dBu(I,J) = product of the 4 surrounding cell masks
        mbu = {
            "sw": cm((gi - 1, gj - 1)) * cm((gi, gj - 1)) * cm((gi - 1, gj)) * cm((gi, gj)),
            "se": cm((gi, gj - 1)) * cm((gi + 1, gj - 1)) * cm((gi, gj)) * cm((gi + 1, gj)),
            "nw": cm((gi - 1, gj)) * cm((gi, gj)) * cm((gi - 1, gj + 1)) * cm((gi, gj + 1)),
            "ne": cm((gi, gj)) * cm((gi + 1, gj)) * cm((gi, gj + 1)) * cm((gi + 1, gj + 1)),
        }
        for cn in ("sw", "se", "nw", "ne"):
            n_rv_t += 1
            rv = mbu[cn] * (dvdx[cn] - dudy[cn]) * r[ix[f"IareaBu_{cn}"]]
            if rel_err(rv, r[ix[f"rv_{cn}"]]) > RTOL:
                n_rv += 1
            n_av_t += 1
            if rel_err(r[ix[f"f_{cn}"]] + r[ix[f"rv_{cn}"]], r[ix[f"av_{cn}"]]) > RTOL:
                n_av += 1
            n_q_t += 1
            if rel_err(r[ix[f"av_{cn}"]] * r[ix[f"Ihq_{cn}"]], r[ix[f"q_{cn}"]]) > RTOL:
                n_q += 1
            # direct pin of the dump's dvdx/dudy columns against the u/v stencil
            n_dv_t += 1
            if rel_err(dvdx[cn], r[ix[f"dvdx_{cn}"]]) > RTOL or \
               rel_err(dudy[cn], r[ix[f"dudy_{cn}"]]) > RTOL:
                n_dv += 1
        # Ih_q law at the NE corner
        n_ihq_t += 1
        hArea_q = (r[ix["hAu_e"]] + r[ix["hAu_en"]]) + (r[ix["hAv_n"]] + r[ix["hAv_ne"]])
        ihq = r[ix["Aq_ne"]] / (hArea_q + r[ix["vol_neglect"]])
        if rel_err(ihq, r[ix["Ihq_ne"]]) > RTOL:
            n_ihq += 1
        # hArea_u / hArea_v laws
        for name, a, b in (("hAu_w", "w", "c"), ("hAu_e", "c", "e"),
                           ("hAu_wn", "nw", "n"), ("hAu_en", "n", "ne"),
                           ("hAv_s", "s", "c"), ("hAv_se", "se", "e"),
                           ("hAv_n", "c", "n"), ("hAv_ne", "e", "ne")):
            n_harea_t += 1
            ha = 0.5 * (r[ix[f"Ah_{a}"]] * r[ix[f"h_{a}"]] +
                        r[ix[f"Ah_{b}"]] * r[ix[f"h_{b}"]])
            if rel_err(ha, r[ix[name]]) > RTOL:
                n_harea += 1
        # KE_c (KE_ARAKAWA)
        n_ke_t += 1
        ke = 0.25 * r[ix["IareaT"]] * (
            r[ix["areaCu_w"]] * r[ix["u_w"]] ** 2 +
            r[ix["areaCu_e"]] * r[ix["u_e"]] ** 2 +
            r[ix["areaCv_s"]] * r[ix["v_s"]] ** 2 +
            r[ix["areaCv_n"]] * r[ix["v_n"]] ** 2)
        if rel_err(ke, r[ix["KE_c"]]) > RTOL:
            n_ke += 1
        # CAu_e
        n_ca_t += 1
        raw = 0.25 * ((r[ix["q_ne"]] * (r[ix["vh_ne"]] + r[ix["vh_n"]])) +
                      (r[ix["q_se"]] * (r[ix["vh_s"]] + r[ix["vh_se"]]))) * r[ix["IdxCu_e"]]
        fv = (r[ix["av_ne"]] * r[ix["v_ne"]], r[ix["av_ne"]] * r[ix["v_n"]],
              r[ix["av_se"]] * r[ix["v_se"]], r[ix["av_se"]] * r[ix["v_s"]])
        ca = min(max(raw, min(fv)), max(fv)) - r[ix["KEx_e"]]
        if rel_err(ca, r[ix["CAu_e"]]) > RTOL:
            n_ca += 1
        # CAu_w
        n_ca_t += 1
        raw = 0.25 * ((r[ix["q_nw"]] * (r[ix["vh_n"]] + r[ix["vh_nw"]])) +
                      (r[ix["q_sw"]] * (r[ix["vh_sw"]] + r[ix["vh_s"]]))) * r[ix["IdxCu_w"]]
        fv = (r[ix["av_nw"]] * r[ix["v_n"]], r[ix["av_nw"]] * r[ix["v_nw"]],
              r[ix["av_sw"]] * r[ix["v_s"]], r[ix["av_sw"]] * r[ix["v_sw"]])
        ca = min(max(raw, min(fv)), max(fv)) - r[ix["KEx_w"]]
        if rel_err(ca, r[ix["CAu_w"]]) > RTOL:
            n_ca += 1
        # CAv_n
        n_ca_t += 1
        raw = -0.25 * ((r[ix["q_nw"]] * (r[ix["uh_w"]] + r[ix["uh_wn"]])) +
                       (r[ix["q_ne"]] * (r[ix["uh_e"]] + r[ix["uh_en"]]))) * r[ix["IdyCv_n"]]
        fu = (-r[ix["av_ne"]] * r[ix["u_en"]], -r[ix["av_ne"]] * r[ix["u_e"]],
              -r[ix["av_nw"]] * r[ix["u_wn"]], -r[ix["av_nw"]] * r[ix["u_w"]])
        ca = min(max(raw, min(fu)), max(fu)) - r[ix["KEy_n"]]
        if rel_err(ca, r[ix["CAv_n"]]) > RTOL:
            n_ca += 1
        # CAv_s (v-point (i, J-1): q_sw/q_se corners, uh at the 4 u-points
        # (I-1,J-1),(I-1,J),(I,J-1),(I,J) = uh_ws, uh_w, uh_es, uh_e)
        n_ca_t += 1
        raw = -0.25 * ((r[ix["q_sw"]] * (r[ix["uh_ws"]] + r[ix["uh_w"]])) +
                       (r[ix["q_se"]] * (r[ix["uh_es"]] + r[ix["uh_e"]]))) * r[ix["IdyCv_s"]]
        fu = (-r[ix["av_se"]] * r[ix["u_e"]], -r[ix["av_se"]] * r[ix["u_es"]],
              -r[ix["av_sw"]] * r[ix["u_w"]], -r[ix["av_sw"]] * r[ix["u_ws"]])
        ca = min(max(raw, min(fu)), max(fu)) - r[ix["KEy_s"]]
        if rel_err(ca, r[ix["CAv_s"]]) > RTOL:
            n_ca += 1
        # grid completeness (note: IareaBu*dxCu*dyCv != 1 on the curvilinear
        # tripole grid - the q-cell areas are computed from the corner
        # coordinates - so only the reciprocal metric pairs are checked)
        n_g_t += 1
        if rel_err(r[ix["areaT"]] * r[ix["IareaT"]], 1.0) > RTOL:
            n_g += 1
        n_g_t += 4
        for a, b in (("IdxCu_w", "dxCu_w"), ("IdxCu_e", "dxCu_e"),
                     ("IdyCv_s", "dyCv_s"), ("IdyCv_n", "dyCv_n")):
            if rel_err(r[ix[a]] * r[ix[b]], 1.0) > RTOL:
                n_g += 1
    check(n_rv == 0, f"coradv rel_vort: {n_rv}/{n_rv_t} corners differ", failures)
    check(n_av == 0, f"coradv abs_vort: {n_av}/{n_av_t} corners differ", failures)
    check(n_q == 0, f"coradv q: {n_q}/{n_q_t} corners differ", failures)
    check(n_dv == 0, f"coradv dvdx/dudy: {n_dv}/{n_dv_t} corners differ", failures)
    check(n_ihq == 0, f"coradv Ih_q: {n_ihq}/{n_ihq_t} NE corners differ", failures)
    check(n_harea == 0, f"coradv hArea: {n_harea}/{n_harea_t} faces differ", failures)
    check(n_ke == 0, f"coradv KE_c: {n_ke}/{n_ke_t} rows differ", failures)
    check(n_ca == 0, f"coradv CAu/CAv: {n_ca}/{n_ca_t} faces differ", failures)
    check(n_g == 0, f"coradv grid: {n_g}/{n_g_t} metrics differ", failures)
    print(f"rv {n_rv_t - n_rv}/{n_rv_t}; av {n_av_t - n_av}/{n_av_t}; "
          f"q {n_q_t - n_q}/{n_q_t}; dvdx/dudy {n_dv_t - n_dv}/{n_dv_t}; "
          f"Ihq_ne {n_ihq_t - n_ihq}/{n_ihq_t}; "
          f"hArea {n_harea_t - n_harea}/{n_harea_t}; "
          f"KE_c {n_ke_t - n_ke}/{n_ke_t}; CAu/CAv {n_ca_t - n_ca}/{n_ca_t}; "
          f"grid {n_g_t - n_g}/{n_g_t}")

    # ------------------------------------------------------------------- pgf
    print("=== mom6_pgf.csv")
    ix = {n: col_ix(c_pg, n) for n in c_pg}
    cells_pg = sorted({(int(r[ix["gi"]]), int(r[ix["gj"]])) for r in r_pg})
    steps_pg = sorted({int(r[ix["dyn_step"]]) for r in r_pg})
    check(len(r_pg) == NCELLS * KE * NSTEPS,
          f"pgf rows {len(r_pg)} != {NCELLS * KE * NSTEPS}", failures)
    check(cells_pg == cells, "pgf cells != cont cells", failures)
    check(steps_pg == steps, "pgf steps != cont steps", failures)
    check(all(int(r[ix["call"]]) == 1 for r in r_pg), "pgf call != 1", failures)
    print(f"rows {len(r_pg)}; cells {len(cells_pg)}; steps {steps_pg[0]}..{steps_pg[-1]}")

    groups = {}
    for r in r_pg:
        key = (int(r[ix["gi"]]), int(r[ix["gj"]]), int(r[ix["dyn_step"]]))
        groups.setdefault(key, []).append(r)
    check(all(len(g) == KE for g in groups.values()), "pgf group size != 75", failures)

    n_eb = n_es = n_xl = n_pa1 = n_san = 0
    n_eb_t = n_es_t = n_xl_t = n_pa1_t = n_san_t = 0
    for g in groups.values():
        for k, r in enumerate(g):
            hz = r[ix["H_to_Z"]]
            for cell in ("c", "e", "n"):
                n_eb_t += 1
                if rel_err(r[ix[f"e_bot_{cell}"]],
                           r[ix[f"e_top_{cell}"]] - r[ix[f"h_{cell}"]] * hz) > RTOL:
                    n_eb += 1
            if k == 0:
                n_es_t += 1
                if r[ix["e_top_c"]] != r[ix["e_sfc_c"]]:
                    n_es += 1
                n_pa1_t += 1
                pa = r[ix["GxRho_ref"]] * (r[ix["e_sfc_c"]] - r[ix["Z_ref"]]) + r[ix["p_atm_c"]]
                if rel_err(pa, r[ix["pa_top_c"]]) > RTOL:
                    n_pa1 += 1
            if k > 0:
                # cross-layer continuity: this row's top interface is the
                # previous row's bottom interface (same array element)
                n_xl_t += 2
                if r[ix["e_top_c"]] != g[k - 1][ix["e_bot_c"]]:
                    n_xl += 1
                if r[ix["pa_top_c"]] != g[k - 1][ix["pa_bot_c"]]:
                    n_xl += 1
    # sanity: constants (this run: g_Earth = 9.8 m s-2, RHO_0 = 1035 kg m-3,
    # Z_ref = 0, H_subroundoff = 1e-20*Angstrom_H = 1e-30 m; all H/Z units
    # are m, so H_to_Z = 1)
    for r in r_pg[:NCELLS * KE]:
        n_san_t += 5
        if r[ix["H_to_Z"]] != 1.0:
            n_san += 1
        if abs(r[ix["h_neglect"]] - 1e-30) > 1e-40:
            n_san += 1
        if rel_err(r[ix["GxRho_ref"]], 9.8 * RHO0) > 1e-9:
            n_san += 1
        if rel_err(r[ix["I_Rho0"]], 1.0 / RHO0) > 1e-9:
            n_san += 1
        if r[ix["Z_ref"]] != 0.0:
            n_san += 1
    check(n_eb == 0, f"pgf e_bot = e_top - h*H_to_Z: {n_eb}/{n_eb_t} differ", failures)
    check(n_es == 0, f"pgf e_top(k=1) != e_sfc: {n_es}/{n_es_t}", failures)
    check(n_xl == 0, f"pgf cross-layer e/pa: {n_xl}/{n_xl_t} differ", failures)
    check(n_pa1 == 0, f"pgf pa_top(k=1) law: {n_pa1}/{n_pa1_t} differ", failures)
    check(n_san == 0, f"pgf constants sanity: {n_san}/{n_san_t} bad", failures)
    print(f"e_bot law {n_eb_t - n_eb}/{n_eb_t}; e_sfc {n_es_t - n_es}/{n_es_t}; "
          f"cross-layer {n_xl_t - n_xl}/{n_xl_t}; pa_top(k=1) {n_pa1_t - n_pa1}/{n_pa1_t}; "
          f"constants {n_san_t - n_san}/{n_san_t}")

    print()
    if failures:
        print(f"FAILED: {len(failures)} checks failed")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
