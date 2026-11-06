#!/usr/bin/env python3
"""Validate the MOM6 stage-1 wave-2 boundary dumps (sw_abs, frazil, bflux,
diabatic_tri, kpp_wide).

Usage: python3 validate_mom6_w2.py DUMP_DIR

Checks, on the CSVs written by mom6_dump_mod from the instrumented
C_JRA_t62_mom6_inst run (see stage1/runs/MOM6_LANE.md):

1. CSV hygiene: headers present, no blank lines, every data row has
   exactly #columns fields; row counts consistent with the traced domain
   (12 ocean cells of PET 0 x 13 thermo steps, discovered from the data:
   sw_abs/bflux/tri 156 rows, frazil 312 (two calls per step)).
2. mom6_sw_abs.csv (join with mom6_bflux.csv on cell+tod):
   - surface SW three-way consistency: per band,
     sw_pen_band(n) == rho0*C_p*penSW_bnd(n)/dt and
     netPen(1) == sum_n penSW_bnd(n)/dt  (rho0=1035, C_p=3992 from the
     run's MOM_input; all internal rescale factors are 1.0);
   - per-band exponential decay pen_sw_rem(n,k) ==
     pen_sw_rem(n,k-1)*exp(-opacity(n,k)*h(k)) using the bflux h_after
     (the h passed to absorbRemainingSW), for every layer above the
     min_SW_heat cutoff (pen_sw_rem(n,k-1) > 1e-6 degC m);
   - total absorption: sum_k dT_k*h_k == sum_n (penSW_bnd(n) -
     pen_sw_rem(n,ke))  (absorbed = surface in minus bottom leftover;
     exact in the model arithmetic, including the thin-column
     redistribution path).
3. mom6_frazil.csv (two rows per cell+step, call=1 pre-diabatic,
   call=2 post-diabatic):
   - heat conservation per call:
     frazil_out - frazil_in == -rho0*C_p*sum_k h_k*(T_after_k - T_before_k);
   - freezing-point law: T_freeze_k == TFREEZE_S0_P0 + DTFREEZE_DS*S_k +
     DTFREEZE_DP*p_k with the run's LINEAR coefficients (0.0, -0.054,
     -7.75e-8) - the dump recomputes calculate_TFreeze with the same
     S and p, so this verifies the EOS coefficients used;
   - T_after_k >= T_freeze_k - 1e-9 for every layer (make_frazil never
     leaves a layer below freezing).
4. mom6_bflux.csv (applyBoundaryFluxesInOut totals vs the column
   before/after states):
   - mass: sum_k (h_after-h_before) == netMassInOut;
   - salt: sum_k (h_after*S_after - h_before*S_before) == netSalt;
   - heat: sum_k (h_after*T_after - h_before*T_before) ==
     netHeat + nonpenSW + sum_n penSW_bnd(n) - leftover, where leftover
     is the bottom pen_sw_rem from mom6_sw_abs.csv (0 for deep columns).
5. mom6_diabatic_tri.csv (tracer_vertdiff_Eulerian tridiagonal solve in
   diabatic_ALE; no surface/bottom fluxes in this call):
   - conservation: sum_k h_k*(T_after_k-T_before_k) == 0 (and S);
   - tridiagonal residual per layer:
     h_k*(T'_k-T_k) == ent_k*(T'_{k-1}-T'_k) + ent_{k+1}*(T'_{k+1}-T'_k)
     with ent_1 = ent_{ke+1} = 0;
   - cross-file: bflux h_after == tri h, bflux T_after == tri T_before,
     bflux S_after == tri S_before (only diagnostics run between the two
     calls).
6. mom6_kpp_wide.csv (widened KPP selection on every PE):
   - selection bits: (sel & 1) == (surfBuoyFlux < 0),
     (sel & 2) == (T_sfc > 20.0, the DUMP_WIDE_WARM_T default), sel != 0;
   - both selection kinds present in the trace (unstable and warm rows);
   - Kt/Ks/Kv update rule as in mom6_kpp.csv: Kt_out = Kd_heat where
     Kd_heat /= 0 else Kt_in (and same for salt/vis);
   - NLT structure: NLT_heat == NLT_salt, NLT(1) == 0, NLT == 0 when
     surfBuoyFlux >= 0;
   - geometry: z_inter(1) == 0, interfaces below the surface negative
     and increasing; kOBL in [1,ke]; z_inter(kOBL) >= -OBLdepth >
     z_inter(kOBL+1); uStar >= 0; OBLdepth > 0.
"""

import argparse
import sys

RHO0 = 1035.0   # RHO_0 default (not overridden in MOM_input)
CP = 3992.0     # C_P in the run's MOM_input
RHO0_CP = RHO0 * CP
TFREEZE_S0_P0 = 0.0       # TFREEZE_S0_P0 default (LINEAR form)
DTFREEZE_DS = -0.054      # DTFREEZE_DS default
DTFREEZE_DP = -7.75e-8    # DTFREEZE_DP in the run's MOM_input
WARM_T = 20.0             # DUMP_WIDE_WARM_T default
NCELLS = 12               # PET 0 ocean cells (SW corner block); sanity only
NSTEPS = 13               # thermo steps traced (cpl_step 1 covers two)


def read_csv(path):
    """Return (description, columns, rows) with the header lines stripped.
    Rows are lists of floats; the trailing comma of each row is stripped."""
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


def kcol(cols, name, k, fmt="{:02d}"):
    """Column index of name_k (k 1-based)."""
    return col_ix(cols, f"{name}_{fmt.format(k)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    args = ap.parse_args()
    dd = args.dump_dir
    failures = []

    # ---------------------------------------------------------------- load
    d_sw, c_sw, r_sw = read_csv(f"{dd}/mom6_sw_abs.csv")
    d_fr, c_fr, r_fr = read_csv(f"{dd}/mom6_frazil.csv")
    d_bf, c_bf, r_bf = read_csv(f"{dd}/mom6_bflux.csv")
    d_tr, c_tr, r_tr = read_csv(f"{dd}/mom6_diabatic_tri.csv")
    d_wd, c_wd, r_wd = read_csv(f"{dd}/mom6_kpp_wide.csv")

    # -------------------------------------------------------------- hygiene
    nsw = int(r_sw[0][col_ix(c_sw, "nsw")])
    ke = len([c for c in c_sw if c.startswith("dT_")])
    nz = ke + 1
    check(nsw == 2, f"sw_abs: nsw = {nsw} (expected 2, OHLMANN_03)", failures)
    cells = sorted({(int(r[col_ix(c_sw, "gi")]), int(r[col_ix(c_sw, "gj")]))
                    for r in r_sw})
    steps = sorted({(int(r[col_ix(c_sw, "cpl_step")]), int(r[col_ix(c_sw, "tod")]))
                    for r in r_sw})
    check(len(cells) == NCELLS,
          f"traced cells: {len(cells)} (expected {NCELLS})", failures)
    check(len(steps) == NSTEPS,
          f"traced steps: {len(steps)} (expected {NSTEPS})", failures)
    check(len(r_sw) == NCELLS * NSTEPS,
          f"sw_abs rows {len(r_sw)} != {NCELLS * NSTEPS}", failures)
    check(len(r_bf) == NCELLS * NSTEPS,
          f"bflux rows {len(r_bf)} != {NCELLS * NSTEPS}", failures)
    check(len(r_tr) == NCELLS * NSTEPS,
          f"tri rows {len(r_tr)} != {NCELLS * NSTEPS}", failures)
    check(len(r_fr) == NCELLS * NSTEPS * 2,
          f"frazil rows {len(r_fr)} != {NCELLS * NSTEPS * 2}", failures)
    check(len(r_wd) > 0, "wide: no rows", failures)
    check(ke == 75, f"ke = {ke} (expected 75)", failures)
    print(f"traced domain: {len(cells)} cells; {len(steps)} steps; "
          f"ke = {ke}; nsw = {nsw}")
    print(f"mom6_sw_abs.csv: {len(r_sw)} rows; mom6_bflux.csv: {len(r_bf)} rows; "
          f"mom6_diabatic_tri.csv: {len(r_tr)} rows; mom6_frazil.csv: {len(r_fr)} rows; "
          f"mom6_kpp_wide.csv: {len(r_wd)} rows")

    # Join keys: (gi, gj, cpl_step, tod)
    def key(r, c, i="gi", j="gj"):
        return (int(r[col_ix(c, i)]), int(r[col_ix(c, j)]),
                int(r[col_ix(c, "cpl_step")]), int(r[col_ix(c, "tod")]))

    bf_by_key = {key(r, c_bf): r for r in r_bf}
    sw_by_key = {key(r, c_sw): r for r in r_sw}
    tr_by_key = {key(r, c_tr): r for r in r_tr}

    # ------------------------------------------------------------ sw_abs
    n_sw3 = n_swdec = n_swtot = n_swtest = 0
    max_dec = 0.0
    for r in r_sw:
        kk = key(r, c_sw)
        rb = bf_by_key.get(kk)
        if rb is None:
            check(False, f"sw_abs: no bflux row for {kk}", failures)
            continue
        dt = rb[col_ix(c_bf, "dt")]
        # three-way surface consistency
        pen_sum = sum(r[col_ix(c_sw, f"sw_pen_band_{n}")] for n in range(1, nsw + 1))
        netpen1 = r[col_ix(c_sw, "netPen_01")]
        if abs(pen_sum - RHO0_CP * netpen1) <= 1e-9 * abs(pen_sum):
            n_sw3 += 1
        for n in range(1, nsw + 1):
            pb = rb[kcol(c_bf, "penSW_bnd", n)]
            sp = r[col_ix(c_sw, f"sw_pen_band_{n}")]
            if abs(sp - RHO0_CP * pb / dt) <= 1e-9 * abs(sp):
                n_sw3 += 1
        # exponential decay with the bflux h_after (testable only above the
        # min_SW_heat cutoff, prev > 1e-6 degC m; below it the model clamps
        # or zeroes SW_trans)
        for n in range(1, nsw + 1):
            prev = None
            for k in range(1, ke + 1):
                cur = r[kcol(c_sw, f"pen_sw_rem_{n}", k)]
                if prev is not None and prev > 1e-6:
                    n_swtest += 1
                    op = r[kcol(c_sw, f"opacity_{n}", k)]
                    h = rb[kcol(c_bf, "h_after", k)]
                    pred = prev * __import__("math").exp(-op * h)
                    e = rel_err(cur, pred)
                    max_dec = max(max_dec, e)
                    if e <= 1e-9:
                        n_swdec += 1
                    else:
                        check(False, f"sw_abs decay: {kk} band {n} k {k} "
                                     f"rel {e:.2e}", failures)
                prev = cur
        # total absorption
        lhs = sum(r[kcol(c_sw, "dT", k)] * rb[kcol(c_bf, "h_after", k)]
                  for k in range(1, ke + 1))
        rhs = sum(rb[kcol(c_bf, "penSW_bnd", n)]
                  - r[kcol(c_sw, f"pen_sw_rem_{n}", ke)]
                  for n in range(1, nsw + 1))
        if abs(lhs - rhs) <= 1e-9 * max(abs(rhs), 1e-30):
            n_swtot += 1
    check(n_sw3 == len(r_sw) * (1 + nsw),
          f"sw_abs surface consistency {n_sw3}/{len(r_sw) * (1 + nsw)}", failures)
    check(n_swtest > 0 and n_swdec == n_swtest,
          f"sw_abs decay law {n_swdec}/{n_swtest} testable layers "
          f"(max rel {max_dec:.2e})", failures)
    check(n_swtot == len(r_sw),
          f"sw_abs total absorption {n_swtot}/{len(r_sw)}", failures)
    print(f"sw_abs: surface consistency {n_sw3}/{len(r_sw) * (1 + nsw)}; "
          f"decay law {n_swdec}/{n_swtest} (max rel {max_dec:.2e}); "
          f"total absorption {n_swtot}/{len(r_sw)}")

    # -------------------------------------------------------------- frazil
    n_fr_cons = n_fr_tf = n_fr_tfrz = 0
    for r in r_fr:
        icall = int(r[col_ix(c_fr, "call")])
        check(icall in (1, 2), f"frazil: call {icall} not in (1,2)", failures)
        fin = r[col_ix(c_fr, "frazil_in")]
        fout = r[col_ix(c_fr, "frazil_out")]
        check(fin >= 0.0 and fout >= 0.0,
              f"frazil: negative frazil heat in/out {fin} {fout}", failures)
        h = [r[kcol(c_fr, "h", k)] for k in range(1, ke + 1)]
        tb = [r[kcol(c_fr, "T_before", k)] for k in range(1, ke + 1)]
        ta = [r[kcol(c_fr, "T_after", k)] for k in range(1, ke + 1)]
        s = [r[kcol(c_fr, "S", k)] for k in range(1, ke + 1)]
        p = [r[kcol(c_fr, "p", k)] for k in range(1, ke + 1)]
        tf = [r[kcol(c_fr, "T_freeze", k)] for k in range(1, ke + 1)]
        dE = RHO0_CP * sum(h[k] * (ta[k] - tb[k]) for k in range(ke))
        if abs((fout - fin) + dE) <= 1e-9 * max(abs(fout - fin), abs(dE), 1e-30):
            n_fr_cons += 1
        ok_tf = True
        for k in range(ke):
            law = TFREEZE_S0_P0 + DTFREEZE_DS * s[k] + DTFREEZE_DP * p[k]
            if rel_err(tf[k], law) > 1e-9:
                ok_tf = False
            if ta[k] < tf[k] - 1e-9:
                check(False, f"frazil: T_after_{k+1} < T_freeze "
                             f"({ta[k]} < {tf[k]})", failures)
            else:
                n_fr_tfrz += 1
        if ok_tf:
            n_fr_tf += 1
    check(n_fr_cons == len(r_fr),
          f"frazil heat conservation {n_fr_cons}/{len(r_fr)}", failures)
    check(n_fr_tf == len(r_fr),
          f"frazil T_freeze law {n_fr_tf}/{len(r_fr)}", failures)
    check(n_fr_tfrz == len(r_fr) * ke,
          f"frazil T_after >= T_freeze {n_fr_tfrz}/{len(r_fr) * ke}", failures)
    print(f"frazil: conservation {n_fr_cons}/{len(r_fr)}; "
          f"T_freeze law {n_fr_tf}/{len(r_fr)}; "
          f"T_after >= T_freeze {n_fr_tfrz}/{len(r_fr) * ke}")

    # ---------------------------------------------------------------- bflux
    n_bm = n_bs = n_bh = 0
    max_bh = 0.0
    for r in r_bf:
        kk = key(r, c_bf)
        hb = [r[kcol(c_bf, "h_before", k)] for k in range(1, ke + 1)]
        ha = [r[kcol(c_bf, "h_after", k)] for k in range(1, ke + 1)]
        tb = [r[kcol(c_bf, "T_before", k)] for k in range(1, ke + 1)]
        ta = [r[kcol(c_bf, "T_after", k)] for k in range(1, ke + 1)]
        sb = [r[kcol(c_bf, "S_before", k)] for k in range(1, ke + 1)]
        sa = [r[kcol(c_bf, "S_after", k)] for k in range(1, ke + 1)]
        dm = sum(ha[k] - hb[k] for k in range(ke))
        if abs(dm - r[col_ix(c_bf, "netMassInOut")]) <= 1e-9 * max(abs(dm), 1e-30):
            n_bm += 1
        ds = sum(ha[k] * sa[k] - hb[k] * sb[k] for k in range(ke))
        if abs(ds - r[col_ix(c_bf, "netSalt")]) <= 1e-9 * max(abs(ds), 1e-30):
            n_bs += 1
        dh = sum(ha[k] * ta[k] - hb[k] * tb[k] for k in range(ke))
        # netHeat already includes the non-penetrating SW (nonpenSW is a
        # diagnostic part of it); the penetrating SW enters per band, and
        # any SW that reaches the bottom leaves the column (leftover).
        rhs = (r[col_ix(c_bf, "netHeat")]
               + sum(r[kcol(c_bf, "penSW_bnd", n)] for n in range(1, 11)))
        rs = sw_by_key.get(kk)
        if rs is not None:
            rhs -= sum(rs[kcol(c_sw, f"pen_sw_rem_{n}", ke)]
                       for n in range(1, nsw + 1))
        e = rel_err(dh, rhs)
        max_bh = max(max_bh, e)
        if e <= 1e-9:
            n_bh += 1
    check(n_bm == len(r_bf), f"bflux mass {n_bm}/{len(r_bf)}", failures)
    check(n_bs == len(r_bf), f"bflux salt {n_bs}/{len(r_bf)}", failures)
    check(n_bh == len(r_bf),
          f"bflux heat {n_bh}/{len(r_bf)} (max rel {max_bh:.2e})", failures)
    print(f"bflux: mass {n_bm}/{len(r_bf)}; salt {n_bs}/{len(r_bf)}; "
          f"heat {n_bh}/{len(r_bf)} (max rel {max_bh:.2e})")

    # ----------------------------------------------------------------- tri
    n_tc = n_tres = 0
    n_cross = 0
    max_res = 0.0
    for r in r_tr:
        kk = key(r, c_tr)
        h = [r[kcol(c_tr, "h", k)] for k in range(1, ke + 1)]
        tb = [r[kcol(c_tr, "T_before", k)] for k in range(1, ke + 1)]
        ta = [r[kcol(c_tr, "T_after", k)] for k in range(1, ke + 1)]
        sb = [r[kcol(c_tr, "S_before", k)] for k in range(1, ke + 1)]
        sa = [r[kcol(c_tr, "S_after", k)] for k in range(1, ke + 1)]
        ent = [r[kcol(c_tr, "ent_t", k)] for k in range(1, nz + 1)]
        if abs(sum(h[k] * (ta[k] - tb[k]) for k in range(ke))) <= 1e-9 * \
                max(abs(sum(h[k] * tb[k] for k in range(ke))), 1e-30):
            n_tc += 1
        if abs(sum(h[k] * (sa[k] - sb[k]) for k in range(ke))) <= 1e-9 * \
                max(abs(sum(h[k] * sb[k] for k in range(ke))), 1e-30):
            n_tc += 1
        # tridiagonal residual (ent_1 = ent_{nz} = 0 checked implicitly:
        # the law below uses them at the boundaries).  Convective interfaces
        # carry ent ~ 1e6 m, so e_in/e_out are cancellations of huge terms
        # and the residual scale must include the operand magnitudes
        # |ent*T'|; with that scale the max residual is ~1e-12 (round-off).
        ok = True
        for k in range(ke):
            lhs = h[k] * (ta[k] - tb[k])
            e_in = ent[k] * (ta[k - 1] - ta[k]) if k >= 1 else 0.0
            e_out = ent[k + 1] * (ta[k + 1] - ta[k]) if k + 1 < ke else 0.0
            res = abs(lhs - (e_in + e_out))
            scale = (abs(lhs) + abs(e_in) + abs(e_out)
                     + abs(ent[k] * ta[k]) + abs(ent[k + 1] * ta[k]) + 1e-30)
            max_res = max(max_res, res / scale)
            if res > 1e-9 * scale:
                ok = False
        if ok:
            n_tres += 1
        rb = bf_by_key.get(kk)
        if rb is not None:
            hok = all(abs(h[k] - rb[kcol(c_bf, "h_after", k + 1)]) <= 1e-12
                      for k in range(ke))
            tok = all(abs(tb[k] - rb[kcol(c_bf, "T_after", k + 1)]) <= 1e-12
                      for k in range(ke))
            sok = all(abs(sb[k] - rb[kcol(c_bf, "S_after", k + 1)]) <= 1e-12
                      for k in range(ke))
            if hok and tok and sok:
                n_cross += 1
    check(n_tc == len(r_tr) * 2,
          f"tri conservation T/S {n_tc}/{len(r_tr) * 2}", failures)
    check(n_tres == len(r_tr),
          f"tri residual {n_tres}/{len(r_tr)} (max {max_res:.2e})", failures)
    check(n_cross == len(r_tr),
          f"cross-file bflux->tri (h, T, S) {n_cross}/{len(r_tr)}", failures)
    print(f"tri: conservation {n_tc}/{len(r_tr) * 2}; "
          f"residual {n_tres}/{len(r_tr)} (max {max_res:.2e}); "
          f"bflux->tri cross {n_cross}/{len(r_tr)}")

    # -------------------------------------------------------------- wide
    n_sel = n_upd = n_nlt = n_nlt1 = n_geom = 0
    n_unstable = n_warm = 0
    max_upd = 0.0
    for r in r_wd:
        sel = int(r[col_ix(c_wd, "sel")])
        bf = r[col_ix(c_wd, "surfBuoyFlux")]
        ts = r[col_ix(c_wd, "T_sfc")]
        if sel == ((1 if bf < 0.0 else 0) | (2 if ts > WARM_T else 0)) and sel != 0:
            n_sel += 1
        if bf < 0.0:
            n_unstable += 1
        if ts > WARM_T:
            n_warm += 1
        ok_upd = True
        for k in range(1, nz + 1):
            kd = r[kcol(c_wd, "Kd_heat", k)]
            ki = r[kcol(c_wd, "Kt_in", k)]
            ko = r[kcol(c_wd, "Kt_out", k)]
            want = kd if kd != 0.0 else ki
            if rel_err(ko, want) > 1e-12:
                ok_upd = False
            kd = r[kcol(c_wd, "Kd_salt", k)]
            ki = r[kcol(c_wd, "Ks_in", k)]
            ko = r[kcol(c_wd, "Ks_out", k)]
            if rel_err(ko, (kd if kd != 0.0 else ki)) > 1e-12:
                ok_upd = False
            kd = r[kcol(c_wd, "Kv", k)]
            ki = r[kcol(c_wd, "Kv_in", k)]
            ko = r[kcol(c_wd, "Kv_out", k)]
            if rel_err(ko, (kd if kd != 0.0 else ki)) > 1e-12:
                ok_upd = False
        if ok_upd:
            n_upd += 1
        nlt = [r[kcol(c_wd, "NLT_heat", k)] for k in range(1, nz + 1)]
        nls = [r[kcol(c_wd, "NLT_salt", k)] for k in range(1, nz + 1)]
        nlt_zero = all(v == 0.0 for v in nlt)
        if all(abs(nlt[k] - nls[k]) <= 1e-12 for k in range(nz)) and \
           (nlt_zero == (bf >= 0.0)):
            n_nlt += 1
        if nlt[0] == 0.0:
            n_nlt1 += 1
        zi = [r[kcol(c_wd, "z_inter", k)] for k in range(1, nz + 1)]
        obl = r[col_ix(c_wd, "OBLdepth")]
        kobl = r[col_ix(c_wd, "kOBL")]
        ust = r[col_ix(c_wd, "uStar")]
        kk = int(kobl)
        ok_geom = (zi[0] == 0.0 and all(zi[k] < 0.0 for k in range(1, nz))
                   and all(zi[k] > zi[k + 1] for k in range(nz - 1))
                   and 1 <= kk <= ke and obl > 0.0 and ust >= 0.0
                   and zi[kk - 1] >= -obl and (kk == ke or -obl > zi[kk]))
        if ok_geom:
            n_geom += 1
    check(n_sel == len(r_wd),
          f"wide selection bits {n_sel}/{len(r_wd)}", failures)
    check(n_unstable > 0, "wide: no unstable rows (surfBuoyFlux<0)", failures)
    check(n_warm > 0, "wide: no warm rows (T_sfc>20)", failures)
    check(n_upd == len(r_wd),
          f"wide Kt/Ks/Kv update rule {n_upd}/{len(r_wd)}", failures)
    check(n_nlt == len(r_wd),
          f"wide NLT structure {n_nlt}/{len(r_wd)}", failures)
    check(n_nlt1 == len(r_wd),
          f"wide NLT(1)==0 {n_nlt1}/{len(r_wd)}", failures)
    check(n_geom == len(r_wd),
          f"wide geometry {n_geom}/{len(r_wd)}", failures)
    print(f"wide: selection {n_sel}/{len(r_wd)} "
          f"({n_unstable} unstable, {n_warm} warm rows); "
          f"update rule {n_upd}/{len(r_wd)}; NLT {n_nlt}/{len(r_wd)} "
          f"(sfc-zero {n_nlt1}); geometry {n_geom}/{len(r_wd)}")

    if failures:
        print("\nFAILURES:")
        for f in failures[:40]:
            print("  " + f)
        print(f"\n{len(failures)} check(s) failed.")
        sys.exit(1)
    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
