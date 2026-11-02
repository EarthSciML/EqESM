#!/usr/bin/env python3
"""Validate the MOM6 stage-1 boundary dumps (EOS, CVMix KPP, surface forcing).

Usage: python3 validate_mom6_dumps.py DUMP_DIR [--eos-tol TOL]

Checks, on the CSVs written by mom6_dump_mod from the instrumented
C_JRA_t62_mom6_inst run (see stage1/runs/MOM6_LANE.md):

1. CSV hygiene: `# description`/`# columns` headers present, no blank
   lines, every data row has exactly #columns fields (trailing comma
   stripped), expected row counts (the traced domain - the number of
   traced cells and steps - is discovered from the CSV keys: PET 0's SW
   corner block has 12 ocean cells and the trace covers 12 coupling
   steps).
2. EOS law (mom6_eos.csv): rho(T,S,p), drho/dT, drho/dS reproduced by the
   Wright 1997 full-range law (MOM_EOS_Wright_full.F90 density_elem /
   calculate_density_derivs_elem, constants at :24-38) to round-off; the
   Python reimplementation mirrors the Fortran operation order exactly.
3. Hydrostatic pressure consistency (mom6_eos.csv + mom6_kpp_obl.csv):
   p_{k+1} - p_k = g*H_to_RZ * h_k with a constant ratio per column (the
   Boussinesq pRef accumulation in KPP_compute_BLD).
4. KPP boundary layer (mom6_kpp.csv): with KPP_IS_ADDITIVE=False, the
   final Kt/Ks/Kv equal the CVMix Kdiffusivity/Kviscosity where those are
   non-zero, else the interior values; NLT_heat == NLT_salt; NLT is zero at
   the surface interface and everywhere when surfBuoyFlux >= 0; the OBL
   depth lies within the kOBL cell (z_inter(kbl+1) <= -OBLdepth <=
   z_inter(kbl), kbl = int(kOBL)); kOBL fractional part in [0,1).
5. OBL-depth boundary (mom6_kpp_obl.csv): kOBL >= 1, OBLdepth > 0 and
   within the column; z_inter identical to the same cell/step in
   mom6_kpp.csv (same column geometry in the two calls).
6. Surface forcing (mom6_sfc_forcing.csv): the fluxes/forces are the
   convert_IOB_to_fluxes/forces mappings of the imported fields (all unit
   conversion factors are 1.0 in this run): latent = q_flux*Lv - fprec*Lf
   - frunoff*Lf - frunoff_glc*Lf (Lv = 2.501e6, Lf = 3.337e5 J/kg),
   latent_evap_diag = q_flux*Lv, latent_fprec_diag = -fprec*Lf,
   latent_frunoff_diag = -frunoff*Lf, sw = sw_vis_dir + sw_vis_dif +
   sw_nir_dir + sw_nir_dif (imported components), evap = q_flux, lprec_o =
   lprec, fprec_o = fprec, lrunoff_o = lrunoff, frunoff_o = frunoff,
   salt_flux_o = salt_flux, seaice_melt_heat_o = seaice_melt_heat,
   seaice_melt_o = seaice_melt, ice_fraction_o = ice_fraction,
   p_surf_full = p and p_surf_full_o = p (fluxes and forces), tau_mag_o =
   GUST_CONST + sqrt(u_flux^2 + v_flux^2) (GUST_CONST = 4.0e-4 Pa; the
   AGRID tau_mag is built from the imported u_flux/v_flux, NOT from
   taux/tauy, which are u-grid averages), ustar_o = sqrt((GUST_CONST +
   sqrt(u_flux^2+v_flux^2))/RHO_0) (RHO_0 = 1035 kg m-3), omega_w2x =
   atan2(v_flux, u_flux).

Exit code 0 if all checks pass; prints a summary of maxima and row counts.
"""
import argparse
import math
import sys

# ---------------------------------------------------------------------------
# Wright 1997 full-range EOS constants (MOM_EOS_Wright_full.F90:24-38)
A0 = 7.133718e-4
A1 = 2.724670e-7
A2 = -1.646582e-7
B0 = 5.613770e8
B1 = 3.600337e6
B2 = -3.727194e4
B3 = 1.660557e2
B4 = 6.844158e5
B5 = -8.389457e3
C0 = 1.609893e5
C1 = 8.427815e2
C2 = -6.931554
C3 = 3.869318e-2
C4 = -1.664201e2
C5 = -2.765195

# Run configuration (C_JRA_t62 reference run, MOM_parameter_doc.all)
GUST_CONST = 4.0e-4   # Pa
RHO_0 = 1035.0        # kg m-3
LF = 3.337e5          # J/kg latent heat of fusion
LV = 2.501e6          # J/kg latent heat of vaporization

NLAYERS = 75          # NK of the tx2_3v2 grid (MOM_input NK)
KE = NLAYERS
NZ = KE + 1


def wright_eos(T, S, p):
    """rho, drho_dT, drho_dS from the Wright full-range law, mirroring the
    Fortran operation order in MOM_EOS_Wright_full.F90 exactly."""
    al0 = A0 + (A1 * T + A2 * S)
    p0 = B0 + (B4 * S + T * (B1 + (T * (B2 + B3 * T) + B5 * S)))
    lam = C0 + (C4 * S + T * (C1 + (T * (C2 + C3 * T) + C5 * S)))
    rho = (p + p0) / (lam + al0 * (p + p0))
    d = lam + al0 * (p + p0)
    I_denom2 = 1.0 / (d * d)
    drho_dT = I_denom2 * (lam * (B1 + (T * (2.0 * B2 + 3.0 * B3 * T) + B5 * S)) -
                          (p + p0) * ((p + p0) * A1 + (C1 + (T * (C2 * 2.0 + C3 * 3.0 * T) + C5 * S))))
    drho_dS = I_denom2 * (lam * (B4 + B5 * T) - (p + p0) * ((p + p0) * A2 + (C4 + C5 * T)))
    return rho, drho_dT, drho_dS


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--eos-tol", type=float, default=1e-10)
    ap.add_argument("--kpp-tol", type=float, default=1e-9)
    args = ap.parse_args()
    D = args.dump_dir.rstrip("/") + "/"
    failures = []

    # ------------------------------------------------------------------
    # 1. CSV hygiene + row counts
    # ------------------------------------------------------------------
    eos_desc, eos_cols, eos_rows = read_csv(D + "mom6_eos.csv")
    obl_desc, obl_cols, obl_rows = read_csv(D + "mom6_kpp_obl.csv")
    kpp_desc, kpp_cols, kpp_rows = read_csv(D + "mom6_kpp.csv")
    sfc_desc, sfc_cols, sfc_rows = read_csv(D + "mom6_sfc_forcing.csv")

    # The traced domain is discovered from the data: PET 0's ocean cells
    # (the SW corner block of the tx2_3v2 grid has only 12 ocean cells -
    # the 77 southernmost rows are Antarctic land) and the number of traced
    # coupling steps (dump_first_steps=12).
    def key_cols(rows, cols):
        ix = {c: i for i, c in enumerate(cols)}
        return {(r[ix["cpl_step"]], r[ix["gi"]], r[ix["gj"]]) for r in rows}
    obl_keys = key_cols(obl_rows, obl_cols)
    NCELLS = len({(k[1], k[2]) for k in obl_keys})
    NSTEPS = len({k[0] for k in obl_keys})
    print(f"traced domain: {NCELLS} cells x {NSTEPS} steps (from CSV keys)")

    expect = {
        "mom6_eos.csv": NCELLS * KE * NSTEPS,
        "mom6_kpp_obl.csv": NCELLS * NSTEPS,
        "mom6_kpp.csv": NCELLS * NSTEPS,
        "mom6_sfc_forcing.csv": NCELLS * NSTEPS,
    }
    got = {
        "mom6_eos.csv": len(eos_rows),
        "mom6_kpp_obl.csv": len(obl_rows),
        "mom6_kpp.csv": len(kpp_rows),
        "mom6_sfc_forcing.csv": len(sfc_rows),
    }
    for f in expect:
        check(got[f] == expect[f], f"row count {got[f]} != {expect[f]}", failures)
        print(f"{f}: {got[f]} rows ({expect[f]} expected)")

    # ------------------------------------------------------------------
    # 2. EOS law
    # ------------------------------------------------------------------
    ix = {c: i for i, c in enumerate(eos_cols)}
    mx_rho = mx_dT = mx_dS = 0.0
    n_checked = 0
    for r in eos_rows:
        T, S, p, rho, dT, dS = (r[ix["T"]], r[ix["S"]], r[ix["p"]],
                                r[ix["rho"]], r[ix["drho_dT"]], r[ix["drho_dS"]])
        rho_c, dT_c, dS_c = wright_eos(T, S, p)
        mx_rho = max(mx_rho, rel_err(rho_c, rho))
        mx_dT = max(mx_dT, rel_err(dT_c, dT))
        mx_dS = max(mx_dS, rel_err(dS_c, dS))
        n_checked += 1
    check(mx_rho < args.eos_tol, f"EOS rho max rel err {mx_rho:.3e} >= {args.eos_tol}", failures)
    check(mx_dT < args.eos_tol, f"EOS drho/dT max rel err {mx_dT:.3e} >= {args.eos_tol}", failures)
    check(mx_dS < args.eos_tol, f"EOS drho/dS max rel err {mx_dS:.3e} >= {args.eos_tol}", failures)
    print(f"EOS law: {n_checked} rows; max rel err rho={mx_rho:.3e} drho/dT={mx_dT:.3e} drho/dS={mx_dS:.3e}")

    # Sanity: pressure increases downward within each (cell, step) column.
    prev_key, prev_p = None, None
    n_bad_p = 0
    for r in eos_rows:
        key = (r[ix["cpl_step"]], r[ix["gi"]], r[ix["gj"]])
        if key == prev_key:
            if not (r[ix["p"]] > prev_p):
                n_bad_p += 1
        prev_key, prev_p = key, r[ix["p"]]
    check(n_bad_p == 0, f"{n_bad_p} non-monotonic p columns", failures)

    # ------------------------------------------------------------------
    # 3. Hydrostatic pressure consistency: dp/dh ratio per column
    #    (Boussinesq pRef accumulation: dp = g_Earth*H_to_RZ * h)
    # ------------------------------------------------------------------
    ixo = {c: i for i, c in enumerate(obl_cols)}
    # h_01..h_75 and p_surf from the OBL dump; p_k from the EOS dump.
    h_names = [f"h_{k:02d}" for k in range(1, KE + 1)]
    eos_by_key = {}
    for re_ in eos_rows:
        eos_by_key.setdefault((re_[ix["cpl_step"]], re_[ix["gi"]], re_[ix["gj"]]), {})[
            int(re_[ix["k"]])] = re_[ix["p"]]
    ratios = []
    n_p1 = 0
    for r in obl_rows:
        key = (r[ixo["cpl_step"]], r[ixo["gi"]], r[ixo["gj"]])
        p_surf = r[ixo["p_surf"]]
        p_k = eos_by_key[key]
        if 1 in p_k and p_k[1] == p_surf:
            n_p1 += 1
        for k in range(1, KE):
            if k in p_k and k + 1 in p_k:
                dp = p_k[k + 1] - p_k[k]
                h = r[ixo[h_names[k - 1]]]
                ratios.append(dp / h)
    # The ratio should be g_Earth*H_to_RZ; check constancy, report the value.
    # (The accumulated pRef carries round-off, so a loose spread is expected
    # for deep thin layers; 1e-6 relative is still far tighter than any
    # physics.)
    check(n_p1 == NCELLS * NSTEPS, f"eos p(k=1) == p_surf: {n_p1}", failures)
    if ratios:
        rmin, rmax = min(ratios), max(ratios)
        spread = (rmax - rmin) / max(abs(rmax), 1e-300)
        check(spread < 1e-6, f"dp/dh ratio spread {spread:.3e}", failures)
        print(f"dp/dh ratio: {rmin:.9e} .. {rmax:.9e} "
              f"(g*H_to_RZ ~ {9.8 / RHO_0:.9e} if H_to_RZ = 1/Rho0)")

    # ------------------------------------------------------------------
    # 4. KPP boundary layer (mom6_kpp.csv)
    # ------------------------------------------------------------------
    ixk = {c: i for i, c in enumerate(kpp_cols)}
    def kcol(name, k):
        return ixk[f"{name}_{k:02d}"]

    n_kt = n_ks = n_kv = n_nlt_eq = n_nlt_sfc = n_nlt_stable = n_obl = n_kobl = 0
    for r in kpp_rows:
        OBL = r[ixk["OBLdepth"]]
        kOBL = r[ixk["kOBL"]]
        bfsfc = r[ixk["surfBuoyFlux"]]
        kbl = int(math.floor(kOBL))
        # OBL within the cell bracketing kOBL (CVMix kOBL = kbl - 1 + frac
        # with kbl the cell containing the OBL, so the bracketing interface
        # pair is (kbl, kbl+1) or (kbl-1, kbl) depending on convention;
        # accept either).  When the OBL reaches the bottom clamp, kOBL = ke+1
        # and the depth lies at or below the deepest interface.
        in_cell = False
        for kb in (kbl, kbl + 1):
            if 1 <= kb <= KE and r[kcol("z_inter", kb + 1)] <= -OBL <= r[kcol("z_inter", kb)]:
                in_cell = True
                break
        if not in_cell and kbl == NZ and -OBL <= r[kcol("z_inter", NZ)]:
            in_cell = True
        if in_cell:
            n_obl += 1
        if 0.0 <= kOBL - kbl < 1.0 and 1 <= kbl <= NZ:
            n_kobl += 1
        for k in range(1, NZ + 1):
            Kd_h = r[kcol("Kd_heat", k)]
            Kd_s = r[kcol("Kd_salt", k)]
            Kv = r[kcol("Kv", k)]
            Kt_in, Ks_in, Kv_in = (r[kcol("Kt_in", k)], r[kcol("Ks_in", k)],
                                   r[kcol("Kv_in", k)])
            Kt_out, Ks_out, Kv_out = (r[kcol("Kt_out", k)], r[kcol("Ks_out", k)],
                                      r[kcol("Kv_out", k)])
            if Kd_h != 0.0:
                if rel_err(Kt_out, Kd_h) < args.kpp_tol:
                    n_kt += 1
            else:
                if rel_err(Kt_out, Kt_in) < args.kpp_tol:
                    n_kt += 1
            if Kd_s != 0.0:
                if rel_err(Ks_out, Kd_s) < args.kpp_tol:
                    n_ks += 1
            else:
                if rel_err(Ks_out, Ks_in) < args.kpp_tol:
                    n_ks += 1
            if Kv != 0.0:
                if rel_err(Kv_out, Kv) < args.kpp_tol:
                    n_kv += 1
            else:
                if rel_err(Kv_out, Kv_in) < args.kpp_tol:
                    n_kv += 1
            # NLT equality and structure
            if r[kcol("NLT_heat", k)] == r[kcol("NLT_salt", k)]:
                n_nlt_eq += 1
            if k == 1 and r[kcol("NLT_heat", 1)] == 0.0:
                n_nlt_sfc += 1
            if bfsfc >= 0.0 and r[kcol("NLT_heat", k)] == 0.0:
                n_nlt_stable += 1
    tot_k = NCELLS * NSTEPS * NZ
    n_stable = sum(1 for r in kpp_rows if r[ixk["surfBuoyFlux"]] >= 0.0)
    check(n_kt == tot_k, f"Kt update rule: {n_kt}/{tot_k}", failures)
    check(n_ks == tot_k, f"Ks update rule: {n_ks}/{tot_k}", failures)
    check(n_kv == tot_k, f"Kv update rule: {n_kv}/{tot_k}", failures)
    check(n_nlt_eq == tot_k, f"NLT_heat==NLT_salt: {n_nlt_eq}/{tot_k}", failures)
    check(n_nlt_sfc == NCELLS * NSTEPS,
          f"NLT(k=1)=0: {n_nlt_sfc}/{NCELLS * NSTEPS}", failures)
    check(n_nlt_stable == NZ * n_stable,
          f"NLT=0 where surfBuoyFlux>=0: {n_nlt_stable}/{NZ * n_stable} "
          f"({n_stable} stable columns)", failures)
    check(n_obl == NCELLS * NSTEPS, f"OBL within kOBL cell: {n_obl}", failures)
    check(n_kobl == NCELLS * NSTEPS, f"kOBL fractional part: {n_kobl}", failures)
    print(f"KPP: Kt/Ks/Kv update rule {n_kt}/{n_ks}/{n_kv} of {tot_k} interfaces; "
          f"NLT structure {n_nlt_sfc}/{n_nlt_stable}; "
          f"OBL/kOBL consistency {n_obl}/{n_kobl} of {NCELLS * NSTEPS} columns")

    # ------------------------------------------------------------------
    # 5. OBL dump consistency with the KPP dump (same column geometry)
    # ------------------------------------------------------------------
    ixob = {c: i for i, c in enumerate(obl_cols)}
    kpp_z = {}
    for r in kpp_rows:
        key = (r[ixk["cpl_step"]], r[ixk["gi"]], r[ixk["gj"]])
        kpp_z[key] = [r[kcol("z_inter", k)] for k in range(1, NZ + 1)]
    n_z = 0
    for r in obl_rows:
        key = (r[ixob["cpl_step"]], r[ixob["gi"]], r[ixob["gj"]])
        if key in kpp_z:
            same = all(r[ixob[f"z_inter_{k:02d}"]] == kpp_z[key][k - 1]
                       for k in range(1, NZ + 1))
            if same:
                n_z += 1
    check(n_z == NCELLS * NSTEPS,
          f"z_inter identical obl vs kpp: {n_z}/{NCELLS * NSTEPS}", failures)
    n_obl_ok = sum(1 for r in obl_rows
                   if r[ixob["OBLdepth"]] > 0.0 and r[ixob["kOBL"]] >= 1.0)
    check(n_obl_ok == NCELLS * NSTEPS, f"OBL dump sanity: {n_obl_ok}", failures)
    print(f"OBL dump: z_inter matches kpp dump on {n_z} columns; "
          f"OBLdepth>0, kOBL>=1 on {n_obl_ok}")

    # ------------------------------------------------------------------
    # 6. Surface forcing mapping
    # ------------------------------------------------------------------
    ixs = {c: i for i, c in enumerate(sfc_cols)}
    n_lat = n_latd = n_sw = n_evap = n_salt = n_tau = n_ustar = n_omega = 0
    n_p = n_ice = n_liq = n_melt = 0
    for r in sfc_rows:
        qf, fp, fr = r[ixs["q_flux"]], r[ixs["fprec"]], r[ixs["frunoff"]]
        frg = r[ixs["frunoff_glc"]]
        lat = r[ixs["latent"]]
        if abs(lat - (qf * LV - fp * LF - fr * LF - frg * LF)) < 1e-6 * max(abs(lat), 1.0):
            n_lat += 1
        if abs(r[ixs["latent_evap_diag"]] - qf * LV) < 1e-6 * max(abs(qf * LV), 1.0):
            n_latd += 1
        if abs(r[ixs["latent_fprec_diag"]] + fp * LF) < 1e-6 * max(abs(fp * LF), 1.0):
            n_latd += 1
        if abs(r[ixs["latent_frunoff_diag"]] + fr * LF) < 1e-6 * max(abs(fr * LF), 1.0):
            n_latd += 1
        sw = r[ixs["sw"]]
        sw_sum = (r[ixs["sw_vis_dir"]] + r[ixs["sw_vis_dif"]] +
                  r[ixs["sw_nir_dir"]] + r[ixs["sw_nir_dif"]])
        if abs(sw - sw_sum) < 1e-6 * max(abs(sw), 1.0):
            n_sw += 1
        if abs(r[ixs["evap"]] - qf) < 1e-9 * max(abs(qf), 1.0):
            n_evap += 1
        if abs(r[ixs["salt_flux_o"]] - r[ixs["salt_flux"]]) < 1e-9 * max(abs(r[ixs["salt_flux"]]), 1.0):
            n_salt += 1
        if (abs(r[ixs["p_surf_full"]] - r[ixs["p"]]) < 1e-6 * max(abs(r[ixs["p"]]), 1.0) and
                abs(r[ixs["p_surf_full_o"]] - r[ixs["p"]]) < 1e-6 * max(abs(r[ixs["p"]]), 1.0)):
            n_p += 1
        if abs(r[ixs["ice_fraction_o"]] - r[ixs["ice_fraction"]]) < 1e-12:
            n_ice += 1
        if (abs(r[ixs["lprec_o"]] - r[ixs["lprec"]]) < 1e-9 * max(abs(r[ixs["lprec"]]), 1.0) and
                abs(r[ixs["fprec_o"]] - r[ixs["fprec"]]) < 1e-9 * max(abs(r[ixs["fprec"]]), 1.0) and
                abs(r[ixs["lrunoff_o"]] - r[ixs["lrunoff"]]) < 1e-9 * max(abs(r[ixs["lrunoff"]]), 1.0) and
                abs(r[ixs["frunoff_o"]] - r[ixs["frunoff"]]) < 1e-9 * max(abs(r[ixs["frunoff"]]), 1.0)):
            n_liq += 1
        if (abs(r[ixs["seaice_melt_o"]] - r[ixs["seaice_melt"]]) < 1e-9 * max(abs(r[ixs["seaice_melt"]]), 1.0) and
                abs(r[ixs["seaice_melt_heat_o"]] - r[ixs["seaice_melt_heat"]]) <
                1e-6 * max(abs(r[ixs["seaice_melt_heat"]]), 1.0)):
            n_melt += 1
        # AGRID: tau_mag/ustar/omega_w2x are built from the imported
        # u_flux/v_flux (mom_surface_forcing_nuopc.F90:956-959), not from
        # forces%taux/tauy (which live on the u grid as averages).
        uf, vf = r[ixs["u_flux"]], r[ixs["v_flux"]]
        tm = math.sqrt(uf * uf + vf * vf)
        if abs(r[ixs["tau_mag_o"]] - (GUST_CONST + tm)) < 1e-6 * max(abs(tm), 1.0):
            n_tau += 1
        us = math.sqrt((GUST_CONST + tm) / RHO_0)
        if abs(r[ixs["ustar_o"]] - us) < 1e-9 * max(abs(us), 1e-12):
            n_ustar += 1
        if abs(r[ixs["omega_w2x"]] - math.atan2(vf, uf)) < 1e-9:
            n_omega += 1
    nr = NCELLS * NSTEPS
    check(n_lat == nr, f"latent law: {n_lat}/{nr}", failures)
    check(n_latd == 3 * nr, f"latent diag laws: {n_latd}/{3 * nr}", failures)
    check(n_sw == nr, f"sw sum: {n_sw}/{nr}", failures)
    check(n_evap == nr, f"evap=q_flux: {n_evap}/{nr}", failures)
    check(n_salt == nr, f"salt_flux: {n_salt}/{nr}", failures)
    check(n_p == 2 * nr, f"p_surf_full=p: {n_p}/{2 * nr}", failures)
    check(n_ice == nr, f"ice_fraction: {n_ice}/{nr}", failures)
    check(n_liq == 4 * nr, f"lprec/fprec/lrunoff/frunoff: {n_liq}/{4 * nr}", failures)
    check(n_melt == 2 * nr, f"seaice_melt(+heat): {n_melt}/{2 * nr}", failures)
    check(n_tau == nr, f"tau_mag law: {n_tau}/{nr}", failures)
    check(n_ustar == nr, f"ustar law: {n_ustar}/{nr}", failures)
    check(n_omega == nr, f"omega_w2x law: {n_omega}/{nr}", failures)
    print(f"Forcing: latent {n_lat}/{nr}, latent diags {n_latd}/{3*nr}, sw {n_sw}/{nr}, "
          f"evap {n_evap}/{nr}, salt {n_salt}/{nr}, p {n_p}/{2*nr}, ice {n_ice}/{nr}, "
          f"liq {n_liq}/{4*nr}, melt {n_melt}/{2*nr}, tau {n_tau}/{nr}, ustar {n_ustar}/{nr}, "
          f"omega {n_omega}/{nr}")

    if failures:
        print("\nFAILURES:")
        for f_ in failures:
            print(" -", f_)
        sys.exit(1)
    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
