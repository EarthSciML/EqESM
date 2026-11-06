#!/usr/bin/env python3
"""Python reference for the MOM6 CVMix KPP column, transcribed from
MOM_CVMix_KPP.F90 (KPP_compute_BLD, KPP_calculate) and
pkg/CVMix-src/src/shared/cvmix_kpp.F90 (coeffs_kpp_low, OBL_depth_low,
kOBL_depth, bulk_Richardson, turbulent_scales_1d_OBL, unresolved_shear,
compute_phi_inv), verified against the stage-1 dumps mom6_kpp_obl.csv and
mom6_kpp.csv to ~1e-12 before transcription into stage2/mom6_kpp*.esm.

Conventions: all arrays 0-based in Python; the Fortran's 1-based layer k is
z_cell[k-1], interface K is z_inter[K-1] (z_inter[0] = 0 = surface).

Run config (C_JRA_t62_mom6_inst):
  RI_CRIT=0.3, VON_KARMAN=0.4, SURF_LAYER_EXTENT=0.1, MINIMUM_OBL_DEPTH=5.0,
  MINIMUM_VT2=1e-10, INTERP_TYPE=quadratic, INTERP_TYPE2=LMD94,
  MATCH_TECHNIQUE=MatchGradient, ENHANCE_DIFFUSION=False, NLT_SHAPE=CVMix,
  N_SMOOTH=3, COMPUTE_MONIN_OBUKHOV=True, COMPUTE_EKMAN=False, LT_K/LT_VT2 off,
  KPP_IS_ADDITIVE=False, KPP_SHORTWAVE_METHOD=MXL_SW, CVMix_ZERO_H_WORK_AROUND=0,
  all unit rescale 1.0.
"""
import os
import sys
import math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dumpcsv import load

# repo root = stage1/analysis/../..
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MAIN = os.environ.get("EQESM_MAIN", ROOT)

# ---------------- constants (run config + CVMix init defaults) ----------------
KAPPA = 0.4           # VON_KARMAN
RI_CRIT = 0.3
SURF_LAYER_EXT = 0.1
MIN_OBL_DEPTH = 5.0
MIN_VTSQR = 1e-10
ZETA_M = -0.2
ZETA_S = -1.0
GORHO = 9.8 / 1035.0  # g/rho0, Boussinesq
NLEV = 75

# derived phi-inv coefficients (cvmix_kpp.F90:369-376)
A_M = (1.0 - 16.0 * ZETA_M) ** (-0.25) * (1.0 - 4.0 * ZETA_M)
C_M = (1.0 - 16.0 * ZETA_M) ** (-0.25) * 12.0
A_S = math.sqrt(1.0 - 16.0 * ZETA_S) * (1.0 + 8.0 * ZETA_S)
C_S = 24.0 * math.sqrt(1.0 - 16.0 * ZETA_S)

# Wright EOS coefficients (MOM_EOS_Wright_full.F90:23-37)
W = dict(a0=7.133718e-4, a1=2.724670e-7, a2=-1.646582e-7,
         b0=5.613770e8, b1=3.600337e6, b2=-3.727194e4, b3=1.660557e2,
         b4=6.844158e5, b5=-8.389457e3,
         c0=1.609893e5, c1=8.427815e2, c2=-6.931554, c3=3.869318e-2,
         c4=-1.664201e2, c5=-2.765195)


def wright_rho(T, S, p):
    """rho = (p + p0)/(lambda + al0*(p+p0)), operation order as in the Fortran."""
    al0 = W["a0"] + (W["a1"] * T + W["a2"] * S)
    p0 = W["b0"] + (W["b4"] * S + T * (W["b1"] + (T * (W["b2"] + W["b3"] * T) + W["b5"] * S)))
    lam = W["c0"] + (W["c4"] * S + T * (W["c1"] + (T * (W["c2"] + W["c3"] * T) + W["c5"] * S)))
    num = p + p0
    return num / (lam + al0 * num)


def phi_inv(zeta, lm):
    """compute_phi_inv (cvmix_kpp.F90:2910-3005), not lStokesMOST."""
    if lm:
        if zeta >= 0.0:
            return 1.0 / (1.0 + 5.0 * zeta)
        if zeta >= ZETA_M:
            return (1.0 - 16.0 * zeta) ** 0.25
        return (A_M - C_M * zeta) ** (1.0 / 3.0)
    else:
        if zeta >= 0.0:
            return 1.0 / (1.0 + 5.0 * zeta)
        if zeta >= ZETA_S:
            return (1.0 - 16.0 * zeta) ** 0.5
        return (A_S - C_S * zeta) ** (1.0 / 3.0)


def turbulent_scale(sigma, OBL_depth, B0, ustar, lm):
    """One value of w_m or w_s (turbulent_scales_1d_OBL, not lStokesMOST,
    ustar != 0).  sigma_loc = min(surf_layer_ext, sigma) (l_LMD_ws=false)."""
    sigma_loc = min(SURF_LAYER_EXT, sigma)
    zeta = sigma_loc * OBL_depth * B0 * KAPPA / (ustar ** 3)
    return phi_inv(zeta, lm) * KAPPA * ustar


def poly_interp_quad(x1, y1, x2, y2, slope):
    """cvmix_math_poly_interp CVMIX_MATH_INTERP_QUAD (cvmix_math.F90:103-138)
    with rhs(3) = slope at x1.  Returns (a, b, c) for y = a + b x + c x^2,
    operation order as in the Fortran (Minv matrix application)."""
    det = -((x2 - x1) ** 2)
    rhs = (y2, y1, slope)
    Minv = np.zeros((3, 3))
    Minv[0, 0] = -1.0 / det; Minv[0, 1] = 1.0 / det; Minv[0, 2] = -1.0 / (x2 - x1)
    Minv[1, 0] = 2.0 * x1 / det; Minv[1, 1] = -2.0 * x1 / det; Minv[1, 2] = (x2 + x1) / (x2 - x1)
    Minv[2, 0] = -(x1 ** 2) / det; Minv[2, 1] = x2 * (2.0 * x1 - x2) / det; Minv[2, 2] = -x2 * x1 / (x2 - x1)
    coeffs = np.zeros(4)
    for k in range(3):
        for k2 in range(3):
            coeffs[k2] = coeffs[k2] + Minv[3 - k2 - 1, k] * rhs[k]
    return tuple(coeffs[:3])  # (a, b, c)


def newton_root(a, b, c, x0):
    """cvmix_math_cubic_root_find (cvmix_math.F90:260-294): Newton on
    f(x) = c x^2 + b x + a, TOL=1e-12 (absolute |f|), 100 iters."""
    root = x0
    f = c * (root ** 2) + b * root + a  # coeffs(4)=0, coeffs(3)=c, coeffs(2)=b, coeffs(1)=a
    for _ in range(100):
        if abs(f) < 1e-12:
            break
        slope = 2.0 * c * root + b
        root = root - f / slope
        f = c * (root ** 2) + b * root + a
    return root


def closed_form_root(a, b, c, x1, x2):
    """Closed-form root of a + b x + c x^2 = 0 that lies in [x1, x2]
    (negative depths), else nearest to the midpoint."""
    if c == 0.0:
        return -a / b if b != 0 else 0.5 * (x1 + x2)
    disc = b * b - 4.0 * a * c
    if disc < 0:
        return 0.5 * (x1 + x2)
    sq = math.sqrt(disc)
    if b >= 0:
        r1 = (-b - sq) / (2.0 * c)
        r2 = (2.0 * a) / (-b - sq)
    else:
        r1 = (-b + sq) / (2.0 * c)
        r2 = (2.0 * a) / (-b + sq)
    lo, hi = min(x1, x2), max(x1, x2)
    cand = [r for r in (r1, r2) if lo <= r <= hi]
    if cand:
        return min(cand, key=lambda r: abs(r - 0.5 * (x1 + x2)))
    return min((r1, r2), key=lambda r: abs(r - 0.5 * (x1 + x2)))


def compute_column(r, r2=None):
    """Full KPP column for one traced (cell, step).  Returns dict of arrays.
    0-based arrays; Fortran 1-based layer k <-> index k-1, interface K <-> K-1.
    r2 (optional) is the mom6_kpp.csv row, supplying the smoothed
    OBLdepth/kOBL used by the coefficient step."""
    T = np.array([r[f"T_{k:02d}"] for k in range(1, NLEV + 1)])
    S = np.array([r[f"S_{k:02d}"] for k in range(1, NLEV + 1)])
    UH = np.array([r[f"U_H_{k:02d}"] for k in range(1, NLEV + 1)])
    VH = np.array([r[f"V_H_{k:02d}"] for k in range(1, NLEV + 1)])
    h = np.array([r[f"h_{k:02d}"] for k in range(1, NLEV + 1)])
    dz = np.array([r[f"dz_{k:02d}"] for k in range(1, NLEV + 1)])
    bf = np.array([r[f"buoyFlux_{k:02d}"] for k in range(1, NLEV + 2)])  # 76
    ustar = r["uStar"]
    p_surf = r["p_surf"]

    # ---- geometry (min_thickness = 0: no inflation) ----
    z_inter = np.zeros(NLEV + 1)  # [0] = surface = 0
    z_cell = np.zeros(NLEV)
    for k in range(NLEV):
        z_inter[k + 1] = z_inter[k] - dz[k]
        z_cell[k] = z_inter[k] - 0.5 * dz[k]

    # ---- surface-layer averages, per layer k (F90:1342-1405) ----
    # SLdepth_0d(k) = 0.1 * max(max(-z_cell(k), -z_inter(2)), 5)
    SLdepth = 0.1 * np.maximum(np.maximum(-z_cell, -z_inter[1]), MIN_OBL_DEPTH)
    # delH(k,kp) = min(max(SLdepth(k) + z_inter(kp), 0), dz(kp)); z_inter(kp)
    # is the TOP interface of 1-based layer kp (0-based index kp-1)
    delH = np.minimum(np.maximum(SLdepth[:, None] + z_inter[:NLEV][None, :], 0.0), dz[None, :])
    hTot = delH.sum(axis=1)
    surfT = (delH @ T) / hTot
    surfS = (delH @ S) / hTot
    surfU = (delH @ UH) / hTot
    surfV = (delH @ VH) / hTot

    # ---- shear and buoyancy (F90:1407-1457) ----
    Uk = UH - surfU
    Vk = VH - surfV
    deltaU2 = Uk ** 2 + Vk ** 2

    # EOS pressure for layer k: pRef_k = p_surf + sum_{k'<k} g*rho0*h(k')
    # (F90:1437-1440; g_Earth*H_to_RZ = 9.8*1035 in MKS), accumulated as a
    # RUNNING sum exactly like the Fortran: pRef_{k+1} = pRef_k + c*h_k.
    G_RHO0 = 9.8 * 1035.0
    pRef = np.zeros(NLEV)
    pRef[0] = p_surf
    for k in range(1, NLEV):
        pRef[k] = pRef[k - 1] + G_RHO0 * h[k - 1]
    rho_surf = np.array([wright_rho(surfT[k], surfS[k], pRef[k]) for k in range(NLEV)])
    rho_k = np.array([wright_rho(T[k], S[k], pRef[k]) for k in range(NLEV)])
    rho_km1 = np.array([wright_rho(T[max(k - 1, 0)], S[max(k - 1, 0)], pRef[k]) for k in range(NLEV)])
    deltaBuoy = GORHO * (rho_k - rho_surf)
    N2 = np.zeros(NLEV + 1)
    for k in range(NLEV):
        N2[k] = GORHO * (rho_k[k] - rho_km1[k]) / (0.5 * (h[max(k - 1, 0)] + h[k]))
    N2[NLEV] = 0.0
    N = np.sqrt(np.maximum(N2, 0.0))

    # ---- w_s, unresolved shear, bulk Ri (F90:1509-1550) ----
    B0 = bf[0] - bf[1:]  # B0(k) = buoyFlux(1) - buoyFlux(k+1), 0-based k
    ws = np.array([turbulent_scale(1.0, -z_cell[k], B0[k], ustar, lm=False) for k in range(NLEV)])
    Vtc = math.sqrt(0.2 / (C_S * SURF_LAYER_EXT)) / (KAPPA ** 2)
    Vt2 = np.zeros(NLEV)
    for k in range(NLEV):
        Nc = N[k + 1]  # N_cntr(kt) = N_iface(kt+1), 1-based kt -> 0-based index kt = k+1
        Cv = 2.1 - 200.0 * Nc if Nc < 0.002 else 1.7
        vt = -Cv * Vtc * z_cell[k] * Nc * ws[k] / RI_CRIT
        Vt2[k] = max(vt, MIN_VTSQR)
    scaling = 1.0 - 0.5 * SURF_LAYER_EXT
    BulkRi = scaling * (-z_cell) * deltaBuoy / (deltaU2 + Vt2)

    # ---- OBL depth (cvmix_kpp.F90:1585-1814, not lStokesMOST) ----
    OBL_limit = abs(z_cell[NLEV - 1])
    if B0[NLEV - 1] > 0.0:  # lstable
        MoninObukhov = ustar ** 3 / (B0[NLEV - 1] * KAPPA)
        OBL_limit = min(OBL_limit, MoninObukhov)
    k = None
    for kk in range(NLEV):  # 0-based: first BulkRi(kk) > RI_CRIT (F90 k+1 1-based)
        if BulkRi[kk] > RI_CRIT:
            k = kk
            break
    if k is None:
        OBL_depth = OBL_limit
    elif k == 0:
        OBL_depth = abs(z_cell[0])
    else:
        x1, x2 = z_cell[k - 1], z_cell[k]
        y1, y2 = BulkRi[k - 1], BulkRi[k]
        slope = 0.0 if k == 1 else (y1 - BulkRi[k - 2]) / (x1 - z_cell[k - 2])
        a, b, c = poly_interp_quad(x1, y1, x2, y2, slope)
        a = a - RI_CRIT
        root_n = newton_root(a, b, c, 0.5 * (x1 + x2))
        root_c = closed_form_root(a, b, c, x1, x2)
        OBL_depth = -root_c  # closed form is what the .esm spells
        if k == 1:
            OBL_depth = max(OBL_depth, -z_cell[0])
        OBL_depth = min(OBL_depth, OBL_limit)
    OBL_depth = max(OBL_depth, MIN_OBL_DEPTH)
    OBL_depth = max(OBL_depth, -z_inter[1])     # MOM6: no shallower than top layer
    OBL_depth = min(OBL_depth, -z_inter[NLEV])  # MOM6: no deeper than bottom iface
    # kOBL (cvmix_kpp.F90:1821-1867): kw=1..75, compare with iface kw+1 and zt(kw)
    kOBL = float(NLEV) + 0.75
    for kw in range(1, NLEV + 1):  # 1-based
        if OBL_depth < abs(z_inter[kw]):  # iface kw+1 (0-based index kw)
            if OBL_depth < abs(z_cell[kw - 1]):
                kOBL = float(kw) + 0.25
            else:
                kOBL = float(kw) + 0.75
            break

    # ---- coefficients (F90:904-931; cvmix coeffs_kpp_low) ----
    # NOTE: uses the SMOOTHED OBL_depth/kOBL (mom6_kpp.csv), not the
    # pre-smoothing values; passed in by the caller for the reference check.
    OBL_depth_s = r2["OBLdepth"]
    kOBL_s = r2["kOBL"]
    B0_MXL = bf[0] - bf[int(kOBL_s)]  # buoyFlux(int(kOBL)+1), 1-based -> 0-based int(kOBL)
    kwup = int(math.floor(kOBL_s))
    ktup = int(round(kOBL_s)) - 1
    sigma = -z_inter / OBL_depth_s  # interfaces 1..76, sigma[K-1] for iface K
    Kd = np.zeros(NLEV + 1)
    Ks = np.zeros(NLEV + 1)
    Kv = np.zeros(NLEV + 1)
    for kw in range(2, kwup + 1):  # 1-based interfaces 2..kwup (kw loop to kwup)
        s = sigma[kw - 1]
        # G via the shape coefficients a=(0,1,-2,1), cvmix_math_evaluate_cubic
        # operation order: (a0 + a1*x) + a2*x^2 + a3*x^3 (cvmix_math.F90:325-332)
        # with x^2 = x*x, x^3 = x*x*x (gfortran expands integer powers)
        G = (s + (-2.0) * (s * s)) + (s * s * s)
        ws_kw = turbulent_scale(s, OBL_depth_s, B0_MXL, ustar, lm=False)
        wm_kw = turbulent_scale(s, OBL_depth_s, B0_MXL, ustar, lm=True)
        Kd[kw - 1] = OBL_depth_s * ws_kw * G
        Ks[kw - 1] = OBL_depth_s * ws_kw * G
        Kv[kw - 1] = OBL_depth_s * wm_kw * G
    # out(2:ktup+1) = OBL(2:ktup+1); ENHANCE_DIFFUSION=False, MatchBoth=False
    # (old values zero): OBL(ktup+1) = 0 when kwup == ktup (kOBL = X.75); when
    # ktup == kwup-1 (kOBL = X.25) the kw loop already computed OBL(kwup).
    # stable columns -> NLT = 0.

    return dict(z_cell=z_cell, z_inter=z_inter, deltaU2=deltaU2,
                deltaBuoy=deltaBuoy, N2=N2, N=N, ws=ws, Vt2=Vt2, BulkRi=BulkRi,
                OBLdepth=OBL_depth, kOBL=kOBL, Kd=Kd, Ks=Ks, Kv=Kv,
                B0_MXL=B0_MXL, root_n=root_n, root_c=root_c, kfirst=k)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-300)


def main(dumps_dir=None):
    if dumps_dir is None:
        dumps_dir = f"{MAIN}/stage1/dumps/mom6"
    meta, rows = load(f"{dumps_dir}/mom6_kpp_obl.csv")
    _, rows2 = load(f"{dumps_dir}/mom6_kpp.csv")

    worst = {}
    worst_obl = 0.0
    worst_kobl = 0
    worst_rootdiff = 0.0
    kfirst_hist = {}
    where = {}
    for r, r2 in zip(rows, rows2):
        c = compute_column(r, r2)
        for k in range(1, NLEV + 1):
            for col, name in [("BulkRi", "BulkRi"), ("deltaU2", "deltaU2"), ("N2", "N2")]:
                got = c[col][k - 1]
                want = r[f"{name}_{k:02d}"]
                d = rel(got, want)
                if d > worst.get(col, 0.0):
                    worst[col] = d
                    where[col] = (r["i"], r["j"], r["tod"], k, got, want)
            for col in ["z_cell", "z_inter"]:
                worst[col] = max(worst.get(col, 0.0), rel(c[col][k - 1], r[f"{col}_{k:02d}"]))
        worst["OBLdepth"] = max(worst.get("OBLdepth", 0.0), rel(c["OBLdepth"], r["OBLdepth"]))
        worst_obl = max(worst_obl, abs(c["OBLdepth"] - r["OBLdepth"]))
        if c["kOBL"] != r["kOBL"]:
            worst_kobl += 1
        worst_rootdiff = max(worst_rootdiff, rel(c["root_n"], c["root_c"]))
        kfirst_hist[c["kfirst"]] = kfirst_hist.get(c["kfirst"], 0) + 1
        for k in range(1, NLEV + 1):
            for col, name in [("Kd", "Kd_heat"), ("Ks", "Kd_salt"), ("Kv", "Kv")]:
                got = c[col][k - 1]
                want = r2[f"{name}_{k:02d}"]
                d = rel(got, want)
                if d > worst.get(col, 0.0):
                    worst[col] = d
                    where[col] = (r["i"], r["j"], r["tod"], k, got, want)
        d = rel(c["B0_MXL"], r2["surfBuoyFlux"])
        if d > worst.get("surfBuoyFlux", 0.0):
            worst["surfBuoyFlux"] = d
            where["surfBuoyFlux"] = (r["i"], r["j"], r["tod"], 0, c["B0_MXL"], r2["surfBuoyFlux"])

    print("rows:", len(rows))
    print("first-crossing k histogram (0-based):", dict(sorted(kfirst_hist.items())))
    print("worst rel diffs (vs obl dump):")
    for col in ["BulkRi", "deltaU2", "N2", "z_cell", "z_inter", "OBLdepth"]:
        print(f"  {col:10s} {worst.get(col, 0):.3e}  at {where.get(col)}")
    print("worst abs OBLdepth diff:", worst_obl)
    print("kOBL mismatches (pre-smoothing):", worst_kobl)
    print("worst |root_newton - root_closed| rel:", worst_rootdiff)
    print("worst rel diffs (vs kpp dump):")
    for col in ["Kd", "Ks", "Kv", "surfBuoyFlux"]:
        print(f"  {col:10s} {worst.get(col, 0):.3e}  at {where.get(col)}")

    nlt_nonzero = 0
    for r2 in rows2:
        for k in range(1, NLEV + 1):
            if r2[f"NLT_heat_{k:02d}"] != 0.0 or r2[f"NLT_salt_{k:02d}"] != 0.0:
                nlt_nonzero += 1
    print("NLT nonzero entries:", nlt_nonzero)


if __name__ == "__main__":
    main(*sys.argv[1:2])
