#!/usr/bin/env python3
"""Identity checks on cice_ridge.csv and cice_itd.csv (CICE ridge + ITD dumps).

Each row is one (cell, step) trace from the D_cice_inst run (f19 grid,
32 traced cells x 12 steps).  The dumps record the icepack ridge_ice /
ridge_shift and linear_itd internals; the checks below recompute every
recorded quantity from the other recorded quantities using the exact
Fortran formulas (icepack_mechred.F90 and icepack_therm_itd.F90), so a
pass pins the dumps to the code.  Config constants are the D-class
namelist values: krdg_partic = 1, krdg_redist = 1, mu_rdg = 4,
Cs = 0.25, astar = 0.05, maxraft = 1, fsnowrdg = 0.5, hi_min = 0.01
(icepack_parameters p01), dt = 3600 s, puny = 1e-11.

cice_ridge.csv checks (per row):
  R1  sum(apartic_0..5) = 1 (normalized participation function)
  R2  aksum = apartic_0 + sum_n apartic_n*(1 - 1/krdg_n)      (ridge_itd)
  R3  closing_net/opning/divu_adv = ridge_prep outputs:
        divu_adv = (1 - asum_entry)/dt, asum_entry = aice0_init + sum aicen_init
        closing_net = Cs*rdg_shear + rdg_conv, raised to max(., -divu_adv)
                      if divu_adv < 0
        opning = closing_net + divu_adv
  R4  hrmin_n = min(2*hi, hi + maxraft), hrexp_n = mu_rdg*sqrt(hi),
      krdg_n = (hrmin_n + hrexp_n)/hi, hi = max(vicen_init_n/aicen_init_n,
      puny) for aicen_init_n > puny (else krdg = 1, hrmin = hrexp = 0)
  R5  transfers: ardg2n_n = ardg1n_n/krdg_n, virdgn_n =
      vicen_init_n*ardg1n_n/aicen_init_n (where ardg1n_n > 0)
  R6  totals: ardg1 = sum ardg1n, ardg2 = sum ardg2n, virdg = sum virdgn
      (niter = 1 in this run, so no iteration accumulation)
  R7  aredistn_nr = farea_nr*ardg2n_1, vredistn_nr = fvol_nr*virdgn_1
      (only when category 1 ridged; farea/fvol = exponential ridge ITD
      integrals over [max(hrmin_1, hin_max(nr-1)), hin_max(nr)], tail
      category expL/(hi1+hexp) forms; sum farea = 1)
  R8  rafting: araftn_n = mraftn_n*ardg2n_n, vraftn_n = mraftn_n*virdgn_n
      with the ridge_itd rafting mask mraftn_n = 1 iff
      |(2*hi + hrexp_n)/hi - krdg_n| < puny
  R9  conservation: sum eicen_final = sum eicen_init, sum sicen/vbrin
      likewise, sum esno_final + esnow_mlt = sum esno_init, sum vicen
      final = init, sum vsnon_final + msnow_mlt/rhos = sum vsnon_init
      (msnow_mlt < 0 snow thrown to ocean is added back)
  R10 aice0_final = 1 - sum aicen_final (asum = aice0 + sum aicen = 1
      at convergence, within puny), asum ~ 1

cice_itd.csv checks (per row):
  I1  hicen_init_n = vicen_init_n/aicen_init_n, hicen_n = vicen_n/aicen_n,
      dhicen_n = hicen_n - hicen_init_n (zero for empty categories)
  I2  hbnew_0 = hin_max_0; hbnew_n (n = 1..4) by the linear interpolation
      between adjacent growth rates (or the single-category/empty
      branches); hbnew_5 = max(3*hicen_5 - 2*hbnew_4, hin_max_4)
  I3  fit_line: g0_n/g1_n/hL_n/hR_n satisfy
        int_{hL}^{hR} g dh = aicen_fit_n      (area)
        int_{hL}^{hR} h g dh = aicen_fit_n*hice_fit_n   (first moment)
      with g(h) = g0 + g1*(h - hL), aicen_fit_1 = aicen_1 - da0_melt,
      hice_fit_1 = vicen_1/aicen_fit_1 (the dh0<0 melt ran first),
      aicen_fit_n = aicen_n, hice_fit_n = hicen_n for n > 1; hL/hR are
      the fit bounds or the h13/h23-shifted values
      (hL = 3*hice - 2*hR or hR = 3*hice - 2*hL)
  I4  boundary transfers: donor_n/daice_n/dvice_n from the integration
      of the donor category's g over the [etamin, etamax] eta range
      (etamin = max(hin_max_n, hL)-hL, etamax = min(hbnew_n, hR)-hL for
      donor n; etamin = 0, etamax = min(hin_max_n, hR_{n+1})-hL_{n+1}
      for donor n+1), then the puny small-transfer and full-category
      clamps against aicen_donor/vicen_donor
      (aicen_1 reduced by da0_melt when the donor is category 1)
  I5  post state: aicen_post_n/vicen_post_n = entry state minus
      da0_melt (n=1) plus the boundary transfers, then the hi_min
      closure on category 1 (da0_himin = aicen*(1 - hicen/hi_min) when
      aicen_post_1 > puny and hicen_1 < hi_min = 0.01); aice_post =
      sum aicen_post, aice0_post = 1 - aice_post
  I6  exit state: aice_exit = sum aicen_exit, aice0_exit = 1 - aice_exit
      (cleanup_itd may zap tiny categories to the ocean, so no
      post/exit equality is asserted)
  I7  conservation totals (linear_itd's own thresholds): |vice_init -
      vice_final| <= puny, |vsno| <= puny, |eice| <= puny*Lfresh*rhoi,
      |esno| <= puny*Lfresh*rhos, |sice| <= puny, |vbri| <= puny*10

Tolerances: rtol = atol = 1e-12 for the algebraic identities (the same
operation order as the Fortran), puny-scaled for the conservation
totals.  Coverage info (niter distribution, remap_flag, nonzero
transfers) is reported at the end.
"""
import csv
import math
import re
import sys

DENORM = re.compile(r"([0-9])([+-])([0-9]{3})$")

def tofloat(s):
    s = s.strip()
    m = DENORM.search(s)
    if m:
        s = s[: m.start(1) + 1] + "E" + m.group(2) + m.group(3)
    return float(s)

PATH_R = sys.argv[1] if len(sys.argv) > 1 else \
    "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cice/cice_ridge.csv"
PATH_I = sys.argv[2] if len(sys.argv) > 2 else \
    "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/cice/cice_itd.csv"

RTOL, ATOL = 1e-12, 1e-12

# D-class configuration
PUNY = 1e-11
CS = 0.25          # ridge_prep: closing_net = Cs*rdg_shear + rdg_conv
ASTAR = 0.05       # ridge_itd participation e-folding
ASTARI = 1.0 / ASTAR
MU_RDG = 4.0       # ridge e-folding scale
MAXRAFT = 1.0      # max rafting thickness
FSNOWRDG = 0.5     # fraction of snow kept on ridges
RHOS = 330.0       # snow density
RHOI = 917.0       # ice density
LFRESH = 3.34e5    # latent heat of fusion
HI_MIN = 0.01      # min ice thickness (icepack p01)
DT = 3600.0        # model timestep
NCAT = 5

def close(a, b):
    return abs(a - b) <= RTOL * abs(b) + ATOL

def load(path):
    with open(path) as f:
        f.readline()  # # description
        cols = f.readline().strip().split(" = ", 1)[1].split(",")
        return list(csv.DictReader(f, fieldnames=cols))

def main():
    r = load(PATH_R)
    i = load(PATH_I)
    if len(r) != 384 or len(i) != 384:
        print(f"FAIL: expected 384 rows each, got {len(r)} ridge, {len(i)} itd")
        sys.exit(1)

    fails = []
    ncheck = [0]

    def check(cond, label, ctx=""):
        ncheck[0] += 1
        if not cond:
            fails.append(label + (f" [{ctx}]" if ctx else ""))

    def ctx(x):
        return f"step {x['step'].strip()} cell ({x['i'].strip()},{x['j'].strip()})"

    # ---------------------------------------------------------------
    # cice_ridge.csv
    # ---------------------------------------------------------------
    niter_dist = {}
    for x in r:
        c = ctx(x)
        ap = [tofloat(x[f"apartic_{n}"]) for n in range(NCAT + 1)]
        krdg = [tofloat(x[f"krdg_{n}"]) for n in range(1, NCAT + 1)]
        hrmin = [tofloat(x[f"hrmin_{n}"]) for n in range(1, NCAT + 1)]
        hrexp = [tofloat(x[f"hrexp_{n}"]) for n in range(1, NCAT + 1)]
        ai = [tofloat(x[f"aicen_init_{n}"]) for n in range(1, NCAT + 1)]
        vi = [tofloat(x[f"vicen_init_{n}"]) for n in range(1, NCAT + 1)]
        a1n = [tofloat(x[f"ardg1n_{n}"]) for n in range(1, NCAT + 1)]
        a2n = [tofloat(x[f"ardg2n_{n}"]) for n in range(1, NCAT + 1)]
        vdn = [tofloat(x[f"virdgn_{n}"]) for n in range(1, NCAT + 1)]
        ardn = [tofloat(x[f"aredistn_{n}"]) for n in range(1, NCAT + 1)]
        vrdn = [tofloat(x[f"vredistn_{n}"]) for n in range(1, NCAT + 1)]
        araft = [tofloat(x[f"araftn_{n}"]) for n in range(1, NCAT + 1)]
        vraft = [tofloat(x[f"vraftn_{n}"]) for n in range(1, NCAT + 1)]
        aice0i = tofloat(x["aice0_init"])
        aice0f = tofloat(x["aice0_final"])
        niter = int(round(tofloat(x["niter"])))
        niter_dist[niter] = niter_dist.get(niter, 0) + 1

        # R1: participation sums to 1
        check(abs(sum(ap) - 1.0) <= 1e-12, "R1 sum apartic", c)

        # R2: aksum identity (ap[0] is the open-water term, added once)
        aks = tofloat(x["aksum"])
        akc = ap[0] + sum(ap[n + 1] * (1.0 - 1.0 / krdg[n]) for n in range(NCAT))
        check(close(aks, akc), "R2 aksum", c)

        # R3: ridge_prep rates
        asum_entry = aice0i + sum(ai)
        divu = (1.0 - asum_entry) / DT
        cl = CS * tofloat(x["rdg_shear"]) + tofloat(x["rdg_conv"])
        if divu < 0.0:
            cl = max(cl, -divu)
        check(close(tofloat(x["closing_net"]), cl), "R3 closing_net", c)
        check(close(tofloat(x["divu_adv"]), divu), "R3 divu_adv", c)
        check(close(tofloat(x["opning"]), cl + divu), "R3 opning", c)

        # R4: ridge ITD parameters from entry state
        for n in range(NCAT):
            if ai[n] > PUNY:
                hi = max(vi[n] / ai[n], PUNY)
                check(close(hrmin[n], min(2.0 * hi, hi + MAXRAFT)),
                      f"R4 hrmin_{n+1}", c)
                check(close(hrexp[n], MU_RDG * math.sqrt(hi)),
                      f"R4 hrexp_{n+1}", c)
                check(close(krdg[n], (hrmin[n] + hrexp[n]) / hi),
                      f"R4 krdg_{n+1}", c)
            else:
                check(krdg[n] == 1.0 and hrmin[n] == 0.0 and hrexp[n] == 0.0,
                      f"R4 empty cat {n+1}", c)

        # R5/R6: transfers and totals
        check(close(tofloat(x["ardg1"]), sum(a1n)), "R6 ardg1", c)
        check(close(tofloat(x["ardg2"]), sum(a2n)), "R6 ardg2", c)
        check(close(tofloat(x["virdg"]), sum(vdn)), "R6 virdg", c)
        for n in range(NCAT):
            if a1n[n] > 0.0:
                check(close(a2n[n], a1n[n] / krdg[n]), f"R5 ardg2n_{n+1}", c)
                check(close(vdn[n], vi[n] * a1n[n] / ai[n]),
                      f"R5 virdgn_{n+1}", c)
            else:
                check(a2n[n] == 0.0 and vdn[n] == 0.0, f"R5 no-transfer {n+1}", c)

        # R7: n = 1 redistribution fractions
        hi1, hexp = hrmin[0], hrexp[0]
        if a1n[0] > 0.0:
            hin = [tofloat(x[f"hin_max_{n}"]) for n in range(NCAT)]
            farea_sum, fvol_sum = 0.0, 0.0
            for nr in range(NCAT):
                if nr < NCAT - 1:
                    if hi1 >= hin[nr + 1]:
                        fa = fv = 0.0  # ridge ITD starts above this bin
                    else:
                        hL = max(hi1, hin[nr])
                        hR = hin[nr + 1]
                        eL = math.exp(-min((hL - hi1) / hexp, 100.0))
                        eR = math.exp(-min((hR - hi1) / hexp, 100.0))
                        fa = eL - eR
                        fv = ((hL + hexp) * eL - (hR + hexp) * eR) / (hi1 + hexp)
                else:
                    hL = max(hi1, hin[NCAT - 1])
                    eL = math.exp(-min((hL - hi1) / hexp, 100.0))
                    fa = eL
                    fv = (hL + hexp) * eL / (hi1 + hexp)
                farea_sum += fa
                fvol_sum += fv
                check(close(ardn[nr], fa * a2n[0]), f"R7 aredistn_{nr+1}", c)
                check(close(vrdn[nr], fv * vdn[0]), f"R7 vredistn_{nr+1}", c)
            check(close(farea_sum, 1.0) and close(fvol_sum, 1.0),
                  "R7 farea/fvol sum", c)

        # R8: rafting mask
        for n in range(NCAT):
            if ai[n] > PUNY:
                hi = max(vi[n] / ai[n], PUNY)
                mraft = 1.0 if abs((2.0 * hi + hrexp[n]) / hi - krdg[n]) < PUNY \
                    else 0.0
            else:
                mraft = 0.0
            check(close(araft[n], mraft * a2n[n]), f"R8 araftn_{n+1}", c)
            check(close(vraft[n], mraft * vdn[n]), f"R8 vraftn_{n+1}", c)

        # R9: conservation (energies dumped as per-category sums)
        def colsum(x, base, suffix):
            # columns are e.g. eicen_init_1 (init) and eicen_1 (final)
            return sum(tofloat(x[f"{base}{suffix}_{n}"]) for n in range(1, NCAT + 1))
        for base, label in (("eicen", "eice"), ("sicen", "sice"),
                            ("vbrin", "vbri")):
            check(close(colsum(x, base, "_init"), colsum(x, base, "")),
                  f"R9 {label}", c)
        check(close(colsum(x, "esno", "_init"),
                    colsum(x, "esno", "") + tofloat(x["esnow_mlt"])),
              "R9 esno", c)
        check(close(sum(vi), sum(tofloat(x[f"vicen_{n}"]) for n in range(1, NCAT + 1))),
              "R9 vice", c)
        vsi = colsum(x, "vsnon", "_init")
        vsf = colsum(x, "vsnon", "")
        check(close(vsf + tofloat(x["msnow_mlt"]) / RHOS, vsi), "R9 vsno", c)

        # R10: final area balance
        asum = aice0f + sum(tofloat(x[f"aicen_{n}"]) for n in range(1, NCAT + 1))
        check(abs(asum - 1.0) <= 1e-10, "R10 asum", c)
        check(close(aice0f, 1.0 - sum(tofloat(x[f"aicen_{n}"])
                                      for n in range(1, NCAT + 1))),
              "R10 aice0", c)
        check(tofloat(x["msnow_mlt"]) >= 0.0, "R10 msnow_mlt >= 0", c)
        check(tofloat(x["esnow_mlt"]) <= 0.0, "R10 esnow_mlt <= 0", c)

    # ---------------------------------------------------------------
    # cice_itd.csv
    # ---------------------------------------------------------------
    remap_dist = {}
    for x in i:
        c = ctx(x)
        hin = [tofloat(x[f"hin_max_{n}"]) for n in range(NCAT + 1)]
        aci = [tofloat(x[f"aicen_init_{n}"]) for n in range(1, NCAT + 1)]
        vci = [tofloat(x[f"vicen_init_{n}"]) for n in range(1, NCAT + 1)]
        ac = [tofloat(x[f"aicen_{n}"]) for n in range(1, NCAT + 1)]
        vc = [tofloat(x[f"vicen_{n}"]) for n in range(1, NCAT + 1)]
        vsc = [tofloat(x[f"vsnon_{n}"]) for n in range(1, NCAT + 1)]
        hci = [tofloat(x[f"hicen_init_{n}"]) for n in range(1, NCAT + 1)]
        hc = [tofloat(x[f"hicen_{n}"]) for n in range(1, NCAT + 1)]
        dhc = [tofloat(x[f"dhicen_{n}"]) for n in range(1, NCAT + 1)]
        hb = [tofloat(x[f"hbnew_{n}"]) for n in range(NCAT + 1)]
        g0 = [tofloat(x[f"g0_{n}"]) for n in range(1, NCAT + 1)]
        g1 = [tofloat(x[f"g1_{n}"]) for n in range(1, NCAT + 1)]
        hL = [tofloat(x[f"hL_{n}"]) for n in range(1, NCAT + 1)]
        hR = [tofloat(x[f"hR_{n}"]) for n in range(1, NCAT + 1)]
        donor = [int(round(tofloat(x[f"donor_{n}"]))) for n in range(1, NCAT)]
        daice = [tofloat(x[f"daice_{n}"]) for n in range(1, NCAT)]
        dvice = [tofloat(x[f"dvice_{n}"]) for n in range(1, NCAT)]
        da0m = tofloat(x["da0_melt"])
        da0h = tofloat(x["da0_himin"])
        remap = int(round(tofloat(x["remap_flag"])))
        remap_dist[remap] = remap_dist.get(remap, 0) + 1

        # I1: initial thicknesses and thickness changes (dhicen is
        # computed from the pre-shift entry thicknesses)
        for n in range(NCAT):
            if aci[n] > PUNY:
                check(close(hci[n], vci[n] / aci[n]), f"I1 hicen_init_{n+1}", c)
            else:
                check(hci[n] == 0.0, f"I1 hicen_init empty {n+1}", c)
            if ac[n] > PUNY:
                check(close(dhc[n], vc[n] / ac[n] - hci[n]),
                      f"I1 dhicen_{n+1}", c)
            else:
                check(dhc[n] == 0.0, f"I1 dhicen empty {n+1}", c)

        # I2: new boundaries (computed from the pre-shift entry hicen)
        if dhc[0] < 0.0:  # melting: left bound stays at hin_max(0)
            check(close(hb[0], hin[0]), "I2 hbnew_0", c)
        else:             # freezing: hbnew(0) = min(dh0, hin_max(1))
            check(close(hb[0], min(dhc[0], hin[1])), "I2 hbnew_0 freeze", c)
        for n in range(NCAT - 1):
            if hci[n] > PUNY and hci[n + 1] > PUNY:
                slope = (dhc[n + 1] - dhc[n]) / (hci[n + 1] - hci[n])
                exp_hb = hin[n + 1] + dhc[n] + slope * (hin[n + 1] - hci[n])
            elif hci[n] > PUNY:
                exp_hb = hin[n + 1] + dhc[n]
            elif hci[n + 1] > PUNY:
                exp_hb = hin[n + 1] + dhc[n + 1]
            else:
                exp_hb = hin[n + 1]
            check(close(hb[n + 1], exp_hb), f"I2 hbnew_{n+1}", c)
        if ac[NCAT - 1] > PUNY:
            check(close(hb[NCAT], max(3.0 * vc[NCAT - 1] / ac[NCAT - 1]
                                      - 2.0 * hb[NCAT - 1],
                                      hin[NCAT - 1])),
                  "I2 hbnew_5", c)
        else:
            check(close(hb[NCAT], hin[NCAT]), "I2 hbnew_5 empty", c)

        # entry state after the dh0<0 melt (aicen_1 reduced, volume kept)
        acm = list(ac)
        acm[0] -= da0m
        vcm = list(vc)

        # I3: fit_line identities
        for n in range(NCAT):
            afit = acm[n]
            hfit = vcm[n] / afit if afit > PUNY else 0.0
            hbL, hbR = hb[n], hb[n + 1]
            if afit > PUNY and hbR - hbL > PUNY:
                dhr = 1.0 / (hR[n] - hL[n])
                wk2 = (hfit - hL[n]) * dhr
                area = g0[n] * (hR[n] - hL[n]) + g1[n] * (hR[n] - hL[n]) ** 2 / 2.0
                moment = (g1[n] * (hR[n] - hL[n]) ** 3 / 3.0
                          + g0[n] * (hR[n] - hL[n]) ** 2 / 2.0
                          + hL[n] * area)
                check(close(area, afit), f"I3 area_{n+1}", c)
                check(close(moment, afit * hfit), f"I3 moment_{n+1}", c)
                # hL/hR are the input bounds or the h13/h23-shifted values
                # (the shift puts hice in the central third of the FINAL
                # range, so test against the final hL/hR)
                hLok = close(hL[n], hbL) or close(hL[n], 3.0 * hfit - 2.0 * hR[n])
                hRok = close(hR[n], hbR) or close(hR[n], 3.0 * hfit - 2.0 * hL[n])
                check(hLok, f"I3 hL_{n+1}", c)
                check(hRok, f"I3 hR_{n+1}", c)
                h13 = (2.0 * hL[n] + hR[n]) / 3.0
                h23 = (hL[n] + 2.0 * hR[n]) / 3.0
                check(h13 - 1e-12 <= hfit <= h23 + 1e-12,
                      f"I3 central third {n+1}", c)
            else:
                check(g0[n] == 0.0 and g1[n] == 0.0 and hL[n] == 0.0
                      and hR[n] == 0.0, f"I3 empty fit {n+1}", c)

        # I4: boundary transfers (donor is 1-based Fortran category)
        for n in range(NCAT - 1):
            nd = donor[n]
            if nd == n + 2:  # transfer from n+1 to n (donor = upper cat)
                etamin = 0.0
                etamax = min(hin[n + 1], hR[nd - 1]) - hL[nd - 1]
                ok_dir = hb[n + 1] <= hin[n + 1]
            elif nd == n + 1:  # transfer from n to n+1 (donor = lower cat)
                etamin = max(hin[n + 1], hL[n]) - hL[n]
                etamax = min(hb[n + 1], hR[n]) - hL[n]
                ok_dir = hb[n + 1] > hin[n + 1]
            else:
                etamin = etamax = 0.0
                ok_dir = True
            check(ok_dir, f"I4 donor dir {n+1}", c)
            exp_da = exp_dv = 0.0
            if etamax > etamin:
                x1 = etamax - etamin
                x2 = (etamax ** 2 - etamin ** 2) / 2.0
                x3 = (etamax ** 3 - etamin ** 3) / 3.0
                exp_da = g1[nd - 1] * x2 + g0[nd - 1] * x1
                exp_dv = g1[nd - 1] * x3 + g0[nd - 1] * x2 + exp_da * hL[nd - 1]
            adon = acm[nd - 1] if nd == 1 else ac[nd - 1]
            vdon = vcm[nd - 1]
            if nd and (exp_da < adon * PUNY or exp_dv < vdon * PUNY):
                exp_da = exp_dv = 0.0
            if exp_da > adon * (1.0 - PUNY) or exp_dv > vdon * (1.0 - PUNY):
                exp_da, exp_dv = adon, vdon
            check(close(daice[n], exp_da), f"I4 daice_{n+1}", c)
            check(close(dvice[n], exp_dv), f"I4 dvice_{n+1}", c)
            if exp_da == 0.0 and exp_dv == 0.0:
                # zeroed by the puny clamps (donor = 0), unless the donor
                # category is empty (clamps cannot fire, donor retained)
                check(nd == 0 or adon <= PUNY, f"I4 donor zeroed {n+1}", c)

        # I5: post state
        acpost = list(acm)
        vcpost = list(vcm)
        for n in range(NCAT - 1):
            nd = donor[n]
            if nd == n + 2:      # n+1 loses, n gains (donor = upper cat)
                vcpost[n] += dvice[n]
                vcpost[n + 1] -= dvice[n]
                acpost[n] += daice[n]
                acpost[n + 1] -= daice[n]
            elif nd == n + 1:    # n loses, n+1 gains (donor = lower cat)
                vcpost[n] -= dvice[n]
                vcpost[n + 1] += dvice[n]
                acpost[n] -= daice[n]
                acpost[n + 1] += daice[n]
        # hi_min closure on category 1
        if acpost[0] > PUNY and vcpost[0] / acpost[0] < HI_MIN:
            exp_da0h = acpost[0] * (1.0 - vcpost[0] / acpost[0] / HI_MIN)
            acpost[0] -= da0h
            check(close(da0h, exp_da0h), "I5 da0_himin", c)
        else:
            check(da0h == 0.0, "I5 da0_himin zero", c)
        for n in range(NCAT):
            check(close(tofloat(x[f"aicen_post_{n+1}"]), acpost[n]),
                  f"I5 aicen_post_{n+1}", c)
            check(close(tofloat(x[f"vicen_post_{n+1}"]), vcpost[n]),
                  f"I5 vicen_post_{n+1}", c)
        aice_post = sum(acpost)
        check(close(tofloat(x["aice_post"]), aice_post), "I5 aice_post", c)
        check(close(tofloat(x["aice0_post"]), 1.0 - aice_post), "I5 aice0_post", c)

        # I1b: the dumped hicen is the post-shift (post-closure) thickness
        for n in range(NCAT):
            if acpost[n] > PUNY:
                check(close(hc[n], vcpost[n] / acpost[n]), f"I1b hicen_{n+1}", c)
            else:
                check(hc[n] == 0.0, f"I1b hicen empty {n+1}", c)

        # I6: exit state
        aice_exit = sum(tofloat(x[f"aicen_exit_{n}"]) for n in range(1, NCAT + 1))
        check(close(tofloat(x["aice_exit"]), aice_exit), "I6 aice_exit", c)
        check(close(tofloat(x["aice0_exit"]), 1.0 - aice_exit), "I6 aice0_exit", c)

        # I7: conservation totals (linear_itd's own thresholds)
        for base, thresh in (("vice", PUNY), ("vsno", PUNY),
                             ("eice", PUNY * LFRESH * RHOI),
                             ("esno", PUNY * LFRESH * RHOS),
                             ("sice", PUNY), ("vbri", PUNY * 10.0)):
            check(abs(tofloat(x[f"{base}_init"]) - tofloat(x[f"{base}_final"]))
                  <= thresh, f"I7 {base}", c)

    print(f"rows: {len(r)} ridge, {len(i)} itd")
    print(f"ridge niter distribution: {dict(sorted(niter_dist.items()))}")
    print(f"itd remap_flag distribution: {dict(sorted(remap_dist.items()))}")
    if fails:
        print(f"FAIL: {len(fails)} of {ncheck[0]} checks failed, e.g.:")
        for f in fails[:10]:
            print("  ", f)
        sys.exit(1)
    print(f"PASS: {ncheck[0]}/{ncheck[0]} identity checks "
          f"({ncheck[0] // (len(r) + len(i))} per row x {len(r) + len(i)} rows)")

if __name__ == "__main__":
    main()
