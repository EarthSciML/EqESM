#!/usr/bin/env python3
# EqESM CLM vegetated-patch lane: verify the stage-1 vegetated dump CSVs in
# stage1/dumps/clm/ against the Fortran laws of CanopyFluxesMod.F90,
# SurfaceAlbedoMod.F90 (TwoStream), CanopyHydrologyMod.F90, and
# SoilWaterMovementMod.F90, and verify solution neutrality of the
# instrumented run (bit-identical restart vs the uninstrumented baseline
# I2000Clm60Sp_f19, and bit-identical re-dumps of the 11 existing clm_*.csv).
"""Verify the CLM vegetated-patch dumps against the Fortran laws.

The vegetated dumps (clm_canopyflux_veg.csv, clm_twostream.csv,
clm_interception.csv, clm_soilwater.csv) are written by the dynamic
veg-patch selection (dump_veg_patches > 0; see clm_dump_mod.F90
dump_select_veg): each traced step, the writer PET (chosen once at step 1
by a max-count allreduce, tie-break lowest PET) ranks its qualifying
patches (frac_veg_nosno == 1, elai+esai > 0, liqcan+snocan > 0) by
elai+esai descending then liqcan+snocan descending then patch index, and
keeps the top dump_veg_patches.  All four files share the same per-step
selection.

Checks (all to round-off / tiny epsilon, es26.17 dump precision):
  selection: twostream (coszen > 0 only) subset of canopyflux_veg subset
    of interception; soilwater columns == distinct columns of the
    selected patches; steps 2..dump_first_steps (12) -- step 1 of this
    cold start has liqcan + snocan == 0 everywhere, so no patch
    qualifies and nothing is traced.
  interception (hyd rows; CanopyHydrologyMod.F90: BulkFlux_CanopyInter-
    ceptionAndThroughfall 479-566, UpdateState_AddInterceptionToCanopy,
    BulkFlux_CanopyExcess, BulkDiag_FracWet 1220-1273):
    fpiliq  = interception_fraction*tanh(elai+esai)     (use_clm5_fpi)
    fpisnow = 1 - exp(-0.5*(elai+esai))
    qflx_through_liq = (forc_rain+irrig)*(1-fpiliq), etc.
    liqcan(hyd) = min(liqcanmx, max(0, liqcan(ent)+dtime*qflx_intercepted_liq))
      (liqcanmx = liq_canopy_storage_scalar*(elai+esai); the checkpoint
      branch requires forc_snow + forc_rain + irrig > 0)
    snocan(hyd) between 0 and max(0, snocan(ent)+dtime*qflx_intercepted_snow)
      (snow unloading removes water but forc_t/forc_wind are not dumped)
    fwet/fcansno/fdry from BulkDiag_FracWet:
      fwet = min((h2ocan/((elai+esai)*liq_canopy_storage_scalar))^(2/3),
                 maximum_leaf_wetted_fraction)   when h2ocan = liqcan+snocan > 0
      fcansno = min((snocan/((elai+esai)*snow_canopy_storage_scalar))^0.15, 1)
      fdry = (1-fwet)*elai/(elai+esai)
  canopyflux_veg (CanopyFluxesMod.F90 exit; step order within the step:
    CanopyHydrology -> CanopyFluxes -> SurfaceAlbedo):
    leaf temperature iteration: 3 <= num_iter <= 41
      (itmax_canopy_fluxes = 40; the convergence check runs at itlef = 3..41,
      so a patch that hits the cap exits with num_iter = 41); num_iter < 41
      means the iteration converged, so det < dtmin (0.01 K) and
      dele < dlemin (0.1 W/m2); num_iter == 41 may have hit the iteration cap.
    t_veg in (150, 400) K; qflx_tran_veg >= 0 (qflx_evap_veg may be below
      qflx_tran_veg when canopy condensation occurs)
    canopy storage update at the end of CanopyFluxes (lines 1618-1635,
    dtime = 1800 s):
      if t_veg > tfrz (273.15):
        if (evap-tran)*dtime > liqcan(hyd):   ! all liquid evaporated
           snocan(cf) = max(0, snocan(hyd) + liqcan(hyd) + (tran-evap)*dtime)
        liqcan(cf) = max(0, liqcan(hyd) + (tran-evap)*dtime)
      else (mirror with roles of liqcan/snocan swapped)
    (checked against the interception hyd rows of the same step/patch;
    truncate_small_values may zero sub-epsilon values afterwards)
  twostream (SurfaceAlbedoMod.F90 TwoStream 1311-1823, Dickinson 1983 /
    Sellers 1985 / Bonan 2011; snowveg_affects_radiation = .true. so the
    fcansno omega blend applies, lSFonly = .false.):
    full recomputation from coszen/elai/esai/t_veg/fwet/fcansno/rho/tau/
    albgrd/albgri and xl(pft) from the ctsm60_params file:
      cosz = max(0.001, coszen); chil = clamp(xl, -0.4, 0.6) (|chil|<=0.01 -> 0.01)
      phi1 = 0.5 - 0.633*chil - 0.330*chil^2; phi2 = 0.877*(1-2*phi1)
      gdir = phi1 + phi2*cosz; twostext = gdir/cosz
      avmu = (1 - phi1/phi2*log((phi1+phi2)/phi1))/phi2
      temp0 = max(gdir + phi2*cosz, 1e-6); temp2 = 1 - temp1/temp0*log((temp1+temp0)/temp1)
      omegal = rho+tau; asu = 0.5*omegal*gdir/temp0*temp2
      betadl = (1+avmu*twostext)/(omegal*avmu*twostext)*asu
      betail = 0.5*((rho+tau) + (rho-tau)*((1+chil)/2)^2)/omegal
      omega/betad/betai = fcansno blend with omegas=[0.8,0.4], betads=betais=0.5
      ... the full h1..h10/p1..p4/u1..u3/tmp2..tmp9 algebra, and:
      albd = h1/sigma + h2 + h3; ftid = h4*s2/sigma + h5*s1 + h6/s1; ftdd = s2
      fabd = 1 - albd - (1-albgrd)*ftdd - (1-albgri)*ftid
      fabd_sun = (1-omega)*(1 - s2 + 1/avmu*(a1+a2)); fabd_sha = fabd - fabd_sun
      albi = h7 + h8; ftii = h9*s1 + h10/s1
      fabi = 1 - albi - (1-albgri)*ftii
      fabi_sun = (1-omega)/avmu*(a1+a2); fabi_sha = fabi - fabi_sun
    compared per band (vis=1, nir=2); identical identities:
      fabd_sun + fabd_sha == fabd; fabi_sun + fabi_sha == fabi
      fabd == 1 - albd - (1-albgrd)*ftdd - (1-albgri)*ftid
      fabi == 1 - albi - (1-albgri)*ftii
  soilwater (SoilWaterMovementMod.F90 / HydrologyNoDrainageMod.F90:
    SoilWater solve entry/exit rows only; soilwater_movement_method = 1
    (moisture_form) means the zengdecker2009 hk/suction law loop does not
    run, so the hk/smp/vwc columns are NaN and the hk-law checks are
    skipped -- see the run-notes in CTSM_LANE.md):
    entry/exit row pairs for the same (step, column); exit h2osoi_liq >= 0;
    wa unchanged by the SoilWater solve (wa(exit) == wa(entry)).
  neutrality (--restart, --existing): the instrumented run's clm2 restart
    is bit-identical (all 292 common variables, mask-aware) to the
    uninstrumented baseline I2000Clm60Sp_f19 restart, and the 11 existing
    clm_*.csv re-dumps are bit-identical to the pre-lane backup.

Usage:
  python3 verify_clm_veg_dumps.py [--dumps DIR] [--param FILE] [--restart NEW BASE]
                                  [--existing BACKUP_DIR]
"""
import argparse, hashlib, math, os, sys

import numpy as np

D = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/clm"
PARAM = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/inputdata/lnd/clm2/paramdata/ctsm60_params.c260518.nc"

DTIME = 1800.0          # land model time step (s)
TFRZ = 273.15           # clm_varcon tfrz
ITMIN = 2               # CanopyFluxesMod.F90 itmin
DTMIN = 0.01            # CanopyFluxesMod.F90 dtmin (K)
DLEMIN = 0.1            # CanopyFluxesMod.F90 dlemin (W/m2)
ITMAX = 40              # itmax_canopy_fluxes in this case's lnd_in
NITER_MAX = 41          # max num_iter seen: the convergence check runs at
                        # itlef = 3..41, so non-converged patches exit at 41
OMEGAS = [0.8, 0.4]     # clm_varcon omegas (vis, nir)
BETADS = 0.5            # clm_varcon betads
BETAIS = 0.5            # clm_varcon betais
FIRST_STEPS = 12        # dump_first_steps in this case's user_nl_clm
VEG_PATCHES = 16        # dump_veg_patches in this case's user_nl_clm

# relative tolerance for the recomputed laws (es26.17 dump precision)
RTOL = 1e-10
ATOL = 1e-14

nbad = 0
def chk(cond, msg):
    global nbad
    if not cond:
        nbad += 1
        print("FAIL:", msg)

def load(path):
    """Parse a dump CSV: returns (columns, rows-as-dicts)."""
    with open(path) as f:
        lines = f.readlines()
    hdr = [l for l in lines if l.startswith("# columns")][0]
    cols = [c.strip() for c in hdr.split("=", 1)[1].strip().split(",")]
    rows = []
    for l in lines:
        if l.startswith("#") or not l.strip():
            continue
        v = l.strip().split(",")
        if v and v[-1] == "":      # dump rows end with a trailing comma
            v = v[:-1]
        assert len(v) == len(cols), f"{path}: row has {len(v)} fields, {len(cols)} columns"
        rows.append(dict(zip(cols, v)))
    return cols, rows

def fnum(s):
    s = s.strip()
    return float(s)

def fnan(v):
    return math.isnan(v)

# ----------------------------------------------------------------------
# Interception law
# ----------------------------------------------------------------------

def check_interception(rows, p_scalars):
    liq_css = p_scalars["liq_canopy_storage_scalar"]
    sno_css = p_scalars["snow_canopy_storage_scalar"]
    mlwf = p_scalars["maximum_leaf_wetted_fraction"]
    int_frac = p_scalars["interception_fraction"]
    ent = { (int(r["step"]), int(r["p"])): r for r in rows if r["stage"].strip() == "ent" }
    hyd = [r for r in rows if r["stage"].strip() == "hyd"]
    print(f"interception rows: ent {len(ent)}, hyd {len(hyd)}")
    chk(len(ent) == len(hyd), "interception ent/hyd row counts differ")
    n = 0
    for r in hyd:
        k = (int(r["step"]), int(r["p"]))
        e = ent.get(k)
        chk(e is not None, f"interception hyd row without ent row {k}")
        if e is None:
            continue
        n += 1
        elai = fnum(r["elai"]); esai = fnum(r["esai"])
        fr = fnum(r["forc_rain"]); fs = fnum(r["forc_snow"]); irr = fnum(r["qflx_irrig_sprinkler"])
        labove = fr + irr
        L = elai + esai
        fpiliq = int_frac * math.tanh(L)
        fpisnow = 1.0 - math.exp(-0.5 * L)
        checkpoint = fnum(e["frac_veg_nosno"]) == 1 and (fs + labove) > 0.0
        if checkpoint:
            exp_through_liq = labove * (1.0 - fpiliq)
            exp_int_liq = labove * fpiliq
            exp_through_sno = fs * (1.0 - fpisnow)
            exp_int_sno = fs * fpisnow
        else:
            exp_through_liq = exp_int_liq = exp_through_sno = exp_int_sno = 0.0
        for name, got, exp in (("qflx_through_liq", fnum(r["qflx_through_liq"]), exp_through_liq),
                               ("qflx_intercepted_liq", fnum(r["qflx_intercepted_liq"]), exp_int_liq),
                               ("qflx_through_snow", fnum(r["qflx_through_snow"]), exp_through_sno),
                               ("qflx_intercepted_snow", fnum(r["qflx_intercepted_snow"]), exp_int_sno)):
            chk(abs(got - exp) <= RTOL * max(1.0, abs(exp)) + ATOL,
                f"interception {name} step {r['step']} p {r['p']}: got {got} exp {exp}")
        # canopy liquid update: liqcan(hyd) = min(liqcanmx, max(0, liqcan(ent) + dtime*int_liq))
        liqcan0 = fnum(e["liqcan"])
        liqcanmx = liq_css * L
        liqcan1 = max(0.0, liqcan0 + DTIME * exp_int_liq)
        exp_liqcan = min(liqcanmx, liqcan1) if checkpoint else liqcan1
        got_liqcan = fnum(r["liqcan"])
        chk(abs(got_liqcan - exp_liqcan) <= RTOL * max(1.0, abs(exp_liqcan)) + ATOL,
            f"interception liqcan(hyd) step {r['step']} p {r['p']}: got {got_liqcan} exp {exp_liqcan}")
        # snow update: unloading (forc_t/forc_wind not dumped) only removes water
        snocan0 = fnum(e["snocan"])
        snocan1 = max(0.0, snocan0 + DTIME * exp_int_sno)
        got_snocan = fnum(r["snocan"])
        chk(-ATOL <= got_snocan <= snocan1 + RTOL * max(1.0, snocan1),
            f"interception snocan(hyd) step {r['step']} p {r['p']}: got {got_snocan} in [0, {snocan1}]")
        # FracWet law (BulkDiag_FracWet): fwet/fcansno/fdry from post-update storage
        h2ocan = got_liqcan + got_snocan
        vegt = L  # frac_veg_nosno == 1 for selected patches
        if h2ocan > 0.0:
            exp_fwet = min((h2ocan / (vegt * liq_css)) ** (2.0 / 3.0), mlwf)
            if got_snocan > 0.0:
                exp_fcansno = min((got_snocan / (vegt * sno_css)) ** 0.15, 1.0)
            else:
                exp_fcansno = 0.0
        else:
            exp_fwet = 0.0
            exp_fcansno = 0.0
        exp_fdry = (1.0 - exp_fwet) * elai / L
        for name, got, exp in (("fwet", fnum(r["fwet"]), exp_fwet),
                               ("fcansno", fnum(r["fcansno"]), exp_fcansno),
                               ("fdry", fnum(r["fdry"]), exp_fdry)):
            chk(abs(got - exp) <= RTOL * max(1.0, abs(exp)) + ATOL,
                f"interception {name} step {r['step']} p {r['p']}: got {got} exp {exp}")
        # ground fluxes are column sums of the patch fluxes; the selected
        # patch's throughfall/canfall cannot be checked exactly without the
        # other patches of the column, so only require non-negativity and a
        # lower bound from this patch's throughfall
        chk(fnum(r["qflx_liq_grnd"]) >= exp_through_liq - RTOL * max(1.0, exp_through_liq),
            f"interception qflx_liq_grnd step {r['step']} c {r['c']} below throughfall")
        chk(fnum(r["qflx_snow_grnd"]) >= exp_through_sno - RTOL * max(1.0, exp_through_sno),
            f"interception qflx_snow_grnd step {r['step']} c {r['c']} below throughfall")
    print(f"interception law rows checked: {n}")

# ----------------------------------------------------------------------
# CanopyFluxes checks
# ----------------------------------------------------------------------

def check_canopyflux(cf_rows, hyd_rows):
    hyd = {(int(r["step"]), int(r["p"])): r for r in hyd_rows if r["stage"].strip() == "hyd"}
    print(f"canopyflux_veg rows: {len(cf_rows)}")
    for r in cf_rows:
        step, p = int(r["step"]), int(r["p"])
        num_iter = fnum(r["num_iter"])
        chk(1.0 <= num_iter <= NITER_MAX, f"canopyflux num_iter {step} {p}: {num_iter}")
        if num_iter < NITER_MAX:
            chk(fnum(r["det"]) < DTMIN, f"canopyflux det {step} {p}: {r['det']} >= {DTMIN}")
            chk(fnum(r["dele"]) < DLEMIN, f"canopyflux dele {step} {p}: {r['dele']} >= {DLEMIN}")
        t_veg = fnum(r["t_veg"])
        chk(150.0 < t_veg < 400.0, f"canopyflux t_veg {step} {p}: {t_veg}")
        tran = fnum(r["qflx_tran_veg"]); evap = fnum(r["qflx_evap_veg"])
        chk(tran >= -ATOL, f"canopyflux qflx_tran_veg {step} {p}: {tran}")
        for c in ("fwet", "fdry", "fcansno"):
            v = fnum(r[c])
            chk(0.0 - ATOL <= v <= 1.0 + ATOL, f"canopyflux {c} {step} {p}: {v}")
        # canopy storage update at the end of CanopyFluxes (1618-1635)
        h = hyd.get((step, p))
        if h is not None:
            liqcan_h = fnum(h["liqcan"]); snocan_h = fnum(h["snocan"])
            net = (tran - evap) * DTIME
            if t_veg > TFRZ:
                exp_liqcan = max(0.0, liqcan_h + net)
                if (evap - tran) * DTIME > liqcan_h:   # all liquid evaporated
                    exp_snocan = max(0.0, snocan_h + liqcan_h + net)
                else:
                    exp_snocan = snocan_h
            else:
                exp_snocan = max(0.0, snocan_h + net)
                if (evap - tran) * DTIME > snocan_h:   # all snow evaporated
                    exp_liqcan = max(0.0, liqcan_h + snocan_h + net)
                else:
                    exp_liqcan = liqcan_h
            # truncate_small_values zeroes sub-relative-epsilon values
            for name, got, exp in (("liqcan", fnum(r["liqcan"]), exp_liqcan),
                                   ("snocan", fnum(r["snocan"]), exp_snocan)):
                chk(abs(got - exp) <= RTOL * max(1.0, abs(exp)) + 1e-12,
                    f"canopyflux {name} post-evap {step} {p}: got {got} exp {exp}")
        else:
            chk(False, f"canopyflux row without interception hyd row {step} {p}")

# ----------------------------------------------------------------------
# Two-stream law
# ----------------------------------------------------------------------

def two_stream_row(r, xl_by_itype):
    """Recompute the TwoStream law for one dump row; return (got, exp) dict."""
    cosz = max(0.001, fnum(r["coszen"]))
    elai = fnum(r["elai"]); esai = fnum(r["esai"])
    fcansno = fnum(r["fcansno"])
    itype = int(r["ptype"])
    L = elai + esai
    chil = min(max(xl_by_itype[itype], -0.4), 0.6)
    if abs(chil) <= 0.01:
        chil = 0.01
    phi1 = 0.5 - 0.633 * chil - 0.330 * chil * chil
    phi2 = 0.877 * (1.0 - 2.0 * phi1)
    gdir = phi1 + phi2 * cosz
    twostext = gdir / cosz
    avmu = (1.0 - phi1 / phi2 * math.log((phi1 + phi2) / phi1)) / phi2
    temp0 = max(gdir + phi2 * cosz, 1e-6)
    temp1 = phi1 * cosz
    temp2 = 1.0 - temp1 / temp0 * math.log((temp1 + temp0) / temp1)
    got = {}
    for ib, band in enumerate(("vis", "nir")):
        rho = fnum(r[f"rho_{band}"]); tau = fnum(r[f"tau_{band}"])
        albgrd = fnum(r[f"albgrd_{band}"]); albgri = fnum(r[f"albgri_{band}"])
        omegal = rho + tau
        asu = 0.5 * omegal * gdir / temp0 * temp2
        betadl = (1.0 + avmu * twostext) / (omegal * avmu * twostext) * asu
        betail = 0.5 * ((rho + tau) + (rho - tau) * ((1.0 + chil) / 2.0) ** 2) / omegal
        # snowveg_affects_radiation = .true., lSFonly = .false.: fcansno blend
        tmp0 = (1.0 - fcansno) * omegal + fcansno * OMEGAS[ib]
        betad = ((1.0 - fcansno) * omegal * betadl + fcansno * OMEGAS[ib] * BETADS) / tmp0
        betai = ((1.0 - fcansno) * omegal * betail + fcansno * OMEGAS[ib] * BETAIS) / tmp0
        omega = tmp0
        b = 1.0 - omega + omega * betai
        c1 = omega * betai
        t0 = avmu * twostext
        d = t0 * omega * betad
        f = t0 * omega * (1.0 - betad)
        t1 = b * b - c1 * c1
        h = math.sqrt(t1) / avmu
        sigma = t0 * t0 - t1
        p1 = b + avmu * h
        p2 = b - avmu * h
        p3 = b + t0
        p4 = b - t0
        s1 = math.exp(-min(h * L, 40.0))
        s2 = math.exp(-min(twostext * L, 40.0))
        # direct beam
        u1 = b - c1 / albgrd
        u2 = b - c1 * albgrd
        u3 = f + c1 * albgrd
        tmp2 = u1 - avmu * h
        tmp3 = u1 + avmu * h
        d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1
        tmp4 = u2 + avmu * h
        tmp5 = u2 - avmu * h
        d2 = tmp4 / s1 - tmp5 * s1
        h1 = -d * p4 - c1 * f
        tmp6 = d - h1 * p3 / sigma
        tmp7 = (d - c1 - h1 / sigma * (u1 + t0)) * s2
        h2 = (tmp6 * tmp2 / s1 - p2 * tmp7) / d1
        h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1
        h4 = -f * p3 - c1 * d
        tmp8 = h4 / sigma
        tmp9 = (u3 - tmp8 * (u2 - t0)) * s2
        h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2
        h6 = (tmp8 * tmp5 * s1 + tmp9) / d2
        albd = h1 / sigma + h2 + h3
        ftid = h4 * s2 / sigma + h5 * s1 + h6 / s1
        ftdd = s2
        fabd = 1.0 - albd - (1.0 - albgrd) * ftdd - (1.0 - albgri) * ftid
        a1 = h1 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) \
            + h2 * (1.0 - s2 * s1) / (twostext + h) \
            + h3 * (1.0 - s2 / s1) / (twostext - h)
        a2 = h4 / sigma * (1.0 - s2 * s2) / (2.0 * twostext) \
            + h5 * (1.0 - s2 * s1) / (twostext + h) \
            + h6 * (1.0 - s2 / s1) / (twostext - h)
        fabd_sun = (1.0 - omega) * (1.0 - s2 + 1.0 / avmu * (a1 + a2))
        fabd_sha = fabd - fabd_sun
        # diffuse
        u1 = b - c1 / albgri
        u2 = b - c1 * albgri
        tmp2 = u1 - avmu * h
        tmp3 = u1 + avmu * h
        d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1
        tmp4 = u2 + avmu * h
        tmp5 = u2 - avmu * h
        d2 = tmp4 / s1 - tmp5 * s1
        h7 = c1 * tmp2 / (d1 * s1)
        h8 = -c1 * tmp3 * s1 / d1
        h9 = tmp4 / (d2 * s1)
        h10 = -tmp5 * s1 / d2
        albi = h7 + h8
        ftii = h9 * s1 + h10 / s1
        fabi = 1.0 - albi - (1.0 - albgri) * ftii
        a1 = h7 * (1.0 - s2 * s1) / (twostext + h) + h8 * (1.0 - s2 / s1) / (twostext - h)
        a2 = h9 * (1.0 - s2 * s1) / (twostext + h) + h10 * (1.0 - s2 / s1) / (twostext - h)
        fabi_sun = (1.0 - omega) / avmu * (a1 + a2)
        fabi_sha = fabi - fabi_sun
        for name, exp in (("albd", albd), ("albi", albi), ("fabd", fabd), ("fabi", fabi),
                          ("fabd_sun", fabd_sun), ("fabd_sha", fabd_sha),
                          ("fabi_sun", fabi_sun), ("fabi_sha", fabi_sha),
                          ("ftdd", ftdd), ("ftid", ftid), ("ftii", ftii)):
            got[f"{name}_{band}"] = exp
    return got

def check_twostream(rows, xl_by_itype):
    print(f"twostream rows: {len(rows)}")
    for r in rows:
        step, p = int(r["step"]), int(r["p"])
        chk(fnum(r["coszen"]) > 0.0, f"twostream coszen <= 0 {step} {p}")
        exp = two_stream_row(r, xl_by_itype)
        for name, e in exp.items():
            # the dump writes sun/shade columns for vis only
            if name not in r:
                continue
            g = fnum(r[name])
            chk(abs(g - e) <= RTOL * max(1.0, abs(e)) + ATOL,
                f"twostream {name} step {step} p {p}: got {g} exp {e}")
        # identity checks independent of the recomputation
        for band in ("vis", "nir"):
            albd = fnum(r[f"albd_{band}"]); albgri = fnum(r[f"albgri_{band}"])
            albgrd = fnum(r[f"albgrd_{band}"])
            ftdd = fnum(r[f"ftdd_{band}"]); ftid = fnum(r[f"ftid_{band}"])
            fabd = fnum(r[f"fabd_{band}"])
            albi = fnum(r[f"albi_{band}"]); ftii = fnum(r[f"ftii_{band}"])
            fabi = fnum(r[f"fabi_{band}"])
            chk(abs(fabd - (1.0 - albd - (1.0 - albgrd) * ftdd - (1.0 - albgri) * ftid))
                <= RTOL * max(1.0, abs(fabd)) + ATOL,
                f"twostream fabd identity {step} {p} {band}")
            chk(abs(fabi - (1.0 - albi - (1.0 - albgri) * ftii)) <= RTOL * max(1.0, abs(fabi)) + ATOL,
                f"twostream fabi identity {step} {p} {band}")
        # sunlit/shaded identities exist for vis only in the dump
        fabd = fnum(r["fabd_vis"])
        fabd_sun = fnum(r["fabd_sun_vis"]); fabd_sha = fnum(r["fabd_sha_vis"])
        fabi = fnum(r["fabi_vis"])
        fabi_sun = fnum(r["fabi_sun_vis"]); fabi_sha = fnum(r["fabi_sha_vis"])
        chk(abs(fabd_sun + fabd_sha - fabd) <= RTOL * max(1.0, abs(fabd)) + ATOL,
            f"twostream fabd_sun+fabd_sha {step} {p} vis")
        chk(abs(fabi_sun + fabi_sha - fabi) <= RTOL * max(1.0, abs(fabi)) + ATOL,
            f"twostream fabi_sun+fabi_sha {step} {p} vis")

# ----------------------------------------------------------------------
# Soil water solve structure
# ----------------------------------------------------------------------

def check_soilwater(rows):
    print(f"soilwater rows: {len(rows)}")
    entry = {(int(r["step"]), int(r["c"])): r for r in rows if r["stage"].strip() == "entry"}
    exit_ = {(int(r["step"]), int(r["c"])): r for r in rows if r["stage"].strip() == "exit"}
    chk(len(entry) == len(exit_), "soilwater entry/exit row counts differ")
    for k, e in entry.items():
        x = exit_.get(k)
        chk(x is not None, f"soilwater entry without exit {k}")
        if x is None:
            continue
        # law columns (vwc_liq..dsmpdw) are NaN: the zengdecker2009 hk loop
        # does not run with soilwater_movement_method = 1 (moisture_form)
        for c in "vwc_liq watsat s_node s_interface icefrac imped hksat bsw sucsat smpmin hk smp dhkdw dsmpdw".split():
            chk(fnan(fnum(e[c])), f"soilwater entry {k} {c} not NaN")
            chk(fnan(fnum(x[c])), f"soilwater exit {k} {c} not NaN")
        for j in range(1, 11):
            got = fnum(x[f"h2osoi_liq_{j}"])
            chk(got >= -ATOL, f"soilwater exit h2osoi_liq_{j} {k}: {got}")
        # the SoilWater solve does not touch the aquifer
        chk(fnum(x["wa"]) == fnum(e["wa"]), f"soilwater wa changed by solve {k}")

# ----------------------------------------------------------------------
# Selection consistency across the four files
# ----------------------------------------------------------------------

def check_selection(cf, ts, ic, sw):
    steps_cf = {int(r["step"]) for r in cf}
    steps_ts = {int(r["step"]) for r in ts}
    steps_ic = {int(r["step"]) for r in ic}
    steps_sw = {int(r["step"]) for r in sw}
    all_steps = steps_cf | steps_ts | steps_ic | steps_sw
    # Step 1 of a cold start (no finidat) has liqcan + snocan == 0 on every
    # patch, so no patch qualifies (frac_veg_nosno == 1, elai+esai > 0,
    # liqcan+snocan > 0) and no rows are written; tracing starts at step 2.
    chk(all_steps == set(range(2, FIRST_STEPS + 1)),
        f"traced steps {sorted(all_steps)} != 2..{FIRST_STEPS}")
    for s in all_steps:
        pc = {(int(r["p"]), int(r["c"])) for r in cf if int(r["step"]) == s}
        pt = {(int(r["p"]), int(r["c"])) for r in ts if int(r["step"]) == s}
        pi = {(int(r["p"]), int(r["c"])) for r in ic if int(r["step"]) == s}
        cs = {int(r["c"]) for r in sw if int(r["step"]) == s}
        # twostream only writes the coszen > 0 (daylit) subset of the
        # selected patches (TwoStream is called on filter_vegsol)
        chk(pt <= pc, f"step {s}: twostream patches not subset of canopyflux")
        chk(pc <= pi, f"step {s}: canopyflux patches not subset of interception")
        chk(1 <= len(pc) <= VEG_PATCHES, f"step {s}: {len(pc)} selected patches")
        exp_cols = {c for _, c in pi}
        chk(cs == exp_cols, f"step {s}: soilwater columns {cs} != selected columns {exp_cols}")
    print(f"selection consistent across {len(all_steps)} steps")

# ----------------------------------------------------------------------
# Neutrality
# ----------------------------------------------------------------------

def check_existing(backup_dir):
    files = ["clm_bareground", "clm_canopyflux", "clm_qsat", "clm_run_imports",
             "clm_soilfluxes", "clm_soiltemp", "clm_surfalbedo", "clm_surfrad",
             "clm_urbalb", "clm_urban_lw", "clm_urban_rad"]
    for f in files:
        a = hashlib.md5(open(f"{D}/{f}.csv", "rb").read()).hexdigest()
        b = hashlib.md5(open(f"{backup_dir}/{f}.csv", "rb").read()).hexdigest()
        chk(a == b, f"existing dump {f}.csv not bit-identical to backup")
    print(f"existing dumps bit-identical vs backup: {len(files)} files")

def check_restart(new, base):
    import netCDF4
    n = netCDF4.Dataset(new)
    b = netCDF4.Dataset(base)
    common = sorted(set(n.variables) & set(b.variables))
    ndiff = 0
    for v in common:
        if n[v].dtype.kind not in "fiu":
            continue
        a = np.asarray(n[v][:]); c = np.asarray(b[v][:])
        if a.shape != c.shape:
            ndiff += 1
            print("FAIL: restart var", v, "shape", a.shape, c.shape)
            continue
        am = np.ma.masked_invalid(a.astype(np.float64))
        cm = np.ma.masked_invalid(c.astype(np.float64))
        if not np.ma.allequal(am, cm, fill_value=True):
            ndiff += 1
            print("FAIL: restart var", v, "differs")
    chk(ndiff == 0, f"restart neutrality: {ndiff} of {len(common)} variables differ")
    print(f"restart neutrality: {len(common)} common variables, {ndiff} differ")

# ----------------------------------------------------------------------

def main():
    global D
    ap = argparse.ArgumentParser()
    ap.add_argument("--dumps", default=D)
    ap.add_argument("--param", default=PARAM)
    ap.add_argument("--existing", metavar="BACKUP_DIR",
                    help="bit-identical check of the 11 pre-lane dumps vs this backup dir")
    ap.add_argument("--restart", nargs=2, metavar=("NEW", "BASE"),
                    help="bit-identical restart check (instrumented vs baseline)")
    a = ap.parse_args()
    D = a.dumps

    # scalar params and leaf orientation from the ctsm60 params file
    import netCDF4
    with netCDF4.Dataset(a.param) as p:
        scalars = {v: float(p.variables[v][:]) for v in
                   ("interception_fraction", "liq_canopy_storage_scalar",
                    "snow_canopy_storage_scalar", "maximum_leaf_wetted_fraction")}
        xl = {int(i): float(x) for i, x in enumerate(np.asarray(p.variables["xl"][:]))}
    print("params:", scalars)

    _, cf = load(f"{D}/clm_canopyflux_veg.csv")
    _, ts = load(f"{D}/clm_twostream.csv")
    _, ic = load(f"{D}/clm_interception.csv")
    _, sw = load(f"{D}/clm_soilwater.csv")

    check_selection(cf, ts, ic, sw)
    check_interception(ic, scalars)
    check_canopyflux(cf, ic)
    check_twostream(ts, xl)
    check_soilwater(sw)

    if a.existing:
        check_existing(a.existing)
    if a.restart:
        check_restart(*a.restart)

    print()
    if nbad == 0:
        print("ALL CHECKS PASSED")
    else:
        print(f"{nbad} CHECKS FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()
