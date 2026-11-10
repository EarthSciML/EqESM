#!/usr/bin/env python3
# EqESM CLM canopy-veg lane: extract numeric test tuples for the stage-2
# .esm tests of CanopyHydrologyMod.F90 (canopy interception/throughfall and
# the wet/dry/snow-covered fractions) and of the end-of-CanopyFluxes canopy
# storage update (CanopyFluxesMod.F90:1618-1635), from the vegetated-patch
# dumps clm_interception.csv (ent/hyd rows) and clm_canopyflux_veg.csv.
#
# The .esm tests blocks are hand-authored (see stage2/clm_canopy.esm and
# stage2/clm_canopy_storage.esm); this script only extracts the numeric
# tuples for a fixed list of (step, patch) keys and verifies the Fortran
# laws on those rows (mirroring stage1/analysis/verify_clm_veg_dumps.py).
#
# Usage:
#   python3 clm_canopy_veg_test_tuples.py [--dumps DIR]
#
# Prints, per selected (step, patch):
#   - interception: the parameter overrides and expected values for the
#     instantaneous (time_span {0,0}) and one-step ODE (time_span {0,1800})
#     pins of clm_canopy.esm (model CanopyHydrology);
#   - storage update: the parameter overrides and expected outputs for the
#     clm_canopy_storage.esm pins (model CanopyStorageUpdate).
#
# The laws (validated against all 176 rows by verify_clm_veg_dumps.py):
#   fpiliq  = interception_fraction*tanh(elai+esai)   (use_clm5_fpi = .true.)
#   fpisnow = 1 - exp(-0.5*(elai+esai))
#   liqcan(hyd) = min(liqcanmx, max(0, liqcan(ent) + dtime*qflx_intercepted_liq))
#   liqcanmx = liq_canopy_storage_scalar*(elai+esai)  (scalar 0.1)
#   fwet  = min((h2ocan/(vegt*liq_css))^(2/3), 0.05)  (h2ocan = liqcan+snocan > 0)
#   fdry  = (1-fwet)*elai/(elai+esai)
#   storage update (t_veg > tfrz only; all 176 veg rows are above freezing):
#     liqcan(cf) = max(0, liqcan(hyd) + (tran-evap)*dtime)
#     snocan(cf) = max(0, snocan(hyd)+liqcan(hyd)+(tran-evap)*dtime)
#                  when (evap-tran)*dtime > liqcan(hyd), else snocan(hyd)
import argparse, math, os, sys

D = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/clm"
DTIME = 1800.0
TFRZ = 273.15
LIQ_CSS = 0.1       # ctsm60_params.c260518.nc liq_canopy_storage_scalar
SNO_CSS = 6.0       # snow_canopy_storage_scalar
MLWF = 0.05         # maximum_leaf_wetted_fraction
INT_FRAC = 1.0      # interception_fraction
REL_EPS = 1e-10     # truncate_small_values custom_rel_epsilon at the call site
FWET_EXP = 0.666666666666  # the Fortran literal 0.666666666666_r8 (BulkDiag_FracWet)
FCANSNO_EXP = 0.15         # the Fortran literal 0.15_r8 (BulkDiag_FracWet)

# (step, patch) keys pinned by the stage-2 tests (all present in the dumps)
KEYS = [
    (2, 50), (2, 104), (2, 135), (2, 198), (2, 250),
    (3, 155), (3, 223),
    (5, 198), (7, 714), (12, 264),
]

def load(path):
    with open(path) as f:
        lines = f.readlines()
    hdr = [l for l in lines if l.startswith("# columns")][0]
    cols = [c.strip() for c in hdr.split("=", 1)[1].strip().split(",")]
    rows = []
    for l in lines:
        if l.startswith("#") or not l.strip():
            continue
        v = l.strip().split(",")
        if v and v[-1] == "":
            v = v[:-1]
        rows.append(dict(zip(cols, v)))
    return cols, rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dumps", default=D)
    a = ap.parse_args()
    _, ic = load(f"{a.dumps}/clm_interception.csv")
    _, cf = load(f"{a.dumps}/clm_canopyflux_veg.csv")
    ent = {(int(r["step"]), int(r["p"])): r for r in ic if r["stage"].strip() == "ent"}
    hyd = {(int(r["step"]), int(r["p"])): r for r in ic if r["stage"].strip() == "hyd"}
    cfr = {(int(r["step"]), int(r["p"])): r for r in cf}

    nbad = 0
    for k in KEYS:
        e, h, c = ent[k], hyd[k], cfr[k]
        step, p = k
        elai = float(e["elai"]); esai = float(e["esai"])
        L = elai + esai
        rain = float(e["forc_rain"]); fs = float(e["forc_snow"])
        irrig = float(e["qflx_irrig_sprinkler"])
        labove = rain + irrig
        fpiliq = INT_FRAC * math.tanh(L)
        fpisnow = 1.0 - math.exp(-0.5 * L)
        checkpoint = float(e["frac_veg_nosno"]) == 1 and (fs + labove) > 0.0
        if checkpoint:
            exp_int_liq = labove * fpiliq
            exp_thr_liq = labove * (1.0 - fpiliq)
            exp_int_sno = fs * fpisnow
            exp_thr_sno = fs * (1.0 - fpisnow)
        else:
            exp_int_liq = exp_thr_liq = exp_int_sno = exp_thr_sno = 0.0
        liqcan0 = float(e["liqcan"]); snocan0 = float(e["snocan"])
        liqcanmx = LIQ_CSS * L
        snocanmx = SNO_CSS * L
        liqcan_h = float(h["liqcan"]); snocan_h = float(h["snocan"])
        liqcan1 = max(0.0, liqcan0 + DTIME * exp_int_liq)
        exp_liqcan = min(liqcanmx, liqcan1) if checkpoint else liqcan1
        # laws checks (mirror of the validator)
        def chk(cond, msg):
            nonlocal nbad
            if not cond:
                nbad += 1
                print(f"  FAIL: {msg}")
        chk(abs(exp_int_liq - float(h["qflx_intercepted_liq"])) <= 1e-10 * max(1.0, abs(exp_int_liq)) + 1e-14,
            f"intercepted_liq {k}")
        chk(abs(exp_thr_liq - float(h["qflx_through_liq"])) <= 1e-10 * max(1.0, abs(exp_thr_liq)) + 1e-14,
            f"through_liq {k}")
        chk(abs(exp_int_sno - float(h["qflx_intercepted_snow"])) <= 1e-10 * max(1.0, abs(exp_int_sno)) + 1e-14,
            f"intercepted_snow {k}")
        chk(abs(exp_thr_sno - float(h["qflx_through_snow"])) <= 1e-10 * max(1.0, abs(exp_thr_sno)) + 1e-14,
            f"through_snow {k}")
        chk(abs(exp_liqcan - liqcan_h) <= 1e-10 * max(1.0, abs(exp_liqcan)) + 1e-14, f"liqcan(hyd) {k}")
        chk(0.0 - 1e-14 <= snocan_h <= max(0.0, snocan0 + DTIME * exp_int_sno) + 1e-10, f"snocan(hyd) {k}")
        h2ocan = liqcan_h + snocan_h
        vegt = L  # frac_veg_nosno == 1 on all selected patches
        if h2ocan > 0.0:
            exp_fwet = min((h2ocan / (vegt * LIQ_CSS)) ** FWET_EXP, MLWF)
            exp_fcansno = min((snocan_h / (vegt * SNO_CSS)) ** FCANSNO_EXP, 1.0) if snocan_h > 0.0 else 0.0
        else:
            exp_fwet = exp_fcansno = 0.0
        exp_fdry = (1.0 - exp_fwet) * elai / L
        chk(abs(exp_fwet - float(h["fwet"])) <= 1e-10 * max(1.0, abs(exp_fwet)) + 1e-14, f"fwet {k}")
        chk(abs(exp_fcansno - float(h["fcansno"])) <= 1e-10 * max(1.0, abs(exp_fcansno)) + 1e-14, f"fcansno {k}")
        chk(abs(exp_fdry - float(h["fdry"])) <= 1e-10 * max(1.0, abs(exp_fdry)) + 1e-14, f"fdry {k}")
        # storage update (end of CanopyFluxes, all rows above freezing)
        tran = float(c["qflx_tran_veg"]); evap = float(c["qflx_evap_veg"])
        t_veg = float(c["t_veg"])
        net = (tran - evap) * DTIME
        evapnet = (evap - tran) * DTIME
        assert t_veg > TFRZ, f"row {k} below freezing (mirror branch unpinned)"
        exp_liq_cf = max(0.0, liqcan_h + net)
        if evapnet > liqcan_h:
            exp_sno_cf = max(0.0, snocan_h + liqcan_h + net)
        else:
            exp_sno_cf = snocan_h
        # truncate_small_values on snocan (custom_rel_epsilon = 1e-10)
        if abs(exp_sno_cf) < REL_EPS * abs(snocan_h):
            exp_sno_cf = 0.0
        got_liq_cf = float(c["liqcan"]); got_sno_cf = float(c["snocan"])
        chk(abs(got_liq_cf - exp_liq_cf) <= 1e-10 * max(1.0, abs(exp_liq_cf)) + 1e-12, f"liqcan(cf) {k}")
        chk(abs(got_sno_cf - exp_sno_cf) <= 1e-10 * max(1.0, abs(exp_sno_cf)) + 1e-12, f"snocan(cf) {k}")

        print(f"=== ({step}, {p}) ===")
        print(f"  interception: rain={rain:.17e} fs={fs:.17e} irrig={irrig:.17e} "
              f"elai={elai:.17e} esai={esai:.17e} L={L:.17e} checkpoint={checkpoint}")
        print(f"    fpiliq={fpiliq:.17e} fpisnow={fpisnow:.17e} liqcanmx={liqcanmx:.17e} snocanmx={snocanmx:.17e}")
        print(f"    liqcan0={liqcan0:.17e} liqcan_h={liqcan_h:.17e} exp_liqcan={exp_liqcan:.17e}")
        print(f"    int_liq={exp_int_liq:.17e} thr_liq={exp_thr_liq:.17e} int_sno={exp_int_sno:.17e} thr_sno={exp_thr_sno:.17e}")
        print(f"    fwet={exp_fwet:.17e} fdry={exp_fdry:.17e} fcansno={exp_fcansno:.17e} h2ocan={h2ocan:.17e} vegt={vegt:.17e}")
        print(f"    dump: hyd liqcan={liqcan_h:.17e} snocan={snocan_h:.17e} fwet={float(h['fwet']):.17e} "
              f"fdry={float(h['fdry']):.17e} fcansno={float(h['fcansno']):.17e}")
        print(f"  storage: t_veg={t_veg:.17e} tran={tran:.17e} evap={evap:.17e} net={net:.17e} evapnet={evapnet:.17e}")
        print(f"    exp_liq_cf={exp_liq_cf:.17e} exp_sno_cf={exp_sno_cf:.17e} (dump {got_liq_cf:.17e} {got_sno_cf:.17e}) "
              f"branch=all-evap:{evapnet > liqcan_h}")
    print()
    if nbad:
        print(f"{nbad} CHECKS FAILED")
        sys.exit(1)
    print(f"ALL {len(KEYS)} ROWS MATCH THE LAWS")

if __name__ == "__main__":
    main()
