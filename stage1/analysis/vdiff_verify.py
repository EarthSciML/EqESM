#!/usr/bin/env python3
"""Independent re-verification of the CAM6 vertical-diffusion dump laws.

Reads the instrumented vdiff solver-context dumps
(``stage1/dumps/cam/cam_vdiff_prof.csv`` and ``cam_vdiff_q.csv``) and
re-derives, purely from the dump columns, the laws transcribed into
``stage2/cam_vertical_diffusion.esm``, printing the maximum relative
errors against the dumped outputs:

* the wind-damping law tau_damp(pver) = -g*(ksrftms +
  max(tau/ws, ksrfmin))/pdel (diffusion_solver.F90:232-239, 274-285)
  against the dumped tau_damp_rate column, with tau = sqrt(taux^2 +
  tauy^2), ws = max(sqrt(u0^2 + v0^2), wsmin), u0/v0 recovered as
  u0 = u1 - ptend_u*dt (the dumped ptend is (u1-u0)/dt by
  construction, F90:1051-1054), ksrfmin = 1e-4, wsmin = 1;
* the surface-layer source identity ptend_u(pver) =
  tau_damp(pver)*u1(pver) + g*tauresx*ramda/pdel(pver) (the pver
  injection of F90:491-499 with ramda = dt/timeres = 1800/7200,
  F90:495), mirrored for v;
* dtk = 0.5(u0^2+v0^2) - 0.5(u1^2+v1^2) (F90:598-609) against the
  dumped dtk column;
* ptend_s = dtk/dt (the do_beljaars branch, F90:610-611, and the
  tendency of F90:1051 with kvh = 0) against the dumped ptend_s
  column;
* the residual-stress memory tauresx = tauresx_enter - tauimpx with
  tauimpx = sum_k (u1-u0)_k pdel_k / (g dt) (F90:536-589, taux =
  tautmsx = taubljx = 0) against the dumped tauresx/tauresy columns,
  where the step-1 entering residual is the observed pbuf initial
  (1.4478837360691117e-04 Pa) and later steps use the previous step's
  dumped updated value;
* the q surface source ptend_q(pver) = g*cflx*dryfac/pdel(pver) =
  g*cflx/pdeldry(pver) (F90:974-975, 1056-1061) against the dumped
  ptend_q of cam_vdiff_q.csv for every constituent row (rows whose
  dumped ptend_q is exactly 0 are reported separately: the law holds to
  rounding only there).

All laws hold to the rounding of the dumped (u1-u0)/dt cancellation;
the errors reported here are the numbers quoted in README-cam.md and in
the test tolerances of cam_vertical_diffusion.esm. The dump window runs
all kvh/kvm/kvq = 0, so no dumped row exercises the diffusion operator;
the operator itself is pinned by the compiled-Fortran driver in
``stage1/analysis/cam_vdiff_driver/``.

Usage: python3 vdiff_verify.py [dump-dir]
(dump-dir defaults to stage1/dumps/cam relative to the repo root).
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dumpcsv import load, select  # noqa: E402

GRAVIT = 9.80616          # gravit, m/s^2
KSRFMIN = 1.0e-4          # diffusion_solver.F90:44, kg/(m2 s)
WSMIN = 1.0               # diffusion_solver.F90:42, m/s
DT = 1800.0               # ztodt, s
TIMERES = 7200.0          # timeres namelist, s
RAMDA = DT / TIMERES      # F90:495, = 0.25


def relerr(actual, expected):
    if expected == 0.0:
        return 0.0 if actual == 0.0 else float("inf")
    return abs(actual - expected) / abs(expected)


def main(dump_dir):
    prof_path = os.path.join(dump_dir, "cam_vdiff_prof.csv")
    q_path = os.path.join(dump_dir, "cam_vdiff_q.csv")
    _, prof = load(prof_path)
    _, qrows = load(q_path)

    steps = sorted({int(r["step"]) for r in prof})
    pver = max(int(r["lev"]) for r in prof)

    print(f"# {os.path.basename(prof_path)}: {len(prof)} rows, steps {steps}, pver {pver}")

    # per-step scalars from the lev-pver rows
    prof_pver = {int(r["step"]): r for r in prof if int(r["lev"]) == pver}

    max_err = {"tau_damp": 0.0, "ptend_u": 0.0, "ptend_v": 0.0, "dtk": 0.0,
               "ptend_s": 0.0, "tauresx": 0.0, "tauresy": 0.0, "ptend_q": 0.0}
    worst = {}

    tauresx_enter = 1.4478837360691117e-4  # observed step-1 pbuf initial, Pa
    tauresy_enter = 5.519501195026867e-4

    for step in steps:
        rows = select(prof, step=step)
        by_lev = {int(r["lev"]): r for r in rows}
        p32 = prof_pver[step]

        # --- wind damping law (esd lev 1 = CAM pver) ---
        u0 = float(p32["u1"]) - float(p32["ptend_u"]) * DT
        v0 = float(p32["v1"]) - float(p32["ptend_v"]) * DT
        tau = math.sqrt(float(p32["taux"]) ** 2 + float(p32["tauy"]) ** 2)
        ws = max(math.sqrt(u0 * u0 + v0 * v0), WSMIN)
        ksrfturb = max(tau / ws, KSRFMIN)
        ksrf = ksrfturb + float(p32["ksrftms"])
        tau_damp = -GRAVIT * ksrf / float(p32["pdel"])
        e = relerr(tau_damp, float(p32["tau_damp_rate"]))
        if e >= max_err["tau_damp"]:
            max_err["tau_damp"] = e
            worst["tau_damp"] = (step, tau_damp, p32["tau_damp_rate"])

        # --- surface-layer source identity (the tau_damp*u + src row) ---
        src_u = GRAVIT * tauresx_enter * RAMDA / float(p32["pdel"])
        e = relerr(tau_damp * float(p32["u1"]) + src_u, float(p32["ptend_u"]))
        if e >= max_err["ptend_u"]:
            max_err["ptend_u"] = e
            worst["ptend_u"] = (step, tau_damp * float(p32["u1"]) + src_u, p32["ptend_u"])
        src_v = GRAVIT * tauresy_enter * RAMDA / float(p32["pdel"])
        e = relerr(tau_damp * float(p32["v1"]) + src_v, float(p32["ptend_v"]))
        if e >= max_err["ptend_v"]:
            max_err["ptend_v"] = e
            worst["ptend_v"] = (step, tau_damp * float(p32["v1"]) + src_v, p32["ptend_v"])

        # --- dtk and ptend_s at every level ---
        for lev in sorted(by_lev):
            r = by_lev[lev]
            u0k = float(r["u1"]) - float(r["ptend_u"]) * DT
            v0k = float(r["v1"]) - float(r["ptend_v"]) * DT
            dtk = 0.5 * (u0k * u0k + v0k * v0k) - 0.5 * (float(r["u1"]) ** 2 + float(r["v1"]) ** 2)
            e = relerr(dtk, float(r["dtk"]))
            if e >= max_err["dtk"]:
                max_err["dtk"] = e
                worst["dtk"] = (step, lev, dtk, r["dtk"])
            e = relerr(dtk / DT, float(r["ptend_s"]))
            if e >= max_err["ptend_s"]:
                max_err["ptend_s"] = e
                worst["ptend_s"] = (step, lev, dtk / DT, r["ptend_s"])

        # --- residual stress memory ---
        # u0 = u1 - ptend*dt, so (u1-u0) = ptend*dt (F90:1051-1054)
        tauimpx = sum(float(by_lev[k]["ptend_u"]) * DT * float(by_lev[k]["pdel"])
                      for k in by_lev) / (GRAVIT * DT)
        tauimpy = sum(float(by_lev[k]["ptend_v"]) * DT * float(by_lev[k]["pdel"])
                      for k in by_lev) / (GRAVIT * DT)
        tauresx = tauresx_enter - tauimpx
        tauresy = tauresy_enter - tauimpy
        e = relerr(tauresx, float(p32["tauresx"]))
        if e >= max_err["tauresx"]:
            max_err["tauresx"] = e
            worst["tauresx"] = (step, tauresx, p32["tauresx"])
        e = relerr(tauresy, float(p32["tauresy"]))
        if e >= max_err["tauresy"]:
            max_err["tauresy"] = e
            worst["tauresy"] = (step, tauresy, p32["tauresy"])
        tauresx_enter = float(p32["tauresx"])  # pbuf memory for the next step
        tauresy_enter = float(p32["tauresy"])

    # --- q surface source (cam_vdiff_q.csv) ---
    # ptend_q(pver) = g*cflx*dryfac/pdel(pver) with dryfac = pdel/pdeldry
    # (F90:1056-1061), i.e. = g*cflx/pdeldry with the pdel factor cancelling;
    # pdeldry comes from the prof dump (cam_vdiff_q.csv has no pdel columns).
    q_zero_abs = 0.0  # max |qsrc| on rows whose dumped ptend_q is exactly 0
    for row in qrows:
        if int(row["lev"]) != pver:
            continue
        step = int(row["step"])
        pdeldry = float(prof_pver[step]["pdeldry"])
        qsrc = GRAVIT * float(row["cflx"]) / pdeldry
        dumped = float(row["ptend_q"])
        if dumped == 0.0:  # zero-flux constituent: law holds to rounding only
            q_zero_abs = max(q_zero_abs, abs(qsrc))
            continue
        e = relerr(qsrc, dumped)
        if e >= max_err["ptend_q"]:
            max_err["ptend_q"] = e
            worst["ptend_q"] = (step, row["cnst_name"], qsrc, row["ptend_q"])

    print("# max relative errors (law vs dumped column):")
    for k in ["tau_damp", "ptend_u", "ptend_v", "dtk", "ptend_s",
              "tauresx", "tauresy", "ptend_q"]:
        print(f"  {k:10s} {max_err[k]:.4e}   worst: {worst[k]}")
    if q_zero_abs:
        print(f"  (dumped-zero ptend_q rows: |qsrc| <= {q_zero_abs:.2e} kg/kg/s, rounding only)")


if __name__ == "__main__":
    dump_dir = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dumps", "cam")
    main(dump_dir)
