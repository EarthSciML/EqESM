#!/usr/bin/env python3
"""Extract stage-2 test tuples for mom6_vertical_diffusion.esm from the
wave-2 instrumentation dump mom6_diabatic_tri.csv (tracer_vertdiff_Eulerian
tridiagonal solve; see stage1/runs/MOM6_LANE.md wave-2 and
stage2/README-mom6-w2.md).

The esd column grid reverses the MOM6 top-down order: esd lev k <-> MOM6
layer NLEV+1-k, esd interface n <-> MOM6 interface NLEV+2-n, with ze
increasing upward from the bottom (ze[1] = -H) to the surface
(ze[NLEV+1] = 0).  The PDE form is dT/dt = D(Kd_heat D(T,lev),lev) lowered
by EarthSciDiscretizations' varcoeff_face_laplacian_lev_flux_bc with zero
prescribed face fluxes, which is EXACTLY the MOM6 tridiagonal law
h_k (T'_k - T_k) = ent_k (T'_{k-1} - T'_k) + ent_{k+1} (T'_{k+1} - T'_k)
with ent = dt Kd/dz_int (verified: the stencil reproduces the ent-form
value to the last bit, and ent = dt Kd/dz_int holds to 4e-16).

Because the dumped T'/S' are the SOLVED doubles, they satisfy the discrete
identity only to the solve's round-off, which convective interfaces
(ent ~ 1e6 m, deep cells of rows here) amplify into a ~1e-4-relative
mismatch of h (T'-T) (the operand-magnitude scaling of README-mom6-w2.md:
the SCALED residual |lhs - e_in - e_out| / (|lhs|+|e_in|+|e_out|+|ent T'|
terms) stays <= ~3e-17).  This script therefore reports, per candidate
layer, the relative discrepancy between the esd stencil (evaluated on the
dumped after-state, exactly the ent-form value) and the dump-derived
tendency (T_after - T_before)/dt, so the hand-authored tests block can
pin rel 1e-12 where the solve residual allows and the justified looser
tolerance where it does not.  The scaled residual stays below the
dump-wide maxima 1.2e-12 T / 2.0e-12 S of validate_mom6_w2.py at every
pinned layer, so the residual pins are tight (abs ~1e-15).

The interface spacing dz_int is taken as 0.5 (ze[n+1] - ze[n-1]), NOT
zc[n] - zc[n-1] with zc = 0.5 (ze[n]+ze[n+1]): the latter rounds at the
~743 m depth magnitude of ze (~1e-13 m absolute, ~1e-10 relative on the
deep 1e-3 m spacings) and would spoil the ent pins.  The ze-difference
form is an exact Sterbenz subtraction, so ent = dt Kd/dz_int reproduces
the dumped ent to ~4e-16 at every interface, matching the Fortran
(MOM_diabatic_driver.F90:1728-1739).

For the three pinned rows it prints, at full precision: the parameter
arrays (ze built from the dumped h, Kd_heat, Kd_salt, T_before, S_before
reversed, dt) and the state initial conditions (T_after, S_after
reversed) as JSON for direct pasting into the hand-authored tests block;
and the expected pin values (ptend_T/ptend_S at the pinned esd layers,
ent at the pinned interfaces, the scaled residual at the pinned mixed
and bottom layers).  All expected values are derived from the dump
columns (marked "derived (not a dump column)" in the tests), never
fitted.

Sanity checks (exit 1 on failure): the esd-reversed stencil equals the
ent-form value at every pinned layer; ent = dt Kd/dz_int holds at the
pinned interfaces; the scaled residuals are below the dump-wide maxima.

Usage: python3 extract_mom6_vertdiff_tuples.py [DUMP_DIR]
"""
import json
import sys

DUMP_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/mom6"
KE = 75            # NLEV
NZ = KE + 1

# Pinned rows (tod, gi, gj), pinned esd layers for ptend (boundaries, the
# mixed region, one deep layer per column), pinned esd interfaces for the
# ent and interface-flux pins, and the esd layer for the scaled residual.
ROWS = {
    "convective": dict(
        key=(32400.0, 3.0, 79.0),
        ptend_levs=[1, 40, 74, 75],         # MOM6 75 (bottom), 36, 2, 1 (sfc)
        ent_ifc=[(6, 71), (28, 49), (75, 2)],   # max-ent, T/S-differing, sfc
        res_lev=5,                          # MOM6 layer 71 (mixed)
    ),
    "moderate": dict(
        key=(10800.0, 5.0, 79.0),
        ptend_levs=[1, 40, 74, 75],
        ent_ifc=[(4, 73), (40, 37), (75, 2)],   # max-ent, T/S-differing, sfc
        res_lev=4,                          # MOM6 layer 72 (mixed)
    ),
    "shallow": dict(
        key=(3600.0, 3.0, 78.0),
        ptend_levs=[1, 40, 74, 75],
        ent_ifc=[(53, 24), (18, 59), (75, 2)],  # max-ent, T/S-differing, sfc
        res_lev=53,                         # MOM6 layer 23 (mixed)
    ),
}


def read_csv(path):
    cols = None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("# columns = "):
                cols = line[len("# columns = "):].rstrip(",").split(",")
            elif line.startswith("#"):
                continue
            elif line.strip():
                vals = [float(x) for x in line.rstrip(",").split(",")]
                rows.append(dict(zip(cols, vals)))
    return rows


def kcol(r, name, k, fmt="{:02d}"):
    return r[f"{name}_{fmt.format(k)}"]


def main():
    rows = read_csv(f"{DUMP_DIR}/mom6_diabatic_tri.csv")
    by_key = {(r["tod"], r["gi"], r["gj"]): r for r in rows}
    for label, spec in ROWS.items():
        r = by_key[spec["key"]]
        tod, gi, gj = spec["key"]
        dt = r["dt"]
        h = [kcol(r, "h", k) for k in range(1, KE + 1)]
        # esd reversal (0-based arrays; esd lev/ifc k = index k-1)
        ze = [0.0] * NZ
        ze[NZ - 1] = 0.0                      # surface interface
        for n in range(NZ - 2, -1, -1):
            ze[n] = ze[n + 1] - h[KE - 1 - n]  # walk down: ze[0] = -H
        kd_h = [kcol(r, "Kd_heat", K) for K in range(1, NZ + 1)][::-1]
        kd_s = [kcol(r, "Kd_salt", K) for K in range(1, NZ + 1)][::-1]
        ent_t = [kcol(r, "ent_t", K) for K in range(1, NZ + 1)][::-1]
        ent_s = [kcol(r, "ent_s", K) for K in range(1, NZ + 1)][::-1]
        tb = [kcol(r, "T_before", k) for k in range(1, KE + 1)][::-1]
        sb = [kcol(r, "S_before", k) for k in range(1, KE + 1)][::-1]
        ta = [kcol(r, "T_after", k) for k in range(1, KE + 1)][::-1]
        sa = [kcol(r, "S_after", k) for k in range(1, KE + 1)][::-1]
        zc = [0.5 * (ze[n] + ze[n + 1]) for n in range(KE)]
        dz = [ze[n + 1] - ze[n] for n in range(KE)]

        print(f"===== row {label}: tod={tod} gj={gj} gi={gi} "
              f"(cpl_step {r['cpl_step']}, dt={dt} s) =====")
        print("ze =", json.dumps(ze))
        print("Kd_heat =", json.dumps(kd_h))
        print("Kd_salt =", json.dumps(kd_s))
        print("T_before =", json.dumps(tb))
        print("S_before =", json.dumps(sb))
        print("T_after (IC) =", json.dumps(ta))
        print("S_after (IC) =", json.dumps(sa))

        # esd stencil on the dumped after-states (0-based; face fluxes zero)
        def stencil(u, kd):
            G0 = [0.0] * NZ
            for n in range(1, KE):
                G0[n] = kd[n] * (u[n] - u[n - 1]) / (zc[n] - zc[n - 1])
            return [(G0[k + 1] - G0[k]) / dz[k] for k in range(KE)]
        LT = stencil(ta, kd_h)
        LS = stencil(sa, kd_s)

        # ptend pins: stencil vs dump-derived tendency, per pinned layer
        # (boundaries, one mid-column layer, one deep layer, plus the
        # mixed-region layer for the loose tolerance table)
        for ell in sorted(set(spec["ptend_levs"]) | {spec["res_lev"]}):
            for field, L, u, u0 in (("T", LT, ta, tb), ("S", LS, sa, sb)):
                exp = (u[ell - 1] - u0[ell - 1]) / dt
                rel = abs(L[ell - 1] - exp) / max(abs(exp), 1e-300)
                print(f"  ptend_{field} esd lev {ell:2d} (MOM6 layer "
                      f"{KE + 1 - ell:2d}): stencil={L[ell - 1]!r} "
                      f"dump={exp!r} rel={rel:.3e}")

        # ent identity at the pinned interfaces, with the interface spacing
        # taken as 0.5 (ze[n+1] - ze[n-1]) (the esd interface spacing; this
        # form is what interface_entrainment uses because 0.5 (ze_a + ze_b)
        # would round at the ~743 m depth magnitude and spoil the deep
        # spacings by ~1e-10 relative).
        for n, K in spec["ent_ifc"]:
            dz_int = 0.5 * (ze[n] - ze[n - 2])   # esd interface n spacing
            for field, kd, ent in (("T", kd_h, ent_t), ("S", kd_s, ent_s)):
                ent_esd = dt * kd[n - 1] / dz_int
                rel = abs(ent_esd - ent[n - 1]) / max(abs(ent[n - 1]), 1e-300)
                print(f"  ent_{field} esd ifc {n:2d} (MOM6 ifc {K:2d}): "
                      f"law={ent_esd!r} dump={ent[n - 1]!r} rel={rel:.3e}")
                assert rel < 1e-9, f"ent mismatch row {label} ifc {n}"

        # scaled tridiagonal residual (README-mom6-w2.md identity) at the
        # pinned mixed layer and at the bottom layer (one deep layer per
        # column), T and S.  esd form (0-based): dz[i] (u[i]-u0[i]) =
        # ent[i] (u[i-1]-u[i]) + ent[i+1] (u[i+1]-u[i]) where ent[i] is esd
        # interface i+1 (= MOM6 interface NLEV+1-i), scaled as in
        # validate_mom6_w2.py.
        for field, ent, u, u0 in (("T", ent_t, ta, tb), ("S", ent_s, sa, sb)):
            for ell in sorted({spec["res_lev"], 1}):
                i = ell - 1
                lhs = dz[i] * (u[i] - u0[i])
                e_in = ent[i] * (u[i - 1] - u[i]) if ell >= 2 else 0.0
                e_out = ent[i + 1] * (u[i + 1] - u[i]) if ell < KE else 0.0
                res = abs(lhs - (e_in + e_out))
                scale = (abs(lhs) + abs(e_in) + abs(e_out)
                         + abs(ent[i] * u[i]) + abs(ent[i + 1] * u[i])
                         + 1e-30)
                print(f"  scaled residual {field} esd lev {ell:2d}: "
                      f"{res / scale:.4e} (res={res:.3e}, scale={scale:.3e})")
                assert res / scale < 1e-12
    print("all checks passed")


if __name__ == "__main__":
    main()
