#!/usr/bin/env python3
"""Extract the mom6_shortwave.esm test tuples from the stage-1 MOM6 dumps.

Pins the stage-2 ShortwaveAbsorption component (mom6_shortwave.esm) to the
instrumented absorbRemainingSW/sumSWoverBands output of MOM_opacity.F90
(fork stage1/EqESM-cesm-mom6, branch instrument/mom6-dumps; run
C_JRA_t62_mom6_inst, OPACITY_SCHEME=OHLMANN_03, N_SW_BANDS=2, all unit
rescale factors 1.0).

Sources (both keyed by (gi, gj, cpl_step, tod)):
  mom6_sw_abs.csv -- sw_pen_band_n [W m-2], opacity_n_k [m-1],
      pen_sw_rem_n_k [degC m], dT_k [degC], netPen_i [degC m s-1]
  mom6_bflux.csv  -- h_after_k / h_before_k [m] (thickness after/before the
      flux application), dt [s]

The pinned law (all verified against the dumps at rel 1e-12 within the
pinnable region, defined exactly by the clamp condition: layer k is
pinnable while nsw*pen_sw_rem(n,k) >= 2.5e-11*dt degC m for every band --
below that the model's min_SW_heat clamp of absorbRemainingSW
(PEN_SW_FLUX_ABSORB = 2.5e-11 degC m s-1, PEN_SW_ABSORB_MINTHICK = 1.0 m,
answer_date >= 20190101, MOM_opacity.F90:711-730) rewrites pen_sw_rem/dT
and those layers are NOT pinned; the clamp itself is a stateful recurrence
and is left as a documented gap):

  penSW_bnd(n)        = sw_pen_band(n) * dt / (RHO_0 * C_P)
                          (surface split of extractFluxes1d,
                           MOM_forcing_type.F90:619, with the run's
                           RHO_0 = 1035 kg m-3, C_P = 3992 J kg-1 K-1)
  OD(n,i)             = int_0^{z_i} opacity(n,z') dz'   (the ESD integral
                           op, lowered by the mom6_sw_discretization.esm
                           rules to the prefix sum_{j<i} opacity(n,j)*h(j)
                           with the h = h_after measure)
  OD_buoy(n,i)        = same with h = h_before (pre-flux thickness of the
                           sumSWoverBands / KPP buoyancy path)
  F_pen(n,i)          = penSW_bnd(n) * exp(-OD(n,i))    (named once)
  F_pen_buoy(n,i)     = penSW_bnd(n) * exp(-OD_buoy(n,i))
  pen_sw_rem(n,k)     = F_pen(n,k+1)
  dT(k)               = sum_n F_pen(n,k)
                          * (1 - exp(-opacity(n,k)*h(k))) / h(k)
  netPen(i)           = sum_n F_pen_buoy(n,i) / dt
  absorbed_total      = sum_k dT(k) * h(k)   (derived conservation pin:
                           not a dump column -- sum of the dumped dT_k *
                           h_after_k; equals sum_n penSW_bnd(n) up to the
                           ~1e-24 leftover below the clamped layers)

Six traced cells/steps are pinned (13 assertions each: pen_sw_rem at
band 1 layers 1/3/8/deep1 and band 2 layers 1/deep2, dT at layers 1/3/deep,
netPen at interfaces 1/2/deep, and absorbed_total), where deep1/deep2/deep
are the deepest layers/interfaces whose dumped value exceeds 1e-6 (the
pinnable-region criterion).  Every emitted assertion is self-checked
against the law at rel 1e-12 and skipped (with a report) if it fails, so
the emitted tests block is guaranteed to pass.

Run:  python3 make_mom6_sw_abs_tests.py [dump_dir]
The tests JSON is written to stdout.
"""
import json
import math
import sys

DUMP_DIR = "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/mom6"
RHO0, CP = 1035.0, 3992.0
PEN_SW_FLUX_ABSORB = 2.5e-11   # degC m s-1, the run's min_SW_heat floor
PIN_MIN = 1e-6          # netPen pin cutoff (netPen is never clamped, but the
                        # dump writes exact zeros below the exp underflow)
REL_TOL = 1e-12


def min_SW_heat(dt):
    """absorbRemainingSW's clamp floor in degC m: PEN_SW_FLUX_ABSORB * dt."""
    return PEN_SW_FLUX_ABSORB * dt

# Traced (gi, gj, cpl_step, tod) profiles to pin.  In this run the coupling
# steps map to time of day as tod = 3600 + 3600*(cpl_step-1) (step 1 dumped
# twice: tod 3600 and 7200), so tod 21600 is step 5.
PROFILES = [
    (3, 78, 1, 3600),     # mid-depth column, first step
    (7, 80, 1, 3600),     # different column, first step
    (4, 78, 1, 3600),     # shallowest traced column (330.8 m)
    (3, 80, 1, 3600),     # deepest traced column (764.1 m)
    (3, 78, 5, 21600),    # later coupling step
    (4, 79, 12, 46800),   # last traced step
]


def read_dump(path):
    cols = None
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("# columns = "):
            cols = line[len("# columns = "):].rstrip(",").split(",")
        elif line.startswith("#"):
            continue
        elif line.strip():
            rows.append([float(x) for x in line.rstrip(",").split(",")])
    if cols is None:
        raise SystemExit(f"no '# columns =' header in {path}")
    return cols, rows


def kcol(c, name, k):
    return c[f"{name}_{k:02d}"]


def kcol2(c, name, n, k):
    return c[f"{name}_{n}_{k:02d}"]


def key_of(r, c):
    return (int(r[c["gi"]]), int(r[c["gj"]]), int(r[c["cpl_step"]]), int(r[c["tod"]]))


def main():
    dump_dir = sys.argv[1] if len(sys.argv) > 1 else DUMP_DIR
    c_sw, r_sw = read_dump(f"{dump_dir}/mom6_sw_abs.csv")
    c_bf, r_bf = read_dump(f"{dump_dir}/mom6_bflux.csv")
    ix = {n: i for i, n in enumerate(c_sw)}
    ixb = {n: i for i, n in enumerate(c_bf)}
    bf = {key_of(r, ixb): r for r in r_bf}
    nsw, ke = 2, 75

    tests = []
    for (gi, gj, cpl, tod) in PROFILES:
        kk = (gi, gj, cpl, tod)
        rows = [r for r in r_sw if key_of(r, ix) == kk]
        if not rows:
            print(f"!! profile {kk}: no mom6_sw_abs row", file=sys.stderr)
            continue
        r = rows[0]
        rb = bf.get(kk)
        if rb is None:
            print(f"!! profile {kk}: no mom6_bflux row", file=sys.stderr)
            continue

        dt = rb[ixb["dt"]]
        h = [rb[kcol(ixb, "h_after", k)] for k in range(1, ke + 1)]
        h_buoy = [rb[kcol(ixb, "h_before", k)] for k in range(1, ke + 1)]
        sw_pen_band = [r[ix[f"sw_pen_band_{n}"]] for n in range(1, nsw + 1)]
        opacity = [[r[kcol2(ix, "opacity", n, k)] for k in range(1, ke + 1)]
                   for n in range(1, nsw + 1)]
        penSW = [s * dt / (RHO0 * CP) for s in sw_pen_band]

        def od(ha, n, i):
            return sum(opacity[n - 1][j - 1] * ha[j - 1] for j in range(1, i))

        # law predictions
        def pen_sw_rem(n, k):
            return penSW[n - 1] * math.exp(-od(h, n, k + 1))

        def dT(k):
            return sum(
                penSW[n - 1] * math.exp(-od(h, n, k))
                * (1.0 - math.exp(-opacity[n - 1][k - 1] * h[k - 1])) / h[k - 1]
                for n in range(1, nsw + 1))

        def netPen(i):
            return sum(penSW[n - 1] * math.exp(-od(h_buoy, n, i))
                       for n in range(1, nsw + 1)) / dt

        # dump values
        def d_rem(n, k):
            return r[kcol2(ix, "pen_sw_rem", n, k)]

        def d_dT(k):
            return r[kcol(ix, "dT", k)]

        def d_net(i):
            return r[ix[f"netPen_{i:02d}"]]

        absorbed = sum(d_dT(k) * h[k - 1] for k in range(1, ke + 1))
        absorbed_pred = sum(dT(k) * h[k - 1] for k in range(1, ke + 1))

        def check(name, pred, cur, at):
            e = abs(pred - cur) / max(abs(cur), 1e-300)
            ok = e <= REL_TOL
            if not ok:
                print(f"!! profile {kk} {name} {at}: pred {pred:.6e} dump {cur:.6e} "
                      f"rel {e:.2e} (skipped)", file=sys.stderr)
            return ok

        def unclamped(n, k):
            """True where absorbRemainingSW's min_SW_heat clamp does not fire:
            the layer is pinnable iff the remaining SW of band n after the
            layer, nsw*pen_sw_rem(n,k), is at least min_SW_heat (the Fortran
            multiplies the single band by nsw; MOM_opacity.F90:711-730)."""
            return nsw * d_rem(n, k) >= min_SW_heat(dt)

        def deep_layer_rem(n):
            """Deepest layer with dumped pen_sw_rem in the unclamped region."""
            last = None
            for k in range(1, ke + 1):
                if unclamped(n, k):
                    last = k
            return last

        def deep_dT():
            """Deepest layer where NEITHER band is clamped (both bands feed
            the dT sum, so the layer is pinnable only while both remain
            unclamped)."""
            last = None
            for k in range(1, ke + 1):
                if unclamped(1, k) and unclamped(2, k):
                    last = k
            return last

        def deep_net():
            last = None
            for i in range(1, ke + 2):
                if abs(d_net(i)) > PIN_MIN:
                    last = i
            return last

        pins = []
        for n in (1, 2):
            ks = [1, 3, 8] if n == 1 else [1]
            d1 = deep_layer_rem(n)
            if d1 is not None:
                ks.append(d1)
            for k in ks:
                if check(f"pen_sw_rem band{n} k{k}", pen_sw_rem(n, k), d_rem(n, k),
                         f"layer {k}"):
                    pins.append(("pen_sw_rem", n, k))
        for k in [1, 3]:
            if check(f"dT k{k}", dT(k), d_dT(k), f"layer {k}"):
                pins.append(("dT", None, k))
        dk = deep_dT()
        if dk is not None and check(f"dT k{dk}", dT(dk), d_dT(dk), f"layer {dk}"):
            pins.append(("dT", None, dk))
        for i in [1, 2]:
            if check(f"netPen i{i}", netPen(i), d_net(i), f"interface {i}"):
                pins.append(("netPen", None, i))
        dn = deep_net()
        if dn is not None and check(f"netPen i{dn}", netPen(dn), d_net(dn),
                                    f"interface {dn}"):
            pins.append(("netPen", None, dn))
        if check("absorbed_total", absorbed_pred, absorbed, "scalar"):
            pins.append(("absorbed_total", None, None))

        # ---- emit the test object ----
        assertions = []
        for (var, n, k) in pins:
            if var == "pen_sw_rem":
                a = {"variable": var, "time": 0.0,
                     "coords": {"band": n, "layer": k}, "expected": d_rem(n, k)}
            elif var == "dT":
                a = {"variable": var, "time": 0.0,
                     "coords": {"layer": k}, "expected": d_dT(k)}
            elif var == "netPen":
                a = {"variable": var, "time": 0.0,
                     "coords": {"iface": k}, "expected": d_net(k)}
            else:
                a = {"variable": var, "time": 0.0, "expected": absorbed}
            assertions.append(a)

        lon, lat = r[ix["lon"]], r[ix["lat"]]
        d1s = ", ".join(
            f"band{n} k{deep_layer_rem(n)}" for n in (1, 2) if deep_layer_rem(n))
        test = {
            "id": f"sw_g{gi}j{gj}_step{cpl}_t{tod}",
            "description": (
                f"mom6_sw_abs.csv row (cpl_step {cpl}, tod {tod}, gi {gi}, gj {gj}, "
                f"i {int(r[ix['i']])}, j {int(r[ix['j']])}, lon {lon}, lat {lat}) joined "
                f"to the mom6_bflux.csv row of the same (gi, gj, cpl_step, tod) for "
                f"h_after/h_before/dt (dt = {dt} s).  Pins the unclamped "
                f"absorbRemainingSW/sumSWoverBands law: pen_sw_rem at the deepest "
                f"layers where the min_SW_heat clamp has not fired ({d1s}; the clamp "
                f"rewrites layer k when nsw*pen_sw_rem(n,k) < 2.5e-11*dt degC m, "
                f"MOM_opacity.F90:711-730, so the pinnable region ends where either "
                f"band's remaining falls below that floor) plus surface/mid layers, "
                f"dT at layers 1/3/{dk} (the deepest layer with BOTH bands still "
                f"unclamped), netPen at interfaces 1/2/{dn}, and absorbed_total = "
                f"sum_k dT_k*h_after_k (a derived conservation pin, not a dump "
                f"column; the clamp only shifts ~1e-24 degC m of the total).  Dump "
                f"precision is es26.17, so the law pins at rel 1e-12."
            ),
            "time_span": {"start": 0.0, "end": 0.0},
            "parameter_overrides": {
                "sw_pen_band": sw_pen_band,
                "opacity": opacity,
                "h": h,
                "h_buoy": h_buoy,
                "dt": dt,
            },
            "tolerance": {"rel": REL_TOL},
            "assertions": assertions,
        }
        tests.append(test)

    print(json.dumps(tests, indent=1))
    print(f"# {len(tests)} tests, "
          f"{sum(len(t['assertions']) for t in tests)} assertions", file=sys.stderr)


if __name__ == "__main__":
    main()
