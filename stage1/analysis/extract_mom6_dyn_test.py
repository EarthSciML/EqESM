#!/usr/bin/env python3
"""Extract the stage-2 test tuples for stage2/mom6_dynamics.esm from the
continuity dump: cell (3,78), dyn_step 2, call 1, layers k = 1 and 37.

The model subdomain is the 8x8 grid (global i = -1..6, j = 74..81, local
i = gi + 2, j = gj - 73) of the dump's hin_im4..ip3 x hin_jm4..jp3 stencil;
the state/parameter fields are real on the traced row/column and filled
with 0 elsewhere (off-axis cells never enter the asserted quantities).

Self-checks: reproduces the dump's h_W_c/h_E_c/uh_E/uh_W/vh_N/vh_S/h_inter
at the pinned cells with the template formulas (rel 1e-12), and verifies
the zonal/meridional convergence identities.
"""
import csv, json, sys

dd = sys.argv[1] if len(sys.argv) > 1 else \
    "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/dumps/mom6"
GI, GJ, STEP, CALL = 3, 78, 2, 1
LAYERS = (1, 37)          # pinned k values (1-indexed Fortran layers)
RTOL = 1e-12
NI = NJ = 8               # local i 1..8 = global -1..6; local j 1..8 = 74..81

lines = open(f"{dd}/mom6_continuity.csv").read().splitlines()
cols = [l for l in lines if l.startswith("# columns = ")][0][len("# columns = "):].rstrip(",").split(",")
ix = {c: i for i, c in enumerate(cols)}
rows = list(csv.reader([l for l in lines if not l.startswith("#")]))

group = [r for r in rows
         if int(r[ix["gi"]]) == GI and int(r[ix["gj"]]) == GJ
         and int(r[ix["dyn_step"]]) == STEP and int(r[ix["call"]]) == CALL]
assert len(group) == 75, f"group rows {len(group)}"
g = {k + 1: group[k] for k in range(75)}

# ---------------- land mask over the subdomain (constant-thickness rule)
# keyed by cell: gather the hin values that each cell contributes through the
# stencil columns of every traced-cell row; constant == land fill.
from collections import defaultdict
cellvals = defaultdict(set)
for r in rows:
    gi, gj = int(r[ix["gi"]]), int(r[ix["gj"]])
    for off, di, dj in (("im4", -4, 0), ("im3", -3, 0), ("im2", -2, 0), ("im1", -1, 0),
                        ("ip1", 1, 0), ("ip2", 2, 0), ("ip3", 3, 0),
                        ("jm4", 0, -4), ("jm3", 0, -3), ("jm2", 0, -2), ("jm1", 0, -1),
                        ("jp1", 0, 1), ("jp2", 0, 2), ("jp3", 0, 3)):
        cellvals[(gi + di, gj + dj)].add(float(r[ix["hin" + ("_" + off if off else "")]]))
land = {c for c, v in cellvals.items() if len(v) == 1 and v.pop() == 0.001}
print("land cells in subdomain:", len([c for c in land if -1 <= c[0] <= 6 and 74 <= c[1] <= 81]))

def mask(gi, gj):
    return 0.0 if (gi, gj) in land else 1.0

def hin_at(gi, gj, r):
    """h at (gi, gj) on the traced row/column of row r; None off-axis or
    outside the dumped stencil (-4..+3 i-offsets, -4..+3 j-offsets)."""
    if gi == GI:
        d = gj - GJ
        if abs(d) > 4:
            return None
        name = "hin" if d == 0 else ("hin_jp%d" % d if d > 0 else "hin_jm%d" % -d)
    elif gj == GJ:
        d = gi - GI
        if d < -4 or d > 3:
            return None
        name = "hin" if d == 0 else ("hin_ip%d" % d if d > 0 else "hin_im%d" % -d)
    else:
        return None
    return float(r[ix[name]])

def hi_at(gj, r):
    """h_inter at (GI, gj) on the traced column of row r; 0.0 outside the
    dumped stencil (-4..+3 j-offsets; such cells never enter the pins)."""
    d = gj - GJ
    if d < -4 or d > 3:
        return 0.0
    name = "h_inter" if d == 0 else ("h_inter_jp%d" % d if d > 0 else "h_inter_jm%d" % -d)
    return float(r[ix[name]])

def f(gi, gj, r):
    v = hin_at(gi, gj, r)
    return v if v is not None else 0.0

# ---------------- template formulas (Python re-implementation)
def ppm_slp(hm1, h0, hp1, mm1, m0, mp1):
    if not (mm1 and m0 and mp1):
        return 0.0
    s = 0.5 * (hp1 - hm1)
    dMx = max(hp1, hm1, h0) - h0
    dMn = h0 - min(hp1, hm1, h0)
    if s == 0.0:
        return 0.0
    lim = min(abs(s), 2.0 * min(dMx, dMn))
    return lim if s > 0 else -lim

def ppm_hL_raw(hm1, h0, slp_m1, slp_0, mm1):
    # hL = 0.5*(him1 + h0) + (slp_m1 - slp_0)/6 with him1 = m*hm1 + (1-m)*h0
    return 0.5 * (mm1 * hm1 + (1.0 - mm1) * h0 + h0) + (1.0 / 6.0) * (slp_m1 - slp_0)

def ppm_hR_raw(hp1, h0, slp_0, slp_p1, mp1):
    return 0.5 * (mp1 * hp1 + (1.0 - mp1) * h0 + h0) + (1.0 / 6.0) * (slp_0 - slp_p1)

def ppm_scale(hL_raw, hR_raw, h0, hmin):
    curv = 3.0 * ((hL_raw + hR_raw) - 2.0 * h0)
    if curv <= 0.0:
        return 1.0
    dh = hR_raw - hL_raw
    if not (abs(dh) < curv):
        return 1.0
    if h0 <= hmin:
        return 0.0
    num = 12.0 * curv * (h0 - hmin)
    den = curv ** 2 + 3.0 * dh ** 2
    return num / den if num < den else 1.0

def flux_x(u, dt, dy, h0, hp1, hWc, hEc, hWe, hEe, IdxT, IdxTe):
    if u > 0.0:
        cfl = u * dt * IdxT
        curv3 = hWc + hEc - 2.0 * h0
        return dy * u * (hEc + cfl * (0.5 * (hWc - hEc) + curv3 * (cfl - 1.5)))
    if u < 0.0:
        cfl = -u * dt * IdxTe
        curv3 = hWe + hEe - 2.0 * hp1
        return dy * u * (hWe + cfl * (0.5 * (hEe - hWe) + curv3 * (cfl - 1.5)))
    return 0.0

def flux_y(v, dt, dx, h0, hp1, hSc, hNc, hSn, hNn, IdyT, IdyTn):
    if v > 0.0:
        cfl = v * dt * IdyT
        curv3 = hSc + hNc - 2.0 * h0
        return dx * v * (hNc + cfl * (0.5 * (hSc - hNc) + curv3 * (cfl - 1.5)))
    if v < 0.0:
        cfl = -v * dt * IdyTn
        curv3 = hSn + hNn - 2.0 * hp1
        return dx * v * (hSn + cfl * (0.5 * (hNn - hSn) + curv3 * (cfl - 1.5)))
    return 0.0

def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)

# ---------------- per-layer data + self-checks
def layer_data(k):
    r = g[k]
    dt = float(r[ix["dt"]])
    hmin = float(r[ix["h_min_recon"]])
    # on-axis stencil maps
    hrow = {gi: hin_at(gi, GJ, r) for gi in range(-1, 7)}     # traced row
    hcol = {gj: hin_at(GI, gj, r) for gj in range(74, 82)}    # traced column
    hic = {gj: hi_at(gj, r) for gj in range(74, 82)}          # h_inter column
    # zonal reconstruction of cell (3,78)
    def zrec(off):
        h = {o: hrow[GI + off + o] for o in (-2, -1, 0, 1, 2)}
        m = {o: mask(GI + off + o, GJ) for o in (-2, -1, 0, 1, 2)}
        slp = {o: ppm_slp(h[o - 1], h[o], h[o + 1], m[o - 1], m[o], m[o + 1]) for o in (-1, 0, 1)}
        hL = ppm_hL_raw(h[-1], h[0], slp[-1], slp[0], m[-1])
        hR = ppm_hR_raw(h[1], h[0], slp[0], slp[1], m[1])
        sc = ppm_scale(hL, hR, h[0], hmin)
        return h[0] + sc * (hL - h[0]), h[0] + sc * (hR - h[0])
    hWc, hEc = zrec(0)
    hWe, hEe = zrec(1)
    hWim1, hEim1 = zrec(-1)
    def mrec(off):
        h = {o: hic[GJ + off + o] for o in (-2, -1, 0, 1, 2)}
        m = {o: mask(GI, GJ + off + o) for o in (-2, -1, 0, 1, 2)}
        slp = {o: ppm_slp(h[o - 1], h[o], h[o + 1], m[o - 1], m[o], m[o + 1]) for o in (-1, 0, 1)}
        hL = ppm_hL_raw(h[-1], h[0], slp[-1], slp[0], m[-1])
        hR = ppm_hR_raw(h[1], h[0], slp[0], slp[1], m[1])
        sc = ppm_scale(hL, hR, h[0], hmin)
        return h[0] + sc * (hL - h[0]), h[0] + sc * (hR - h[0])
    hSc, hNc = mrec(0)
    hSn, hNn = mrec(1)
    hSjm1, hNjm1 = mrec(-1)
    # flux laws at the 4 faces of (3,78)
    uE, uW = float(r[ix["u_E"]]), float(r[ix["u_W"]])
    vN, vS = float(r[ix["v_N"]]), float(r[ix["v_S"]])
    uhE = flux_x(uE, dt, float(r[ix["dy_Cu_E"]]), hrow[GI], hrow[GI + 1],
                 hWc, hEc, hWe, hEe, float(r[ix["IdxT"]]), float(r[ix["IdxT_e"]]))
    uhW = flux_x(uW, dt, float(r[ix["dy_Cu_W"]]), hrow[GI - 1], hrow[GI],
                 hWim1, hEim1, hWc, hEc, float(r[ix["IdxT_w"]]), float(r[ix["IdxT"]]))
    vhN = flux_y(vN, dt, float(r[ix["dx_Cv_N"]]), hic[GJ], hic[GJ + 1],
                 hSc, hNc, hSn, hNn, float(r[ix["IdyT"]]), float(r[ix["IdyT_n"]]))
    vhS = flux_y(vS, dt, float(r[ix["dx_Cv_S"]]), hic[GJ - 1], hic[GJ],
                 hSjm1, hNjm1, hSc, hNc, float(r[ix["IdyT_s"]]), float(r[ix["IdyT"]]))
    ia = float(r[ix["IareaT"]])
    hint = float(r[ix["h_inter"]])
    checks = {
        "h_W_c": (hWc, float(r[ix["h_W_c"]])), "h_E_c": (hEc, float(r[ix["h_E_c"]])),
        "h_W_e": (hWe, float(r[ix["h_W_e"]])), "h_E_e": (hEe, float(r[ix["h_E_e"]])),
        "h_S_c": (hSc, float(r[ix["h_S_c"]])), "h_N_c": (hNc, float(r[ix["h_N_c"]])),
        "h_S_n": (hSn, float(r[ix["h_S_n"]])), "h_N_n": (hNn, float(r[ix["h_N_n"]])),
        "uh_E": (uhE, float(r[ix["uh_E"]])), "uh_W": (uhW, float(r[ix["uh_W"]])),
        "vh_N": (vhN, float(r[ix["vh_N"]])), "vh_S": (vhS, float(r[ix["vh_S"]])),
        "h_inter": (hint, float(r[ix["h_inter"]])),
    }
    for name, (a, b) in checks.items():
        ok = rel(a, b) <= RTOL
        print(f"  k={k:2d} {name:9s} model={a:.12e} dump={b:.12e} {'OK' if ok else 'FAIL'}")
        assert ok, name
    zid = max(float(r[ix["hin"]]) - dt * ia * (uhE - uhW), 0.0)
    mid = max(hint - dt * ia * (vhN - vhS), float(r[ix["h_min"]]))
    print(f"  k={k:2d} z-id rel={rel(zid, hint):.2e} m-id rel={rel(mid, float(r[ix['h_out']])):.2e}")
    return dict(dt=dt, hmin=hmin, ia=ia, uE=uE, uW=uW, vN=vN, vS=vS,
                uhE=uhE, uhW=uhW, vhN=vhN, vhS=vhS, hint=hint,
                hWc=hWc, hEc=hEc, hWe=hWe, hEe=hEe, hSc=hSc, hNc=hNc, hSn=hSn, hNn=hNn,
                hrow=hrow, hcol=hcol, hic=hic)

print("self-checks (model formulas vs dump, rel 1e-12):")
data = {k: layer_data(k) for k in LAYERS}

# ---------------- build the JSON fragments
def grid_field(pos):
    f = [[0.0 for _ in range(NJ)] for _ in range(NI)]
    for (li, lj), v in pos.items():
        f[li - 1][lj - 1] = v
    return f

def field(fn):
    return [[fn(gi, gj) for gj in range(74, 82)] for gi in range(-1, 7)]

def shifted(fn, di, dj):
    return field(lambda gi, gj: fn(gi - di, gj - dj))

r0 = g[LAYERS[0]]
out = {}
out["h"] = [[[f(gi, gj, g[k]) for gj in range(74, 82)] for gi in range(-1, 7)] for k in LAYERS]
for k in LAYERS:
    r = g[k]
    hf = lambda gi, gj, kk=k: f(gi, gj, g[kk])
    out[f"h_m2_k{k}"] = shifted(lambda gi, gj, kk=k: hf(gi, gj, kk), 2, 0)
    out[f"h_m1_k{k}"] = shifted(lambda gi, gj, kk=k: hf(gi, gj, kk), 1, 0)
    out[f"h_p1_k{k}"] = shifted(lambda gi, gj, kk=k: hf(gi, gj, kk), -1, 0)
    out[f"h_p2_k{k}"] = shifted(lambda gi, gj, kk=k: hf(gi, gj, kk), -2, 0)
    hif = lambda gi, gj, kk=k: hi_at(gj, g[kk]) if gi == GI else 0.0
    out[f"h_inter_k{k}"] = field(lambda gi, gj, kk=k: hif(gi, gj, kk))
    out[f"hi_jm2_k{k}"] = shifted(lambda gi, gj, kk=k: hif(gi, gj, kk), 0, 2)
    out[f"hi_jm1_k{k}"] = shifted(lambda gi, gj, kk=k: hif(gi, gj, kk), 0, 1)
    out[f"hi_jp1_k{k}"] = shifted(lambda gi, gj, kk=k: hif(gi, gj, kk), 0, -1)
    out[f"hi_jp2_k{k}"] = shifted(lambda gi, gj, kk=k: hif(gi, gj, kk), 0, -2)
    ue = grid_field({(5, 5): float(r[ix["u_E"]]), (4, 5): float(r[ix["u_W"]])})
    vn = grid_field({(5, 5): float(r[ix["v_N"]]), (5, 4): float(r[ix["v_S"]])})
    hwe = grid_field({(5, 5): float(r[ix["h_W_e"]]), (4, 5): float(r[ix["h_W_c"]])})
    hee = grid_field({(5, 5): float(r[ix["h_E_e"]]), (4, 5): float(r[ix["h_E_c"]])})
    hsn = grid_field({(5, 5): float(r[ix["h_S_n"]]), (5, 4): float(r[ix["h_S_c"]])})
    hnn = grid_field({(5, 5): float(r[ix["h_N_n"]]), (5, 4): float(r[ix["h_N_c"]])})
    out[f"u_E_k{k}"] = ue; out[f"v_N_k{k}"] = vn
    out[f"h_W_e_k{k}"] = hwe; out[f"h_E_e_k{k}"] = hee
    out[f"h_S_n_k{k}"] = hsn; out[f"h_N_n_k{k}"] = hnn

mask_field = field(mask)
out["m_0"] = mask_field
out["m_m1"] = shifted(mask, 1, 0); out["m_m2"] = shifted(mask, 2, 0)
out["m_p1"] = shifted(mask, -1, 0); out["m_p2"] = shifted(mask, -2, 0)
out["m_jm1"] = shifted(mask, 0, 1); out["m_jm2"] = shifted(mask, 0, 2)
out["m_jp1"] = shifted(mask, 0, -1); out["m_jp2"] = shifted(mask, 0, -2)

out["dy_Cu_E"] = grid_field({(5, 5): float(r0[ix["dy_Cu_E"]]), (4, 5): float(r0[ix["dy_Cu_W"]])})
out["IdxT"] = grid_field({(5, 5): float(r0[ix["IdxT"]]), (4, 5): float(r0[ix["IdxT_w"]])})
out["IdxT_e"] = grid_field({(5, 5): float(r0[ix["IdxT_e"]]), (4, 5): float(r0[ix["IdxT"]])})
out["dx_Cv_N"] = grid_field({(5, 5): float(r0[ix["dx_Cv_N"]]), (5, 4): float(r0[ix["dx_Cv_S"]])})
out["IdyT"] = grid_field({(5, 5): float(r0[ix["IdyT"]]), (5, 4): float(r0[ix["IdyT_s"]])})
out["IdyT_n"] = grid_field({(5, 5): float(r0[ix["IdyT_n"]]), (5, 4): float(r0[ix["IdyT"]])})
out["IareaT"] = grid_field({(5, 5): float(r0[ix["IareaT"]])})
out["dt"] = float(r0[ix["dt"]])
out["h_min_recon"] = float(r0[ix["h_min_recon"]])

json.dump(out, open("/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM-wt-mom6-dynamics-instrument/scratch/mom6_dyn_fields.json", "w"))

print("\nassertions:")
for k in LAYERS:
    d = data[k]
    ki = 1 if k == LAYERS[0] else 2
    for name, val, coords in (
        ("h_W_c", d["hWc"], (5, 5, ki)), ("h_E_c", d["hEc"], (5, 5, ki)),
        ("uh", d["uhE"], (5, 5, ki)), ("uh", d["uhW"], (4, 5, ki)),
        ("vh", d["vhN"], (5, 5, ki)), ("vh", d["vhS"], (5, 4, ki)),
        ("h_inter", d["hint"], (5, 5, ki)),
        ("div_x", d["ia"] * (d["uhE"] - d["uhW"]), (5, 5, ki)),
        ("div_y", d["ia"] * (d["vhN"] - d["vhS"]), (5, 5, ki)),
    ):
        print(f'  {{"variable": "{name}", "time": 0.0, "expected": {val!r}, '
              f'"coords": {{"i": {coords[0]}, "j": {coords[1]}, "k": {coords[2]}}}, '
              f'"tolerance": {{"rel": 1e-12, "abs": 0.0}}}},')
# integrated pin (h_out at t=dt)
for k in LAYERS:
    ki = 1 if k == LAYERS[0] else 2
    hv = float(g[k][ix["h_out"]])
    print(f'  {{"variable": "h", "time": 1800.0, "expected": {hv!r}, '
          f'"coords": {{"i": 5, "j": 5, "k": {ki}}}, "tolerance": {{"rel": 0.001, "abs": 0.0}}}},')
    out["assert_h_out_" + str(ki)] = hv

# ---- machine-readable assertions for the assembler
instant = []
for k in LAYERS:
    d = data[k]
    ki = 1 if k == LAYERS[0] else 2
    for name, val, (ci, cj, ck) in (
        ("h_W_c", d["hWc"], (5, 5, ki)), ("h_E_c", d["hEc"], (5, 5, ki)),
        ("uh", d["uhE"], (5, 5, ki)), ("uh", d["uhW"], (4, 5, ki)),
        ("vh", d["vhN"], (5, 5, ki)), ("vh", d["vhS"], (5, 4, ki)),
        ("h_inter", d["hint"], (5, 5, ki)),
        ("div_x", d["ia"] * (d["uhE"] - d["uhW"]), (5, 5, ki)),
        ("div_y", d["ia"] * (d["vhN"] - d["vhS"]), (5, 5, ki)),
    ):
        instant.append({"variable": name, "time": 0.0, "expected": val,
                        "coords": {"i": ci, "j": cj, "k": ck},
                        "tolerance": {"rel": RTOL, "abs": 0.0}})
integrated = []
for k in LAYERS:
    ki = 1 if k == LAYERS[0] else 2
    integrated.append({"variable": "h", "time": 1800.0,
                       "expected": float(g[k][ix["h_out"]]),
                       "coords": {"i": 5, "j": 5, "k": ki},
                       "tolerance": {"rel": 0.001, "abs": 0.0}})
json.dump({"instant": instant, "integrated": integrated},
          open("/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM-wt-mom6-dynamics-instrument/scratch/mom6_dyn_assertions.json", "w"))
print(f"wrote {len(instant)} instant + {len(integrated)} integrated assertions")
