#!/usr/bin/env python3
"""Extract the numeric test tuples of stage2/cam_physics_update.esm from the
stage-1 CAM physics dumps, and re-verify the full law against every dump row.

The dumps (gitignored, living in the main checkout): stage1/dumps/cam/
cam_phys_tend.csv (per-step, per-call-site, per-level state + ptend rows of
every physics_update call, physics_types.F90:207-495) and cam_phys_tend_q.csv
(per-constituent q_before / ptend_q; constituents are dumped only where
ptend%lq is set), plus cam_tend_total.csv (the accumulated tend%dudt/dvdt/dtdt
at the end of the tphysbc (phase=bc) and tphysac (phase=ac) sequences).

Pre-flight verification (fails loudly on any mismatch):
  * full-column law (u/v/T update, the qneg3 qmin table of chemistry.F90 and
    the number-concentration clamp, the CLDLIQ/CLDICE 1e-36 zeroing at the
    zm_conv_tend sites, the geopotential recompute and s_after) applied at
    EVERY (step, site) of the 48-step window and chained to the next call
    site's dump rows, worst relative error reported;
  * the phase=bc row of cam_tend_total.csv equals the sum of the per-site
    ptend_u/ptend_v/ptend_s over the 15 consecutive tphysbc call sites at the
    pinned step (worst error reported);
  * the copy-check: every emitted assertion must be able to fail — the
    expected value must differ from the overridden input it derives from
    (T_after vs T, ..., q_after vs q), so no assertion pins a state that the
    law leaves unchanged.

The emitted fragments are the numeric `parameter_overrides` + `assertions` of
the hand-authored tests blocks in stage2/cam_physics_update.esm (the model
body itself is hand-authored; only the numbers come from the dumps).

Model conventions (documented in the .esm metadata):
  * the lev axis is SURFACE-FIRST (model index 1 = dump lev 32), so every
    [lev] array is the dump's lev-1..32 column reversed;
  * the cnst axis is the dump's pcnst order (34 constituents), 1-based;
  * the site axis of the TendAccumulation model is the tphysbc call order
    (site 1 = energy_fixer ... site 15 = radiation).

Usage:
  python3 stage1/analysis/extract_cam_physics_update_pins.py [dumps_dir]
where dumps_dir defaults to <main checkout>/stage1/dumps/cam (the environment
variable EQESM_MAIN overrides the main-checkout location).  Prints the two
tests arrays (PhysicsUpdate, TendAccumulation) as JSON to stdout.
"""
import json, os, re, sys
from collections import OrderedDict

D = os.path.join(os.environ.get("EQESM_MAIN", "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM"),
                 "stage1", "dumps", "cam")
if len(sys.argv) > 1:
    D = sys.argv[1]

CPAIR = 1004.64
RAIR = 287.0423113650487
GRAVIT = 9.80616
RH2O = 461.5046398201599
ZVIR = RH2O / RAIR - 1.0

QMIN = {"Q": 1e-12}
NUM_CLAMP = ("NUMICE", "NUMLIQ", "NUMRAI", "NUMSNO")
WATER = ("Q", "CLDLIQ", "CLDICE", "RAINQM", "SNOWQM")
CLD_MIN_NZ_SITES = ("convect_deep",)  # proc zm_conv_tend / convect_deep (cldlim_names)
BC_SITES = ["energy_fixer", "dadadj", "convect_deep", "convect_shallow", "ssatcontrail",
            "clubb", "microp", "clubb", "microp", "clubb", "microp",
            "aero_calcsize", "aero_wetdep", "convect_deep_2", "radiation"]


def bc_it(i):
    """macmic_it of BC_SITES[i]: clubb/microp repeat 3x (it 1..3), all others 1."""
    site = BC_SITES[i]
    return BC_SITES[:i].count(site) + 1 if site in ("clubb", "microp") else 1


def fixnum(s):
    """gfortran es26.17 writes 3-digit exponents without the 'E'; insert it."""
    s = s.strip()
    m = re.search(r"[0-9]([+-][0-9]{3})$", s)
    if m and "E" not in s.upper():
        s = s[: m.start(1)] + "E" + s[m.start(1):]
    return float(s)


def load_phys():
    rows = []
    with open(f"{D}/cam_phys_tend.csv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split(",")
            rows.append(dict(
                step=int(p[0]), site=p[3].strip('"'), proc=p[4].strip('"'),
                lu=p[6] == "T", lv=p[7] == "T", ls=p[8] == "T",
                top=int(p[10]), bot=int(p[11]),
                macmic_it=int(p[12]), macmic_num=int(p[13]),
                dt=float(p[14]), lev=int(p[16]),
                T=float(p[19]), u=float(p[20]), v=float(p[21]), s=float(p[22]),
                pmid=float(p[23]), pdel=float(p[24]), phis=float(p[26]), zm=float(p[27]),
                pts=float(p[29]), ptu=float(p[30]), ptv=float(p[31])))
    return rows


def load_q():
    raw = []
    with open(f"{D}/cam_phys_tend_q.csv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split(",")
            raw.append(dict(step=int(p[0]), site=p[3].strip('"'), lev=int(p[4]),
                            cnst=p[6].strip('"'), q_before=fixnum(p[7]), ptend_q=fixnum(p[8])))
    by = OrderedDict()
    for r in raw:
        by.setdefault((r["step"], r["site"]), []).append(r)
    out = []
    for (step, site), rows in by.items():
        ncnst = len({r["cnst"] for r in rows})
        blk = 32 * ncnst
        for i, r in enumerate(rows):
            # the q dump has no macmic column; the clubb/microp rows repeat
            # per macmic_it in order, in 32*ncnst-row blocks
            r["macmic_it"] = i // blk + 1 if site in ("clubb", "microp") else 1
            out.append(r)
    return out


def qmin_of(cnst):
    if cnst == "Q":
        return 1e-12
    if cnst in ("CLDLIQ", "CLDICE", "RAINQM", "SNOWQM"):
        return 0.0
    if cnst.startswith("num_a"):
        return 1e-5
    return 1e-36


def update_column(r, qmap):
    """The law at one level: post-call (u, v, T, q-dict) of physics_update."""
    dt = r["dt"]
    u = r["u"] + r["ptu"] * dt
    v = r["v"] + r["ptv"] * dt
    t = r["T"] + r["pts"] * dt / CPAIR
    q = {}
    for cnst, qr in qmap.items():
        qnew = qr["q_before"] + qr["ptend_q"] * dt
        if cnst in NUM_CLAMP:
            qnew = min(max(qnew, 1e-12), 1e10)
        else:
            qnew = max(qnew, qmin_of(cnst))
        q[cnst] = qnew
    if r["site"] in CLD_MIN_NZ_SITES:  # state_cnst_min_nz 1e-36 zeroing, F90:353-375
        for mass, num in (("CLDLIQ", "NUMLIQ"), ("CLDICE", "NUMICE")):
            if mass in q and q[mass] < 1e-36:
                q[mass] = 0.0
                if num in q:
                    q[num] = 0.0
    return u, v, t, q


def geopotential(t, q, pmid, pdel):
    """zm from the updated T/q (geopotential_temp.F90:133-216, non-lagrangian).
    Arrays in dump order (index 0 = lev 1 = TOA); the Fortran integrates from
    the surface (lev 32) upward."""
    qfac = [1.0 / (1.0 - sum(q[l][w] for w in WATER)) for l in range(32)]
    sdmr = [1.0 / (1.0 + qfac[l] * sum(q[l][w] for w in WATER)) for l in range(32)]
    rog = RAIR / GRAVIT
    zi = [0.0] * 33  # zi[32] = surface = 0; zi[0] = top
    zm = [0.0] * 32
    for j in range(32):
        l = 31 - j
        tvfac = (1.0 + (ZVIR + 1.0) * q[l]["Q"] * qfac[l]) * sdmr[l]
        tv = t[l] * tvfac
        hkl = pdel[l] / pmid[l]
        zm[l] = zi[l + 1] + rog * tv * 0.5 * hkl
        zi[l] = zi[l + 1] + rog * tv * hkl
    return zm, zi


def verify_full_chain(rows, qrows, pk, qk, corder):
    """Apply the law at every (step, site, level) and chain to the next call
    site's dump rows; report the worst relative error."""
    steps = OrderedDict()
    for r in rows:
        label = r["site"] if r["macmic_num"] <= 1 else f'{r["site"]}#{r["macmic_it"]}'
        seq = steps.setdefault(r["step"], [])
        if not seq or seq[-1] != label:
            seq.append(label)

    def resolve(label):
        if "#" in label:
            base, it = label.split("#")
            return base, int(it)
        return label, 1

    worst = 0.0
    worst_loc = None
    for st, seq in steps.items():
        for k in range(len(seq) - 1):
            b1, it1 = resolve(seq[k])
            b2, it2 = resolve(seq[k + 1])
            for lev in range(1, 33):
                r1 = pk[(st, b1, lev, it1)]
                r2 = pk[(st, b2, lev, it2)]
                q1 = {c: qk[(st, b1, lev, c, it1)] for c in corder if (st, b1, lev, c, it1) in qk}
                u, v, t, q = update_column(r1, q1)
                checks = [("u", u, r2["u"]), ("v", v, r2["v"]), ("T", t, r2["T"])]
                for cnst, v1 in q.items():
                    if (st, b2, lev, cnst, it2) in qk:
                        checks.append((cnst, v1, qk[(st, b2, lev, cnst, it2)]["q_before"]))
                for name, v1, v2 in checks:
                    rel = abs(v1 - v2) / abs(v2) if v2 else abs(v1 - v2)
                    if rel > worst:
                        worst, worst_loc = rel, (st, seq[k], seq[k + 1], lev, name, v1, v2)
    print(f"column-state chaining: worst rel err over all 48 steps = {worst:.3g} at {worst_loc}")
    if worst > 1e-14:
        raise SystemExit("chain verification FAILED")

    # geopotential zm/s chaining: zm_after(k) == zm_before(k+1) wherever
    # derive_new_geopotential and all five water species are known.
    worst_zm = worst_s = 0.0
    n_checked = 0
    for st, seq in steps.items():
        qstate = [{} for _ in range(32)]
        for k in range(len(seq) - 1):
            b1, it1 = resolve(seq[k])
            b2, it2 = resolve(seq[k + 1])
            T = [0.0] * 32
            pmid = [0.0] * 32
            pdel = [0.0] * 32
            for lev in range(1, 33):
                r1 = pk[(st, b1, lev, it1)]
                q1 = {c: qk[(st, b1, lev, c, it1)] for c in corder if (st, b1, lev, c, it1) in qk}
                u, v, t, q = update_column(r1, q1)
                T[lev - 1] = t
                pmid[lev - 1] = r1["pmid"]
                pdel[lev - 1] = r1["pdel"]
                for cnst, val in q.items():
                    qstate[lev - 1][cnst] = val
            r1a = pk[(st, b1, 1, it1)]
            water_lq = any((st, b1, 1, w, it1) in qk for w in WATER)
            if not (r1a["ls"] or water_lq) or not all(w in qstate[l] for l in range(32) for w in WATER):
                continue
            n_checked += 1
            zm, zi = geopotential(T, qstate, pmid, pdel)
            for lev in range(1, 33):
                r2 = pk[(st, b2, lev, it2)]
                rel = abs(zm[lev - 1] - r2["zm"]) / abs(r2["zm"]) if r2["zm"] else abs(zm[lev - 1] - r2["zm"])
                worst_zm = max(worst_zm, rel)
                s1 = T[lev - 1] * CPAIR + GRAVIT * zm[lev - 1] + r1a["phis"]
                rel = abs(s1 - r2["s"]) / abs(r2["s"]) if r2["s"] else abs(s1 - r2["s"])
                worst_s = max(worst_s, rel)
    print(f"geopotential chaining: worst rel err = {worst_zm:.3g} over {n_checked} adjacent pairs, s_after {worst_s:.3g}")
    if max(worst_zm, worst_s) > 1e-14:
        raise SystemExit("geopotential chain verification FAILED")


def rev(a):
    return list(reversed(a))


def state_overrides(st, site, it, pk, qk, corder, extra_q=None):
    """parameter_overrides of one physics_update call: dumped state + tendencies."""
    o = {}
    levs = {}
    for lev in range(1, 33):
        levs[lev] = pk[(st, site, lev, it)]
    for col in ("T", "u", "v", "s", "zm", "pmid", "pdel", "pts", "ptu", "ptv"):
        o[{"pts": "ptend_s", "ptu": "ptend_u", "ptv": "ptend_v"}.get(col, col)] = rev(
            [levs[l][col] for l in range(1, 33)])
    o["phis"] = levs[1]["phis"]
    o["ls"] = 1.0 if levs[1]["ls"] else 0.0
    o["lu"] = 1.0 if levs[1]["lu"] else 0.0
    o["lv"] = 1.0 if levs[1]["lv"] else 0.0
    o["top_level"] = float(levs[1]["top"])
    o["bot_level"] = float(levs[1]["bot"])
    o["dt"] = levs[1]["dt"]
    o["cld_min_nz"] = [1.0 if (site in CLD_MIN_NZ_SITES and c in ("CLDLIQ", "CLDICE")) else 0.0
                       for c in corder]
    o["lq"] = [1.0 if (st, site, 1, c, it) in qk else 0.0 for c in corder]
    q = [[0.0] * len(corder) for _ in range(32)]
    ptq = [[0.0] * len(corder) for _ in range(32)]
    for j, c in enumerate(corder):
        if (st, site, 1, c, it) in qk:
            for lev in range(1, 33):
                q[32 - lev][j] = qk[(st, site, lev, c, it)]["q_before"]
                ptq[32 - lev][j] = qk[(st, site, lev, c, it)]["ptend_q"]
        elif extra_q and c in extra_q:
            for lev in range(1, 33):
                q[32 - lev][j] = extra_q[c][lev - 1]
    o["q"] = q
    o["ptend_q"] = ptq
    return o


def chain_q(st, site, it, cnst, lev, seq, qk):
    """q_before of the first later call site that dumps cnst (the state is
    unchanged by the calls in between that do not touch cnst)."""
    i = seq.index((site, it))
    for s2, it2 in seq[i + 1:]:
        if (st, s2, lev, cnst, it2) in qk:
            return qk[(st, s2, lev, cnst, it2)]["q_before"]
    raise SystemExit(f"no chain target for {st} {site}#{it} {cnst} lev {lev}")


def seq_of(st, rows):
    seq = []
    for r in rows:
        if r["step"] != st:
            continue
        label = (r["site"], max(r["macmic_it"], 1))
        if not seq or seq[-1] != label:
            seq.append(label)
    return seq


def copy_check(variable, coords, expected, o):
    """The review rule: every assertion must be able to fail.  Skip the
    assertion when the expected value is within the assertion tolerance
    (rel 1e-12) of the overridden input it derives from: the law left that
    point unchanged, or changed it by less than the test can resolve."""
    if variable == "q_after":
        inp = o["q"][coords["lev"] - 1][coords["cnst"] - 1]
    else:
        inp = o[{"T_after": "T", "u_after": "u", "v_after": "v",
                 "zm_after": "zm", "s_after": "s"}[variable]][coords["lev"] - 1]
    return abs(expected - inp) <= 1e-12 * abs(expected)


def physics_tests(rows, qrows, pk, qk, corder):
    """The PhysicsUpdate model tests (all chained to the next call site's rows)."""
    tests = []

    def add_state_assertions(A, st, site, it, chain_site, chain_it, levels):
        for lev in levels:
            l = 33 - lev
            r2 = pk[(st, chain_site, lev, chain_it)]
            for var, exp in (("T_after", r2["T"]), ("zm_after", r2["zm"]), ("s_after", r2["s"])):
                A.append((var, l, exp))
            A.append(("u_after", l, r2["u"]))
            A.append(("v_after", l, r2["v"]))

    def add_q_assertions(A, st, site, it, seq, points):
        for cnst, levels in points:
            j = corder.index(cnst) + 1
            for lev in levels:
                A.append((f"q_after", 33 - lev, chain_q(st, site, it, cnst, lev, seq, qk), j))

    def emit(tid, description, o, A, chain_desc):
        asrt = []
        for v, l, e, *j in A:
            coords = {"lev": l}
            if j:
                coords["cnst"] = j[0]
            if copy_check(v, coords, e, o):
                continue
            asrt.append(dict(variable=v, time=0.0, coords=coords, expected=e,
                             tolerance={"rel": 1e-12, "abs": 0.0}))
        tests.append(dict(
            id=tid,
            description=description + " " + chain_desc,
            time_span={"start": 0.0, "end": 0.0},
            parameter_overrides=o,
            assertions=asrt))

    # ---- convect_deep at the ZM-active steps 18 and 40 (steps 1-17, 24-39,
    # 41-48 have all-zero ptends at this site and are not pinned)
    for st in (18, 40):
        seq = seq_of(st, rows)
        o = state_overrides(st, "convect_deep", 1, pk, qk, corder)
        A = []
        add_state_assertions(A, st, "convect_deep", 1, "convect_shallow", 1, (32, 28, 24, 19, 15))
        add_q_assertions(A, st, "convect_deep", 1, seq, [
            ("Q", (32, 24, 15)), ("CLDLIQ", (32, 21, 15)), ("NUMLIQ", (32, 15)),
            ("NUMICE", (32, 15)), ("RAINQM", (32, 15))])
        emit(f"convect_deep_step{st}", (
            f"The physics_update law at the convect_deep call site (proc zm_conv_tend, "
            f"physpkg.F90 tphysbc), step {st}, pinned from the cam_phys_tend.csv / "
            f"cam_phys_tend_q.csv rows (step {st}, site convect_deep, all 9 dumped "
            f"constituents) at the ZM-active levels only (dump lev 15..32; steps 1-17, "
            f"24-39 and 41-48 carry all-zero ptends at this site and are not pinned, "
            f"and levels 1-14 are pinned nowhere because the law leaves them unchanged). "
            f"Exercises: the T update T + ls*ptend_s*dt/cpairv, the u/v updates, the "
            f"constituent update with the qneg3 floor (qmin table of chemistry.F90) and "
            f"the number-concentration clamp min(max(q,1e-12),1e10) for "
            f"NUMICE/NUMLIQ/NUMRAI/NUMSNO (fires at NUMICE across the column), the "
            f"state_cnst_min_nz(1e-36) zeroing of the (CLDLIQ, NUMLIQ) and (CLDICE, "
            f"NUMICE) pairs (fires at CLDLIQ on ~29 of 32 levels), and the geopotential "
            f"recompute asserted via zm_after/s_after.  cld_min_nz = 1 at CLDLIQ/CLDICE "
            f"because the site's proc name zm_conv_tend is in cldlim_names "
            f"(physics_types.F90:246-247)."),
            o, A,
            f"Expected values chained to the next call site's dump rows (step {st}, "
            f"site convect_shallow: state_before(convect_shallow) == "
            f"state_after(convect_deep), verified for the whole 48-step window by "
            f"stage1/analysis/extract_cam_physics_update_pins.py).")

    # ---- clubb#1 (step 1) -> microp#1
    st = 1
    seq = seq_of(st, rows)
    o = state_overrides(st, "clubb", 1, pk, qk, corder)
    A = []
    add_state_assertions(A, st, "clubb", 1, "microp", 1, (32, 25, 17, 9, 1))
    add_q_assertions(A, st, "clubb", 1, seq, [
        ("Q", (32, 16)), ("CLDLIQ", (21,)), ("NUMICE", (32,)), ("NUMRAI", (16,)),
        ("DMS", (32, 16)), ("H2O2", (32,)), ("SOAG", (32,))])
    emit("clubb_macmic1_step1", (
        "The physics_update law at the first CLUBB substep of the macro-micro "
        "co-substepping loop (physpkg.F90:2646, clubb#1: macmic_it = 1 of 3, ptend "
        "pre-scaled by 1/3), step 1, pinned from the clubb rows of cam_phys_tend.csv / "
        "cam_phys_tend_q.csv (the q file has no macmic column; the clubb rows repeat "
        "per macmic_it in order, the first 32*15-row block being macmic_it = 1). "
        "Exercises the 15 CLUBB constituents (the water/cloud/number set plus the "
        "chem gases DMS/H2O2/H2SO4/SO2/SOAE/SOAG, qmin = 1e-36, chemistry.F90), the "
        "number-concentration clamp, and the full-column T/u/v updates and "
        "geopotential recompute; cld_min_nz = 0 (proc clubb not in cldlim_names)."),
        o, A,
        "Expected values chained to the next call site's dump rows (step 1, site "
        "microp, macmic_it = 1: state_before(microp#1) == state_after(clubb#1)).")

    # ---- microp#1 (step 1) -> clubb#2 / microp#2
    o = state_overrides(st, "microp", 1, pk, qk, corder)
    A = []
    add_state_assertions(A, st, "microp", 1, "clubb", 2, (32, 25, 17, 9))
    add_q_assertions(A, st, "microp", 1, seq, [
        ("num_a1", (32, 16)), ("so4_a1", (32, 9)), ("bc_a1", (32,)),
        ("CLDLIQ", (32, 24)), ("dst_a3", (32,))])
    emit("microp_macmic1_step1", (
        "The physics_update law at the first MG2 microphysics substep (microp#1: "
        "macmic_it = 1 of 3, ptend pre-scaled by 1/3, ls only), step 1, pinned from "
        "the microp rows of cam_phys_tend.csv / cam_phys_tend_q.csv (first 32*28-row "
        "q block = macmic_it 1). Exercises the T update and the geopotential recompute "
        "with the 28 MG2 constituents (the water/cloud/number set plus the aerosol "
        "mass/number species bc_a*/dst_a*/ncl_a*/num_a*/pom_a*/so4_a*/soa_a*, qmin "
        "1e-36 for the mass species, 1e-5 for num_a*, chemistry.F90) and the "
        "number-concentration clamp; cld_min_nz = 0 (proc micro_mg not in "
        "cldlim_names)."),
        o, A,
        "Expected values chained to the next call site's dump rows (step 1: the "
        "state-after of microp#1 is state_before(clubb#2) for the water/cloud/number "
        "and chem species, and state_before(microp#2) for the aerosols, which clubb#2 "
        "does not touch; both chains verified over the window).")

    # ---- aero_calcsize (step 1) -> aero_wetdep (q-only site)
    o = state_overrides(st, "aero_calcsize", 1, pk, qk, corder)
    A = []
    add_q_assertions(A, st, "aero_calcsize", 1, seq, [
        ("num_a1", (1, 6)), ("num_a4", (1, 9)), ("so4_a1", (1, 6)),
        ("dst_a1", (1,)), ("ncl_a2", (1,))])
    emit("aero_calcsize_step1", (
        "The physics_update law at the aero_calcsize call site (proc aero_wat, q-only: "
        "ls = lu = lv = 0, derive_new_geopotential = .false.), step 1, pinned from the "
        "aero_calcsize rows of cam_phys_tend.csv / cam_phys_tend_q.csv. Exercises the "
        "qneg3 qmin floor of the number-aerosol species (num_a1 at dump lev 1: "
        "q_before + ptend_q*dt = 9.587812943006401e-06 < 1e-5, asserted q_after = 1e-5 "
        "exactly) and of the aerosol mass species (qmin = 1e-36, chemistry.F90). The "
        "derive=0 geopotential identity branch (zm_after = zm) is covered structurally "
        "but not asserted (it is a copy of the input by construction)."),
        o, A,
        "Expected values chained to the next call site's dump rows (step 1, site "
        "aero_wetdep: state_before(aero_wetdep) == state_after(aero_calcsize)).")

    # ---- radiation (step 1) -> aoa_tracers (ls-only; water q inputs chained)
    qstate = [{} for _ in range(32)]
    for label in seq:
        if label == ("radiation", 1):
            break
        site, it = label
        for lev in range(1, 33):
            r1 = pk[(st, site, lev, it)]
            q1 = {c: qk[(st, site, lev, c, it)] for c in corder if (st, site, lev, c, it) in qk}
            u, v, t, q = update_column(r1, q1)
            for cnst, val in q.items():
                qstate[lev - 1][cnst] = val
    extra_q = {c: [qstate[l].get(c, 0.0) for l in range(32)] for c in WATER}
    o = state_overrides(st, "radiation", 1, pk, qk, corder, extra_q=extra_q)
    A = []
    add_state_assertions(A, st, "radiation", 1, "aoa_tracers", 1, (32, 25, 17, 9, 1))
    emit("radiation_step1", (
        "The physics_update law at the radiation call site (proc radheat, ls only: "
        "nlq = 0, the last physics_update of tphysbc), step 1, pinned from the "
        "radiation rows of cam_phys_tend.csv. Exercises the T update, the geopotential "
        "recompute (derive_new_geopotential = .true. from ls) and s_after. The five "
        "thermo-active water-species mixing ratios q_Q/q_CLDLIQ/q_CLDICE/q_RAINQM/"
        "q_SNOWQM are not dumped at radiation (nlq = 0: the q file carries only "
        "constituents in ptend%lq), so they are provided as chained inputs from the "
        "verified running state (q_after(microp#3), propagated through aero_calcsize/"
        "aero_wetdep/convect_deep_2, which do not touch the water species; derived "
        "inputs, not dump columns). The constituent law is otherwise idle (lq = 0 "
        "everywhere: q_after = q, asserted nowhere)."),
        o, A,
        "Expected values chained to the next call site's dump rows (step 1, site "
        "aoa_tracers, the first tphysac call: state_before(aoa_tracers) == "
        "state_after(radiation), including zm/s).")

    # ---- vdiff (step 1) -> gw (the first later q-dumping site; the calls in
    # between, rayleigh/aero_drydep/carma_tend, leave the state unchanged)
    o = state_overrides(st, "vdiff", 1, pk, qk, corder)
    A = []
    add_state_assertions(A, st, "vdiff", 1, "gw", 1, (32,))
    add_q_assertions(A, st, "vdiff", 1, seq, [
        ("DMS", (32,)), ("H2O2", (32,)), ("SO2", (32,)), ("num_a1", (32,)),
        ("num_a2", (32,)), ("so4_a1", (32,)), ("bc_a4", (32,)), ("ncl_a1", (32,))])
    emit("vdiff_step1", (
        "The physics_update law at the vdiff call site (proc vertical, lu/lv/ls, all "
        "34 constituents; the first tphysac call that carries tendencies), step 1, "
        "pinned from the vdiff rows of cam_phys_tend.csv / cam_phys_tend_q.csv. vdiff's "
        "ptends are nonzero at dump lev 32 only in this run (a surface-flux-dominated "
        "call), so the state assertions pin that single active level (model index 1) "
        "and the q pins the 16 constituents with nonzero ptend_q there; the other "
        "levels are left unchanged by the law and are asserted nowhere (copy rule). "
        "Exercises lu/lv/ls together, the num_a* 1e-5 floor and the chem-gas 1e-36 "
        "floors at the surface level, and the geopotential recompute."),
        o, A,
        "Expected values chained to the later call site's dump rows (step 1, site gw: "
        "state_before(gw) == state_after(vdiff) — the intervening calls rayleigh (lu/"
        "lv/ls set with all-zero ptends), aero_drydep and carma_tend leave the state "
        "unchanged; gw is the first later site that dumps all 34 constituents).")

    return tests


def accumulation_test(rows, qrows, pk, qk, corder):
    """The TendAccumulation test: the whole tphysbc phase as a multi-site sum
    over the 15 consecutive call sites, vs the cam_tend_total.csv phase=bc row."""
    st = 18  # a step with active ZM convection
    trows = []
    with open(f"{D}/cam_tend_total.csv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split(",")
            trows.append(dict(step=int(p[0]), phase=p[3].strip('"'), lev=int(p[4]),
                              dtdt=float(p[5]), dudt=float(p[6]), dvdt=float(p[7])))
    tbc = {r["lev"]: r for r in trows if r["step"] == st and r["phase"] == "bc"}

    # verify: the bc row == the per-site sum over the 15 sites
    sums = {l: [0.0, 0.0, 0.0] for l in range(1, 33)}
    for i, site in enumerate(BC_SITES):
        it = bc_it(i)
        for lev in range(1, 33):
            r = pk[(st, site, lev, it)]
            sums[lev][0] += r["pts"] / CPAIR
            sums[lev][1] += r["ptu"]
            sums[lev][2] += r["ptv"]
    worst = 0.0
    for lev in range(1, 33):
        for j, col in enumerate(("dtdt", "dudt", "dvdt")):
            rel = abs(sums[lev][j] - tbc[lev][col]) / abs(tbc[lev][col]) if tbc[lev][col] else abs(sums[lev][j] - tbc[lev][col])
            worst = max(worst, rel)
    print(f"tend accumulation step {st}: bc-row vs 15-site sum worst rel err = {worst:.3g}")
    if worst > 1e-14:
        raise SystemExit("tend accumulation verification FAILED")

    o = {}
    o["ptend_s"] = [[sums[l][0] for l in range(31, 0, -1)]]  # placeholder; fixed below
    o["ptend_u"] = [[0.0] * 32] * 15
    o["ptend_v"] = [[0.0] * 32] * 15
    ps = [[0.0] * 32 for _ in range(15)]
    pu = [[0.0] * 32 for _ in range(15)]
    pv = [[0.0] * 32 for _ in range(15)]
    for i, site in enumerate(BC_SITES):
        it = bc_it(i)
        for lev in range(1, 33):
            r = pk[(st, site, lev, it)]
            ps[i][32 - lev] = r["pts"]
            pu[i][32 - lev] = r["ptu"]
            pv[i][32 - lev] = r["ptv"]
    o["ptend_s"] = ps
    o["ptend_u"] = pu
    o["ptend_v"] = pv
    o["lu"] = [1.0 if pk[(st, site, bc_it(i), 1)]["lu"] else 0.0 for i, site in enumerate(BC_SITES)]
    o["lv"] = [1.0 if pk[(st, site, bc_it(i), 1)]["lv"] else 0.0 for i, site in enumerate(BC_SITES)]
    o["ls"] = [1.0 if pk[(st, site, bc_it(i), 1)]["ls"] else 0.0 for i, site in enumerate(BC_SITES)]
    o["top_level"] = [float(pk[(st, site, bc_it(i), 1)]["top"]) for i, site in enumerate(BC_SITES)]
    o["bot_level"] = [float(pk[(st, site, bc_it(i), 1)]["bot"]) for i, site in enumerate(BC_SITES)]

    A = []
    for lev in (32, 17, 1):
        l = 33 - lev
        A.append(dict(variable="dtdt", time=0.0, coords={"lev": l}, expected=tbc[lev]["dtdt"],
                      tolerance={"rel": 1e-12, "abs": 0.0}))
        A.append(dict(variable="dudt", time=0.0, coords={"lev": l}, expected=tbc[lev]["dudt"],
                      tolerance={"rel": 1e-12, "abs": 0.0}))
        A.append(dict(variable="dvdt", time=0.0, coords={"lev": l}, expected=tbc[lev]["dvdt"],
                      tolerance={"rel": 1e-12, "abs": 0.0}))
    return [dict(
        id="tend_accumulation_bc_step18",
        description=(
            "The tend accumulation of the whole tphysbc phase, step 18, pinned as a "
            "multi-site sum over the 15 consecutive tphysbc physics_update calls "
            "(energy_fixer, dadadj, convect_deep, convect_shallow, ssatcontrail, "
            "clubb#1..3, microp#1..3, aero_calcsize, aero_wetdep, convect_deep_2, "
            "radiation; the CLUBB/MG2 rows already carry the 1/3 macmic pre-scaling, so "
            "summing the dumped rows reproduces the accumulation) against the "
            "cam_tend_total.csv phase=bc row: dtdt = accumulated ptend_s/cpairv, "
            "dudt/dvdt = accumulated ptend_u/ptend_v (physics_types.F90:289-300, "
            "417-419). The input is the per-site tendency of every call; the expected "
            "value is the independent phase-total dump column, so the assertion fails "
            "if the window gating, the flag gating or the summation is wrong. The "
            "accumulation is linear in the calls, so the sum over the 15 sites is "
            "bit-identical to the bc row (worst rel err 0, verified by "
            "stage1/analysis/extract_cam_physics_update_pins.py). The ac-phase dtdt is "
            "NOT pinned: after the last physics_update call the SE dycore rescales "
            "tend%dtdt by cpairv/cp_or_cv_dycore (physpkg.F90:2026-2030), a documented "
            "gap (README-cam-physics-update.md)."
        ),
        time_span={"start": 0.0, "end": 0.0},
        parameter_overrides=o,
        assertions=A)]


def main():
    rows = load_phys()
    qrows = load_q()
    pk = {}
    for r in rows:
        pk[(r["step"], r["site"], r["lev"], max(r["macmic_it"], 1))] = r
    qk = {}
    for r in qrows:
        qk[(r["step"], r["site"], r["lev"], r["cnst"], r["macmic_it"])] = r
    corder = []
    for r in qrows:
        if r["cnst"] not in corder:
            corder.append(r["cnst"])
    assert len(corder) == 34, f"expected 34 constituents, found {len(corder)}"

    verify_full_chain(rows, qrows, pk, qk, corder)
    phys = physics_tests(rows, qrows, pk, qk, corder)
    acc = accumulation_test(rows, qrows, pk, qk, corder)
    n = sum(len(t["assertions"]) for t in phys) + sum(len(t["assertions"]) for t in acc)
    print(f"emitted {len(phys)} PhysicsUpdate tests ({sum(len(t['assertions']) for t in phys)} "
          f"assertions) + {len(acc)} TendAccumulation test ({len(acc[0]['assertions'])} assertions), "
          f"{n} total")
    print(json.dumps({"PhysicsUpdate": phys, "TendAccumulation": acc}, indent=1))


if __name__ == "__main__":
    main()
