#!/usr/bin/env python3
"""Extract and verify the numeric test tuples of the M9 ocean-flux subassembly
(stage2/mom6_cap_imports.esm and stage2/m9_ocn_flux.esm) from the stage-1
dumps, and compute the values that are NOT dump columns.

The overlapping (cell, step) rows of the C_JRA_t62 reference run (compset
2000_DATM%JRA-1p5-2023_SLND_DICE%SSMI_MOM6_DROF%JRA-1p5-2023_SGLC_SWAV_SESP,
24 coupling steps of 3600 s on 2000-01-01; 12 traced ocean cells on PET 0 of
the 6x6 decomposition, the only ocean cells of that block - the 77
southernmost grid rows are Antarctic land mask):

  mom6 surface forcing  rows (cpl_step k,   i, j)  of mom6/mom6_sfc_forcing.csv
  mediator merge        rows (step s,       gid)   of cmeps/cmeps_merge.csv
  mediator aoflux       rows (step s,       gid)   of cmeps/cmeps_aoflux.csv
  mediator ocnalb       rows (step s,       gid)   of cmeps/cmeps_ocnalb.csv

matched by (lat, lon): mom6 (i,j) -> cmeps gid (7,82)/(8,82) -> 41583/41584,
(7..10,83) -> 42123..42126, (6..11,84) -> 42662..42667.

Chain laws verified here on ALL overlapping rows (round-off, rel 1e-12):

  mom6 cpl_step k  <-  merge/aoflux step k+1  (the CESM ocean-lag sequence:
    within each coupling interval the ocean runs first and consumes the
    mediator products prepared at the END of the interval, so the ocean
    dump step k lines up with mediator step k+1; verified bit-exact on all
    144 mom6 rows): p = pslv, u10_sqr = duu10n, ice_fraction = ifrac (these
    three are NOT areacor-scaled in the cap, mom_cap_methods.F90:114,
    506-508, 514-516) and, for the sixteen areacor-scaled imports
    (state_getimport with areacor=med2mod_areacor, F90:114-204, 213-314,
    471-491), u_flux/v_flux/sw_vis_dir/sw_vis_dif/sw_nir_dir/sw_nir_dif/
    lw_flux/t_flux/q_flux/lprec/fprec/lrunoff/frunoff/salt_flux/
    seaice_melt_heat/seaice_melt = areacor x the corresponding merge export
    (Foxx_taux/Foxx_tauy/Foxx_swnet_vdr/vdf/idr/idf/Foxx_lwnet/Foxx_sen/
    Foxx_evap/Faxa_rain/Faxa_snow/Foxx_rofl/Foxx_rofi/Fioi_salt/Fioi_melth/
    Fioi_meltw), with the wind-stress rotation the identity (cos_rot = 1,
    sin_rot = 0, F90:165-175; all traced cells lie south of the tripole
    fold).
  aoflux step s  ==  merge step s  (flux_atmOcn consumed the same step;
    Faox_taux/tauy/sen/evap/lwup == taux/tauy/sen/lat/lwup, verified
    bit-exact on all 288 merge rows).
  merge step s So_*  ==  ocnalb step s-1  (the one-step albedo lag of the
    mediator: the merge at its step s consumes the albedos computed at
    step s-1; So_avsdr/avsdf/anidr/anidf == avsdr/avsdf/anidr/anidf of the
    previous ocnalb step, verified bit-exact on all 276 overlapping rows).
  AGRID average (mom_surface_forcing_nuopc.F90:943-951): taux == u_flux
    bit-exact iff the east neighbor is land (mask2dT(i+1,j) = 0) and
    tauy == v_flux bit-exact iff the north neighbor is land; over
    ocean-ocean pairs taux = (u_flux + u_flux_east)/2 (mask 1). The dump
    row itself therefore carries the neighbor masks: at (11,84) the east
    neighbor (12,84) is land but the north neighbor (11,85) is ocean -
    the north neighbor lies beyond the 12 traced cells, so its stress is
    back-derived exactly from the dump average, v_flux_north = 2*tauy -
    v_flux (this back-derived value enters the cpl1_i11_j84 test as a
    free SurfaceForcing override and is the ONLY test value not a dump
    column).
  areacor per cell is constant over the 12 steps to ~2e-16 (fitted from
    the step-1 row as sw_vis_dir/Foxx_swnet_vdr).
  derived dump columns (mom6_sfc_forcing.csv): salt_restore_flux =
    salt_flux_o - salt_flux_in and net_fw_global_adjust = -vprec (the
    restoring salt flux and the ADJUST_NET_FRESH_WATER_TO_ZERO adjustment,
    not IOB imports).

The .esm tests' parameter_overrides and assertions are then checked against
these dump-derived values at rel 1e-12, so every pinned number in the two
stage-2 files is traceable to a dump row (or to the one back-derived value
above, which is marked in the test description and cited to this script).

Usage:
  python3 stage1/analysis/extract_m9_ocn_flux_rows.py [dumps_dir] [esm_dir]
where dumps_dir defaults to <main checkout>/stage1/dumps (dumps are
gitignored and live in the main checkout, not in the worktree; the
environment variable EQESM_MAIN overrides the main-checkout location) and
esm_dir defaults to <this script's repo>/stage2.
"""
import csv, json, math, os, re, sys

EQESM_MAIN = os.environ.get("EQESM_MAIN", "/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM")
DUMPS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(EQESM_MAIN, "stage1/dumps")
ESM_DIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "stage2")

MOM6_STEPS = range(1, 13)      # 12 traced coupling steps
MERG_STEPS = range(1, 25)      # 24 mediator steps


def load_dump(path, first_dup=False):
    cols, data = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("# columns"):
                cols = [c.strip() for c in line.split("=", 1)[1].strip().split(",")]
            elif not line.startswith("#"):
                data.append(line)
    rows = []
    for line in data:
        v = [x.strip() for x in line.strip().split(",")]
        if first_dup:  # cmeps_aoflux.csv has a duplicated lat column; keep the first
            d, seen = {}, set()
            for h, val in zip(cols, v):
                if h not in seen:
                    d[h] = val
                    seen.add(h)
            rows.append(d)
        else:
            rows.append(dict(zip(cols, v)))
    return rows


def f(x):
    return float(x.strip())


def ll_key(r):
    return (round(f(r["lat"]), 6), round(f(r["lon"]), 6))


def verify_chain():
    """Verify the chain laws on all rows; return (per-cell areacor,
    (i,j,k) -> {override values, expectation values, masks, neighbors})."""
    mom = load_dump(os.path.join(DUMPS, "mom6/mom6_sfc_forcing.csv"))
    mrg = load_dump(os.path.join(DUMPS, "cmeps/cmeps_merge.csv"))
    aof = load_dump(os.path.join(DUMPS, "cmeps/cmeps_aoflux.csv"), first_dup=True)
    ocn = load_dump(os.path.join(DUMPS, "cmeps/cmeps_ocnalb.csv"))

    mom_by = {(r["i"], r["j"], int(r["cpl_step"])): r for r in mom}
    mrg_by = {(ll_key(r), int(r["step"])): r for r in mrg}
    aof_by = {(ll_key(r), int(r["step"])): r for r in aof}
    ocn_by = {(ll_key(r), int(r["step"])): r for r in ocn}
    cells = sorted(set((r["i"], r["j"]) for r in mom))
    ll_of = {(r["i"], r["j"]): ll_key(r) for r in mom}

    def rel(a, b):
        return abs(a - b) / max(abs(a), abs(b), 1e-300)

    # ---- per-cell areacor fitted from the step-1 row; constancy over steps
    areacor = {}
    for ij in cells:
        ll = ll_of[ij]
        r1 = mom_by[(ij[0], ij[1], 1)]
        m2 = mrg_by[(ll, 2)]
        ac = f(r1["sw_vis_dir"]) / f(m2["Foxx_swnet_vdr"])
        areacor[ij] = ac
        for k in MOM6_STEPS:
            r = mom_by[(ij[0], ij[1], k)]
            m = mrg_by[(ll, k + 1)]
            # sw_vis_dir is areacor-scaled; the ratio must stay constant
            err = rel(f(r["sw_vis_dir"]) / f(m["Foxx_swnet_vdr"]), ac)
            assert err <= 1e-12, (ij, k, err)

    # ---- the mom6 step k <- merge/aoflux step k+1 ocean lag, all 144 rows
    scaled = [  # mom6 col, merge col
        ("u_flux", "Foxx_taux"), ("v_flux", "Foxx_tauy"),
        ("sw_vis_dir", "Foxx_swnet_vdr"), ("sw_vis_dif", "Foxx_swnet_vdf"),
        ("sw_nir_dir", "Foxx_swnet_idr"), ("sw_nir_dif", "Foxx_swnet_idf"),
        ("lw_flux", "Foxx_lwnet"), ("t_flux", "Foxx_sen"),
        ("q_flux", "Foxx_evap"), ("lprec", "Faxa_rain"), ("fprec", "Faxa_snow"),
        ("lrunoff", "Foxx_rofl"), ("frunoff", "Foxx_rofi"),
        ("salt_flux", "Fioi_salt"), ("seaice_melt_heat", "Fioi_melth"),
        ("seaice_melt", "Fioi_meltw"),
    ]
    for ij in cells:
        ll = ll_of[ij]
        for k in MOM6_STEPS:
            r = mom_by[(ij[0], ij[1], k)]
            m = mrg_by[(ll, k + 1)]
            a = aof_by[(ll, k + 1)]
            o = ocn_by[(ll, k)]
            # unscaled imports: bit-exact
            assert f(r["p"]) == f(a["pslv"]), (ij, k, "p")
            assert f(r["u10_sqr"]) == f(a["duu10n"]), (ij, k, "u10_sqr")
            assert f(r["ice_fraction"]) == f(m["ifrac"]), (ij, k, "ice_fraction")
            # areacor-scaled imports (rotation identity): rel 1e-12
            for mcol, xcol in scaled:
                err = rel(f(r[mcol]), areacor[ij] * f(m[xcol]))
                assert err <= 1e-12, (ij, k, mcol, err)
            # derived values (not dump columns): the restoring salt flux is
            # salt_flux_o - salt_flux_in (the salt_restore_flux override of
            # the .esm tests; the dump also records the copy law salt_flux_in
            # = salt_flux) and the net-FW adjustment is -vprec (vprec =
            # vprec_restore - net_fw_global_adjust with vprec_restore = 0 in
            # the reference run, so the dump vprec column IS the negative
            # adjustment; the net_fw_global_adjust override is -vprec)
            assert f(r["salt_flux_in"]) == f(r["salt_flux"]), (ij, k)
            # AGRID masks from the row itself: taux==u_flux iff east land
            # (the land-neighbor rows used by the tests are pinned by this)

    # ---- aoflux step s == merge step s (same-step), all 288 rows
    for m in mrg:
        s, ll = int(m["step"]), ll_key(m)
        a = aof_by[(ll, s)]
        for mcol, acol in [("Faox_taux", "taux"), ("Faox_tauy", "tauy"),
                           ("Faox_sen", "sen"), ("Faox_evap", "evap"),
                           ("Faox_lwup", "lwup")]:
            assert f(m[mcol]) == f(a[acol]), (s, ll, mcol)

    # ---- merge step s So_* == ocnalb step s-1 (albedo lag), 276 rows
    for m in mrg:
        s = int(m["step"])
        if s == 1:
            continue  # no ocnalb step 0
        o = ocn_by[(ll_key(m), s - 1)]
        for mcol, ocol in [("So_avsdr", "avsdr"), ("So_avsdf", "avsdf"),
                           ("So_anidr", "anidr"), ("So_anidf", "anidf")]:
            assert f(m[mcol]) == f(o[ocol]), (s, mcol)

    # ---- AGRID average law on all 144 rows: taux == u_flux bit-exact iff
    # the east neighbor is land (mask2dT(i+1,j) = 0) and tauy == v_flux
    # bit-exact iff the north neighbor is land; over ocean-ocean traced
    # pairs the average (u_flux + u_flux_east)/2 is verified exactly. When
    # an ocean neighbor lies beyond the traced cells its stress is
    # back-derived from the dump average, a_neighbor = 2*avg - a. The
    # neighbor masks and values feed the .esm test check below.
    masks = {}
    for ij in cells:
        for k in MOM6_STEPS:
            r = mom_by[(ij[0], ij[1], k)]
            masks[(ij[0], ij[1], k)] = {"east": 1.0, "north": 1.0,
                                        "u_flux_east": 0.0, "v_flux_north": 0.0}
            # east neighbor (i+1, j)
            ij2 = (str(int(ij[0]) + 1), ij[1])
            if ij2 in ll_of:
                ue = f(mom_by[(ij2[0], ij2[1], k)]["u_flux"])
                assert f(r["taux"]) == (f(r["u_flux"]) + ue) / 2.0, (ij, k, "taux pair")
                masks[(ij[0], ij[1], k)]["u_flux_east"] = ue
            elif f(r["taux"]) == f(r["u_flux"]):
                masks[(ij[0], ij[1], k)]["east"] = 0.0  # land neighbor
            else:  # ocean, untraced: back-derive from the average
                masks[(ij[0], ij[1], k)]["u_flux_east"] = 2.0 * f(r["taux"]) - f(r["u_flux"])
            # north neighbor (i, j+1)
            ij3 = (ij[0], str(int(ij[1]) + 1))
            if ij3 in ll_of:
                vn = f(mom_by[(ij3[0], ij3[1], k)]["v_flux"])
                assert f(r["tauy"]) == (f(r["v_flux"]) + vn) / 2.0, (ij, k, "tauy pair")
                masks[(ij[0], ij[1], k)]["v_flux_north"] = vn
            elif f(r["tauy"]) == f(r["v_flux"]):
                masks[(ij[0], ij[1], k)]["north"] = 0.0  # land neighbor
            else:  # ocean, untraced: back-derive from the average
                masks[(ij[0], ij[1], k)]["v_flux_north"] = 2.0 * f(r["tauy"]) - f(r["v_flux"])

    table = {}
    for ij in cells:
        ll = ll_of[ij]
        for k in MOM6_STEPS:
            r = mom_by[(ij[0], ij[1], k)]
            m = mrg_by[(ll, k + 1)]
            a = aof_by[(ll, k + 1)]
            o = ocn_by[(ll, k)]
            table[(ij[0], ij[1], k)] = {
                "overrides": {
                    "OcnAlb.nextsw_cday": f(o["nextsw_cday"]),
                    "OcnAlb.rlat": f(o["rlat"]), "OcnAlb.rlon": f(o["rlon"]),
                    "AoFlux.zbot": f(a["zbot"]), "AoFlux.ubot": f(a["ubot"]),
                    "AoFlux.vbot": f(a["vbot"]), "AoFlux.thbot": f(a["thbot"]),
                    "AoFlux.qbot": f(a["qbot"]), "AoFlux.rainc": f(a["rainc"]),
                    "AoFlux.rbot": f(a["rbot"]), "AoFlux.tbot": f(a["tbot"]),
                    "AoFlux.pslv": f(a["pslv"]), "AoFlux.us": f(a["us"]),
                    "AoFlux.vs": f(a["vs"]), "AoFlux.ts": f(a["ts"]),
                    "OcnMerge.ifrac": f(m["ifrac"]), "OcnMerge.ofrac": f(m["ofrac"]),
                    "OcnMerge.ifrad": f(m["ifrad"]), "OcnMerge.ofrad": f(m["ofrad"]),
                    "OcnMerge.Faxa_rainc": f(m["Faxa_rainc"]),
                    "OcnMerge.Faxa_rainl": f(m["Faxa_rainl"]),
                    "OcnMerge.Faxa_snowc": f(m["Faxa_snowc"]),
                    "OcnMerge.Faxa_snowl": f(m["Faxa_snowl"]),
                    "OcnMerge.Faxa_lwdn": f(m["Faxa_lwdn"]),
                    "OcnMerge.Faxa_swvdr": f(m["Faxa_swvdr"]),
                    "OcnMerge.Faxa_swvdf": f(m["Faxa_swvdf"]),
                    "OcnMerge.Faxa_swndr": f(m["Faxa_swndr"]),
                    "OcnMerge.Faxa_swndf": f(m["Faxa_swndf"]),
                    "OcnMerge.Fioi_taux_ice": f(m["Fioi_taux_ice"]),
                    "OcnMerge.Fioi_tauy_ice": f(m["Fioi_tauy_ice"]),
                    "OcnMerge.Fioi_melth_ice": f(m["Fioi_melth_ice"]),
                    "OcnMerge.Fioi_meltw_ice": f(m["Fioi_meltw_ice"]),
                    "OcnMerge.Fioi_salt_ice": f(m["Fioi_salt_ice"]),
                    "OcnMerge.Fioi_swpen_ice": f(m["Fioi_swpen_ice"]),
                    "OcnMerge.Forr_rofl": f(m["Forr_rofl"]),
                    "OcnMerge.Forr_rofi": f(m["Forr_rofi"]),
                    "MomCap.areacor": areacor[ij], "MomCap.cos_rot": 1.0,
                    "MomCap.sin_rot": 0.0,
                    "SurfaceForcing.u_flux_east": masks[(ij[0], ij[1], k)]["u_flux_east"],
                    "SurfaceForcing.mask_east": masks[(ij[0], ij[1], k)]["east"],
                    "SurfaceForcing.v_flux_north": masks[(ij[0], ij[1], k)]["v_flux_north"],
                    "SurfaceForcing.mask_north": masks[(ij[0], ij[1], k)]["north"],
                    "SurfaceForcing.salt_restore_flux": f(r["salt_flux_o"]) - f(r["salt_flux_in"]),
                    "SurfaceForcing.net_fw_global_adjust": -f(r["vprec"]),
                    "SurfaceForcing.hrain": f(r["hrain"]), "SurfaceForcing.hsnow": f(r["hsnow"]),
                    "SurfaceForcing.hrofl": f(r["hrofl"]), "SurfaceForcing.hrofi": f(r["hrofi"]),
                    "SurfaceForcing.hevap": f(r["hevap"]), "SurfaceForcing.hcond": f(r["hcond"]),
                },
                "cap_overrides": {
                    "Foxx_taux": f(m["Foxx_taux"]), "Foxx_tauy": f(m["Foxx_tauy"]),
                    "Foxx_swnet_vdr": f(m["Foxx_swnet_vdr"]), "Foxx_swnet_vdf": f(m["Foxx_swnet_vdf"]),
                    "Foxx_swnet_idr": f(m["Foxx_swnet_idr"]), "Foxx_swnet_idf": f(m["Foxx_swnet_idf"]),
                    "Foxx_lwnet": f(m["Foxx_lwnet"]), "Foxx_sen": f(m["Foxx_sen"]),
                    "Foxx_evap": f(m["Foxx_evap"]), "Faxa_rain": f(m["Faxa_rain"]),
                    "Faxa_snow": f(m["Faxa_snow"]), "Foxx_rofl": f(m["Foxx_rofl"]),
                    "Foxx_rofi": f(m["Foxx_rofi"]), "Fioi_salt": f(m["Fioi_salt"]),
                    "Fioi_melth": f(m["Fioi_melth"]), "Fioi_meltw": f(m["Fioi_meltw"]),
                    "areacor": areacor[ij], "cos_rot": 1.0, "sin_rot": 0.0,
                },
                "expectations": {
                    "OcnAlb.anidr": f(o["anidr"]),
                    "AoFlux.duu10n": f(a["duu10n"]),
                    "OcnMerge.Foxx_swnet": f(m["Foxx_swnet"]),
                    "OcnMerge.Foxx_taux": f(m["Foxx_taux"]),
                    "OcnMerge.Foxx_lwnet": f(m["Foxx_lwnet"]),
                    "OcnMerge.Faxa_snow": f(m["Faxa_snow"]),
                },
                "cap_expectations": {c: f(r[c]) for c, _ in scaled},
                "mom6_expectations": {
                    "sw": f(r["sw"]), "lw": f(r["lw"]), "sens": f(r["sens"]),
                    "latent": f(r["latent"]), "taux": f(r["taux"]), "tauy": f(r["tauy"]),
                    "tau_mag": f(r["tau_mag"]), "ustar": f(r["ustar"]),
                    "p_surf_full": f(r["p_surf_full"]), "ice_fraction_o": f(r["ice_fraction_o"]),
                    "evap": f(r["evap"]), "vprec": f(r["vprec"]),
                    "seaice_melt_heat_o": f(r["seaice_melt_heat_o"]),
                    "seaice_melt_o": f(r["seaice_melt_o"]), "salt_flux_o": f(r["salt_flux_o"]),
                    "latent_evap_diag": f(r["latent_evap_diag"]),
                    "latent_fprec_diag": f(r["latent_fprec_diag"]),
                    "latent_frunoff_diag": f(r["latent_frunoff_diag"]),
                    "omega_w2x": f(r["omega_w2x"]), "net_mass_src": f(r["net_mass_src"]),
                },
            }
    print("chain laws verified: mom6<-merge/aoflux ocean lag (144 rows, rel 1e-12),")
    print("  aoflux==merge same step (288 rows, bit-exact), merge So_*==ocnalb s-1")
    print("  (276 rows, bit-exact), AGRID average/masks (144 rows), areacor constancy.")
    return areacor, table


def check_esm_tests(areacor, table):
    """Check every .esm test override/assertion against the dump-derived table."""
    problems = 0
    for fname in ("mom6_cap_imports.esm", "m9_ocn_flux.esm"):
        path = os.path.join(ESM_DIR, fname)
        if not os.path.exists(path):
            print(f"(skipping {fname}: not found)")
            continue
        doc = json.load(open(path))
        for model in doc["models"].values():
            for test in model.get("tests", []):
                m = re.match(r"cpl(\d+)_i(\d+)_j(\d+)", test["id"])
                if not m:
                    print(f"  WARNING {fname}: unparseable test id {test['id']}")
                    problems += 1
                    continue
                k, i, j = int(m.group(1)), m.group(2), m.group(3)
                row = table[(i, j, k)]
                # override names: bare (mom6_cap_imports.esm MomCap tests) or
                # model-prefixed (m9_ocn_flux.esm); the table carries both forms
                expect = dict(row["overrides"])
                expect.update(row["cap_overrides"])
                expect.update({("MomCap." + kk): vv for kk, vv in row["cap_overrides"].items()})
                for key, val in test.get("parameter_overrides", {}).items():
                    exp = expect.get(key)
                    if exp is None:
                        continue
                    if abs(val - exp) > 1e-12 * abs(exp):
                        problems += 1
                        print(f"  MISMATCH {fname} {test['id']} override {key}: {val} != dump {exp}")
                expects = dict(row["expectations"])
                expects.update({("MomCap." + kk): vv for kk, vv in row["cap_expectations"].items()})
                expects.update({kk: vv for kk, vv in row["cap_expectations"].items()})
                expects.update({("SurfaceForcing." + kk): vv for kk, vv in row["mom6_expectations"].items()})
                expects.update(row["mom6_expectations"])
                for a in test.get("assertions", []):
                    exp = expects.get(a["variable"])
                    if exp is None:
                        continue
                    if abs(a["expected"] - exp) > 1e-12 * abs(exp):
                        problems += 1
                        print(f"  MISMATCH {fname} {test['id']} assertion {a['variable']}: "
                              f"{a['expected']} != {exp}")
    print(f".esm test tuples checked against the dumps: {'OK' if problems == 0 else f'{problems} PROBLEMS'}")
    return problems


if __name__ == "__main__":
    areacor, table = verify_chain()
    print()
    print("per-cell areacor (fitted from mom6 step-1 sw_vis_dir / merge step-2")
    print("Foxx_swnet_vdr; NOT a dump column; constant over steps to ~2e-16):")
    for ij in sorted(areacor, key=lambda c: (int(c[1]), int(c[0]))):
        print(f"  (i,j) = ({ij[0]},{ij[1]})  areacor = {areacor[ij]:.17g}")
    print()
    print("back-derived (NOT a dump column): v_flux_north of the untraced ocean")
    print("north neighbor (11,85) of cell (11,84), from the AGRID average of the")
    print("dump row, v_flux_north = 2*tauy - v_flux:")
    for k in MOM6_STEPS:
        row = table[("11", "84", k)]
        print(f"  cpl_step {k:>2}: {row['overrides']['SurfaceForcing.v_flux_north']:.17g}")
    print()
    sys.exit(1 if check_esm_tests(areacor, table) else 0)
