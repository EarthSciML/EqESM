#!/usr/bin/env python3
"""Generate the curated tests block for stage2/cice_shortwave.esm from the
cice_dedd dump (compute_dEdd_3bd + solution_dEdd, icepack_shortwave.F90).

The dump has one row per (cell, category, surface type, spectral band); the
three rows of a (cell, category, srftyp) triple carry the per-band quantities,
with the accumulated outputs final on the band-3 row.  Two triples are used:

  snow: step 1, i=2, j=4, ncat=2, srftyp=1 (hs = 0.15364 m, hi = 1.2089 m,
        frsnw = 290.36 um, aerosols in all four slots); all three bands run
        the full snow/ice column; band 3 truncates at the interface below
        the snow SSL (trntdr drops below trmin = 1e-3).
  bare: same cell, ncat=1, srftyp=0 (hs = 0, hi = 0.34181 m); band 3
        truncates at the interface below the sea-ice SSL.

Index mapping (model -> dump):
  layer k  (1..13)        = Fortran layer k-1            -> tau_{k-1}, rdir_{k-1}, ...
  iface k  (1..14)        = Fortran iface k-1            -> trndir_{k-1}, dfdir_{k-1}, ...
  upward rupdir/rupdif    model iface k = Fortran iface 14-k -> rupdir_{14-k}
  snow_layer sl (1..3)    = Fortran snow index k=sl (Qs_0 = Qs_1 duplicate) -> Qs_{sl}
  dzk_{k-1}, frsnw_{sl}, Sabs_{sl}, Iabs_{il}, fthrul_{pl} as dumped (band-3 row
  for the running accumulations).

Curation (189 assertions: 113 snow + 76 bare; every assertion can fail):
  dzk layers {1,2,5,6,13} (geometry: snow SSL, snow interior, ice SSL,
    ice delta layer, ice interior; the bare-group snow layers are
    identically zero and skipped); frsnw all 3 snow layers; the snow IOPs
    Qs/ws_snw/gs_snw/ks_snw at snow layers {1,3} x all bands (snow group
    only - the dump records zeros for srftyp = 0 rows); the mixed tau/w0/g
    at layer 6 x all bands (snow) / layer 6 x all bands (bare; the bare
    layer-1 tau/w0/g are identically zero and skipped); the guard-zeroed
    layer props m_* at one active layer per band (band1 layer 6, band2
    layer 1, band3 layer 1 for snow; band1/2 layer 6 and band3 layer 5 for
    bare - the last active layer before the band-3 cutoff); the downward
    trndir/trntdr/rdndif/trndif at the top and bottom interfaces of band 1
    plus one guarded-interface pin per truncated band (snow band2 iface 5
    and band3 iface 2, bare band3 iface 6 - the dump value 0 exercises the
    trmin guard branch; the iface-1 rdndif boundary value 0 is skipped);
    the upward rupdir/rupdif at the top of the column (iface 14) of band 1
    and band 3 (the truncated-column reflectivity; the ocean-interface
    iface-1 values equal the band albodr constant and are skipped);
    dfdir/dfdif at ifaces {1,14} of band 1; the 11 scalar accumulations
    avdr/avdf/aidr/aidf/fsfc/fint/fthru/fthruvdr/fthruvdf/fthruidr/fthruidf
    (band-3 row); Sabs snow layers {1,2,3} (snow group; Sabs is identically
    zero for bare ice, so it is skipped there), Iabs ice layers {1,2,8} and
    fthrul levels {1,5,9} (band-3 row).
  The unfailable pins are not emitted: the extractor skips them with
  tolerance-aware copy/zero checks (|value - reference| <= rel*|value|,
  matching the assertion tolerance) instead of asserting them - rupdir/
  rupdif at iface 1 (a copy of the band ocean albedo albodr), rdndif at
  iface 1 (the boundary value 0) and the identically-zero dzk layers 1-2
  and tau/w0/g at layer 1 of the bare group (the hs = 0 zero region).

Tolerances: rel 1e-12 for the layer/table/accumulated quantities; rel 1e-9
for the interface quantities trndir/trntdr/rdndif/trndif/rupdir/rupdif/
dfdir/dfdif (round-off amplification in the near-pole snow band-3 algebra -
the dump build evaluates mu*rdr*gwt as mu*(rdr*gwt) and e*e*mu*mu as
(e*e)*(mu*mu); the ~2e-10 relative excursions appear only in the band-3
guarded-interface chain).

Output: the replacement text for the `"tests": []` line of
stage2/cice_shortwave.esm (same 1-space-per-level JSON layout, base indent
3), written to scratch/dedd_tests_block.txt; the tuple values are spliced
into the hand-authored tests block by hand (file tools only).
"""
import csv
import json
import re
import sys
from pathlib import Path

# gfortran prints denormal reals without the 'E' (e.g. 5.58...-322); fix that
DENORM = re.compile(r"([0-9])([+-])([0-9]{3})$")


def tofloat(s):
    s = s.strip()
    m = DENORM.search(s)
    if m:
        s = s[: m.start(1) + 1] + "E" + m.group(2) + m.group(3)
    return float(s)


def is_copy(v, ref, rel=1e-9):
    """True when the dump value v equals the reference constant ref within
    the assertion tolerance (|v - ref| <= rel*|v|): an assertion of such a
    value could not fail (it is a copy of a constant), so the pin is not
    emitted."""
    return abs(v - ref) <= rel * abs(v) if v != 0.0 else v == ref


def is_zero(v):
    """True when the dump value v is identically zero (|v| <= rel*|v| only
    when v == 0.0): an assertion of a constant-zero region could not fail,
    so the pin is not emitted."""
    return v == 0.0


PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path("/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM"
         "/stage1/dumps/cice/cice_dedd.csv")
OUT = Path("scratch/dedd_tests_block.txt")

# scalar (per-cell) input columns, dumped identically on all three band rows
SCALAR_IN = ["hs", "hi", "hp", "fnidr", "coszen", "swvdr", "swvdf", "swidr",
             "swidf", "kalg"]
# effective sea-ice/ponded-ice IOPs, dumped per band row (the Fortran
# ki_ssl(ns) etc. arrays, icepack_shortwave_data.F90:235-257)
BAND_IN = ["ki_ssl", "wi_ssl", "gi_ssl", "ki_dl", "wi_dl", "gi_dl", "ki_int",
           "wi_int", "gi_int", "ki_p_ssl", "wi_p_ssl", "gi_p_ssl",
           "ki_p_int", "wi_p_int", "gi_p_int"]
SNOW_IN = ["rsnw", "rhosnw"]           # per (snow_layer) columns rsnw_1..3
AERO_IN = [f"aero_mp_{i}" for i in range(1, 13)]

# per-assertion tolerance for the interface flux quantities (rel 1e-9, see
# module docstring); everything else uses the test-level rel 1e-12
IFACE_VARS = ("trndir", "trntdr", "rdndif", "trndif", "rupdir", "rupdif",
              "dfdir", "dfdif")

GROUPS = [
    {
        "id": "cice_shortwave_snow",
        "ncat": 2,
        "srftyp": 1,
        "desc": ("Snow-covered category (step 1, i=2, j=4, ncat 2, srftyp 1): "
                 "hs = 0.15364 m, hi = 1.2089 m, frsnw = 290.36 um and "
                 "aerosols in all four slots, so every band runs the full "
                 "snow/ice column with the SNICAR snow IOPs, the chlorophyll "
                 "absorption and the bulk aerosol mixing; band 3 truncates "
                 "at the interface below the snow SSL (trntdr < trmin = 1e-3 "
                 "there), so the layer props of layers 2..13 are guard-zeroed "
                 "and feed the upward recurrence as zeros. The pins cover "
                 "the layer geometry dzk and the grain-radius adjustment "
                 "frsnw, the snow-table IOPs at the SSL and the bottom snow "
                 "layer, the mixed tau/w0/g at the sea-ice delta layer "
                 "(the kfrsnl fresnel layer), the "
                 "delta-Eddington layer props at one active layer per band "
                 "(including the band-3 truncated-column top layer), the "
                 "downward interface chain at the top and bottom of band 1 "
                 "plus the first guarded interface of bands 2 and 3 "
                 "(expected 0, exercising the trmin guard), the upward "
                 "chain at the top of the column (the ocean-interface "
                 "values equal the band albodr constant and are not "
                 "pinned), the dfdir/dfdif flux fractions, the "
                 "band-weighted scalar accumulations and the Sabs/Iabs/"
                 "fthrul layer absorptions. Tolerances: rel 1e-9 on the "
                 "interface quantities (round-off amplification of the "
                 "near-pole band-3 algebra, see README-cice-shortwave.md), "
                 "rel 1e-12 elsewhere."),
        "dzk_layers": [1, 2, 5, 6, 13],
        "tau_layers": [6],
        "m_pins": {1: 6, 2: 1, 3: 1},          # band -> layer
        "down_pins": {1: [1, 14], 2: [5, 6], 3: [2, 3]},   # band -> ifaces
        # the band-2/3 pins straddle the trmin cutoff: the last pass-through
        # value (band2 iface 5 = 6.1e-4, band3 iface 2 = 1.9e-4, both below
        # trmin but still emitted because their top interface is active) and
        # the first guard zero one interface deeper
        "pass_through_pins": [(2, 5), (3, 2)],
        "guard_pins": [(2, 6), (3, 3)],  # (band, iface) below the cutoff
        "up_pins": {1: [1, 14], 3: [1, 14]},
        "df_pins": [1, 14],
        "snow_iops": True,
        "Sabs": True,
    },
    {
        "id": "cice_shortwave_bare",
        "ncat": 1,
        "srftyp": 0,
        "desc": ("Bare-ice category (same cell, ncat 1, srftyp 0): hs = 0 so "
                 "the snow branch is inactive and the column is the sea-ice "
                 "SSL/delta-layer/interior stack with the unconditional ice "
                 "aerosol mixing and the ocean albedo; band 3 truncates at "
                 "the interface below the sea-ice SSL. The pins cover the "
                 "layer geometry (the snow-layer dzk values are identically "
                 "zero and not pinned), frsnw, the mixed tau/w0/g at the "
                 "delta layer (the SSL values are identically zero and not "
                 "pinned), the delta-Eddington layer props at one active "
                 "layer per band (band 3 at the last layer above the "
                 "cutoff), the downward chain at the top and bottom of "
                 "band 1 plus the first guarded interface of band 3 "
                 "(expected 0), the upward chain at the top of the column, "
                 "the dfdir/dfdif fractions, the scalar accumulations and "
                 "the Iabs/fthrul layer absorptions. The snow IOPs "
                 "Qs/ws_snw/gs_snw/ks_snw and Sabs are not pinned: the dump "
                 "records zeros for srftyp = 0 rows and Sabs is identically "
                 "zero for bare ice. Tolerances: rel 1e-9 on the interface "
                 "quantities, rel 1e-12 elsewhere."),
        "dzk_layers": [1, 2, 5, 6, 13],
        "tau_layers": [1, 6],
        "m_pins": {1: 6, 2: 6, 3: 5},
        "down_pins": {1: [1, 14], 3: [6, 7]},
        "pass_through_pins": [(3, 6)],
        "guard_pins": [(3, 7)],
        "up_pins": {1: [1, 14], 3: [1, 14]},
        "df_pins": [1, 14],
        "snow_iops": False,
        "Sabs": False,
    },
]

SCALARS = ["avdr", "avdf", "aidr", "aidf", "fsfc", "fint", "fthru",
           "fthruvdr", "fthruvdf", "fthruidr", "fthruidf"]


def main():
    with open(PATH) as f:
        f.readline()
        cols = f.readline().strip().split(" = ", 1)[1].split(",")
        rows = list(csv.DictReader(f, fieldnames=cols))

    tests = []
    for g in GROUPS:
        pick = [r for r in rows
                if (r["step"].strip() == "1" and r["i"].strip() == "2"
                    and r["j"].strip() == "4"
                    and r["ncat"].strip() == str(g["ncat"])
                    and r["srftyp"].strip() == str(g["srftyp"]))]
        assert len(pick) == 3, (g["id"], len(pick))
        bands = {int(r["ns"].strip()): r for r in pick}
        b1, b2, b3 = bands[1], bands[2], bands[3]
        # the band rows of a triple must carry identical scalar inputs
        for c in SCALAR_IN + AERO_IN:
            vals = {tofloat(bands[ns][c]) for ns in (1, 2, 3)}
            assert len(vals) == 1, (g["id"], c)
        for c in SNOW_IN:
            for i in (1, 2, 3):
                vals = {tofloat(bands[ns][f"{c}_{i}"]) for ns in (1, 2, 3)}
                assert len(vals) == 1, (g["id"], c, i)

        overrides = {"srftyp": g["srftyp"]}
        for c in SCALAR_IN:
            overrides[c] = tofloat(b1[c])
        for c in BAND_IN:
            overrides[c] = [tofloat(bands[ns][c]) for ns in (1, 2, 3)]
        for c in SNOW_IN:
            overrides[c] = [tofloat(b1[f"{c}_{i}"]) for i in (1, 2, 3)]
        # dump layout is species-major: aero_mp_{(sp-1)*4+1..4} = species sp
        # in the snow SSL, snow interior, ice SSL and ice interior slots
        # (the Fortran loop `do na = 1, 4*n_aero, 4` with slots na+0..na+3,
        # F90:1476-1485, 2451-2547, 2605-2717); the model axis is
        # [aero (species), aero_slot]
        overrides["aero_mp"] = [[tofloat(b1[f"aero_mp_{(sp - 1) * 4 + slot}"])
                                 for slot in range(1, 5)] for sp in (1, 2, 3)]

        A = []

        def add(var, val, coords=None, rel=None):
            a = {"variable": var, "time": 0.0, "expected": val}
            if coords:
                a["coords"] = coords
            if rel is not None:
                a["tolerance"] = {"rel": rel, "abs": 0.0}
            A.append(a)

        def band_row(ns):
            return bands[ns]

        # layer geometry and grain radii; skip identically-zero layers (the
        # bare-group snow layers with hs = 0, which could not fail)
        for k in g["dzk_layers"]:
            v = tofloat(b1[f"dzk_{k - 1}"])
            if is_zero(v):
                continue
            add("dzk", v, {"layer": k})
        for sl in (1, 2, 3):
            add("frsnw", tofloat(b1[f"frsnw_{sl}"]), {"snow_layer": sl})

        # snow-table IOPs (snow group only)
        if g["snow_iops"]:
            for sl in (1, 3):
                for ns, b in ((1, b1), (2, b2), (3, b3)):
                    for var, col in (("Qs", "Qs"), ("ws_snw", "ws"),
                                     ("gs_snw", "gs"), ("ks_snw", "ks")):
                        add(var, tofloat(b[f"{col}_{sl}"]),
                            {"band": ns, "snow_layer": sl})

        # mixed layer IOPs (after the aerosol mixing); skip identically-zero
        # values (the bare-group layer 1 with hs = 0, which could not fail)
        for k in g["tau_layers"]:
            for ns, b in ((1, b1), (2, b2), (3, b3)):
                for var in ("tau", "w0", "g"):
                    v = tofloat(b[f"{var}_{k - 1}"])
                    if is_zero(v):
                        continue
                    add(var, v, {"band": ns, "layer": k})

        # guard-zeroed layer props feeding the upward recurrence
        m_vars = ("m_rdir", "m_tdir", "m_rdif_a", "m_rdif_b",
                  "m_tdif_a", "m_tdif_b", "m_trnlay")
        m_cols = ("rdir", "tdir", "rdif_a", "rdif_b",
                  "tdif_a", "tdif_b", "trnlay")
        for ns, k in sorted(g["m_pins"].items()):
            b = band_row(ns)
            for var, col in zip(m_vars, m_cols):
                add(var, tofloat(b[f"{col}_{k - 1}"]),
                    {"band": ns, "layer": k})

        # downward interface chain (rel 1e-9); skip the iface-1 rdndif
        # boundary value 0 (identically zero, could not fail)
        for ns, ks in sorted(g["down_pins"].items()):
            b = band_row(ns)
            for k in ks:
                for var in ("trndir", "trntdr", "rdndif", "trndif"):
                    v = tofloat(b[f"{var}_{k - 1}"])
                    if var == "rdndif" and k == 1 and is_zero(v):
                        continue
                    add(var, v, {"band": ns, "iface": k}, rel=1e-9)

        # upward chain: model iface k = Fortran iface 14-k; skip the
        # ocean-interface (iface 1) values, a copy of the band ocean albedo
        # albodr (0.01 for band 1, 0 otherwise), which could not fail
        for ns, ks in sorted(g["up_pins"].items()):
            b = band_row(ns)
            albodr = 0.01 if ns == 1 else 0.0
            for k in ks:
                for var in ("rupdir", "rupdif"):
                    v = tofloat(b[f"{var}_{14 - k}"])
                    if k == 1 and is_copy(v, albodr):
                        continue
                    add(var, v, {"band": ns, "iface": k}, rel=1e-9)

        # down-up flux fractions
        for k in g["df_pins"]:
            for var in ("dfdir", "dfdif"):
                add(var, tofloat(b1[f"{var}_{k - 1}"]),
                    {"band": 1, "iface": k}, rel=1e-9)

        # scalar accumulations (final on the band-3 row)
        for var in SCALARS:
            add(var, tofloat(b3[var]))

        if g["Sabs"]:
            for sl in (1, 2, 3):
                add("Sabs", tofloat(b3[f"Sabs_{sl}"]), {"snow_layer": sl})
        for il in (1, 2, 8):
            add("Iabs", tofloat(b3[f"Iabs_{il}"]), {"ice_layer": il})
        for pl in (1, 5, 9):
            add("fthrul", tofloat(b3[f"fthrul_{pl}"]), {"pen_level": pl})

        # curation sanity: guard pins must be exactly zero, active pins
        # nonzero, and every asserted variable must differ between the two
        # groups where that is expected (the two columns are distinct)
        for ns, k in g["guard_pins"]:
            v = tofloat(band_row(ns)[f"trntdr_{k - 1}"])
            assert v == 0.0, (g["id"], "guard pin not zero", ns, k, v)
        for ns, k in g["pass_through_pins"]:
            v = tofloat(band_row(ns)[f"trntdr_{k - 1}"])
            assert 0.0 < v < 0.001, (g["id"], "pass-through pin", ns, k, v)
        for var, col in zip(m_vars, m_cols):
            for ns, k in g["m_pins"].items():
                v = tofloat(band_row(ns)[f"{col}_{k - 1}"])
                assert v != 0.0, (g["id"], "active m_ pin is zero", var, ns, k)
        for ns, ks in g["up_pins"].items():
            albodr = 0.01 if ns == 1 else 0.0
            for k in ks:
                for var in ("rupdir", "rupdif"):
                    v = tofloat(band_row(ns)[f"{var}_{14 - k}"])
                    if k == 1:     # ocean interface: albodr per band, so
                        assert is_copy(v, albodr), \
                            (g["id"], "ocean up pin is not albodr",
                             var, ns, k, v)
                    else:          # top of the column: the actual albedo
                        assert v != 0.0, (g["id"], "top up pin is zero",
                                          var, ns, k)

        tests.append({
            "id": g["id"],
            "description": g["desc"],
            "time_span": {"start": 0.0, "end": 0.0},
            "tolerance": {"rel": 1e-12, "abs": 0.0},
            "parameter_overrides": overrides,
            "assertions": A,
        })
        print(f"{g['id']}: {len(A)} assertions")

    # emit the replacement text for the `   "tests": []` line, matching the
    # file's 1-space-per-level JSON layout with base indent 3
    lines = json.dumps(tests, indent=1).split("\n")
    block = ["   \"tests\": ["]
    block += ["    " + ln for ln in lines[1:-1]]
    block.append("   ]")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(block) + "\n")
    print(f"wrote {OUT} ({sum(len(t['assertions']) for t in tests)} "
          f"assertions total)")


if __name__ == "__main__":
    main()
