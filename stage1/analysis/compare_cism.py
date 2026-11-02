#!/usr/bin/env python3
"""Compare CISM output files between the reference and instrumented runs.

Usage: python3 compare_cism.py REF_DIR INST_DIR [PREFIX]

Bit-compares every matching file in the two directories (default prefix
T1850Gg_cism_inst.cism.gris.*, i.e. the cism.h0i history files and the
cism restart file) between the dump-off reference run (run_ref_dumpoff)
and the dump-on instrumented run (run).  Global and variable attributes
that are expected to differ (case name, timestamps, build metadata) are
listed in SKIP_GLOBAL / SKIP_VAR_ATTR and ignored; any *other* differing
attribute is reported so the list can be extended consciously.
"""
import glob
import os
import sys

import netCDF4
import numpy as np

SKIP_GLOBAL = {
    "case_title", "caseid", "source", "host", "username", "version",
    "inithist", "Version", "revision_id", "build_revision",
    "history_software_version", "history_creation_date", "history_script",
    "history_command", "model_doi_url", "git_version", "input_data_origin",
    # Case-name/timestamp attrs (differ by construction between runs)
    "title", "date_created", "history",
}

# Variable attributes that can legitimately differ between runs
SKIP_VAR_ATTR = {"time_origin"}


def compare_file(ref_path, inst_path, label):
    nbad = 0
    ref, inst = netCDF4.Dataset(ref_path), netCDF4.Dataset(inst_path)

    ra, ia = ref.ncattrs(), inst.ncattrs()
    for a in sorted(set(ra) | set(ia)):
        if a in SKIP_GLOBAL:
            continue
        rv = ref.getncattr(a) if a in ra else "<absent>"
        iv = inst.getncattr(a) if a in ia else "<absent>"
        if rv != iv:
            print(f"  [{label}] GLOBAL ATTR DIFF: {a}: {rv!r} != {iv!r}")
            nbad += 1

    rv = set(ref.variables)
    iv = set(inst.variables)
    for v in sorted(rv ^ iv):
        print(f"  [{label}] VARIABLE SET DIFF: {v} only in "
              + ("ref" if v in rv else "inst"))
        nbad += 1
    for v in sorted(rv & iv):
        r, i = ref.variables[v], inst.variables[v]
        for a in sorted(set(r.ncattrs()) | set(i.ncattrs())):
            if a in SKIP_VAR_ATTR:
                continue
            rva = r.getncattr(a) if a in r.ncattrs() else "<absent>"
            iva = i.getncattr(a) if a in i.ncattrs() else "<absent>"
            if rva != iva:
                print(f"  [{label}] var {v} attr {a}: {rva!r} != {iva!r}")
                nbad += 1
        rd, id_ = r[...], i[...]
        if rd.dtype != id_.dtype or not np.array_equal(rd, id_):
            mx = np.nanmax(np.abs(rd.astype(float) - id_.astype(float))) \
                if rd.size else 0.0
            print(f"  [{label}] var {v}: DATA DIFF "
                  f"(dtypes {rd.dtype} vs {id_.dtype}, "
                  f"equal={np.array_equal(rd, id_)}, max|d|={mx})")
            nbad += 1
    return nbad


def main():
    ref_dir, inst_dir = sys.argv[1], sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) > 3 else "T1850Gg_cism_inst.cism.gris.*"

    ref_files = sorted(glob.glob(os.path.join(ref_dir, prefix)))
    inst_files = sorted(glob.glob(os.path.join(inst_dir, prefix)))
    if not ref_files:
        print(f"no reference files matching {os.path.join(ref_dir, prefix)}")
        return 2
    ref_names = {os.path.basename(p) for p in ref_files}
    inst_names = {os.path.basename(p) for p in inst_files}
    nbad = 0
    for n in sorted(ref_names - inst_names):
        print(f"FILE SET DIFF: {n} only in reference")
        nbad += 1
    for n in sorted(inst_names - ref_names):
        print(f"FILE SET DIFF: {n} only in instrumented")
        nbad += 1
    for n in sorted(ref_names & inst_names):
        rf, inf = os.path.join(ref_dir, n), os.path.join(inst_dir, n)
        nbad += compare_file(rf, inf, n)

    msg = ("PASS: bit-identical (modulo expected case metadata)"
           if nbad == 0 else f"FAIL: {nbad} differences")
    print(f"\n{msg}")
    return 0 if nbad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
