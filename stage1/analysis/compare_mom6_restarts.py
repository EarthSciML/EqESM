#!/usr/bin/env python3
"""Bit-compare MOM6 run outputs between the reference and instrumented runs.

Usage: python3 compare_mom6_restarts.py REF_RUNDIR INST_RUNDIR [--files F1 F2 ...]

Compares, file by file, all netCDF variables (exact bit equality of data)
and all attributes (global and per-variable), for the run-output files that
must be solution-neutral: mom6.r, mom6.r_stoch, cpl.r, datm.r, dice.r,
drof.r, ocean.stats.nc.  Global attributes that legitimately differ
(case name, timestamps, build metadata) are skipped consciously; any other
difference is reported.

By default the two run dirs are compared on the intersection of their file
lists, or the user can pass an explicit file list (names relative to the
run dirs).
"""
import argparse
import os
import sys

import netCDF4

SKIP_GLOBAL = {
    "case_title", "caseid", "source", "host", "username", "version",
    "inithist", "Version", "revision_id", "build_revision",
    "history_software_version", "history_creation_date", "history_script",
    "history_command", "model_doi_url", "git_version", "input_data_origin",
    "title", "date_created", "history", "filename", "experiment",
    "ice_present", "ocn_present", "atm_present", "lnd_present", "rof_present",
    "glc_present", "wav_present", "esp_present",
}


def compare_file(ref_path, inst_path, nbad):
    ref, inst = netCDF4.Dataset(ref_path), netCDF4.Dataset(inst_path)
    for a in sorted(set(ref.ncattrs()) | set(inst.ncattrs())):
        if a in SKIP_GLOBAL:
            continue
        rv = ref.getncattr(a) if a in ref.ncattrs() else "<absent>"
        iv = inst.getncattr(a) if a in inst.ncattrs() else "<absent>"
        if rv != iv:
            print(f"  {os.path.basename(ref_path)} GLOBAL ATTR {a}: {rv!r} != {iv!r}")
            nbad[0] += 1
    rv = set(ref.variables)
    iv = set(inst.variables)
    for v in sorted(rv ^ iv):
        print(f"  {os.path.basename(ref_path)} VARIABLE SET: {v} only in "
              + ("ref" if v in rv else "inst"))
        nbad[0] += 1
    for v in sorted(rv & iv):
        r, i = ref.variables[v], inst.variables[v]
        for a in sorted(set(r.ncattrs()) | set(i.ncattrs())):
            rva = r.getncattr(a) if a in r.ncattrs() else "<absent>"
            iva = i.getncattr(a) if a in i.ncattrs() else "<absent>"
            if rva != iva:
                print(f"  {os.path.basename(ref_path)} var {v} attr {a}: {rva!r} != {iva!r}")
                nbad[0] += 1
        rd, id_ = r[...], i[...]
        if rd.dtype != id_.dtype or not (rd == id_).all():
            print(f"  {os.path.basename(ref_path)} var {v}: DATA DIFF "
                  f"(dtypes {rd.dtype} vs {id_.dtype})")
            nbad[0] += 1
    ref.close()
    inst.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref_dir")
    ap.add_argument("inst_dir")
    ap.add_argument("--files", nargs="+", default=None,
                    help="file names (relative to the run dirs); default: "
                         "intersection of the two dirs")
    args = ap.parse_args()

    if args.files is None:
        common = sorted(set(os.listdir(args.ref_dir)) & set(os.listdir(args.inst_dir)))
        files = [f for f in common if f.endswith(".nc") or f == "ocean.stats"]
    else:
        files = args.files

    nbad = [0]
    nfiles = 0
    for f in files:
        rp, ip = os.path.join(args.ref_dir, f), os.path.join(args.inst_dir, f)
        if not (os.path.exists(rp) and os.path.exists(ip)):
            print(f"SKIP {f}: missing in one dir")
            continue
        print(f"comparing {f}")
        compare_file(rp, ip, nbad)
        nfiles += 1

    print(f"\n{nfiles} files compared, {nbad[0]} differences")
    sys.exit(1 if nbad[0] else 0)


if __name__ == "__main__":
    main()
