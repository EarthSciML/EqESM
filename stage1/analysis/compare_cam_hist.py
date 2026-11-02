#!/usr/bin/env python3
"""Compare two CAM history files for solution neutrality.

Usage: python3 compare_cam_hist.py REF.nc INST.nc [more.nc ...]

Bit-compares all variables (data + attributes) between the reference run
(FSCAMARM97_ne3) and the instrumented run (FSCAMARM97_ne3_cam_inst).
Global attributes that legitimately differ (case name, timestamps, build
metadata) are listed in SKIP_GLOBAL; any *other* differing global
attribute is reported so the list can be extended consciously.

Pass as many instrumented files as needed; each is compared against the
reference file in the same position.  A directory argument is expanded
to its *.nc files (sorted); then every file is compared pairwise by
sorted filename (the two cases use the same history-file naming).
"""
import sys
import glob
import os
import numpy as np
import netCDF4

# Global attrs that legitimately differ between the two cases.
SKIP_GLOBAL = {
    "case_title", "caseid", "source", "host", "username", "version",
    "inithist", "Version", "revision_id", "build_revision",
    "history_software_version", "history_creation_date", "history_script",
    "history_command", "model_doi_url", "git_version", "input_data_origin",
    "title", "date_created", "history", "logname",
    # CESM/NUOPC provenance attrs
    "CASEROOT", "CASE", "SRCROOT", "OBJROOT", "BATCH_SYSTEM", "COMPILER",
    "DEBUG", "MPILIB", "PROJECT", "CHARGE_ACCOUNT", "JOBID", "USER",
    "EXEROOT", "RUNDIR", "DOUT_S", "DOUT_S_ROOT", "DOUT_S_SAVE_TIMING",
    "DOUT_S_SAVE_COUPLER_RESTART", "DOUT_S_SAVE_ATMO_RESTART",
    "DOUT_S_SAVE_ICE_RESTART", "DOUT_S_SAVE_LND_RESTART",
    "DOUT_S_SAVE_OCN_RESTART", "DOUT_S_SAVE_ROF_RESTART",
    "DOUT_S_SAVE_GLC_RESTART", "DOUT_S_SAVE_WAV_RESTART",
    "DOUT_S_SAVE_ESP_RESTART", "DOUT_S_SAVE_PIO",
    "DOUT_S_SAVE_HISTORY", "DOUT_L_MS", "DOUT_S", "DOUT_S_IGNORE_ERRORS",
    "DOUT_S_FCNTL", "DOUT_S_FCNTL_LND", "DOUT_S_FCNTL_ICE",
    "DOUT_S_FCNTL_OCN", "DOUT_S_FCNTL_ROF", "DOUT_S_FCNTL_GLC",
    "DOUT_S_FCNTL_WAV", "DOUT_S_FCNTL_CPL", "DOUT_S_FCNTL_ESP",
    "DOUT_S_FCNTL_ATM", "DOUT_S_FCNTL_HISTORY", "DOUT_S_FCNTL_PIO",
    "DOUT_S_ROOT", "DOUT_S_SAVE_CPL_RESTART", "DOUT_S_SAVE_TIMING",
    "DOUT_S_SAVE_COUPLER_RESTART",
}


def compare_one(ref_path, inst_path):
    """Compare one pair; return number of differences found."""
    nbad = 0
    try:
        ref = netCDF4.Dataset(ref_path)
        inst = netCDF4.Dataset(inst_path)
    except OSError as e:
        print(f"OPEN FAIL: {e}")
        return 1

    ra, ia = ref.ncattrs(), inst.ncattrs()
    for a in sorted(set(ra) | set(ia)):
        if a in SKIP_GLOBAL:
            continue
        rv = ref.getncattr(a) if a in ra else "<absent>"
        iv = inst.getncattr(a) if a in ia else "<absent>"
        if rv != iv:
            print(f"  GLOBAL ATTR DIFF: {a}: {rv!r} != {iv!r}")
            nbad += 1

    rv = set(ref.variables)
    iv = set(inst.variables)
    for v in sorted(rv ^ iv):
        print(f"  VARIABLE SET DIFF: {v} only in "
              + ("ref" if v in rv else "inst"))
        nbad += 1
    for v in sorted(rv & iv):
        r, i = ref.variables[v], inst.variables[v]
        if r.dimensions != i.dimensions:
            print(f"  var {v}: DIMENSION DIFF {r.dimensions} vs {i.dimensions}")
            nbad += 1
        for a in sorted(set(r.ncattrs()) | set(i.ncattrs())):
            rva = r.getncattr(a) if a in r.ncattrs() else "<absent>"
            iva = i.getncattr(a) if a in i.ncattrs() else "<absent>"
            if rva != iva:
                print(f"  var {v} attr {a}: {rva!r} != {iva!r}")
                nbad += 1
        rd, id_ = r[...], i[...]
        if rd.shape != id_.shape:
            print(f"  var {v}: SHAPE DIFF {rd.shape} vs {id_.shape}")
            nbad += 1
        elif rd.dtype != id_.dtype or not np.array_equal(rd, id_):
            if rd.dtype.kind in "fc" and id_.dtype.kind in "fc":
                # numeric: report max abs difference
                d = np.abs(rd.astype(float) - id_.astype(float))
                print(f"  var {v}: DATA DIFF (dtypes {rd.dtype} vs {id_.dtype}, "
                      f"equal={np.array_equal(rd, id_)}, max|d|="
                      f"{np.nanmax(d)})")
            else:
                # char/other: report a few differing elements
                neq = np.count_nonzero(rd != id_)
                print(f"  var {v}: DATA DIFF (dtypes {rd.dtype} vs {id_.dtype}, "
                      f"equal={np.array_equal(rd, id_)}, neq={neq})")
            nbad += 1
    return nbad


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return 1

    # Expand directories to their *.nc files, sorted.
    refs, insts = [], []
    for a in args:
        if os.path.isdir(a):
            files = sorted(glob.glob(os.path.join(a, "*.nc")))
            (refs if len(refs) <= len(insts) else insts).extend(files)
        elif len(refs) <= len(insts):
            refs.append(a)
        else:
            insts.append(a)

    if len(refs) != len(insts):
        print(f"Need equal numbers of ref/inst files: {len(refs)} vs {len(insts)}")
        return 1

    total = 0
    for rp, ip in zip(refs, insts):
        if os.path.basename(rp) != os.path.basename(ip):
            print(f"NOTE: filename mismatch {os.path.basename(rp)} vs "
                  f"{os.path.basename(ip)} - comparing anyway")
        print(f"comparing {os.path.basename(rp)} vs {os.path.basename(ip)}")
        total += compare_one(rp, ip)

    msg = ("PASS: bit-identical (modulo expected case metadata)"
           if total == 0 else f"FAIL: {total} differences")
    print(f"\n{msg}")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
