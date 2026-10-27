#!/usr/bin/env python3
"""Compare two CICE history files for solution neutrality.

Usage: python3 compare_iceh.py REF.nc INST.nc

Bit-compares all variables (data + attributes) between the reference run
(D_f19g17) and the instrumented run (D_cice_inst).  Global attributes that
are expected to differ (case name, timestamps, build metadata) are listed
in SKIP_GLOBAL and ignored; any *other* differing global attribute is
reported so the list can be extended consciously.
"""
import sys
import numpy as np
import netCDF4

# Global attrs that legitimately differ between the two cases.
SKIP_GLOBAL = {
    "case_title", "caseid", "source", "host", "username", "version",
    "inithist", "Version", "revision_id", "build_revision",
    "history_software_version", "history_creation_date", "history_script",
    "history_command", "model_doi_url", "git_version", "input_data_origin",
    # Case-name/timestamp attrs (differ by construction between runs)
    "title", "date_created", "history",
}

def main():
    ref_path, inst_path = sys.argv[1], sys.argv[2]
    ref, inst = netCDF4.Dataset(ref_path), netCDF4.Dataset(inst_path)

    nbad = 0
    # Global attributes
    ra, ia = ref.ncattrs(), inst.ncattrs()
    for a in sorted(set(ra) | set(ia)):
        if a in SKIP_GLOBAL:
            continue
        rv = ref.getncattr(a) if a in ra else "<absent>"
        iv = inst.getncattr(a) if a in ia else "<absent>"
        if rv != iv:
            print(f"GLOBAL ATTR DIFF: {a}: {rv!r} != {iv!r}")
            nbad += 1

    # Variables
    rv = set(ref.variables)
    iv = set(inst.variables)
    for v in sorted(rv ^ iv):
        print(f"VARIABLE SET DIFF: {v} only in "
              + ("ref" if v in rv else "inst"))
        nbad += 1
    for v in sorted(rv & iv):
        r, i = ref.variables[v], inst.variables[v]
        for a in sorted(set(r.ncattrs()) | set(i.ncattrs())):
            rva = r.getncattr(a) if a in r.ncattrs() else "<absent>"
            iva = i.getncattr(a) if a in i.ncattrs() else "<absent>"
            if rva != iva:
                print(f"  var {v} attr {a}: {rva!r} != {iva!r}")
                nbad += 1
        rd, id_ = r[...], i[...]
        if rd.dtype != id_.dtype or not np.array_equal(rd, id_):
            print(f"  var {v}: DATA DIFF (dtypes {rd.dtype} vs {id_.dtype}, "
                  f"equal={np.array_equal(rd, id_)}, max|d|="
                  f"{np.nanmax(np.abs(rd.astype(float) - id_.astype(float)))})")
            nbad += 1

    msg = ("PASS: bit-identical (modulo expected case metadata)"
           if nbad == 0 else f"FAIL: {nbad} differences")
    print(f"\n{msg}")
    return 0 if nbad == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
