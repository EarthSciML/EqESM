#!/usr/bin/env python3
"""Extract the step-1 vdiff data block for main_vdiff_driver.F90.

Mechanical transcription (AGENTS.md: scripts MAY translate numeric test
tuples from dumps): reads the step-1 rows of the instrumented run's
cam_vdiff_prof.csv and emits the three Fortran DATA statements embedded
in the driver -- pint (interface pressures, CAM order), dpidz_sq, and
u0 = u1 - ptend_u*ztodt (the entering wind, CAM order). Values are
printed with repr() (shortest round-trip), the same precision the .esm
files use.

Usage: python3 extract_step1_data.py cam_vdiff_prof.csv
"""

import sys

ZTODT = 1800.0  # physics timestep, s


def main(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            rows.append(line.strip().split(","))
    r1 = [r for r in rows if int(r[0]) == 1]
    r1.sort(key=lambda r: int(r[3]))  # CAM lev 1 = model top
    assert len(r1) == 32
    pint = [float(r[5]) for r in r1] + [float(r1[31][6])]
    dpidz = [float(r[11]) for r in r1] + [float(r1[31][12])]
    u0 = [float(r[25]) - float(r[29]) * ZTODT for r in r1]

    def lines(name, values):
        # _r8 kind suffix: a bare literal in a DATA statement is default
        # (single) kind and is rounded before conversion to real(r8).
        vals = ", ".join(f"{v!r}_r8" for v in values)
        return f"  data {name} / {vals} /"

    print("! pint (33), dpidz_sq (33), u0 (32) -- step 1, CAM order:")
    for name, values in (("pint", pint), ("dpidz_sq", dpidz), ("u0", u0)):
        print(lines(name, values))


if __name__ == "__main__":
    main(sys.argv[1])
