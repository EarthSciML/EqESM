#!/bin/bash
# Build the standalone vdiff operator reference driver.
#
# -cpp -DUSE_CONTIGUOUS= is the CAM gfortran build's configuration (CAM
# turns the F2008 contiguous attribute off for gfortran,
# cam/doc/ChangeLog:84560): with -DUSE_CONTIGUOUS=contiguous the source
# `real(r8), USE_CONTIGUOUS intent(in)` lacks the comma the attribute
# requires and does not compile.
set -euo pipefail
cd "$(dirname "$0")"
gfortran -cpp -DUSE_CONTIGUOUS= -ffree-line-length-none -O2 -o vdiff_driver \
  shr_kind_mod.F90 shr_log_mod.F90 shr_sys_mod.F90 \
  coords_1d.F90 linear_1d_operators.F90 main_vdiff_driver.F90
./vdiff_driver > reference_output.txt
echo "built and ran vdiff_driver; reference_output.txt updated"
