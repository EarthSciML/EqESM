#!/bin/bash
# Build ESMF 8.9.1 from source into stage1/build/esmf, following the
# netCDF local-deps pattern (build-netcdf.sh).
#
# CMEPS (the CESM 3.0 mediator) requires ESMF via find_package(ESMF REQUIRED)
# and the cmeps buildnml reads $ESMFMKFILE (see
# components/cmeps/cime_config/buildnml:735). NCAR machines load an
# esmf-mpi/8.9.1 module; the cluster has no ESMF module, so we build it.
# Derecho's pin is 8.9.1 — use the same.
#
# Toolchain: gcc 13.3.0 + openmpi 5.0.1 (the sibling-projects' pair), with
# the hand-built netCDF stack in stage1/build/local on the path.
#
# Usage: bash stage1/build/build-esmf.sh [jobs]

set -euo pipefail

JOBS="${1:-8}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PREFIX="$REPO_ROOT/stage1/build/esmf"
SRC_DIR="$REPO_ROOT/stage1/build/local-src"
LOG_DIR="$REPO_ROOT/stage1/build"
mkdir -p "$SRC_DIR"

export ESMF_DIR="$SRC_DIR/esmf-8.9.1"
export ESMF_INSTALL_PREFIX="$PREFIX"
export ESMF_COMPILER="gfortran"
export ESMF_COMM="openmpi"
export ESMF_OS="Linux"
export ESMF_ABI="64"
export ESMF_BOPT="O"
export ESMF_OPTLEVEL="2"
# Static libs are what CIME/CMEPS link against; build only the static lib
# (no shared) to keep the install small and self-contained.
export ESMF_SHARED_LIBRARIES="OFF"
export ESMF_NETCDF="split"
export ESMF_NETCDF_INCLUDE="$REPO_ROOT/stage1/build/local/include"
export ESMF_NETCDF_LIBPATH="$REPO_ROOT/stage1/build/local/lib"
# ESMF_PNETCDF must be UNSET, not "OFF": common.mk does `ifdef ESMF_PNETCDF`
# to add -DESMF_PNETCDF=1, so any defined value (even "OFF") enables it and
# the IO code then #includes <pnetcdf.h>. No pnetcdf on the cluster.
unset ESMF_PNETCDF || true
export ESMF_MPIRUN="mpirun"
export ESMF_MPIMPMDRUN="mpirun"

if [[ ! -d "$ESMF_DIR" ]]; then
  echo "==> Downloading ESMF 8.9.1..."
  cd "$SRC_DIR"
  curl -sL -o esmf-8.9.1.tar.gz \
    "https://github.com/esmf-org/esmf/archive/refs/tags/v8.9.1.tar.gz"
  tar xzf esmf-8.9.1.tar.gz
fi

if [[ ! -f "$PREFIX/lib/libesmf.a" ]]; then
  echo "==> Building ESMF 8.9.1 (this takes a while)..."
  cd "$ESMF_DIR"
  make -j "$JOBS" >"$LOG_DIR/esmf-make.log" 2>&1
  make install >"$LOG_DIR/esmf-install.log" 2>&1
fi

echo "==> ESMF installed under $PREFIX"
ls "$PREFIX/lib" | head -5
ESMFMKFILE_PATH="$(find "$PREFIX/lib" -name esmf.mk | head -1)"
echo "ESMFMKFILE=$ESMFMKFILE_PATH"
