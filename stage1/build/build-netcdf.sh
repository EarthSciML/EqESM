#!/bin/bash
# Build the netCDF stack (hdf5 -> netcdf-c -> netcdf-fortran) from source
# into stage1/build/local, following the EqIce local-deps pattern.
#
# The cluster has no netCDF module and no system netcdf install (see
# stage1/runs/TOOLCHAIN.md). CESM/CIME needs netCDF-C and netCDF-Fortran;
# pnetcdf is optional (PIO2 runs netcdf-only) and skipped for now.
#
# Toolchain: gcc 13.3.0 + openmpi 5.0.1 (the sibling-projects' pair).
# Serial hdf5 (no --enable-parallel) is sufficient for smoke runs; revisit
# if parallel I/O performance becomes a need.
#
# Usage: bash stage1/build/build-netcdf.sh [jobs]

set -euo pipefail

JOBS="${1:-8}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PREFIX="$REPO_ROOT/stage1/build/local"
SRC_DIR="$REPO_ROOT/stage1/build/local-src"
LOG_DIR="$REPO_ROOT/stage1/build"
mkdir -p "$SRC_DIR"

export PATH="$PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:${LD_LIBRARY_PATH:-}"
export CPPFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib"

# --- hdf5 (serial, with Fortran + high-level) ---
HDF5_VER="1.14.6"
HDF5_TGZ="hdf5-$HDF5_VER.tar.gz"
if [[ ! -f "$PREFIX/include/hdf5.h" ]]; then
  echo "==> Building hdf5 $HDF5_VER..."
  cd "$SRC_DIR"
  if [[ ! -d "hdf5-$HDF5_VER" ]]; then
    # support.hdfgroup.org now 404s this path; use the GitHub release asset.
    curl -sL -o "$HDF5_TGZ" \
      "https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/$HDF5_TGZ"
    tar xzf "$HDF5_TGZ"
  fi
  cd "hdf5-$HDF5_VER"
  ./configure --prefix="$PREFIX" --enable-fortran --enable-hl \
      >"$LOG_DIR/hdf5-configure.log" 2>&1
  make -j "$JOBS" >"$LOG_DIR/hdf5-make.log" 2>&1
  make install >"$LOG_DIR/hdf5-install.log" 2>&1
fi

# --- netcdf-c ---
NC_VER="4.9.3"
NC_TGZ="netcdf-c-$NC_VER.tar.gz"
if [[ ! -f "$PREFIX/include/netcdf.h" ]]; then
  echo "==> Building netcdf-c $NC_VER..."
  cd "$SRC_DIR"
  if [[ ! -d "netcdf-c-$NC_VER" ]]; then
    curl -sL -o "$NC_TGZ" \
      "https://downloads.unidata.ucar.edu/netcdf-c/$NC_VER/$NC_TGZ"
    tar xzf "$NC_TGZ"
  fi
  cd "netcdf-c-$NC_VER"
  ./configure --prefix="$PREFIX" --disable-dap --disable-doxygen \
      >"$LOG_DIR/netcdf-c-configure.log" 2>&1
  make -j "$JOBS" >"$LOG_DIR/netcdf-c-make.log" 2>&1
  make install >"$LOG_DIR/netcdf-c-install.log" 2>&1
fi

# --- netcdf-fortran ---
NF_VER="4.6.1"
NF_TGZ="netcdf-fortran-$NF_VER.tar.gz"
if [[ ! -f "$PREFIX/include/netcdf.inc" ]]; then
  echo "==> Building netcdf-fortran $NF_VER..."
  cd "$SRC_DIR"
  if [[ ! -d "netcdf-fortran-$NF_VER" ]]; then
    curl -sL -o "$NF_TGZ" \
      "https://downloads.unidata.ucar.edu/netcdf-fortran/$NF_VER/$NF_TGZ"
    tar xzf "$NF_TGZ"
  fi
  cd "netcdf-fortran-$NF_VER"
  ./configure --prefix="$PREFIX" \
      >"$LOG_DIR/netcdf-fortran-configure.log" 2>&1
  make -j "$JOBS" >"$LOG_DIR/netcdf-fortran-make.log" 2>&1
  make install >"$LOG_DIR/netcdf-fortran-install.log" 2>&1
fi

echo "==> netCDF stack installed under $PREFIX"
"$PREFIX/bin/nf-config" --all | head -12
