# Shared runtime environment for CESM Stage-1 runs.
# Source this from run scripts:  source stage1/runs/env.sh
#
# Provides:
#   PATH/LD_LIBRARY_PATH - gcc 13.3.0 + openmpi 5.0.1 toolchain (the
#                          sibling-projects' pair), local netCDF stack
#   NETCDF_PATH          - for CIME machine configs / builds
#
# See stage1/runs/TOOLCHAIN.md for the probe that motivated this.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

MPI_PREFIX="/sw/apps/mpi/openmpi/5.0.1/gcc/13.3.0"
GCC_PREFIX="/sw/apps/gcc/13.3.0"
LOCAL_PREFIX="$REPO_ROOT/stage1/build/local"

export PATH="$MPI_PREFIX/bin:$GCC_PREFIX/bin:$LOCAL_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_PREFIX/lib:$GCC_PREFIX/lib64:$LOCAL_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$LOCAL_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

export NETCDF_PATH="$LOCAL_PREFIX"
