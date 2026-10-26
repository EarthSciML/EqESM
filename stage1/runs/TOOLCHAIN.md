# M1 toolchain probe — building CESM (CIME) on the Illinois Campus Cluster

Probed 2026-10-25 (interactive session, SLURM job 11812286, node ccc0233).
All probes were read-only. Where a probe errored or returned nothing, it is
recorded as such.

## 1. Summary of the environment

| Item | Finding |
|---|---|
| Cluster | Illinois Campus Cluster ("icc"), SLURM 26.05.5 |
| OS | Red Hat Enterprise Linux 9.8 (Plow), kernel 5.14.0-687.58.1.el9_8.x86_64 |
| Node profile | 40 cores, ~188 GB RAM (176992 MB in cgroup), 8 GPUs per node (current node ccc0233) |
| Module system | Lmod 9.4.2; `MODULEPATH=/sw/Core:/sw/spack:/sw/user/modulefiles:/usr/share/lmod/lmod/modulefiles/Core`. Site modulefiles under `/sw/user/modulefiles`, installs under `/sw/apps/...` (e.g. `/sw/apps/gcc/13.3.0`, `/sw/apps/mpi/openmpi/5.0.1/gcc/13.3.0`). Note: `/sw/spack` in MODULEPATH is an empty dir. |
| Scheduler | SLURM; partitions seen: `scavenger` (1 day), `secondary` (4 h), `ctessum` (3 days, nodes ccc0232–0236), `eng-research-gpu` (2 d), `ic-express` (8 h), `ic-highmem` (3 d), `IllinoisComputes` (3 d), `IllinoisComputes-GPU` (3 d) |
| SLURM account | `ctessum` (assocs: `campusclu`, `ctessum-c…`, `ctessum-ic`, `ctessum`; default QOS `scav_high`, MaxJobs 40). Current job: `--partition=ctessum --account=ctessum --qos=scav_high`. |
| Disk | `/projects`: Lustre "taiga", 22 PB total, 5.8 PB avail (72 % used). `/scratch.local`: 794 GB, 250 GB avail (69 % used); `/scratch.local/ctessum` exists and is owned by us. `/tmp` is RAM-backed — do not store builds there. |
| Network | `https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/` → HTTP 200 (SVN server up; `/trunk/README` 404s, so that exact path is wrong but the server is reachable). `https://web.lcrc.anl.gov/public/e3sm/inputdata/` → HTTP 200. `https://github.com` → HTTP 200. All CESM inputdata sources reachable. **`data.cesm.ucar.edu` does NOT resolve (DNS rc=2)** — point CIME at `svn-ccsm-inputdata` + the LCRC mirror only. **No `svn` client** on the login node — CIME svn-protocol inputdata fetches need `--protocol=https` (or install svn). |
| Existing CESM work | **None found.** Bounded searches for `*cesm*`, `*inputdata*`, `*cime*` under `/projects/illinois/eng/cee` (maxdepth 4, excluding ctessum) and `/scratch.local` (maxdepth 3) returned nothing. `/projects/illinois/eng/cee` holds ~30 group dirs (alipour, barros, ctessum, gardoni, …); no shared inputdata mirror. |

## 2. Probe results

### 2.1 Modules available (relevant subset of `module avail`)

Present:
- Compilers: `gcc/12.4.0`, `gcc/13.3.0` (default), `intel/compiler/2025.0.0`/`2025.0.4` (D), `intel/compiler-intel-llvm/2025.0.4`, `amd/aocc-compiler/5.0.0`
- MPI: `openmpi/5.0.1-gcc-13.3.0`, `mvapich2/2.3.7-gcc-13.3.0`, `intel/mpi/2021.14`
- Math: `intel/mkl/2025.0`, `gsl/2.8`, `fftw/3.3.10`, `gmp/6.3.0`, `mpfr/4.2.1`, `lapack`, `openblas`, `petsc`
- Languages/tools: `python/3.10.16`, `python/3.11.11`, `python/3.13.2` (D), `git/2.46.3`, `git/2.48.1` (D), `anaconda3/2024.10`, `miniforge3/25.9.1`, `texlive/2024`, `cuda/12.8` (D), `java/23`

Absent (checked via `module avail` grep and `module spider`): **netcdf, hdf5, pnetcdf, cmake, perl, nco, ncl, julia, apptainer, spack, singularity** — no modules for any of these.

### 2.2 Versions (recommended module set loaded: `gcc/13.3.0` + `openmpi/5.0.1-gcc-13.3.0`)

- `gcc` 13.3.0, `gfortran` 13.3.0
- `mpirun`/`mpiexec` — Open MPI 5.0.1 (`/sw/apps/mpi/openmpi/5.0.1/gcc/13.3.0/bin/`)
- `cmake` 3.31.8 — system `/usr/bin/cmake` (no module needed)
- `perl` 5.32.1 — system `/usr/bin/perl` (no module needed)
- `python3` 3.9.25 — system `/usr/bin/python3`; `python/3.13.2` module also available (`/sw/apps/python/3.13.2/bin/python3`)
- `apptainer` 1.5.4-1.el9 — system `/usr/bin/apptainer` (present, not a module). Runs **unprivileged** (user-namespace): no `starter-suid` helper under `/usr/libexec/apptainer/bin`, `unshare -U --map-root-user` works (rc=0), so image builds are viable without setuid — the EqWeFiC "setuid mode" note is stale or cluster-state-dependent.
- System python 3.9.25 already has CIME's Python deps: `lxml` 4.6.5, `pytest` 8.4.2, `numpy` 2.0.2, `pyyaml` (netCDF4/h5py are pip wheels — Python-only, not linkable).
- `git` 2.48.1 (module default; system git also exists)
- `ncdump`, `nc-config`, `nf-config`, `ncks`, `ncl`, `spack` — **none found in PATH**

### 2.3 netCDF availability

**No netCDF anywhere on the cluster.** Checked:
- `module avail | grep -i netcdf` → nothing; `module spider netcdf` → "Unable to find"
- `/usr/local` (bin/include/lib) → nothing; `/usr/local/bin` holds only pmix tools
- `/opt` → only mellanox, puppetlabs, slurm, xcat
- `/sw/apps` → has `hdf5/` but it contains only `src/` and `test/` (source, not an install); no netcdf dir anywhere under `/sw`
- Spack install (see §2.4) has no netcdf package built

Neither netCDF-C nor netCDF-Fortran is present. **netCDF must be built from source or obtained via container.**

### 2.4 Spack

- A full Spack source tree exists at `/sw/apps/spack/spack` (dated 2025-02-12, `bin/spack` present) and a modulefile tree at `/sw/apps/spack/2025-02-12/modulefiles/linux-rhel9-x86_64/{Core,aocc,intel,openmpi}` (only Core bootstrap packages: gcc, openmpi, gmake, libtool, …; `opt/spack` holds only `aocc-5.0.0`).
- `/sw` and `/sw/apps/spack` are **not writable** by us, so the shared spack cannot be extended. The MODULEPATH entry `/sw/spack` is an empty directory.
- Conclusion: treat spack as unavailable; hand-build netCDF or use apptainer.

### 2.5 Scheduler / account (for the CIME machine config)

- `sinfo` (SLURM 26.05.5): `ctessum` partition = 5 nodes, 3-day limit, QOS default `scav_high`. Bigger partitions (`scavenger` ~60+ nodes, `secondary` ~130 nodes) are shared; `ctessum` is the group's own partition.
- `sacctmgr show assoc user=ctessum`: accounts `campusclu`, `ctessum-c…`, `ctessum-ic`, `ctessum` (Default QOS `scav_high`, MaxJobs 40).
- Current session runs with `SLURM_JOB_ACCOUNT=ctessum`, `SLURM_JOB_PARTITION=ctessum`, `SLURM_JOB_QOS=scav_high`.
- No `SBATCH_*` defaults in `~/.bashrc`; nothing to inherit. Use `--account=ctessum --partition=ctessum` explicitly.

### 2.6 Existing CESM/inputdata on the cluster

None (see Summary). No NCAR inputdata mirror; the group layout under `/projects/illinois/eng/cee` is one dir per PI/user group. Inputdata will have to be downloaded from NCAR (reachable, §1) or the LCRC E3SM mirror (reachable).

### 2.7 Sibling-project precedent (EqIce, EqWeFiC)

- **EqIce** (`../eqice/stage1/env.sh`): `module purge; module load gcc/13.3.0; module load openmpi/5.0.1-gcc-13.3.0` (+ `gsl/2.8`, `fftw/3.3.10`), used to build PETSc/PISM. This is the group's standard toolchain.
- **EqWeFiC** (`PLAN.md` §0.5, §11): same modules; "No netCDF module. The spack netCDF/OpenMPI install under `libs/spack` is stale (rhel7 paths, gcc 12.2.0)"; apptainer image builds work on this cluster — WRF was built inside `tools/apptainer/wrf-build.def` (Ubuntu 24.04, gfortran 13.3, `libnetcdf-dev libnetcdff-dev netcdf-bin libhdf5-dev libopenmpi-dev`, python3-netCDF4), producing a 334 MB `.sif`; `APPTAINER_CACHEDIR`/`APPTAINER_TMPDIR` must be under `/projects`/`data/`, never `/tmp`.

## 3. Recommended toolchain recipe for the CESM/CIME build

### 3.1 Native build (preferred)

```bash
module purge
module load gcc/13.3.0
module load openmpi/5.0.1-gcc-13.3.0
# python/3.13.2 optional; system python 3.9.25 also satisfies CIME (>=3.7)
```

- **Compiler:** gcc 13.3.0 + gfortran 13.3.0 — matches the EqIce/EqWeFiC toolchain exactly (openmpi 5.0.1 module is built against gcc 13.3.0; keep the ABI consistent by loading both).
- **MPI:** openmpi 5.0.1 (module). CIME's `MPILIB=mpi-serial` is available for smoke tests without MPI; use openmpi for real runs.
- **System tools that are already present (no module needed):** cmake 3.31.8 (`/usr/bin`), perl 5.32.1 (`/usr/bin`), python3 (`/usr/bin` or module), git.
- **netCDF:** **not available — must be built.** Hand-build in a `.gitignore`d directory under the worktree (or under `/projects/.../ctessum/ctessum/code/.cache`, per the XDG_CACHE_HOME convention), NOT in `/tmp` (RAM-backed) and NOT in the shared `ctessum/` root:
  1. `hdf5` (1.14.x, `--enable-fortran` not needed; parallel HDF5 optional — skip for now; netCDF-4 classic works without it)
  2. `netcdf-c` (4.9.x, `--prefix=<prefix>` with `CPPFLAGS/LDFLAGS` pointing at hdf5)
  3. `netcdf-fortran` (4.6.x, needs netcdf-c installed first; produces `netcdf.mod`/`nf-config`)
  Then point CIME at it with `NETCDF_PATH=<prefix>` (and `HDF5_PATH` if PIO is built with HDF5). netcdf-c can be configured `--disable-netcdf-4 --disable-dap` for a minimal hdf5-free build if desired, but CESM's PIO2 prefers netCDF-4; build hdf5 + netcdf-c + netcdf-fortran (all three are small, quick builds).
- **pnetcdf:** not available and not required for a first build — PIO2 works netcdf-only (add later if parallel I/O performance matters).
- **Account/partition for the CIME machine config** (`config_machines.xml`): `project="ctessum"`, queue/partition `ctessum` (3-day limit; use `secondary` for short 4 h test jobs). Default QOS is `scav_high`. CIME's `case.submit` will pass `--account=ctessum`; the `ctessum` partition is small (5 nodes) — for larger runs use `scavenger`.

### 3.2 Apptainer fallback

- `apptainer 1.5.4` is present at `/usr/bin/apptainer` and image builds are proven on this cluster (EqWeFiC built its 334 MB WRF toolchain image here in setuid mode; recipe pattern in `EqWeFiC/tools/apptainer/wrf-build.def`).
- A CESM/CIME build *inside* a container is viable but heavier: CIME drives configure/build/test via perl + python + the compilers, all of which must be inside the image (or bind-mounted), and MPI runs inside a container add the host/container MPI-matching question on top. **Recommend the native hand-built netCDF first; keep apptainer as the fallback** if the netCDF build fights back (e.g. an Ubuntu-based image with `libnetcdf-dev libnetcdff-dev` like EqWeFiC's, bind-mounting `/projects` and the worktree).
- Container builds must set `APPTAINER_CACHEDIR`/`APPTAINER_TMPDIR` under `/projects` (e.g. the worktree or `ctessum/code/.cache`), never `/tmp`.

## 4. Open questions / risks

1. **No netCDF (C or Fortran) anywhere** — the single biggest gap. Must hand-build (recommended) or containerize before CIME will configure. Watch for CIME's requirement of both `nc-config` (C) and `nf-config` (Fortran) plus `netcdf.mod`/`typesizes.mod`.
2. **No pnetcdf, no parallel HDF5, no hdf5 module** — fine for single-node runs; plan for `PIO_NETCDF_FORMAT=64bit_data` and revisit if multi-node I/O throughput matters.
3. **No NCO/NCL modules** — irrelevant to CIME builds (post-processing only); note if we later want `ncks` for diagnostics (could come from conda/miniforge module if needed).
4. **Spack is not usable** (not writable, nothing prebuilt, stale installs) — do not plan around it.
5. **Inputdata has to be downloaded from NCAR/LCRC** — no local mirror. Both endpoints are reachable now; CESM compsets pull multi-GB inputdata, so budget disk under `/projects` (Lustre, 5.8 PB free) and mind the 40-job MaxJobs / scav_high QOS limits. The SVN inputdata server answers 200 on `/trunk/` (the `README` path 404s — path is wrong, server is fine).
6. **`ctessum` partition is small (5 nodes, 3-day limit)** — fine for F-case testing; large production runs would need `scavenger`/`secondary` and a `--partition` override in `config_batch.xml`/`config_machines.xml`.
7. **Group directory is shared** — keep all build artifacts (netCDF prefix, caches, spack-like trees) in a `.gitignore`d subdir of the EqESM worktree or under `ctessum/code/.cache`, not in the shared `ctessum/` root and not in `/tmp` (RAM).
8. **Version pinning:** gcc 13.3.0 + openmpi 5.0.1 are the sibling-projects' pinned pair; intel 2025.0.4 + intel/mpi 2021.14 exist as an alternative if a Fortran-compiler-specific CESM bug appears (not recommended to start there). Note: `openmpi/5.0.8` exists under `/sw/apps` but has **no module** — stick with 5.0.1.
9. **hdf5 download URL:** `support.hdfgroup.org` 404s the classic release path (site reorganized) — use the GitHub release asset (`https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz`).
