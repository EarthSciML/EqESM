# CIME machine configuration for the Illinois Campus Cluster ("illinois")

This directory holds a draft CIME machine definition for the Illinois Campus
Cluster (ICC) so that CESM 3.0 can be built and run on the `ctessum`
partition (nodes ccc0232–ccc0236, SLURM, 40 cores / ~188 GB per node).

| File | Purpose |
|---|---|
| `config_machines.xml` | the `<machine MACH="illinois">` entry (+ top-level `NODENAME_REGEX`) |
| `config_batch.xml` | machine-specific SLURM settings: `sbatch`, `-p/--account/--qos` submit args, queues |

Facts were taken from `stage1/runs/TOOLCHAIN.md` (probe of 2026-10-25) and
verified again on node `ccc0233` (2026-10-26): module names `gcc/13.3.0` and
`openmpi/5.0.1-gcc-13.3.0` exist under `/sw/user/modulefiles`; Lmod lives at
`/usr/share/lmod/lmod` (`init/{sh,csh,perl,env_modules_python.py}` and
`libexec/lmod` all present); the hand-built netCDF/HDF5 stack is installed at
`/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/stage1/build/local`
(`netcdf.h`, `netcdf.mod`, `hdf5.h` verified present).

## Where the files go

CIME discovers machine definitions from three places (see
`cime/CIME/XML/machines.py` and `cime/CIME/XML/batch.py` in the CESM 3.0 /
cime 6.7.1 checkout):

1. **`--extra-machines-dir <dir>` flag to `create_newcase`/`create_test`
   (preferred).** Point it at this directory (the one containing this
   README). CIME appends `<dir>/config_machines.xml` to the merged machine
   list **and** `<dir>/config_batch.xml` to the merged batch config, which is
   exactly the pair we need. No modification of the shared CESM checkout.
   ```
   create_newcase --machine illinois --extra-machines-dir $REPO/stage1/cime/machine ...
   ```
2. **`$HOME/.cime/config_machines.xml` + `$HOME/.cime/config_batch.xml`.**
   Copy both files into `$HOME/.cime/` (creating it if needed). Same effect as
   (1), but applies to every CIME invocation for the user. Note that when
   `$HOME/.cime/config_batch.xml` exists, CIME does **not** also read the
   machine-specific `ccs_config/machines/<MACH>/config_batch.xml` — fine for
   us, we ship our own.
3. **A copy into the checkout** at `ccs_config/machines/illinois/`
   (`config_machines.xml` + `config_batch.xml`), the layout used by the
   bundled machines. This modifies the shared read-only checkout, so we do
   *not* do this in stage 1; it is the right move if the machine config is
   ever upstreamed (e.g. to ccs_config_cesm).

Also supported: `$HOME/.cime/illinois/config_machines.xml` (a per-machine
subdirectory of `~/.cime`, cime 3.0 layout). It only covers the machine
entry, not the batch overrides (batch.py reads the flat
`$HOME/.cime/config_batch.xml`), so we do not recommend it alone.

**Machine auto-detection:** with (1) or (2), hostnames matching the
top-level `NODENAME_REGEX` value `^ccc02.*` (the ctessum compute nodes)
resolve to MACH=`illinois`. If a login node does not match (login hostnames
are unverified — TBD), force it with `CIME_MACHINE=illinois` or
`create_newcase --machine illinois`.

**Operational note (2026-10-26):** `--extra-machines-dir` is only honored
by `create_newcase`/`create_test`. The subsequent case commands
(`case.setup`, `case.build`, `check_input_data`, …) re-discover machines
from the persistent paths and will fail with "No machine illinois found"
unless the config is also in `$HOME/.cime/` (option 2). **Use option 2
(`$HOME/.cime/config_machines.xml` + `config_batch.xml`) on every machine
where cases are created**; keep the canonical copies in this directory
(repo) and `cp` them to `$HOME/.cime/` after editing.

## Compiler macros (CMake)

Compiler settings are CMake macros, not `config_compilers.xml`:
`ccs_config/machines/cmake_macros/` (`Macros.cmake` includes universal,
`${COMPILER}.cmake`, `${OS}.cmake`, `${MACH}.cmake`,
`${COMPILER}_${OS}.cmake`, `${COMPILER}_${MACH}.cmake`).

We checked `ccs_config/machines/cmake_macros/gnu.cmake` (read 2026-10-26):
it sets `MPICC=mpicc`, `MPICXX=mpicxx`, `MPIFC=mpif90`, `SCC=gcc`,
`SCXX=g++`, `SFC=gfortran`, `FFLAGS="-fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none"`,
and — via the `CMAKE_Fortran_COMPILER_VERSION >= 10` branch —
`-fallow-argument-mismatch -fallow-invalid-boz` (the gfortran ≥ 10 flags).
**gcc/gfortran 13.3.0 is
fully covered; no machine-specific `.cmake` file is needed.** If we ever need
tweaks, the file to add would be `cmake_macros/gnu_illinois.cmake`
(auto-included for COMPILER=gnu on MACH=illinois), e.g. to point
`NETCDF_PATH`/`PNETCDF_PATH` at the hand-built stack exactly like
`CNL.cmake` does (`set(NETCDF_PATH "$ENV{NETCDF_DIR}")`); for now we set the
paths via `<environment_variables>` instead.

## How netCDF is found

- **CMake build (used by CESM 3.0):** `cime/CMakeLists.txt` calls
  `find_package(NetCDF COMPONENTS C Fortran)`. The machine entry exports
  `NETCDF_ROOT`, `HDF5_ROOT`, and `CMAKE_PREFIX_PATH` pointing at the
  hand-built prefix `<repo>/stage1/build/local` (the same prefix
  `stage1/runs/env.sh` exports as `NETCDF_PATH`).
- **Legacy Makefile path (`cime/CIME/Tools/Makefile`):** uses `NETCDF_PATH`
  (and `NETCDF_C_PATH`/`NETCDF_FORTRAN_PATH`); we export `NETCDF_PATH` too,
  so both build paths work.
- The cluster has **no netCDF/HDF5 modules** (TOOLCHAIN.md §2.3), so the
  stack is the hand-built one (hdf5 1.14.x + netcdf-c + netcdf-fortran,
  installed by `stage1/build/build-netcdf.sh`). The prefix path is
  hardcoded in the XML — if the stack ever moves, update the five `<env>`
  entries in `<environment_variables>`.

## Batch system (SLURM)

`<BATCH_SYSTEM>slurm</BATCH_SYSTEM>` uses the generic slurm entry in
`ccs_config/machines/config_batch.xml` (directives `--job-name`, `--nodes`,
`--ntasks-per-node`, `--output`, `--exclusive`) plus our machine-specific
`config_batch.xml`, which supplies the missing `batch_submit` (`sbatch`) and
submit args. With the machine defaults (`PROJECT=ctessum`,
`JOB_QUEUE=ctessum`) the effective submission is:

```
sbatch --time $JOB_WALLCLOCK_TIME -p ctessum --account ctessum --qos scav_high
```

(`--partition=ctessum --account=ctessum --qos=scav_high` in slurm's
`=`-notation — exactly the directives our interactive probes use).
`PROJECT_REQUIRED=TRUE` plus the machine-level `<PROJECT>ctessum</PROJECT>`
default means `create_newcase` always has an account; override with
`--project` or `./xmlchange PROJECT=...` if the group ever adds accounts.

`<mpirun mpilib="default">` runs `mpirun -np {{ total_tasks }}` (Open MPI
5.0.1). We deliberately omit `-npernode` (removed in Open MPI 5.0): rank
placement follows the SLURM allocation (`--ntasks-per-node` +
`--exclusive`).

## First build

CMEPS requires ESMF: the CMake build does `find_package(ESMF REQUIRED)`
(`components/cmeps/CMakeLists.txt:30`) and the cmeps buildnml reads
`$ESMFMKFILE` (`components/cmeps/cime_config/buildnml:735`). NCAR machines
load a prebuilt `esmf-mpi/8.9.1` module (derecho); the cluster has **no ESMF
module, and ESMF is NOT bundled in the CESM checkout** (the earlier draft
claim that CIME builds it from bundled source was wrong — verified
2026-10-26). We build **ESMF 8.9.1** (derecho's pin) from source with
`stage1/build/build-esmf.sh` (gcc 13.3 + openmpi 5.0.1, netcdf split mode,
static libs) into `stage1/build/esmf`, and the machine entry sets
`ESMFMKFILE=$REPO/stage1/build/esmf/lib/esmf.mk` and `ESMF_DIR`. Expect a
long initial build. Nothing else in this config is ESMF-specific beyond
`ESMF_RUNTIME_PROFILE=ON` (required on all NUOPC systems for CMEPS/ESMF
timing output; CESM 3.0 is NUOPC-only) and
`ESMF_RUNTIME_PROFILE_OUTPUT=SUMMARY` (derecho convention).

## Inputdata

`DIN_LOC_ROOT` points at `<repo>/stage1/inputdata` (directory created; still
empty — populate with `./check_input_data --download` when creating cases).
`data.cesm.ucar.edu` does not resolve on this cluster (TOOLCHAIN.md §1), but
the OSDF endpoint, `ftp.cgd.ucar.edu`, and `svn-ccsm-inputdata` are
reachable; no svn client exists, so keep CIME inputdata on wget/https
protocol (default). `DIN_LOC_ROOT_CLMFORC` defaults to
`$DIN_LOC_ROOT/atm/datm7` per convention.

## TBD items / caveats

1. **`/scratch.local` is node-local disk** (verified: `/dev/mapper/
   local-scratch`, not a shared FS). `CIME_OUTPUT_ROOT=/scratch.local/
   ctessum/eqesm-cime` therefore only works when the SLURM job runs on the
   same node that created/built the case (single-node jobs on the same node).
   *Verify whether ICC keeps jobs on the submitting node for this path
   before relying on it; if not, move `CIME_OUTPUT_ROOT` (and
   `DOUT_S_ROOT`/`BASELINE_ROOT`) to shared `/projects` storage.*
2. **Login-node hostnames unverified.** `^ccc02.*` covers compute nodes
   ccc0232–ccc0236 (FQDN `ccc02NN.campuscluster.illinois.edu` verified).
   If login nodes do not match, use `CIME_MACHINE=illinois`/`--machine`.
3. **`mpirun` inside an sbatch script** (not via `srun`): Open MPI 5.0.1
   discovers the allocation from SLURM environment variables; if it reports
   "unable to find the number of slots", we will need to add `--map-by` /
   `--bind-to` args or run through `srun`. Untested.
4. **`--qos scav_high` is hardcoded** (it is the ctessum account's default
   QOS). If other partitions need a different QOS, parameterize it. Untested.
5. **PIO is netcdf-only** (no pnetcdf/parallel-HDF5 on the cluster): expect
   `PIO_NETCDF_FORMAT=64bit_data` for I/O; fine for single-node runs.
6. **Walltime default** = queue `walltimemax` (72 h on `ctessum`, 24 h
   `scavenger`, 4 h `secondary`); shorten with `./xmlchange
   JOB_WALLCLOCK_TIME=HH:MM:SS`.
7. **XSD note:** the top-level `NODENAME_REGEX` in `config_machines.xml` is
   not expressible in `config_machines.xsd` — but the upstream
   `ccs_config/machines/config_machines.xml` has the identical structure, and
   CIME never schema-validates these files at runtime (no schema entries in
   `CIME/data/config/*/config_files.xml`), so `xmllint --schema` fails on
   upstream files too. The `<machine>` element itself is XSD-clean (verified
   with xmllint against `cime/CIME/data/config/xml_schemas/
   config_machines.xsd`).
