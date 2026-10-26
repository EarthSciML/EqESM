# CESM 3.0 alpha — compset and grid survey for stage-1 reference runs

Survey date: 2026-10-25. Source: local checkout of
`ESCOMP/CESM` at `cesm3.0-alphabranch`, commit `e59f251` (2026-10-16,
"Merge pull request #435 from fischer-ncar/cesm3.0-alphabranch"), with
submodules cime 6.7.1, ccs_config_cesm 1.0.92, cam6_4_212, ctsm5.4.063,
cesm3_cice6_6_3_22, mom mi_260909, mosart 1.1.13, cmeps 1.1.73,
cdeps 1.0.112, mizuRoute v3.1.1, rtm 1_0_89, ww3 main_0.1.1.

All file paths below are relative to `stage1/cesm/` in the main checkout.
All recommendations assume a custom CIME machine on our cluster
(Slurm, gcc 13.3, openmpi 5.0.1) — see §5.

---

## 1. Executive summary (what changed vs CESM2)

1. **Compsets are split into per-component files.** There is no monolithic
   `config/cesm/allactive/config_compsets.xml` anymore. `cime_config/config_compsets.xml`
   holds only the fully coupled CESM3 (CAM7) compsets; each component
   (`components/{cam,clm,cice,mom,cmeps,mosart,mizuRoute,ww3,cism}/cime_config/config_compsets.xml`)
   defines its own. `create_newcase -list` merges them.
2. **CESM3's supported compsets are CAM7-based** (`CAM70%LT` low-top /
   `CAM70%MT` mid-top + `CLM60` + `CICE` + `MOM6` + `MOSART` + `DGLC%NOEVOLVE` + `WW3`).
   Classic CAM6 compsets (`F2000climo`, `FHIST`, …) still exist and are tested, but the
   "C" (CESM3) variants are `FHISTC_*`, `F1850C_*`, `B1850C_*`, `BHISTC_*`.
3. **No `I2000`/`I1850` CICE compsets.** The I-class (active CICE, data atm/ocn)
   is now the **`D` compsets** (`D`/`DIAF`/`DTEST`, defined by CICE). The old I-class
   *land* compsets (`I1850Clm60Sp` etc., active CLM + data atm/ice/ocn + **active MOSART**)
   live in the CLM compset file.
4. **No `T62_g17` grid** (removed; T62 is only paired with MOM6 tripole grids or g37).
   **No `ne4_ne4` grid** — the smallest CAM-SE aliases are `ne3_ne3_mg37` /
   `ne3pg3_ne3pg3_mt232` (~10°, "Spectral Elem 10-deg").
   PLAN.md §4.1's "ne4_ne4" assumption should be updated to these.
5. **Standalone MOSART compset is broken as shipped**: `R2000MOSART`
   (`2000_SATM_DLND%LCPL_SICE_SOCN_MOSART_SGLC_SWAV`) is commented
   "NOT yet functional" and uses `DLND%LCPL`, which is **not a valid DLND mode** in
   this branch (valid: `SCPL`, `GCPL`, `RCPL`). The `RCPL` mode ("river output
   coupling mode to drive mosart with data from a previous run") is tested
   (`components/cdeps/dlnd/.../testlist_dlnd.xml`), so `R2000MOSART` should work
   with the modifier fixed to `DLND%RCPL`. Practical MOSART runs today: any F-class
   or I-class compset (both run MOSART coupled to CLM).
6. **Fully coupled CESM3 runs use MOM6 on 2/3° tripole grids** (`tx2_3v2`/`tx2_3v3`),
   tested at `ne30pg3_t232`. Ultra-coarse MOM6 testmods exist:
   `mom/tx10deg` (36×32 ≈ 10°) and `mom/tx3deg` (108×96 ≈ 3°) — these are the
   smoke-test route for C/G/B-class runs.
7. **Compilers config moved from `config_compilers.xml` to CMake macros**
   (`ccs_config/machines/cmake_macros/gnu.cmake` etc.), and machine configs to
   per-machine directories (`ccs_config/machines/<MACH>/config_machines.xml`).
   A custom machine can be added via `$HOME/.cime/config_machines.xml` without
   touching the shared checkout (see §5).
8. **Inputdata** still uses classic wget/svn servers; no git-lfs mechanism on this
   branch. New highest-priority server is an OSDF endpoint
   `https://osdf-data.gdex.ucar.edu/ncar/gdex/d651077/cesmdata/inputdata/` (§6).
9. `ESMF_RUNTIME_PROFILE=ON` must be set as an environment variable on every
   NUOPC machine for timing output (derecho sets it "required on all systems").

---

## 2. Compset inventory (by definition file)

Compset naming: `TIME_ATM[%phys]_LND[%phys]_ICE[%phys]_OCN[%phys]_ROF[%phys]_GLC[%phys]_WAV[%phys][_ESP][_BGC]`.
Models available: ATM `CAM40/50/60/70` (+`SATM/DATM/XATM`), LND `CLM50/60` (+`SLND/DLND/XLND`),
ICE `CICE` (+`DICE/SICE/XICE`), OCN `MOM6` (+`DOCN/SOCN/XOCN`), ROF `MOSART/RTM/MIZUROUTE` (+`DROF/SROF/XROF`),
GLC `CISM2/DGLC/SGLC`, WAV `WW3/DWAV/SWAV/XWAV`.

### 2.1 `cime_config/config_compsets.xml` — fully coupled CESM3 (B-class, 18 compsets)

| Alias | Longname | Notes |
|---|---|---|
| `B1850C_MTt4s`, `B1850CM` | `1850_CAM70%MT%CT4S2_CLM60%BGC-CROP_CICE_MOM6_MOSART_DGLC%NOEVOLVE_WW3` | CESM3 1850 coupled, mid-top + simplified chem |
| `B1850C_MTso` | `1850_CAM70%MT_CLM60%BGC-CROP_CICE_MOM6_MOSART_DGLC%NOEVOLVE_WW3` | CESM3 1850 coupled, mid-top, no chem |
| `BHISTC_MTso`, `BHISTC_MTt4s`, `BHISTCM` | `HIST_CAM70%MT[%CT4S2]_CLM60%BGC-CROP_CICE_MOM6_MOSART_DGLC%NOEVOLVE_WW3` | CESM3 historical |
| `B1850C_LTso`, `BHISTC_LTso` | `..._CAM70%LT_CLM60%BGC-CROP_...` | low-top variants |
| `BHISTC_LT_DART` | `HIST_CAM70%LT_CLM60%BGC_CICE_MOM6_MOSART_DGLC%NOEVOLVE_SWAV` | DART (multi-instance) variant |
| `B1850E_MTt4s`, `B1850EM`, `BHISTE_MTt4s`, `BHISTEM` | same + `_BGC%BPRP` | emissions-driven (CMIP-style) |
| `ETEST`, `E1850TEST` | `2000/1850_CAM60_CLM50%SP_CICE_DOCN%SOM_MOSART_SGLC_SWAV_TEST` | classic CESM2-style coupled w/ SOM ocean, `_TEST` suffix |

No plain `B1850` alias (CAM6+SOM-ocean B1850 survives only as `E1850TEST`/`ETEST`).

### 2.2 `components/cam/cime_config/config_compsets.xml` — F-class & friends (128 compsets)

- **CAM6 science-supported**: `F2000climo` (`2000_CAM60_CLM50%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV`,
  science_support grid `f09_f09_mg17`), `FHIST` (same w/ `HIST_`).
- **CAM7 F-class**: `FHISTC_LTso`, `FHISTC_MTso` (`HIST_CAM70%LT|MT_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV`),
  `F1850C_LTso`, `F1850C_MTso`, `FLT1850_TESTINGONLY_v0c`, `FMT1850_TESTINGONLY_v0c`, `F2000dev`
  (`2000_CAM70_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV`).
- **CAM6 tested**: `F2010climo` (`2010_CAM60_CLM50%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV`),
  `F1850`, `FHIST_BGC`, `FHIST_C5/C4`, `FHIST_DARTC6`, `FHIST_BDRD`.
- **SCAM (single column, 17 IOP compsets + FCSCAM)**: `FSCAMARM95`, `FSCAMARM97`,
  `FSCAMATEX`, `FSCAMBOMEX`, `FSCAMCGILSS6/11/12`, `FSCAMDYCOMSRF01/02`, `FSCAMGATE3`,
  `FSCAMMPACE`, `FSCAMRICO`, `FSCAMSPARTICUS`, `FSCAMTOGA2`, `FSCAMTWP06`, `FSCAMCAMFRC`,
  `FCSCAM` (CAM-Chem SCAM). All `2000_CAM60%SCAM*_CLM50%SP_CICE%PRES_DOCN%DOM_SROF_SGLC_SWAV`.
  Entries auto-set `PTS_MODE=TRUE`, `PTS_LAT/PTS_LON`, `STOP_OPTION=nhours`, `STOP_N`
  (e.g. ARM97: 695 nhours), `START_TOD`, `RUN_STARTDATE`. **No SCAM tests in testlist_cam.xml**
  (SCAM is regression-tested on NCAR machines out-of-band).
- **Simplified**: `FHS94`, `FKESSLER`, `FTJ16`, `FADIAB`, `FGRAYRAD`, PORT compsets (`PC4`…`PC7`).
- **Aquaplanet**: `QPC4/5/6/7`, `QPX2000`, `QPMOZ`, `QSC4/5/6` (SOM aqua), etc.
- **WACCM**: `FWHIST` (supported), `FW2000climo`, `FWsc*`, `FX*` (WACCM-X); **cam-chem**: `FC2000climo`, `FCHIST`, …
- `FMOZ` is the only F-compset using `RTM` instead of `MOSART`.

### 2.3 `components/clm/cime_config/config_compsets.xml` — I-class (land, 112 compsets)

- `I2000Clm60Sp` = `2000_DATM%GSWP3v1_CLM60%SP_SICE_SOCN_MOSART_DGLC%NOEVOLVE_SWAV`
  (active CLM + **active MOSART**, data atm/ice/ocn) — the workhorse land test compset.
- Same family: `I1850Clm60Sp`, `I2000Clm50Sp`, `I1850Clm50Sp`, `IHistClm50/60*`,
  `ISSP245/370/585*`, BGC/BGC-CROP/FATES variants, `I*NoAnthroRtm` (RTM), `I*Miz*` (mizuRoute).
- Single-point: `I1PtClm50SpRsGs`, `I1PtClm60SpRsGs`, `I1PtClm60BgcRsGs`,
  `I1PtClm60FatesRsGs` (`2000_DATM%1PT_CLM*_SICE_SOCN_SROF_SGLC_SWAV`), plus urban
  `1x1_*` compsets with `STOP_OPTION=nsteps` (158/331 steps).

### 2.4 `components/cice/cime_config/config_compsets.xml` — D-class (sea ice)

- `D` = `2000_DATM%JRA-RYF9091_SLND_CICE_DOCN%SOM_DROF%JRA-RYF9091_SGLC_SWAV`
  (active CICE, data atm + slab ocean + data rof) — the I-class-for-ice.
- `DIAF` = same with `DATM%JRA-1p4-2018`; `DTEST` = `D` + `_TEST` suffix (used in CICE tests).

### 2.5 `components/mom/cime_config/config_compsets.xml` — C/G-class (ocean)

- **C-class (MOM6 alone)**: `C` (`2000_DATM%NYF_SLND_DICE%SSMI_MOM6_DROF%NYF_SGLC_SWAV`),
  `C_IAF`, `C_JRA`, `C_JRA_RYF`, `CW_JRA` (+WW3), `C1850MARBL_JRA`, `C1850ABIOMARBL_JRA`.
- **G-class (MOM6+CICE)**: `G`, `G_IAF`, `G_JRA`, `G_JRA_RYF`, `GW_JRA`,
  `G1850MARBL_JRA`, `G1850ABIOMARBL_JRA`, `GW1850MARBL_JRA`, `GW1850MARBL_CPLHIST`,
  `GW1850MARBL_CPLHIST_PHYS_CYCLE`, `GW1850ABIOMARBL_JRA`.
- No `C2000`/`C1850` aliases; use `C`/`C_JRA`/`C1850MARBL_JRA`.

### 2.6 `components/cmeps/cime_config/config_compsets.xml` — all-data / coupler testbeds

- `X` = `2000_XATM_XLND_XICE_XOCN_XROF_XGLC_XWAV` — all-stub compset; X-components are
  CMEPS test comps (`components/cmeps/med_test_comps/{xatm,xlnd,xice,xocn,xrof,xglc,xwav}`);
  tested on `f19_g17` (`SMS_D`/`ERS_Ln9` in `components/cmeps/cime_config/testdefs/testlist_drv.xml`).
- `A` = `2000_DATM%NYF_SLND_DICE%SSMI_DOCN%DOM_DROF%NYF_SGLC_SWAV` — all-data (real DATM/DICE/DOCN/DROF);
  `AA` (with DLND), `ADSOM`, `ADAQP3`, `AIAF`, `ADWAV`, `ADESP`, `S` (`SATM…_SESP`).

### 2.7 ROF compsets (standalone routing)

- `components/mosart/...`: `R2000MOSART` = `2000_SATM_DLND%LCPL_SICE_SOCN_MOSART_SGLC_SWAV`
  — **marked NOT functional; DLND%LCPL invalid** (see §3).
- `components/mizuRoute/...`: `R2000MizuRoute` (same shape, `MIZUROUTE`).
- `components/cism/...`: `T1850Gg`/`T1850Ga`/`T1850Gag` (CISM2 + DLND%SCPL), `JK1850D`.
- `components/ww3/...`: `WW3test`.

---

## 3. Which ROF: mosart vs drof

| Compset family | ROF model |
|---|---|
| F-class (CAM) | **MOSART** (all except `FMOZ`/some WACCM-4 which use RTM) |
| I-class (CLM) | **MOSART** (majority; RTM in `*NoAnthroRtm`, mizuRoute in `*Miz*`; `SROF` in 1Pt) |
| B-class coupled | **MOSART** |
| C/G-class (MOM6) | **DROF** (data river, `rx1` grid; JRA streams on `JRA025`) |
| D-class (CICE) | **DROF** |
| A-class (all-data) | **DROF** (`rx1`) |
| X-class | **XROF** (stub) |

**MOSART standalone route**: `R2000MOSART` is the intended route but is **not functional
as shipped** — it references `DLND%LCPL`; the DLND modes in this branch are `SCPL`, `GCPL`,
`RCPL` (per `components/cdeps/dlnd/cime_config/config_component.xml`), where RCPL is the
"river output coupling mode to drive mosart with data from a previous run" and is tested
(`2000_SATM_DLND%RCPL_SICE_SOCN_SROF_SGLC_SWAV` on `f09_f09_mg17`). Two options:
- patch the compset lname to `2000_SATM_DLND%RCPL_SICE_SOCN_MOSART_SGLC_SWAV`
  (a `usermods_dirs` entry or a local compset override), or
- run MOSART coupled: F-class (`F2000climo`) or I-class (`I2000Clm60Sp`) both exercise
  the full CLM→MOSART runoff path, which is what stage-1 instrumentation needs anyway.

CIME's `mosartOff` testmod (`MOSART_MODE=null`) exists for MOSART-off comparisons
(`components/mosart/cime_config/testdefs/testmods_dirs/mosart/mosartOff`).

---

## 4. Grids

`ccs_config/config_grids_nuopc.xml` includes `modelgrid_aliases_nuopc.xml` (aliases,
572 model_grid entries) and `component_grids_nuopc.xml` (domain/mesh definitions).
Default component grids (when an alias omits a component): rof `r05` for MOSART/RTM/XROF,
`rx1` for DROF, `r05mz` for mizuRoute; glc `null` for SGLC; wav `null` for SWAV.
Aliases carry `compset=`/`not_compset=` regexes; most F/I/D/A-class grids are
`not_compset="_MOM"`; MOM6-coupled grids use `tx2_*` ocnice.

### 4.1 Coarse aliases for F-class / I-class / D-class (no MOM6)

| Alias | atm | lnd | ocnice | mask | Notes |
|---|---|---|---|---|---|
| `ne3_ne3_mg37` | ne3np4 (~10°) | ne3np4 | ne3np4 | gx3v7 | smallest SE; aquaplanet-tested (`QPC5/6/7`) |
| `ne3pg3_ne3pg3_mt232` | ne3np4.pg3 | ne3np4.pg3 | ne3np4.pg3 | tx2_3v2 | CAM7 smoke grid (FHISTC tests) |
| `ne5_ne5_mg37` | ne5np4 | ne5np4 | ne5np4 | gx3v7 | |
| `f10_f10_mg37` | 10x15 (24×19) | 10x15 | 10x15 | gx3v7 | coarsest FV; F2000climo-tested |
| `f19_f19_mg17` | 1.9x2.5 (144×96) | 1.9x2.5 | 1.9x2.5 | gx1v7 | classic 2°; F2000climo/FW* tests |
| `f09_f09_mg17` | 0.9x1.25 (288×192) | 0.9x1.25 | 0.9x1.25 | gx1v7 | **F2000climo/FHIST science-support grid** |
| `C48_C48_mg17` | C48 (~3.75°) | C48 | C48 | gx1v7 | F2000climo-tested |
| `C96_C96_mg17` | C96 (~1.9°) | C96 | C96 | gx1v7 | F2000climo/FHIST/FHS94-tested |
| `ne30_ne30_mg17` / `ne30pg3_ne30pg3_mg17` | ne30np4(.pg3) (~1°) | same | same | gx1v7 | default-resolution SE tests |
| `f09_g17`, `f19_g17` | 0.9x1.25 / 1.9x2.5 | same | gx1v7 | gx1v7 | D-class, E1850TEST, X-compset tests |
| `g17_g17` | gx1v7 | gx1v7 | gx1v7 | — | DATM interpolation test (D-class) |
| `f10_f10_mg37`, `f19_f19_mg17`, `f09_f09_mg17`, `C96_C96_mt233`, `C96_C96_mg17` | | | | | I-class (CLM) tested grids |

**No `T62_g17` and no `ne4_ne4` in CESM3 alpha.** T62 survives only as
`T62_g37`, `T62_t232/t233/t025`, `T62_oQU*` (MPAS, unused by CESM).

### 4.2 MOM6-coupled aliases (C/G/B-class)

| Alias | atm | ocnice (MOM6/CICE) | Notes |
|---|---|---|---|
| `TL319_t201` | TL319 (~0.5°) | tx2_0v1 (**2°**, 180×128) | coarsest MOM6 tripole; GW_JRA-tested |
| `TL319_t232` | TL319 | tx2_3v2 (2/3°, 540×480) | D/C/G-test grid |
| `TL319_t233` | TL319 | tx2_3v3 (2/3°) | C/G-test grid |
| `T62_t232`, `T62_t233` | T62 | tx2_3v2/v3 | C/G-tested |
| `T62_t025` | T62 | tx0.25v1 (0.25°) | heavy |
| `ne30pg3_t232` / `ne30pg3_t233` | ne30np4.pg3 (1°) | tx2_3v2/v3 | **the tested CESM3 coupled grid** |
| `ne3pg3_t232` / `ne3pg3_t233` | ne3np4.pg3 (~10°) | tx2_3v2/v3 | ultra-coarse coupled (comment: "for low-res B compset testing, e.g. for DART") |
| `ne16pg3_t201` | ne16np4.pg3 (~1°) | tx2_0v1 (2°) | coarse-coupled alternative |
| `C96_t232`/`C96_t233` | C96 | tx2_3v2/v3 | note: `not_compset="_MOM"` — not for B-class |

**Ultra-coarse MOM6 testmods** (`components/mom/cime_config/testdefs/testmods_dirs/mom/`):
`tx10deg` (OCN/ICE 36×32, ~10°, NTASKS=6, `ROF2OCN_*_RMAPNAME=unset`) and
`tx3deg` (108×96, ~3°, NTASKS=32). Apply to any tx2_* grid to make a tiny ocean/ice.
Tested e.g. `SMS_D_Ld2 T62_t232 C mom/tx10deg`, `ERI_cG_C2_Ld8 T62_t232 G mom/tx10deg--mom/leapyear`.

### 4.3 Single-point / SCAM grids

1x1 point grids (`1x1_brazil`, `1x1_smallvilleIA`, `1x1_numaIA`, `1x1_mexicocityMEX`,
`1x1_vancouverCAN`, `1x1_urbanc_alpha`, `CLM_USRDAT`, `5x5_amazon`) exist for
DATM+CLM compsets. SCAM compsets auto-set `PTS_MODE`/`PTS_LAT`/`PTS_LON`; any
F-compatible grid can be used for case creation (the run collapses to one column;
domain file comes from the grid's `PTS_DOMAINFILE`).

---

## 5. Machines / compilers — what a custom machine entry needs

New ccs_config layout: `ccs_config/machines/config_machines.xml` (top level) holds only
`NODENAME_REGEX` (hostname→MACH mapping) and `default_run_suffix`; each machine lives in
its own directory `ccs_config/machines/<MACH>/config_machines.xml`. Compiler settings
are CMake macros in `ccs_config/machines/cmake_macros/*.cmake`
(`Macros.cmake` includes, in order: universal, `${COMPILER}.cmake`, `${OS}.cmake`,
`${MACH}.cmake`, `${COMPILER}_${OS}.cmake`, `${COMPILER}_${MACH}.cmake`).
No `config_compilers.xml` in ccs_config (only a stale template at
`userdefined_laptop_template/config_compilers.xml`).

**Adding our cluster**: per `cime/CIME/XML/machines.py` (v3.0), machines are discovered
by scanning subdirectories of the machines dir (minus `cmake_macros`,
`userdefined_laptop_template`); additionally `$HOME/.cime/config_machines.xml` and
`--extra-machines-dir` are appended (so we can define the machine without editing the
shared checkout). See `cime/doc/source/ccs/model-configuration/support-a-new-machine.rst`.

Required/typical fields of a `<machine MACH="...">` entry (per `config_machines.xsd`,
template `userdefined_laptop_template/config_machines.xml`, derecho entry):

| Field | Our cluster (suggested) |
|---|---|
| `DESC` | free text |
| `OS` | `LINUX` |
| `COMPILERS` | `gnu` |
| `MPILIBS` | `openmpi` (blank ⇒ single) |
| `CIME_OUTPUT_ROOT` | scratch dir for cases (e.g. `$ENV{SCRATCH}/cime`) |
| `DIN_LOC_ROOT` | **required** — inputdata dir (see §6) |
| `DIN_LOC_ROOT_CLMFORC` | optional (datm7 forcing mirror) |
| `DOUT_S_ROOT` | `$CIME_OUTPUT_ROOT/archive/$CASE` |
| `BASELINE_ROOT` | for test baselines |
| `CCSM_CPRNC` | optional (cprnc binary) |
| `GMAKE_J` | parallel compile jobs (compile_threads) |
| `BATCH_SYSTEM` | `slurm` (template already in `ccs_config/machines/config_batch.xml`; sbatch directives `--nodes/--ntasks-per-node/--exclusive`) |
| `SUPPORTED_BY` | us |
| `MAX_TASKS_PER_NODE`, `MAX_MPITASKS_PER_NODE` | cores/node (e.g. 128); `MAX_MEM_PER_NODE`, `MEM_PER_TASK` optional |
| `PROJECT_REQUIRED` | `TRUE` if sbatch needs an account; PROJECT is set via `./xmlchange PROJECT=...` |
| `<mpirun mpilib="default">` | e.g. `mpirun`/`mpiexec` with `<arg name="ntasks"> -np {{ total_tasks }}</arg>` (jinja; see stampede3: `ibrun -n {{ total_tasks }}`) |
| `<module_system type="module">` | init_path/cmd_path per shell + `<modules>` load commands (gcc/13.3, openmpi/5.0.1, netcdf-c/fortran, parallel-netcdf, cmake, python); or `type="none"` + `<environment_variables>` |
| `<environment_variables>` | `OMP_STACKSIZE`, **`ESMF_RUNTIME_PROFILE=ON` (required on all nuopc systems for timing output)**, netcdf hints |
| `NODENAME_REGEX` | add `<value MACH="ourcluster">...</value>` to the top-level file (or to the `$HOME/.cime` copy) |

Compiler side: `gnu.cmake` already sets `MPICC=mpicc`, `MPICXX=mpicxx`, `MPIFC=mpif90`,
`SCC=gcc`, `SFC=gfortran`, `FFLAGS=-fconvert=big-endian …`, and the gfortran ≥10 flags
(`-fallow-argument-mismatch -fallow-invalid-boz`); gcc 13.3 is fully covered.
NetCDF is found by the top-level CMake build via `find_package(NetCDF COMPONENTS C Fortran)`
(`cime/CMakeLists.txt`) — set `NETCDF_ROOT`/`CMAKE_PREFIX_PATH` or load a netcdf module;
the legacy Makefile path (`cime/CIME/Tools/Makefile`) uses `NETCDF_PATH`/`NETCDF_C_PATH`+
`NETCDF_FORTRAN_PATH` — per-machine cmake files can set `NETCDF_PATH`/`PNETCDF_PATH`
exactly like `CNL.cmake` (`set(NETCDF_PATH "$ENV{NETCDF_DIR}")`).
First build will also compile ESMF from CIME's bundled source — expect a long
initial build; derecho loads prebuilt `esmf-mpi/8.9.1` modules, our machine will build it.

---

## 6. Inputdata

`ccs_config/config_inputdata.xml` lists servers in precedence order (wget protocol):

1. `https://osdf-data.gdex.ucar.edu/ncar/gdex/d651077/cesmdata/inputdata/` (new OSDF CDN)
2. `https://ftp.cgd.ucar.edu/cesm/inputdata/` (+ `inputdata_checksum.dat` verification)
3. `https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/` (wget) and (svn)
4. `https://redoak.cs.toronto.edu/twitcher/ows/proxy/thredds/fileServer/datasets/CESM/inputdata/` (mirror)

No git-lfs inputdata mechanism on this branch. `DIN_LOC_ROOT` is a required per-machine
XML variable — convention: NCAR sets `$ENV{CESMDATAROOT}/inputdata`
(derecho: `/glade/campaign/cesm/cesmdata/inputdata`); we need our own mirror directory
(populate with `./check_input_data --download`), and optionally `DIN_LOC_ROOT_CLMFORC`
for datm7 forcing (`$DIN_LOC_ROOT/atm/datm7`).

---

## 7. Shortening runs (CIME / CESM3)

- Run length: `./xmlchange STOP_OPTION=<opt>,STOP_N=<n>` where valid options (CMEPS
  driver, `components/cmeps/cime_config/namelist_definition_drv.xml`) are
  `none, never, nseconds, nminutes, nhours, ndays, monthly, nmonths, nyears, date, end`;
  default is `ndays`/`5`. Defaults produce restart files at the end.
- `RESUBMIT=n` re-submits n times; `RESUBMIT_SETS_CONTINUE_RUN` (FALSE ⇒ each submit
  restarts over).
- Start date: `RUN_STARTDATE` (compset entries set 0001-01-01 for 1850/2000, 1850-01-01
  for HIST, etc.); `STOP_YMD`/`STOP_DATE` for fixed-date stops.
- Coupling: `ATM_NCPL` (default 48; 288 forced for WACCM-X), `RUNDT`; driver namelist
  via `user_nl_cpl` (`stop_option` group `CLOCK_attributes`).
- SCAM compsets auto-set `STOP_OPTION=nhours` + `STOP_N` + `START_TOD` (e.g. ARM95:
  418 nhours); 1Pt urban compsets auto-set `STOP_OPTION=nsteps`.
- Machine/batch: `JOB_WALLCLOCK_TIME`, `PROJECT`, `QUEUE` via `./xmlchange`.
- `--run-unsupported` is required at case creation for compset+grid combos without
  `science_support` (F2000climo on anything but `f09_f09_mg17`; all B/C/D/X combos).
- Test infrastructure: `create_test` with testmods (`allactive/defaultio`,
  `cam/outfrq9s` — 9-second output, `mom/tx10deg`, …).

---

## 8. Recommended smoke-run configurations (stage-1 reference runs)

"Coarse" = coarsest grid the compset is tested on; "default" = the CESM3 tested
config. All F/I/D/B-class runs need `--run-unsupported` off the science grid.

| # | Goal | Compset (alias) | Grid | Run length | Notes |
|---|---|---|---|---|---|
| 1 | Atmosphere+land (CAM6) | `F2000climo` (2000_CAM60_CLM50%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV) | `f19_f19_mg17` (or `C96_C96_mg17`, `f10_f10_mg37`; science grid `f09_f09_mg17`) | `STOP_OPTION=nhours,STOP_N=6` (or ndays,1) | MOSART active — doubles as a coupled MOSART run |
| 1′ | Atmosphere+land (CAM7) | `FHISTC_LTso` (HIST_CAM70%LT_CLM60%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV) | `ne3pg3_ne3pg3_mt232` or `ne30pg3_ne30pg3_mg17` | nhours,6 | CAM7 low-top; also `F1850C_LTso`, `F2000dev` |
| 2 | MOSART (coupled) | `I2000Clm60Sp` (2000_DATM%GSWP3v1_CLM60%SP_SICE_SOCN_MOSART_DGLC%NOEVOLVE_SWAV) | `f10_f10_mg37` or `f19_f19_mg17`/`f09_f09_mg17` | ndays,1 | the tested CLM+MOSART route |
| 2′ | MOSART standalone | `R2000MOSART` (2000_SATM_DLND%LCPL_SICE_SOCN_MOSART_SGLC_SWAV) — **fix `DLND%LCPL`→`DLND%RCPL`** | `f19_f19_mg17` (rof r05) | ndays,1 | broken as shipped; RCPL mode is tested (dlnd tests) |
| 3 | Sea ice (CICE) | `D` (2000_DATM%JRA-RYF9091_SLND_CICE_DOCN%SOM_DROF%JRA-RYF9091_SGLC_SWAV); also `DTEST` | `f19_g17`/`f09_g17` (tested: `TL319_t232`) | ndays,1 | DOCN%SOM slab ocean; CICE on gx1v7 |
| 4 | Ocean (MOM6) | `C_JRA` (2000_DATM%JRA-1p5-2023_SLND_DICE%SSMI_MOM6_DROF%JRA-1p5-2023_SGLC_SWAV) or `C` (DATM%NYF) | `T62_t232` (+ `mom/tx10deg` testmod ⇒ 36×32 ocean) or `TL319_t201` | ndays,1 | coarsest MOM6 tripole = tx2_0v1 (2°); **SUCCESS 2026-10-26** (`C_JRA_t62`, full tx2_3v2, 7:53 wall, MOM6 23 steps, CFL 0.075, mass err 6e-19) — see notes below |
| 5 | Fully coupled (stage 3) | `B1850C_LTso` (1850_CAM70%LT_CLM60%BGC-CROP_CICE_MOM6_MOSART_DGLC%NOEVOLVE_WW3); also `B1850C_MTso`, `BHISTC_LTso` | `ne30pg3_t232` (tested CESM3 grid); coarse: `ne3pg3_t232` (+ `mom/tx10deg`) | ndays,1 | full CESM3 coupled: CAM7+CLM6+CICE+MOM6+MOSART+WW3 |
| 6 | CAM single column | `FSCAMARM97` (or any `FSCAM*`, `FCSCAM`) | **`ne3_ne3_mg37` (required — SCAM only works with the SE dycore**; `cam/bld/configure` dies otherwise; IOP ic files keyed to `hgrid="ne3np4"`; CAM's own tests use `SCT_D_Ln7` + `cam/scm_prep` on `ne3_ne3_mg37`; `scam_mandatory` usermods: `NTASKS=1`, `REST_OPTION=never`, `CLM_FORCE_COLDSTART=on`) | auto (nhours, e.g. 695) | PTS_MODE/lat/lon auto-set by compset; STOP_N per IOP |
| — | Coupler testbed | `X` (2000_XATM_XLND_XICE_XOCN_XROF_XGLC_XWAV) | `f19_g17` | ndays,1 | all-stub; tested `SMS_D`/`ERS_Ln9` |
| — | Land single point | `I1PtClm60SpRsGs` | `1x1_brazil` / `CLM_USRDAT` | auto (ndays) | CLM60 single-point |

Grids shared between runs: `f19_f19_mg17` / `f09_f09_mg17` are usable by F-class
(CAM), I-class (CLM+MOSART), and (as `f19_g17`/`f09_g17`) D-class (CICE) —
one inputdata footprint (`f09`/`f19` + gx1v7 + r05) covers goals 1, 2, 3.
MOSART's `r05` rof grid is the same across F and I compsets.

### C-class run notes (2026-10-26, `C_JRA_t62` @ `T62_t232`, 1 day)

- **`mom/tx10deg` testmod unusable as shipped** — the 6 files it references
  (`ESMF_mesh_tx10deg_20260219.nc`, `ocean_hgrid_tx10deg_20260219.nc`, …)
  exist on no inputdata server and no generator is in the checkout →
  `bugs.md` #2. Ran the full `tx2_3v2` (2°) ocean instead: 72 PEs (2 nodes),
  7:53 wall, MOM6 23 steps, stable (CFL 0.075, mass err 6e-19).
- **Multi-node launch chain (all three walls documented in TOOLCHAIN.md #11):**
  `mpirun` blocked by `pam_slurm_adopt` (ssh daemons) → `srun -n N` blocked by
  PMI1/2 vs PMIx mismatch → working launcher `srun --mpi=pmix -n N` (machine
  config `config_machines.xml`, committed). Submit from the interactive
  session also requires `env -u SLURM_MEM_PER_NODE … ./case.submit`
  (SLURM_MEM_* mutual-exclusion).
- **`AUTO_MASKTABLE=True` FATAL on 36 OCN PEs** at tx2_3v2 ("Couldn't
  auto-eliminate any land blocks") — workaround `AUTO_MASKTABLE = False` in
  `user_nl_mom` (the `tx10deg` testmod does the same). Run crashed at init
  ~10 min in, batch hung until manually cancelled — check `run/ocn.log.*`
  for MOM FATALs, not just `cesm.log`.
- Input data: 619 JRA55 v1.5 files (~775 GB) + MOM ICs — the background
  download needs a retry loop (`check_input_data --download` ×3) to survive
  transient server failures; first unattempted attempt died mid-download.

### Notes for stage-1 instrumentation targets

- MOSART first (PLAN.md): capture from `I2000Clm60Sp`/`f10_f10_mg37` or `F2000climo`/`f19_f19_mg17`.
- CICE (cesm3_cice6_6_3_22) runs on tripole grids `tx2_3v2/v3` and `gx1v7`; D-class on `f19_g17` gives gx1v7.
- MOM6 (mi_260909) with `mom/tx10deg` gives a 36×32 ocean — cheap instrumented runs.
- CMEPS mediator is exercised in every run (NUOPC); the `X` compset isolates it.

---

## 9. Surprises / blockers (also worth a bugs.md entry)

1. **`R2000MOSART` standalone compset is non-functional as shipped** (`DLND%LCPL` is not
   a valid dlnd mode; modes are SCPL/GCPL/RCPL) — compset comment says "NOT yet functional".
2. **No `ne4_ne4`, no `T62_g17`** grid aliases (PLAN.md references ne4_ne4 — needs update
   to `ne3_ne3_mg37`/`ne3pg3_ne3pg3_mt232`).
2a. **SCAM requires the SE dycore** — FSCAM cases must use an SE grid
   (`ne3_ne3_mg37`), not FV (`f19_f19` etc.). Verified in `cam/bld/configure`
   (`if ($scam eq 'ON' and !($dyn_pkg eq 'se'))` → die).
3. **Compset files restructured** into per-component files; `-list` output differs from CESM2.
4. I-class renamed: sea-ice runs are `D` compsets; CLM I-class compsets moved to CLM's file.
5. B-class is MOM6-based only (`B1850C_*`); classic SOM-ocean B1850 survives only as `ETEST`/`E1850TEST`.
6. No SCAM tests in the CAM testlist (SCAM validated out-of-band at NCAR).
7. `mom/tx10deg`/`tx3deg` testmods use `ROF2OCN_*_RMAPNAME=unset` — fine for runs, not for answer comparisons against full-ocean runs.
8. Inputdata moved primary server to OSDF; `ftp.cgd.ucar.edu` still second; no git-lfs mechanism.
9. Machine/compiler config is CMake-based now (`cmake_macros/*.cmake`); `config_compilers.xml` is gone.
