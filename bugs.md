# Bugs and issues found in the CESM Fortran/C++ codes

Centralized list per PLAN.md §8. Each entry records: location, description,
impact, whether the `.esm` reproduces or corrects it, and status. (CESM is
mostly Fortran; "the C++ code" in AGENTS.md means the reference
implementation regardless of language.)

Format template:

## N. <one-line summary>

- **Location:** `<file>:<line>/<function>`
- **Description:** what the code does vs. what it should do.
- **Impact on results:** none / quantitative effect.
- **Reproduced/corrected in .esm:** which `.esm` component matches or
  corrects it (or n/a).
- **Status:** open / patched locally / reported upstream / fixed upstream.

## 1. R2000MOSART standalone compset is non-functional as shipped

- **Location:** `components/mosart/cime_config/config_compsets.xml`
  (`R2000MOSART` alias, lname
  `2000_SATM_DLND%LCPL_SICE_SOCN_MOSART_SGLC_SWAV`), found via the M1
  compset survey (2026-10-25, CESM `cesm3.0-alphabranch` @ `e59f251`,
  mosart 1.1.13).
- **Description:** the compset is commented "NOT yet functional" and
  references `DLND%LCPL`, which is not a valid dlnd mode on this branch
  (valid modes per `components/cdeps/dlnd/cime_config/config_component.xml`:
  `SCPL`, `GCPL`, `RCPL`). `RCPL` ("drive mosart with data from a previous
  run's coupler history") is tested in the dlnd testlist, so the intended
  fix is `DLND%LCPL` → `DLND%RCPL` (or use an F/I-class compset, which run
  MOSART coupled to CLM).
- **Impact on results:** cannot create a standalone-MOSART reference case
  as shipped; workaround: `usermods_dirs`/local compset override with
  `DLND%RCPL`, or couple via F-class/I-class.
- **Reproduced/corrected in .esm:** n/a (run-configuration bug, not a
  physics bug).
- **Status:** open (recorded; workaround documented in
  `stage1/runs/COMPSET_SURVEY.md` §3).

## 2. MOM `tx10deg`/`tx3deg` testmods reference inputdata files that are not published

- **Location:** `components/mom/cime_config/testdefs/testmods_dirs/mom/{tx10deg,tx3deg}/`
  (added `f0f3f07`, 2026-02-21; grid-independence fix `e5edf77`, 2026-07-09).
- **Description:** the testmods' `user_nl_mom`/`shell_commands`/`user_nl_cice`
  reference six files per grid under
  `$DIN_LOC_ROOT/ocn/mom/tx2_3v2/tx10deg/` (e.g. `ocean_hgrid_tx10deg_20260219.nc`,
  `ESMF_mesh_tx10deg_20260219.nc`, `cice_grid_tx10deg_20260219.nc`,
  `ocean_topog_tx10deg_20260219.nc`, `seawifs-clim-1997-2010-tx10deg_20260219.nc`,
  `state_restore_tx10deg_20260219.nc`; same for `tx3deg`), but **no `tx10deg`/`tx3deg`
  directory exists on any inputdata server** — checked 2026-10-26:
  `svn-ccsm-inputdata.cgd.ucar.edu` (directory listing of `ocn/mom/tx2_3v2/` shows no
  `tx10deg`), the OSDF endpoint (`osdf-data.gdex.ucar.edu`, cert-expired anyway), the
  LCRC E3SM mirror (`web.lcrc.anl.gov`), and `ftp.cgd.ucar.edu` (all 404). The
  `check_input_data --download` run for a `C_JRA`/`T62_t232` case with the `tx10deg`
  testmod ends with "Could not find all inputdata on any server" (5 missing refs, all
  the same `ESMF_mesh_tx10deg_20260219.nc` plus the MOM files). No grid-generation
  tool exists in the checkout (`MOM6/` has no coarse-grid generator; the testmod
  commit message says nothing about generating them).
- **Impact on results:** cannot create the ultra-coarse MOM smoke cases
  (`C_JRA`+`tx10deg`, `mom/tx3deg` variants) as shipped; workaround: run the
  C-class case at full `tx2_3v2` resolution (inputdata present) or a coarser
  atmosphere grid without the testmod.
- **Reproduced/corrected in .esm:** n/a (run-configuration/inputdata gap).
- **Status:** open (recorded; may be worth reporting upstream — the files
  were presumably generated once and never uploaded).
