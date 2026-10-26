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
