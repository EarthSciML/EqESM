# Stage-1 reference runs — configurations and pinned source

Reference implementation: CESM **cesm3.0-alphabranch**, checked out at
`stage1/cesm` (gitignored) with `bin/git-fleximod update`. All submodules are
pinned to the tags below (recorded from `git-fleximod status`,
2026-10-25). Version bumps are deliberate updates, recorded here.

## Pinned source

| Repo | Tag / commit |
|---|---|
| CESM (parent) | `cesm3.0-alphabranch` @ `e59f251` |
| ccs_config | `ccs_config_cesm1.0.92` |
| cime | `cime6.7.1` |
| share | `share1.1.22` |
| FMS (libraries/FMS) | `fi_260310` (nested fms `dev/ncar_0.1.2`) |
| cam | `cam6_4_212` |
| └ chem_proc | `chem_proc5_0_07` |
| └ atmos_phys | `atmos_phys0_31_000` (mmm-physics `20251208-CAM-SIMA`, pumas `pumas_cam-release_v1.41`, **rte-rrtmgp `v1.9.2`**, tempo `tempo_v3.1.2`) |
| └ cam-sima | `sima0_07_000` |
| └ fv3 | `fv3int_100726` (atmos_cubed_sphere `FV3-202204-public`) |
| └ geoschem | `14.5.3`; hemco `hemco-cesm2_3_0_hemco3_10_3` (HEMCO `3.10.3`) |
| └ clubb | `clubb_4ncar_20260109_ddf5110` |
| └ cosp2 | `v2.2.1`; rrtmgp-data `v1.8`; mpas `v8.4.0`; carma `carma4_10` |
| clm (CTSM) | `ctsm5.4.063` |
| cice | `cesm3_cice6_6_3_22` |
| mom | `mi_260909` (MARBL `marbl0.50.3`) |
| cism | `cismwrap_3_00_002` (source_cism `cism_main_3.00.002`) |
| cdeps | `cdeps1.0.112` |
| cmeps | `cmeps1.1.73` |
| mosart | `mosart1.1.13` |
| rtm / mizuRoute / ww3 | `rtm1_0_89` / `v3.1.1` / `main_0.1.1` (deferred components) |
| parallelio / mpi-serial | `pio2_6_8` / `MPIserial_2.5.4` |
| tuv-x | `v0.14.3` |
| FTorch | optional, not checked out |

Notes for the component inventory (PLAN.md §3): CAM 6.4 radiation is
**RTE-RRTMGP** (`rte-rrtmgp` v1.9.2, with `rrtmgp-data` v1.8) — the classic
RRTMG is gone; photolysis is **tuv-x** v0.14.3 (and Fast-JX lives in the
CAM-chem stack via the geoschem/hemco path). The CESM `tuv-x` submodule is
the photolysis library CAM uses.

## Reference runs (candidates; finalized per compset survey)

Each run records: exact `create_newcase` command, grid, timestep, coupling
interval, namelist diffs, precision, and regeneration command. Layout:
`stage1/runs/<name>.json` + `stage1/dumps/<component>/`.

- `F2000climo` (or the CESM3.0 F-class equivalent) at a coarse grid — the
  atmosphere+land reference. See `COMPSET_SURVEY.md`.
- All-data compset with `rof = mosart` — standalone MOSART reference.
- I-class compset (active CICE, data atm/ocean) — sea ice.
- C-class compset (active MOM) — ocean.
- B/G-class fully coupled — coupler and stage 3.
- Single-column harnesses: CAM SCAM (FSCAM), Icepack column driver,
  CLM single-point, MOM6 idealized tests, CISM standalone.
