# M9: mediator ocean-export chain -> MOM6 surface forcing — gate result: feasible (3 of 4 links); CICE links not traceable on current dumps

> **Manager note (2026-11-10):** only this README (gate verdicts, chain
> laws, dump spec) and `stage1/analysis/extract_m9_ocn_flux_rows.py` are
> merged. The worker's `m9_ocn_flux.esm` / `mom6_cap_imports.esm` were
> emitted by an assembler script, restated the 86 AoFlux couplings, and
> fitted `areacor` from the dump; they are being rewritten by hand in the
> follow-up job `m9-ocn-flux2` (branch `wt/m9-ocn-flux`).

`stage2/m9_ocn_flux.esm` **was written** (with `stage2/mom6_cap_imports.esm`),
coupling the traced part of the CESM ocean-flux chain end to end: the
mediator atmosphere-ocean bulk fluxes (aoflux) -> the ocean-export merge ->
the MOM6 nuopc-cap import transforms -> the MOM6 surface forcing. The CICE
<-> CMEPS links of the same chain are **not traceable** in the current
stage-1 dumps (the CICE lane ran a different case); this README records the
per-link gate analysis and the exact dump fields/points a stage-1 job would
need to unblock the CICE part.

## Gate question and answer

**Q:** For each CMEPS <-> MOM6 and CICE <-> CMEPS coupling link, do the two
sides of the link (the producing component's exports and the consuming
component's imports) overlap on any (cell, step) in the stage-1 dumps, in
the same run configuration?

**A:** Three of the five links are traceable and are coupled in
`stage2/m9_ocn_flux.esm`; the two CICE links are not.

| link | verdict | evidence (stage-1 dumps, C_JRA_t62 DATM+DICE+MOM6) |
|---|---|---|
| aoflux -> merge (Faox_* exports vs merge inputs) | **FEASIBLE, same step** | `cmeps_aoflux.csv` and `cmeps_merge.csv` share (gid, step): Faox_taux/tauy/sen/evap/lwup == aoflux taux/tauy/sen/lat/lwup **bit-exact on all 288 rows** (12 traced cells x 24 steps) |
| ocnalb -> merge (So_* albedos vs merge So_* inputs) | **FEASIBLE, one-step lag** | merge step s So_avsdr/avsdf/anidr/anidf == ocnalb step s-1 avsdr/avsdf/anidr/anidf **bit-exact on all 276 overlapping rows** (the mediator's albedo lag: the merge consumes the albedos computed in the previous coupling step; `cmeps_merge.esm` documents it) |
| merge -> mom6 imports (Foxx_*/Faxa_*/Fioi_* exports vs mom_import) | **FEASIBLE, one-step ocean lag + areacor + rotation** | mom6 cpl_step k == merge/aoflux step k+1 on all 144 mom6 rows: p = pslv, u10_sqr = duu10n, ice_fraction = ifrac bit-exact (the unscaled imports); the 16 areacor-scaled imports equal the merge exports times the per-cell med2mod_areacor at rel <= 1e-12, with the wind-stress rotation the identity on all traced cells (south of the tripole fold) |
| CICE -> merge (Fioi_*_ice exports vs merge Fioi_*_ice inputs) | **NOT FEASIBLE** | the CICE dumps come from a **different case** (JRA-RYF9091/CICE/DOCN%SOM gx1v7; `stage1/runs/CICE_LANE.md`) than the merge dumps (C_JRA with DICE%SSMI tx2_3v2): no shared (cell, step), different ice model (the merge's Fioi_*_ice values are the DICE data-model dummies), different grids |
| merge -> CICE imports (Fioi_* exports vs cice_import) | **NOT FEASIBLE** | same mismatch; the CICE lane's import dumps are on the CICE gx1v7 grid of its own case, not on the C_JRA merge's grid |

The two C_JRA runs that produced the feasible side are the same
configuration with bit-identical trajectories (see next section), so the
overlap rows above are genuine same-run links.

## The lane runs (same configuration, different traces)

Both traces ran compset
`2000_DATM%JRA-1p5-2023_SLND_DICE%SSMI_MOM6_DROF%JRA-1p5-2023_SGLC_SWAV_SESP`
at `tx2_3v2` (atm/lnd/ice) with the MOM6 `gx1v7`-derived `tx10deg`-family
ocean grid, 2000-01-01, 1 day, 24 coupling steps of 3600 s (`ATM_NCPL=24`,
`OCN_NCPL=$ATM_NCPL`; `stage1/runs/MOM6_LANE.md`, `CMEPS_LANE.md`):

| | CMEPS lane (`C_JRA_t62_inst`) | MOM6 lane (`C_JRA_t62_mom6_inst`) |
|---|---|---|
| traces | `cmeps/cmeps_aoflux.csv`, `cmeps/cmeps_merge.csv`, `cmeps/cmeps_ocnalb.csv` — 12 traced ocean points x 24 steps, PET 0 | `mom6/mom6_sfc_forcing.csv` — 12 traced ocean cells x 12 coupling steps, PET 0 (the first 12 coupling steps, tods 0, 7200..43200) |
| neutrality | restart vs the reference `C_JRA_t62`: bit-identical (CMEPS lane checks) | restart vs the reference `C_JRA_t62`: bit-identical (MOM6 lane checks; the mom6 dumps were regenerated in the wave-3 run `C_JRA_t62_mom6_dyn_inst` run 12282697, the same trajectory) |

The mom6 dump step semantics are the CESM ocean lag: within each coupling
interval the ocean runs **first** and consumes the mediator products
prepared at the **end** of the interval, so mom6 cpl_step k lines up with
merge/aoflux step k+1 (and ocnalb step k). This offset is verified
bit-exact on all 144 mom6 rows by
`stage1/analysis/extract_m9_ocn_flux_rows.py`.

The 12 traced points are the only ocean cells of PET 0's block of the 6x6
decomposition (the 77 southernmost grid rows are Antarctic land mask): the
cells are shared by all four dumps by (lat, lon), mom6 (i,j) -> cmeps gid
(7,82)/(8,82) -> 41583/41584, (7..10,83) -> 42123..42126, (6..11,84) ->
42662..42667.

## The chain as implemented (what the .esm models)

Cited against the stage-1 forks (commit pins in `stage1/runs/*_instrument.json`):

1. **aoflux** (`cesm/flux_atmocn/flux_atmocn_Large.F90`, `flux_atmOcn_large`):
   the Large & Yeager bulk fluxes from the atmosphere state — implemented
   and tested standalone in `stage2/cmeps_aoflux.esm` (`AoFlux`).
2. **ocnalb** (`mediator/med_phases_ocnalb_mod.F90`): the ocean albedos from
   the solar geometry — `stage2/cmeps_ocnalb.esm` (`OcnAlb`).
3. **merge** (`mediator/med_phases_prep_ocn_accum_mod.F90` +
   `med_phases_prep_ocn_mod.F90`): the ocean-export merge — the ice/ocean
   fraction-weighted sums of the aoflux, DICE and river exports plus the
   `med_phases_prep_ocn_custom` shortwave computation —
   `stage2/cmeps_merge.esm` (`OcnMerge`).
4. **mom6 cap imports** (`components/mom/MOM6/config_src/drivers/nuopc_cap/
   mom_cap_methods.F90` `mom_import`): every flux import is scaled by
   `med2mod_areacor = mesh_areas/model_areas` (`mom_cap.F90:1345`,
   `state_getimport` with `areacor`, F90:114-204, 213-314, 471-491) and the
   wind stress is rotated from true zonal/meridional to the local
   coordinates (F90:165-175) — `stage2/mom6_cap_imports.esm` (`MomCap`,
   standalone, 4 tests). The p/u10_sqr/ice_fraction imports are **not**
   scaled (Sa_pslv F90:114, Si_ifrac F90:506-508, So_duu10n F90:514-516).
5. **mom6 surface forcing** (`mom_surface_forcing_nuopc.F90`
   `convert_IOB_to_fluxes/forces`): the import-to-fluxes/forces mapping —
   `stage2/mom6_surface_forcing.esm` (`SurfaceForcing`, standalone,
   merged earlier).

`stage2/m9_ocn_flux.esm` mounts all five by reference and wires them with
`variable_map` couplings at their shared step; the one-step ocean lag and
the one-step albedo lag are realized by the tests (the ocnalb inputs of
dump step k drive `OcnAlb` while `AoFlux`/`OcnMerge` run on step k+1 inputs
and `SurfaceForcing` asserts the step-k mom6 outputs). Because the CLI
cannot mount a file's internal couplings by reference, the 86
AoFlux-internal couplings of `cmeps_aoflux.esm` are restated verbatim in
the subassembly (same limitation and workaround as `m9_datm_clm.esm`).

## What the dumps contain (the reference side)

- `cmeps_aoflux.csv`: 288 rows = 12 cells x 24 steps; the aoflux inputs
  (zbot/ubot/vbot/thbot/qbot/rbot/tbot/pslv/us/vs/ts), the full iteration
  state, and the outputs (taux/tauy/sen/lat/lwup/evap, tref/qref,
  duu10n/u10res).
- `cmeps_merge.csv`: 288 rows; the merge inputs (ifrac/ofrac/ifrad/ofrad,
  the Faxa_* atm fields, the Fioi_*_ice DICE fields, Forr_*), the custom
  shortwave intermediates (fswabsv/fswabsi, the *_scaled fractions) and the
  ocean exports (Foxx_swnet[_vdr/vdf/idr/idf], Foxx_lwnet/sen/taux/tauy/
  evap, Faxa_rain/snow, Foxx_rofl/rofi, Fioi_melth/meltw/salt). The Fioi_*
  values are the **DICE data-model dummies** (see the CICE section).
- `cmeps_ocnalb.csv`: 288 rows; nextsw_cday/rlat/rlon inputs and the
  anidr/avsdr/anidf/avsdf outputs (also the ifrac/ofrac/ifrad/ofrad the
  phase sees).
- `mom6_sfc_forcing.csv`: 144 rows = 12 cells x 12 steps; the imports as
  MOM6 consumes them (after areacor/rotation) and the resulting fluxes and
  forces (including the AGRID taux/tauy, tau_mag/ustar/omega_w2x and
  net_mass_src), plus the neighbor rows needed for the AGRID stencils.

Two values are **not** dump columns and are derived by the committed
extractor `stage1/analysis/extract_m9_ocn_flux_rows.py`:

- the per-cell `med2mod_areacor`, fitted from the mom6 step-1
  sw_vis_dir / merge step-2 Foxx_swnet_vdr ratio (constant over the 12
  steps to ~2e-16; 0.99998068916425..0.99998089684951 by cell);
- the neighbor stress of an **untraced ocean** cell: the dump rows carry
  the AGRID masks themselves (taux == u_flux bit-exact iff the east
  neighbor is land; tauy == v_flux iff the north neighbor is land). At
  (11,84) the east neighbor (12,84) is land but the north neighbor (11,85)
  — which lies beyond the 12 traced cells — is ocean, so its stress is
  back-derived exactly from the dump average, v_flux_north = 2*tauy - v_flux
  (the only test override that is not a dump value; marked as such in the
  cpl1_i11_j84 test description).

## Exactly which dump fields/points a stage-1 job would need for the CICE links

The CICE lane's dumps (`stage1/dumps/cice/*.csv`, from case
`JRA-RYF9091/CICE/DOCN%SOM` gx1v7; `stage1/runs/CICE_LANE.md`) and the CMEPS
dumps (case `C_JRA_t62` with DICE%SSMI) do **not** overlap: different
compsets (real CICE vs the DICE data model), different grids (gx1v7 vs
tx2_3v2), different traced points, no shared (cell, step). To unblock the
two CICE <-> CMEPS links, a stage-1 job in the **C_JRA_t62 configuration
with the CICE data model replaced by real CICE** (compset
`2000_DATM%JRA-1p5-2023_SLND_CICE%..._MOM6_DROF%JRA-1p5-2023_SGLC_SWAV_SESP`,
same atm/ocean/rof as the traced runs so the mom6 side stays comparable)
must produce **two new dumps**:

**A. CICE export side (cice fork — owned by the CICE lane):** at each
coupling step s = 1..24, on the traced ocean points (the CICE grid's
ice-cell columns that map onto the 12 traced ocn cells; the ice grid is
gx1v7, the ocn grid tx2_3v2, so the mapping file `cice2ocn`/`ocn2cice`
weights define the contributing stencil): `Fioi_taux_ice`, `Fioi_tauy_ice`,
`Fioi_melth_ice`, `Fioi_meltw_ice`, `Fioi_salt_ice`, `Fioi_swpen_ice` (the
six fields the merge consumes in `med_phases_prep_ocn_accum`), i.e. the
`x2o_Fioi_*` values as produced by `cice_export` (`cice/cimecice/ice_import_export.F90`).

**B. CICE import side (cmeps side, on the C_JRA run with real CICE):** the
same six fields **as received by CICE** (the `fldptr` values in
`ice_import_export.F90` `state_getimport` after `med2mod_areacor`) and the
CICE imports that come from the merge/ocean side (`Foxx_swnet_vdr/.../idf`,
`Foxx_lwnet/sen/taux/tauy/evap`, `Faxa_rain/snow`, `Forr_rofl/rofi`) on the
same (cell, step) filter — this pins the ice side of the coupled test. The
ocean side of the same run would then also exercise the real-ice merge
(`Fioi_melth/meltw/salt` nonzero where the DICE run has dummies) and the
`Fioi_swpen_ice` branch, which the current merge dump only records as 0.

## Coverage caveats for the current tests

- The CICE links are entirely unpinned (`Fioi_*_ice` are DICE dummies; the
  tests cover the merge fraction-weighted sums with those values).
- The enthalpy imports (Foxx_hrain/hsnow/hrofl/hrofi/hevap/hcond) take the
  same areacor scaling in the cap (F90:245-314) but their mediator exports
  are untraced (no Foxx_h* column in `cmeps_merge.csv`), so they enter the
  tests as free `SurfaceForcing` parameters overridden from the mom6 dump.
- The rotation is pinned only as the identity (cos_rot = 1, sin_rot = 0);
  the traced cells all lie south of the tripole fold, so the nonzero
  rotation branch is unexercised.
- The mom6 dump has no `med2mod_areacor` column (the cap applies it
  internally); the tests use the fitted per-cell value above, and the
  ratio's step-constancy to ~2e-16 is verified by the extractor on all 144
  rows.
- `Fioi_salt` is 0 on every traced row, so the salt-import areacor scaling
  (F90:471-491) is pinned only as 0 = areacor x 0.
- No CESM defect was found during this gate analysis; the missing CICE
  overlap and the missing areacor/enthalpy dump columns are dump-coverage
  gaps, so nothing is appended to `bugs.md`.

## Step 2 status

Implemented and merged on `wt/m9-ocn-flux`:
`stage2/mom6_cap_imports.esm` (MomCap, 4 tests) and `stage2/m9_ocn_flux.esm`
(M9 subassembly, 8 tests; 272 assertions total, all pass at rel 1e-12),
with `stage1/analysis/extract_m9_ocn_flux_rows.py` reproducing every pinned
number from the dumps. The CICE part stays gated on the stage-1 job
specified above.
