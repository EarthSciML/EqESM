# CISM lane (ice sheets, milestone M4)

Scope: CISM as built in CESM (`components/cism` = cismwrap_3_00_002, source
`source_cism` = cism_main_3.00.002) is to be re-implemented by **reusing the
EqIce components** (`../eqice/stage2/`, PISM-derived, 3443+ tests). This lane
delivers (1) the CISM→EqIce process mapping, (2) the cheapest runnable
evolving-CISM CESM configuration, (3) the first instrumented slice of CISM
boundary dumps (velocity solve, temperature step, thickness/SMB evolution,
basal friction law, cap coupling boundary), (4) validation vs the Fortran.

Sibling lanes' patterns: CICE (`stage1/runs/cice_instrument.json`, fork
`stage1/EqESM-cesm-cice`), MOSART, CTSM.

## Ownership / collision guards

- **Ours:** `stage1/runs/CISM_LANE.md`, `stage1/runs/cism_instrument.json`,
  fork `stage1/EqESM-cesm-cism` (branch `instrument/cism-dumps`; cism
  wrapper submodule + source_cism submodule on the same branch name),
  `stage1/dumps/cism/`, `stage1/analysis/*cism*`, stage-2 CISM work (later).
- **EqIce (read-only):** `/projects/.../code/eqice` — do not modify; import
  by reference in stage 2.
- **Coordinate:** `PROGRESS.md` (manager-maintained), the shared
  `stage1/inputdata`, the `ctessum` partition, `bugs.md`.

## 1. The reference configuration (what a CESM evolving-Greenland run does)

Chosen run (cheapest that evolves CISM; see §3): **T1850Gg**
(`1850_SATM_DLND%SCPL_SICE_SOCN_SROF_CISM2%GRIS-EVOLVE_SWAV`) on
**f09_g17_gris20** (land 0.9x1.25 — required by the DLND%SCPL default stream;
glc = gris20, the coarsest Greenland grid, 76×141 cells).

Options in force for GRIS-EVOLVE (from `cime_config/namelist_definition_cism.xml`
and `source_cism/libglide/glide_types.F90` defaults; the CESM namelist exposes
only a subset — the rest are the source defaults):

| Option | Value | Meaning |
|---|---|---|
| `dycore` | 2 | **Glissade** (HO FEM). Glide (SIA) and glam are *not supported* in CESM (`valid_values=2`) |
| `which_ho_approx` | 4 | HO_APPROX_DIVA — depth-integrated viscosity approximation (2-D membrane solve + 3-D recovery) |
| `temperature` | 1 | TEMP_PROGNOSTIC — **cold-ice temperature** (enthalpy option 3 "under construction", unused) |
| `temp_init` | 3 | advective-diffusive initialization |
| `flow_law` | 2 | FLWA_PATERSON_BUDD (Paterson & Budd 1982, with the −10 °C activation-energy switch) |
| `evolution` | 3 | EVOL_INC_REMAP — incremental-remapping transport |
| `which_ho_babc` | 3 | HO_BABC_PSEUDO_PLASTIC_OLD (the CESM-Greenland default basal law) |
| `which_ho_bwat` | 2 | HO_BWAT_LOCAL_TILL (local till water: dbwat/dt = bmlt·ρi/ρw − c_drainage, clamped [0, bwat_till_max]) |
| `which_ho_effecpress` | 4 | HO_EFFECPRESS_BWAT_BVP — Bueler & van Pelt (2015) Eq. 23 (bwat → N) |
| `which_ho_ground` | 2 | subgrid grounding line, "deluxe" (f_ground_cell weighting) |
| `which_ho_ground_bmlt` | 1 | bmlt_ground scaled by floating fraction |
| `which_ho_thermal_timestep` | 0 | thermal solve **before** transport (vertical diffusion + strain heating; advection via transport) |
| `which_ho_ice_age` | 1 | compute ice age (increment by dt each step; advected by remapping) |
| `bmlt_float` | 0 | none for Greenland (no ice shelves in this config) |
| `ocean_data_domain` | 0 | ocean data internal (unused with bmlt_float=0) |
| `isostasy` | 0 | **off** in the CESM default |
| `gthf` | 1 | prescribed 2-D geothermal flux |
| `calving_domain` / `marine_margin` / `calving_init` | 1 / 1 / 1 | calve wherever the criterion is met; `remove_icebergs=.true.`, `remove_ice_caps=.true.` |
| `smb_input` | 0 | acab from the coupler (m/yr ice) |
| `enable_*_anomaly` | all .false. | no SMB/artm/snow corrections |
| `glc_smb` | .true. | SMB comes from the land model via the mediator (never the internal PDD scheme) |

Time stepping: NCPL_BASE_PERIOD=year, GLC_NCPL=1 → CISM is called once per
year; the gris20 icesheet config sets `dt = 0.5 yr`, so each call runs 2
half-year `glissade_tstep`s.

## 2. CISM → EqIce mapping

EqIce components (PISM-derived, `../eqice/stage2/`): `ice_rheology`,
`sia`, `ssa`, `stress_balance`, `ice_energy` (+`ice_energy_3d`), `ice_age`
(+`ice_age_3d`), `basal_resistance`, `hydrology`, `bed_deformation`,
`calving`, `geometry`, `surface_forcing`, `ocean_coupling` (+ their
`*_templates` libraries).

| # | CISM process | CISM code (file:line) | EqIce component | Verdict |
|---|---|---|---|---|
| 1 | Glissade time step (order: age → isostasy → gthf → climate forcing → thermal → bmlt_float → transport → calving → cleanup → thermal(2nd, off) → diagnostic/velocity) | `libglissade/glissade.F90:1172` (`glissade_tstep`) | `stress_balance` + `geometry` + `ice_energy` subassembly (EqIce stage3 `eqice.esm` has the same ordering pattern) | **same structure** — same operator split and order as PISM's `IceModel::step`; the reference runs this order |
| 2 | Velocity solve — HO/DIVA membrane stress balance (Picard FEM, preconditioned CG; DIVA = 2-D depth-integrated solve + 3-D recovery) | `libglissade/glissade_velo.F90:42` (`glissade_velo_driver`), `libglissade/glissade_velo_higher.F90:662` (`glissade_velo_higher_solve`) | `ssa.esm` (membrane balance), `stress_balance.esm` (assembly) | **different numerics, same PDE family**: CISM solves the SSA/BP/L1L2/DIVA membrane equation with a FEM staggered-grid assembly; EqIce's SSA is the SSAFD finite-difference form on the flow-line verification case. The DIVA approximation has **no PISM counterpart** (PISM's SSA is closest). The *SSA* subcase (which_ho_approx=1) is the same PDE as EqIce's |
| 3 | SIA velocity (option −1; not the CESM default) | `libglissade/glissade_velo_sia.F90:83` | `sia.esm` | **same law, different discretization** (PISM: FV diffusivity chain; CISM: local SIA quadrature). Not exercised by the reference run — instrument later if wanted |
| 4 | Flow factor (Paterson-Budd): A(T*) = 1.733e3·exp(−139e3/(8.314·T*)) [T* ≥ 263.15 K], 3.613e-13·exp(−60e3/(8.314·T*)) [T* < 263.15 K], T* = clamp(T − T_pmp, −50, 0) + 273.15, ×E enhancement | `libglissade/glissade_therm.F90:2408` (`glissade_flow_factor`), constants `libglimmer/glimmer_physcon.F90:96-103` | `ice_rheology.esm` (Paterson-Budd softness) | **same law, same constants** as PISM's Paterson-Budd; small variant in the T* construction (CISM clamps T−T_pmp to [−50,0] before adding 273.15) and in where the enhancement factor multiplies |
| 5 | Temperature step (vertical diffusion + strain heating on sigma levels, implicit tridiagonal; surface Dirichlet artm; pmp limit; bmlt_ground from the bed energy balance) | `libglissade/glissade_therm.F90:530` (`glissade_therm_driver`), `libglissade/glissade.F90:1419` (`glissade_thermal_solve`) | `ice_energy.esm` / `ice_energy_3d.esm` (PISM enthalpy column/3-D) | **same physics, different formulation**: CISM solves for temperature T (deg C) with a cold-ice pmp limit and waterfrac; PISM solves the enthalpy equation E with a cold/temperate switch. The vertical diffusion/advection operators are both implicit column solves (CISM tridiagonal vs PISM BOMBPROOF). Energy balance at the bed → bmlt_ground is the same mechanism |
| 6 | Thickness + tracer evolution (incremental remapping advection, then SMB + basal mass balance applied) | `libglissade/glissade.F90:1643` (`glissade_thickness_tracer_solve`), `libglissade/glissade_transport.F90:276` (`glissade_transport_driver`), `libglissade/glissade_mass_balance.F90:819` (`glissade_apply_smb`) | `geometry.esm` (mass continuity dH/dt = −div Q) | **same law (continuity), different advection numerics**: CISM incremental remapping (monotone, 2nd order) vs PISM flux-limited FV upwind; the SMB source term dH/dt += acab is the same |
| 7 | Surface mass balance input (acab, m/yr ice) and surface air temperature (artm) from the coupler; lapse-rate/downscaling only if anomalies on (all off in reference) | `libglissade/glissade_mass_balance.F90:160` (`glissade_prepare_climate_forcing`); cap import `drivers/cpl/nuopc/glc_import_export.F90:295`; conversion qsmb (kg/m²/s) → acab (m/yr ice, ×rhow/rhoi, zeroing rules) `libglad/glad_main.F90:802-823` | `surface_forcing.esm` (PISM PIK atmosphere, Martin 2011) | **different source, same application**: in CESM the SMB is CLM-computed (via DLND%SCPL + mediator renormalization); EqIce's PIK parameterization is a different SMB law. The *application* (acab → dH/dt) is identical. The CESM cap/mediator path (elevation classes, renormalize_smb) is CMEPS territory |
| 8 | Basal friction — pseudo-plastic drag: tau_c = tan(phi(bed))·N (phi linear in bed elevation between phimin@bedmin and phimax@bedmax), beta = tau_c/(u0^q·speed^(1−q)), speed = sqrt(u²+v²+0.01²), tau_b = beta·u | `libglissade/glissade_basal_traction.F90:69` (`glissade_calcbeta`), PSEUDO_PLASTIC_OLD branch at :250 | `basal_resistance.esm` / `basal_templates.esm` (pseudo-plastic drag) | **same drag law** — EqIce's `beta = tauc·(eps²+u²+v²)^((q−1)/2)·u0^(−q)` with eps=0.01 m/yr is exactly CISM's formula. **Different yield stress**: CISM tau_c = tan(phi(bed))·N with the bed-elevation phi ramp; EqIce tauc = c0 + tan(phi)·N_till with phi from data/inversion |
| 9 | Effective pressure: N = N0·(δ·P/N0)^s·10^((e0/Cc)(1−s)), s = bwat/bwat_till_max, capped at overburden, zero for floating; min with the ocean-penetration term | `libglissade/glissade_basal_traction.F90:668` (`glissade_calc_effecpress`), BWAT_BVP branch at :944 | `basal_templates.esm` (Mohr-Coulomb effective pressure) | **same law** — EqIce's `N_interior = N0·(δP/N0)^s·10^((e0/Cc)(1−s))` (Tulaczyk 2000 / Bueler & van Pelt 2015) is the identical functional form with s = W/W_till_max. Same parameters family (N0, δ, e0, Cc); CISM additionally takes min(·, P), zeroes floating cells, and applies the ocean_p min. **CISM quirk**: PSEUDO_PLASTIC_OLD computes phi from `topg` (cell centers) while beta lives on the staggered grid — the code comments "There is a minor bug… should use stagtopg"; kept for backward compatibility (recorded, reproduce in .esm) |
| 10 | Basal water — local till: dbwat/dt = bmlt·ρi/ρw − c_drainage/scyr, clamp [0, bwat_till_max] | `libglissade/glissade_basal_water.F90:146` (`glissade_calcbwat`) | `hydrology.esm` (PISM Routing/Shreve q = −K∇ψ) | **different law (as implemented in EqIce)**: CISM default is the *local till* water balance; EqIce's hydrology component models the Shreve routing model (PISM also has local till — not in EqIce stage2). CISM's flux-routing option (which_ho_bwat=3) is closer to EqIce's routing but is off |
| 11 | Bed deformation / isostasy (elastic lithosphere + Kelvin relaxing asthenosphere; **off** in the CESM default, isostasy=0) | `libglissade/glissade_isostasy.F90:193` (`glissade_isostasy_solve`), `libglissade/glissade_isostasy_elastic.F90`, `glissade_isostasy_kelvin.F90` | `bed_deformation.esm` (PISM PointwiseIsostasy) | **same physics family, off in the reference run** — PISM's iso model (elastic + relaxing asthenosphere) matches CISM's isostasy options; not exercised by T1850Gg |
| 12 | Calving (threshold-based: remove thin/floating marine-margin ice; remove icebergs and ice caps) | `libglissade/glissade_calving.F90:618` (`glissade_calving_solve`) | `calving.esm` (PISM EigenCalving / von Mises / Hayhurst rate laws) | **different law**: CISM's calving is a thickness/margin criterion (with iceberg/cap removal); PISM's are strain-rate-based rate laws. At gris20 the run may or may not trigger calving |
| 13 | Basal melt — grounded (from the bed energy balance in the thermal solve) and floating (bmlt_float=0 in the reference) | `libglissade/glissade_therm.F90` (bmlt_ground), `libglissade/glissade_bmlt_float.F90:843` (`glissade_bmlt_float_solve`) | `ice_energy.esm` (bed energy balance), `ocean_coupling.esm` (PISM PIK sub-shelf melt) | grounded: **same mechanism** (bed energy balance); floating: **different law, off in the reference** (PISM PIK Beckmann-Goosse vs CISM ISMIP6/thermal-forcing options) |
| 14 | Ice age (increment by dt; advected as a tracer by remapping) | `libglissade/glissade.F90:1203-1213`; transport tracers | `ice_age.esm` / `ice_age_3d.esm` | **same equation (dτ/dt = 1 + advection), simpler numerics**: CISM increments in place and advects with remapping; PISM solves the age PDE with its own advection stencil |
| 15 | Cap coupling boundary (imports qsmb/tsfc with the tsfc −273.15 K clamp quirk; exports topo, hflx, rofi, rofl, ice_covered, ice_sheet_grid_mask, frac fields) | `drivers/cpl/nuopc/glc_import_export.F90:295/347`, `source_glc/glc_RunMod.F90:59/150` (glad_gcm), `libglad/glad_main.F90:581` | none (CMEPS/mediator lane) | **missing in EqIce** — CESM-specific coupling conversions (qsmb kg/m²/s → acab m/yr ice with ×rhow/rhoi, zeroing on ocean/thin cells; tsfc → artm with the where(tsfc<−250) = 0 quirk). These become new CISM coupling tests |
| 16 | Ocean coupling fields (salinity/tocn imports; ocean_data_domain=0) | cap + `glad_main.F90` | `ocean_coupling.esm` | **off in the reference run** (Greenland, bmlt_float=0); different scheme family when on |

### Mapping summary

- **Same law** (reproduce the CISM formulas in EqIce components / new CISM tests): flow factor / Paterson-Budd rheology (#4), basal pseudo-plastic drag + Bueler–van Pelt effective pressure (#8, #9 — the strongest "same law" pair), SIA (#3, off in reference), SSA membrane balance (#2 subcase), thickness continuity + SMB application (#6, #7-application), grounded basal melt (#13), ice age (#14), step ordering (#1).
- **Different law / formulation** (new CISM components or parameterized variants): DIVA/HO velocity solve numerics (#2), temperature-vs-enthalpy formulation (#5), incremental-remapping vs FV advection (#6), SMB source (CLM vs PIK, #7), yield-stress construction (#8), local-till vs routing hydrology (#10), threshold vs rate calving (#12), floating melt (off, #13), ocean coupling (off, #16).
- **Missing in EqIce / CESM-only**: cap coupling conversions (#15), the mediator SMB path (elevation classes, renormalization — CMEPS lane), DIVA approximation (#2).

Counts: **same = 9** (rheology, SIA, SSA, effecpress, drag, continuity/SMB-application, grounded bmlt, age, step order) · **different = 7** (velocity numerics, energy formulation, advection numerics, SMB source, yield stress, hydrology, calving) · **missing/off = 4** (cap conversions, mediator SMB, DIVA, isostasy/ocean-melt — off in the reference).

## 3. Chosen run configuration

**T1850Gg / f09_g17_gris20** — the cheapest CESM configuration that evolves
CISM: SATM/SICE/SOCN/SROF/SWAV are stubs, DLND%SCPL replays saved land snow
coupling from a previous CESM2 CMIP6-piControl run (the "saved land SMB"
route), and **CISM2%GRIS-EVOLVE is the only active model**. The f09 land grid
is *required* (the DLND%SCPL default stream is only defined for l%0.9x1.25);
gris20 is the coarsest Greenland grid (76×141, dt=0.5 yr — the CICE lane's
f19_g17 case already pulled the 4 km mesh but the 20 km files are tiny).

- Compset: `1850_SATM_DLND%SCPL_SICE_SOCN_SROF_CISM2%GRIS-EVOLVE_SWAV` (+ SESP, added by CIME).
- Grid: `a%null_l%0.9x1.25_oi%null_r%null_g%gris20_w%null_z%null_m%gx1v7`.
- Run length: STOP_OPTION=nyears, STOP_N=2 (2 coupling years = 4 CISM 0.5-yr steps).
- PEs: 1 task for every component incl. CISM (TOTALPES=1; CISM PET 0 owns the whole 20 km domain — the dump is PET-0-only by construction).
- Input data (all downloaded 2026-11-01 to `stage1/inputdata/`, svn-ccsm-inputdata server):
  - `glc/cism/Greenland/gland20.input_c150415.nc` (348 KB — CISM 20 km mesh+state),
  - `share/meshes/gland_20km_c150511_ESMFmesh.nc` (658 KB — mediator GLC mesh),
  - `lnd/dlnd7/CPLHIST_SNO/b.e21.B1850.f09_g17.CMIP6-piControl.001_c210324/b.e21.B1850.f09_g17.CMIP6-piControl.001.cpl.hl2x1yr_glc.{1971..2000}-01-01.nc` (30×9.97 MB — DLND SCPL stream; model year 1 ↔ data year 1971 by DLND_CPLHIST_YR_ALIGN=1).
- Alternative if gris20 misbehaves: `f09_g17_gris4` (4 km, files already in `stage1/inputdata/glc/cism/Greenland/`).

## 4. First instrumented slice (this lane)

Dump module `cism_dump_mod.F90` in `source_cism/libglissade/` (auto-picked-up
by the CMake `file(GLOB libglissade/*.F90)`), namelist-gated (`dump_nml`
group in `cism_in`, added to the CIME definition + buildnml's
`_CISM_IN_GROUPS`, read in `source_glc/glc_InitMod.F90`), PET-0 only
(`main_task` of the CISM communicator — with TOTALPES=1 it owns the whole
grid), first N traced cells × first M steps, OMP-critical writes, one CSV per
boundary with `# description = …` header (units + F90 file:line).

Boundaries (each gives the named pair a test):

1. **`cism_velo`** (D, velocity solve): entry/exit of `glissade_velo_higher_solve` — (thck, topg, usrf, stagthck, stagtopg, flwa, effecpress, eus, masks) → (uvel, vvel 3-D at the traced vertex). Tests #2 + #4 together.
2. **`cism_thermo`** (I-at-substep, temperature step): entry/exit of `glissade_therm_driver` — (temp 3-D, thck, artm, bheatflx, bfricflx, dissip) → (temp 3-D, flwa, bpmp, bmlt_ground, waterfrac). Tests #5.
3. **`cism_transport`** (I, thickness/tracer evolution): entry/exit of `glissade_thickness_tracer_solve` — (thck, tracers temp/age, uvel, vvel, acab, bmlt) → (thck, tracers). Tests #6 + #7-application.
4. **`cism_basal_traction`** (D, pointwise law): inside `glissade_calcbeta` PSEUDO_PLASTIC_OLD — per traced vertex, per Picard iteration (bounded): (u, v, speed, bed, phi, tanphi, N_stag, tau_c, beta) and inside `glissade_calc_effecpress` BWAT_BVP: (bwat, bwat_till_max, N0, delta, e0, Cc, overburden, relative_bwat, N). Tests #8 + #9.
5. **`cism_cap_imports` / `cism_cap_exports`** (I, coupling boundary): in `source_glc/glc_RunMod.F90` around `glad_gcm` — (tsfc, qsmb per traced cell per coupling year) and (topo, hflx, rofi, rofl, ice_covered, ice_sheet_grid_mask). Tests #15.
6. **`cism_flowlaw`** (D, pointwise): inside `glissade_flow_factor` — (temp, pmptemp, E) → flwa per traced column layer. Direct pin of #4.

Dumps land in `stage1/dumps/cism/`; run config + validation in
`stage1/runs/cism_instrument.json`; solution-neutrality = instrumented vs
dump-off run of the same case, bit-identical restart/history modulo timestamps.

## 4b. Run record (2026-11-02)

Case `T1850Gg_cism_inst` (T1850Gg/f09_g17_gris20, TOTALPES=1, STOP_N=2)
built from the fork `stage1/EqESM-cesm-cism` (branch `instrument/cism-dumps`).
Run protocol: two startup runs of the same case/exe — run 1 dump-off
(`run_ref_dumpoff` snapshot), run 2 dump-on (`dump_cism=.true.` in
user_nl_cism) — compared bit-identical via `stage1/analysis/compare_cism.py`.

- Old-exe pair (first attempt): 11978727 (dump-off) / 11978818 (dump-on),
  exit 0:0 both. Dumps revealed two dump-module defects and one design
  flaw, fixed in `cism_dump_mod.F90`:
  1. the 400-char `cols` header buffers truncated the trailing column
     names of the state traces (thermo lost `btemp_float,bmlt_ground`;
     velo/transport lost their last layer columns) — buffers now 2000;
  2. the `a5` phase format right-justified `exit` with a leading space —
     now unpadded `a`;
  3. traced cells were the first 32 cells of the *static* active-grid
     mask (`ice_sheet_grid_mask`), which are ice-free coastal cells with
     all-zero physics (artm=bpmp=bwat=acab=0, no velocities) — the
     enumeration now uses the initial ice coverage (`ice_covered` =
     thck > 0 from `glad_get_initial_outputs`) at `cism_dump_init`,
     so the traces carry real interior-Greenland ice.
- Intermediate pair (after fixes 1–3): 11978998 (dump-off) /
  11979075 (dump-on), exit 0:0 both.  Old-exe dump-off outputs preserved
  as `run_ref_dumpoff_oldexe` (in the case dir) and
  `stage1/dumps/cism_oldexe_20261102/`.
- Final exe pair: 11979644 (dump-off) / 11980950 (dump-on, final),
  exit 0:0 both; 11979768 and 11980002 were dump-on reruns before
  fixes 5 and 6 respectively.  Fix 4 (an *instrumentation* bug, not a
  CESM bug — no `bugs.md` entry): every model-side value in the first
  two pairs was shifted by (−2,−2).  The model's local arrays are
  halo-indexed (nhalo=2: `model%general%ewn/nsn` = 80×145 local vs
  76×141 owned) while the cap arrays (`qsmb`, `ice_covered`) are
  owned/global-indexed; the trace enumerated cells from the cap's
  global `ice_covered` but read the model arrays at those indices, so
  traced cell (i,j) dumped global cell (i−2,j−2) (e.g. (25,11) showed
  thck=0, topg=282.47 = global (23,9)).  `cism_dump_init` now takes
  the model local dims and computes off = (local−owned)/2 = 2; model
  reads go to (i+off_i, j+off_j) and `dump_cell_ix` translates
  call-site local indices back to the traced global cell.  Also
  removed the trailing comma `dump_reals` emitted after each row,
  which then exposed Fix 5: the cap trace wrote `ice_covered` right
  after the `dump_reals` block with no separator (the old trailing
  comma had doubled as the delimiter), so every cap row had 14 fields
  against the 15-column header — now written with a leading comma and
  verified with a standalone driver before the rerun.
  Fix 6 (instrumentation, no `bugs.md` entry): the transport trace's
  temperature columns were off by one — it wrote `temp(1..upn),temp(0)`
  into the header slots `temp_000..temp_upn`, shifting every level and
  landing the surface value `temp(0)` (level 0) in the bed slot
  `temp_upn`; now writes `temp(0..upn)` in vertical order (level 0 =
  surface, level upn = bed), verified with the standalone driver
  before the rerun.
  Intermediate snapshots preserved as `run_ref_dumpoff_exeB` and the
  (superseded) `dumps/cism_oldexe_20261102/`.
- Effective-pressure trace: legitimately empty for this run — the
  gland20 initial state has bwat = 0 everywhere and the 2-year cold-base
  run melts nothing (bmlt_ground = 0), so the BWAT_BVP branch never
  fires. The law is pinned from the source formula + glide_types
  defaults; the trace machinery is exercised whenever bwat > 0.
- Validation: `.work/cism_validate.py` (law recomputation, structural
  row counts, step-to-step chaining, cross-file identity checks),
  rewritten after the offset analysis to chain by the true solve order
  (thermo → transport → diagnostic velocity: within-step velo-entry ==
  transport-exit, thermo == transport-entry; cross-step transport
  entry(s+1) == exit(s) exactly and velo velocities to rtol 1e-12 for
  the m/yr↔m/s scale round-trips) and to apply the SMB zeroing rules
  (acab = 0 where thck ≤ 0 and acab < 0, or thck ≤ 0 and topg < eus,
  using the same step's transport-entry thck);
  `stage1/runs/cism_instrument.json` records the full protocol.  PASS on
  the final dumps; instrumented vs dump-off run bit-identical modulo
  timestamps (`stage1/analysis/compare_cism.py`).

## 4c. Thermo-extension rerun (2026-11-03, stage-2 lane row 5)

The stage-2 temperature lane needs four more pieces from the thermo trace:
`dttem` (the 0.5-yr solve step), `dissip_001..010` (interior strain heating
divided by rhoi·shci, deg/s), and the grid coordinates `sigma_001..011` /
`stagsigma_001..010`.  A dissip = 0 reproduction of the backward-Euler step
fails the operator identity at 3–5% of the step-1 temperature change, and the
pin needs the exit rows plus the exact geometry.  `cism_dump_thermo` gained the
four column groups (both call sites, `glissade_therm.F90:727,1331`); committed
at all three fork levels — source_cism `b33660b6`, cism wrapper `142a4cb`,
fork superproject `a95fe50` (branch `instrument/cism-dumps`, not pushed).
Case rebuilt and rerun (dump-on only; the extension is solution-neutral by
construction — it only adds trace columns): job 12053984 COMPLETED 22:47.
New `cism_thermo.csv` (256 rows × 78 cols) validated in the stage-2 lane:
sigma/stagsigma match the glide_calc_sigma SIGMA_COMPUTE_GLIDE predictions
exactly, dttem = 1.5768e7 s on every row, dissip entry == exit, and the
operator identity (T_exit − T_entry)/dttem = L·T_exit + dissip holds to
≤6e-9 relative (pinned in `stage2/cism_thermo.esm`).  Old dumps backed up to
`/scratch.local/ctessum/cism_dumps_backup/`.  Note: the stage-1 validator
`.work/cism_validate.py` was removed in the shared-repo cleanup; the
extended-dump validation now lives in the stage-2 test suite
(`stage2/cism_thermo.esm`, `stage2/README-cism.md`).

## 4d. DIVA neighbour-stencil rerun (2026-11-05, stage-2 lane row 6)

The stage-2 DIVA velocity lane (`stage2/cism_velo.esm`) needs vertex and
neighbour quantities the original `cism_velo.csv` trace does not carry: the
4-cell stencil of the NE-corner vertex (usrf/thck/ice_mask/land_mask of the
SW/SE/NW/NE cells) for the hybrid-margin surface gradient (`dusrf_dx/dy`) and
the 4-cell mean staggered thickness, the quadrature-point effective
viscosities (`efvs_q1..q4` per layer) that enter the Goldberg (2011) `omega_k`
integrand, the `compute_integrals_diva` outputs (`omega_k` at the interfaces,
the four cell `omega`s and their vertex mean `stag_omega`, `stag_omega_k`),
and the solve-exit DIVA state (`beta_eff_x/y`, `beta_internal`,
`stag_theta_slope_x/y`, `stag_diva_slope_factor_x/y`, `dusrf_dx/dy`).
`cism_dump_mod.F90` gained the `cism_velo_nbr.csv` trace (122 columns,
described in `stage1/runs/cism_instrument.json`), committed at the fork
superproject `db624ec` with the cism wrapper bump `76d170d` (components/cism =
`cismwrap_3_00_002-4-g2a8316d`, branch `instrument/cism-dumps`, not pushed).
Case rerun (dump-on only; the extension is solution-neutral by construction —
it only adds trace columns): job 12118489 COMPLETED exit 0:0 (2026-11-05,
08:41). 256 rows = 32 cells × 4 steps × 2 phases (128 exit rows).
Validation: `stage1/analysis/extract_cism_velo_test_rows.py` recomputes every
DIVA pointwise law on all 128 exit rows — efvs quad-point arithmetic mean
(exact), `omega_k` cumulative integral (≤4e-16 rel), `omega` trapezoid (exact),
`stag_omega` 4-cell mean (≤3e-16), `theta_x/y` from the plain vertex gradient
of `z_mean = usrf − thck/2` (≤7e-14), the F slope factors (exact),
`beta_eff` (exact), `dusrf_dx/dy` hybrid gradient with the max_slope clamp
(≤7e-13), `taudx/dy` (≤7e-13), the bwat = 0 identity
`effecpress_stag = rhoi·g·stagthck` (≤3e-16), the basal-recovery law
`uvel(nz) = uvel_2d·beta_eff_x/beta_internal` (exact) and the full 3-D
profiles (≤3e-16) — and emits the numeric test tuples pinned in
`stage2/cism_velo.esm` (136 assertions, suite 2004P/0F/0E).  Note: the
beta chain (`beta_raw → beta_internal → beta_eff_x/y → uvel_2d`) is asserted
at rel 1e-5 because the dumped `beta_internal` is the last *accepted*
calcbeta pass of the Picard iteration while the exit velocities come from the
next iteration's accepted `beta_eff` (worst chain deviation 6.0e-7 over the
128 rows); `cism_instrument.json`'s `dump_files` registry gains the
`cism_velo_nbr.csv` entry with this commit.

## 5. Risks

- **gland20.input_c150415.nc is a 2015-era file** — if CISM 3.0 rejects it
  (format/version check), fall back to gris4 (4 km, current file present).
- **Single-PE NUOPC run** — untested in previous lanes (CICE ran 36 PEs); if
  the driver or a component rejects TOTALPES=1, move to 1 PE per component
  (9 PEs) — the dump stays PET-0-only and cell numbering must then use
  CISM's own decomposition (PET 0's local slab).
- **DLND stream**: the SCPL mode needs the exact file pattern
  `…cpl.hl2x1yr_glc.%y-01-01.nc`; wrong RUN_STARTDATE/YR_ALIGN would read
  out-of-range years (cycle mode covers it).
- **CIME namelist trap**: `dump_nml` must be added to `_CISM_IN_GROUPS` in
  `cime_config/buildnml`, or nmlgen validates but never writes it (the CICE
  lane hit the same).
- The shared queue is contended with other lanes' runs; submit early, rerun
  late.
