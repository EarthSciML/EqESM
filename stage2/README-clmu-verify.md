# CLMU urban-leaf verification — EarthSciModels `components/urban_canopy/`

Verification of the CLMU (urban) leaves of
[EarthSciModels](https://github.com/EarthSciML/EarthSciModels)
`components/urban_canopy/` against the stage-1 CLM instrumentation dumps
(`stage1/dumps/clm/`). Companion to `README-clm.md`, which covers the
vegetated CLM leaves (`clm_canopy.esm`); this README covers the urban leaves
only. Work done on EarthSciModels branch `clmu-verify2` (worktree
`EarthSciModels-wt-clmu-verify2`), EqESM branch `wt/clmu-verify`.

## Dump inventory (what the dumps can and cannot pin)

Twelve CLM dumps exist. Only three trace urban (CLMU) columns/landunits:

| Dump | Rows | Urban coverage |
|---|---|---|
| `clm_urbalb.csv` | 880 | UrbanAlbedo exit (UrbanAlbedoMod.F90): per-column albd/albi/albgrd/albgri and per-unit-incident direct-beam absorbed solar `sabs_*_dir` (band 1 = vis), for landunits 252-267, all 5 ctypes (71 roof, 72 sunwall, 73 shadewall, 74 road_perv, 75 road_imperv), 11 steps. Landunit 255 is the daytime landunit (coszen > 0 at steps 2-12); the rest are nighttime (albedo = view factors). |
| `clm_urban_rad.csv` | 960 | UrbanRadiation exit, per patch (16 landunits × 5 ctypes × 12 steps): coszen, forc_solad/forc_solai (gridcell), sabg/sabv/fsa/fsa_u, eflx_lwrad_out/eflx_lwrad_net. |
| `clm_urban_lw.csv` | 192 | UrbanRadiation net_longwave exit, per landunit (16 landunits × 12 steps): canyon_hwr, lwdown, t_roof/t_improad/t_perroad/t_sunwall/t_shadewall, lwnet_*/lwup_* per surface and canyon. |

The other nine dumps trace only the vegetated/bare (ctype = 1) branch:
`clm_soiltemp.csv` (t_soisno/tk/cv/fact/imelt/xmf/h2osoi_liq/h2osoi_ice/
qflx_snomelt_lyr/qflx_snofrz_lyr/t_grnd/tkmg), `clm_soilfluxes.csv`
(eflx_*), `clm_run_imports.csv` (htvp/thv/thm/qg/dqgdT/soilalpha),
`clm_canopyflux.csv`, `clm_bareground.csv`, `clm_surfrad.csv`,
`clm_surfalbedo.csv` and the scalar `clm_qsat.csv` — all ctype = 1 only.
`clm_surfalbedo.csv` additionally has a `kind=urban` section (176 rows,
landunits 252-267, 1 column each), but it is a pass-through of the
UrbanAlbedoMod outputs (albgrd/albgri for urban columns equal the
`clm_urbalb.csv` rows; albsnd/albsni are NaN because urban columns do not
run SNICAR; the soil-albedo inputs albsod/albsoi = 0.2 are the non-urban
defaults), so it pins no urban physics beyond `clm_urbalb.csv`.

**Coverage boundary:** the dumps cover exactly the CLM urban **radiation**
(shortwave albedo/absorption and longwave net/up at each surface). Nothing
is dumped for the urban temperature solve, urban hydrology, the
canyon-atmosphere interface, or the offline forcing path — see
[Coverage gaps](#coverage-gaps).

## What was verified (branch `clmu-verify2`, this round)

The remaining radiation laws of `urban_radiation.esm` are now pinned to the
dumps. Test counts: `urban_radiation.esm` carries 23 tests / 77 assertions,
10 of them CESM-dump-grounded; `./esm test components/urban_canopy/` from the
EarthSciModels worktree: **1078 assertions pass, 0 fail, 0 err**.

| Test (urban_radiation.esm) | Dump row | Pins | Tolerance |
|---|---|---|---|
| `cesm_dump_canyon_shortwave_daytime` | clm_urbalb.csv l=255 step 2 (tod 3600), ctypes 71-75 | S_net_roof, S_net_road_dir_vis, S_net_sunwall_dir_vis, S_net_shdwall_dir_vis, S_refl_sky_dir_vis, α_uc_dir_vis at coszen = 0.91580348106827802 (θ_0 saturates at π/2) | roof rel 1e-12; canyon rel 1e-5 |
| `cesm_dump_canyon_shortwave_low_sun` | clm_urbalb.csv l=255 step 12 (tod 21600) | same six outputs at coszen = 0.39402549836194734 (partial sunlit: θ_0 = 0.31121988851515503) | same |
| `cesm_dump_canyon_longwave_hwr0p8` | clm_urban_lw.csv l=267 step 1 | L_net_sunwall, L_net_shdwall, L_net_road, L_net_uc (abs 0.001); L_up_roof, L_net_roof (rel 1e-12) | abs 0.001 canyon; rel 1e-12 roof |
| `cesm_law_direct_beam_factors`, `cesm_forcing_canyon_shortwave_multi_reflection` | (descriptions updated) | unchanged pins; the "no daytime canyon shortwave dump available" claim is now false — the tests remain closed-form pins of the intermediate direct-beam factors and of the full 4-band solution with generic albedos | unchanged |

The earlier round (commit `e21ee29`) pinned the view factors (nighttime
albedo == view factors at hwr = 4.8), the daytime roof shortwave, the roof
longwave, and the canyon longwave at hwr = 1.4 (l=255) and hwr = 4.8 (l=252).

### Method notes

- **Unit-incident direct beam.** The dumped `sabs_*_dir` are solar absorbed
  per unit incident direct flux, band 1 (vis), with sdir = 1
  (UrbanAlbedoMod.F90:1309-1312, dump description), so the SW tests run the
  .esm with S_atm_dir_vis = 1 and all other incident fluxes 0. Expected
  values are the dumped sabs (roof/sunwall/shadewall rows directly, road =
  wt_imp·sabs_improad + wt_perv·sabs_perroad), and the dumped albd_vis
  combination 0.5·albd_75 + 0.5·albd_74 + H_W·(albd_72 + albd_73) for the
  canyon reflected radiation.
- **Surfdata intrinsic albedos.** The l=255 group's intrinsic albedos come
  from `surfdata_1.9x2.5_hist_2000_78pfts_c251022.nc`: ALB_ROOF_DIR = 0.14,
  ALB_WALL_DIR = 0.2166, ALB_IMPROAD_DIR = 0.13, ALB_PERROAD_DIR = 0.08,
  WTROAD_PERV = 0.5. The pervious/impervious road split reduces *exactly* to
  a single effective road of albedo wt_imp·0.13 + wt_perv·0.08 = 0.105:
  CLM combines the two roads' reflected fluxes linearly before distributing
  to the walls, so the two-road iteration and the single-road closed form
  agree to ~2e-16. The dumped albd_vis roof row (float32(0.14)) and the
  wtroad-weighted lwnet_canyon combination independently confirm ALB_ROOF_DIR
  and WTROAD_PERV. No input is invented: every test input is a dumped value
  or a surfdata value.
- **Longwave at l=267.** All-pervious road (wtroad_perv = 1.0), so the road
  reduces to the pervious road alone (t_road = t_perroad,
  ε_road = float32(0.95)); em_wall = 0.9008 (surfdata, cross-checked on the
  l=255 rows: a single em_wall reproduces all three canyon surfaces to the
  iteration-truncation level); em_roof = float32(0.9), recovered exactly from
  lwup_roof via UrbanRadiationMod.F90:707. This row pins the shallow-canyon
  (H/W = 0.8) re-emission solution, the opposite extreme of the hwr = 4.8
  test.
- **Tolerances.** SW canyon surfaces rel 1e-5: CLM solves the canyon
  reflection iteratively (errcrit = 1e-5, n = 50), so the dumped values carry
  up to ~3.2e-7 of truncation (measured); the .esm is the exact closed form
  of the same linear system. Roof rel 1e-12: sabs_roof = 1 − albd_roof is a
  single pass (no iteration). LW canyon surfaces abs 0.001: CLM's net_longwave
  is a fixed-point iteration with a 0.001 W/m² cutoff (measured truncation
  ~5e-6); the .esm closed form is exact. LW roof rel 1e-12: single-pass
  formula (verified to 5.7e-14 abs against the dump).

## Coverage gaps

The following urban leaves of `components/urban_canopy/` are transcribed but
have **no dump coverage**; they cannot be pinned from the current stage-1
dumps. Nothing was loosened or invented to force a pin.

| Leaf | CLM Fortran | Missing fields |
|---|---|---|
| `temps/` (urban_surface_thermal, soil_thermal, snow_thermal, interface_thermal_conductivity, heat_flux, phase_change_energy, layer_phase_change_energy, total_phase_change_energy, phase_change_adjustment, snow_melt_no_layers, snow_soil_blended_heat_capacity, snow_layer_geometry, uniform_grid, exponential_grid, adjusted_layer_thickness, freezing_point_depression, surface_energy_flux, heating_cooling_flux) | SoilTemperatureMod.F90 urban branch (UrbanLayerTemperature / UrbanTemperatures), ThermalConductivityMod.F90, SnowTemperatureMod.F90 | t_soisno of the 5 urban soil/building layers and snow layers, tk, cv, fact, imelt, xmf, qflx_snomelt_lyr/qflx_snofrz_lyr on urban columns, eflx_soil_grnd/eflx_sh_grnd at canyon surfaces. `clm_soiltemp.csv` dumps exactly these quantities but ctype = 1 only — the urban layer branch is never traced. |
| `temps/building_temperature.esm`, `temps/waste_heat_ac.esm`, `temps/waste_heat_allocation.esm` | UrbBuildTempOleson2015Mod.F90 (building interior temperature T_iB, air-conditioning waste heat H_wasteheat/H_aircond) | T_iB, H_wasteheat, H_aircond, roof length L_roof. No dump traces the building-energy path. |
| `hydro/` (pervious_road_water_balance, impervious_water_balance, surface_runoff_infiltration, soil_water_flux, soil_water_content_calc, soil_water_equilibrium, water_table_depth, aquifer_water_balance, groundwater_drainage, interface_hydraulic_conductivity, soil_hydraulic_properties, surface_layer_update, impervious_runoff, snow_*) | HydrologyNoDrainageMod.F90 urban branch (PerviousRoadFluxes/ImperviousRoadFluxes), SoilHydrologyMod.F90, SnowHydrologyMod.F90 | h2osoi_liq/h2osoi_ice on urban layers, qflx_urban_runoff/qflx_imperv_runoff/qflx_over, water_input, soilalpha on urban columns (clm_run_imports.csv soilalpha is ctype = 1 only), water-table depth z_v, snow water/ice/density on urban surfaces. |
| `clmu_atmosphere.esm` | UrbanFluxesMod.F90 (Oleson Ch. 5): canyon-height atmospheric state | e_atm, rho_atm, z_atm and the derived canyon meteorology (thm/thv/htvp at canyon height, eflx_sh/lh/lwrad at canyon level). clm_run_imports.csv dumps htvp/thv/thm but ctype = 1 only. |
| `offline_mode_forcing.esm` | — (not a CESM code path) | The CLMU offline tool (Oleson 2010 Ch. 1/6) that derives S_atm_dir/dif (vis/nir), L_atm_down, q_rain/q_snow, u/v/θ at canyon height from observed radiation/weather. CESM CLM6 receives forc_* directly from the coupler and never executes this path, so no dump can exist: **unverifiable by construction**, not merely uncovered. |

`urban_surface_thermal.esm` is a trivial identity wrapper (λ_surf/c_surf
passthrough from the surface dataset, Table 1.3 of Oleson 2010); its
prescribed λ/c values are surfdata parameters and are not pinned either.

## Bugs

No CLM bugs found in this round; nothing was appended to `bugs.md`. The θ_0
units fix of commit `e21ee29` was an .esm transcription error (asin already
returns an angle; the trailing `*one_rad` produced angle²), not a CLM defect.

## Testing

- EarthSciModels worktree, branch `clmu-verify2` (HEAD `cd56c1d`):
  `./esm validate components/urban_canopy/urban_radiation.esm` ✓;
  `./esm test components/urban_canopy/` — 1078 assertions, 0 fail, 0 err
  (urban_radiation.esm: 23 tests / 77 assertions).
- EqESM worktree, branch `wt/clmu-verify`: full suite
  `./esm test stage2/` — 1171P/0F/0E (see commit message; README-only
  change, counts unchanged from baseline). Run with
  `export ESD_ROOT=/projects/illinois/eng/cee/ctessum/ctessum/code/EarthSciDiscretizations`.
