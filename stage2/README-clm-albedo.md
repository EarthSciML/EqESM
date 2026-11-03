# CLM surface albedo — ground albedo + bare-patch shortwave partition

Stage-2 implementation of the CLM surface-albedo piece of SurfaceAlbedoMod.F90
(ctsm5.4.063, cesm3.0-alphabranch; checkout
`stage1/cesm/components/clm`): the ground albedo chain — the SoilAlbedo
soil law, the albsat/albdry soil color-class tables, and the snow/soil
ground-albedo blend — plus the bare-patch branch of the SurfaceRadiation
absorbed/reflected shortwave partition. `clm_surfalbedo.esm` is the component
file, with the pointwise laws imported by reference from the shared
`clm_templates.esm` library; the tests pin the stage-1 instrumentation dumps
(`stage1/dumps/clm/clm_surfalbedo.csv` and `clm_surfrad.csv`; the W3
extension carries the SoilAlbedo wetness input and color class as columns).

The traced run is the D-class single-point column experiment of
README-clm.md: 16 patches on 16 columns (bare, deep snowpacks), 12 steps of
dtime = 1800 s. The SurfaceAlbedo dump traces the `nourban` and `urban`
column filters separately, 11 radiation steps (2..12; doalb runs at the end
of the step, clm_driver.F90:1325-1364); the SurfaceRadiation dump traces the
non-urban patches, 12 steps (SurfaceRadiation runs early in the step,
SurfaceRadiationMod.F90:658, so step s consumes the albedo computed at step
s−1's doalb).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `clm_templates.esm` | 4 new pointwise laws: `soil_albedo_wetness` (inc = max(0.11 − 0.40·h2osoi_vol, 0)), `soil_albedo` (albsod = min(albsat + inc, albdry)), `ground_albedo_direct` (albgrd = albsod·(1−frac_sno_albedo) + albsnd·frac_sno_albedo), `ground_albedo_diffuse` (albgri = albsoi·(1−frac_sno_albedo) + albsni·frac_sno_albedo) | `SurfaceAlbedoMod.F90`:1244, 1249, 816-817 |
| `clm_surfalbedo.esm` | `function_tables` albsat_20_vis/nir, albdry_20_vis/nir (the mxsoil_color = 20 branch; the reference run's surfdata has mxsoil_color = 20, so the 8-class branch at F90:204-208 is not transcribed), 1-D nearest over the discrete color class | `SurfaceAlbedoMod.F90`:210-218 |
| | Model `GroundAlbedo` — the SoilAlbedo chain + the snow/soil blend: albsat/albdry table lookups by color class, inc, albsod = min(albsat + inc, albdry), albsoi = albsod (the direct/diffuse soil albedos are the same value, F90:1250), albgrd/albgri blends; 13 equations | `SurfaceAlbedoMod.F90`: SurfaceAlbedo 1118-1128, 1244-1250, 816-817 |
| | Model `SurfaceRadiationBare` — the bare-patch (ptype = 0) branch of the shortwave partition: sabv = 0, sabg/sabg_soil/sabg_snow = Σ_ib forc_rad·(1−albedo), fsa = sabv + sabg, fsr = Σ_ib albedo·forc_rad, parveg = 0; 7 equations | `SurfaceRadiationMod.F90`:688-711, 897-899; `SurfaceAlbedoMod.F90`:1118-1128 |

Both models are the pointwise laws written scalar (one column / one patch of
the stage-2 skeleton), the same convention as the other table-using stage-2
law (cam_wv_saturation.esm): the surface-data color class is a scalar table
query. The stage-3 top-level model lifts them onto the column/patch index
sets. `GroundAlbedo`'s albsnd/albsni (the SNICAR snow albedos) are inputs
here — SnowSnicarMod is a separate component.

## Laws

**Soil albedo (SoilAlbedo, SurfaceAlbedoMod.F90:1241-1251).** For soil/crop
columns with coszen > 0,

```
inc    = max(0.11 − 0.40·h2osoi_vol(c,1), 0)
albsod = min(albsat(isoicol(c),ib) + inc, albdry(isoicol(c),ib))
albsoi = albsod
```

with h2osoi_vol(c,1) the top-layer volumetric soil water [m3/m3] and
isoicol(c) the surface-data soil color class (1..mxsoil_color). The albsat/
albdry tables (F90:210-218, mxsoil_color = 20 branch) give the saturated
(wet) and dry soil albedo of each color class per waveband:

| class | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| albsat vis | .25 | .23 | .21 | .20 | .19 | .18 | .17 | .16 | .15 | .14 | .13 | .12 | .11 | .10 | .09 | .08 | .07 | .06 | .05 | .04 |
| albsat nir | .50 | .46 | .42 | .40 | .38 | .36 | .34 | .32 | .30 | .28 | .26 | .24 | .22 | .20 | .18 | .16 | .14 | .12 | .10 | .08 |
| albdry vis | .36 | .34 | .32 | .31 | .30 | .29 | .28 | .27 | .26 | .25 | .24 | .23 | .22 | .20 | .18 | .16 | .14 | .12 | .10 | .08 |
| albdry nir | .61 | .57 | .53 | .51 | .49 | .48 | .45 | .43 | .41 | .39 | .37 | .35 | .33 | .31 | .29 | .27 | .25 | .23 | .21 | .16 |

**Ground albedo (SurfaceAlbedoMod.F90:816-817).** The snow/soil blend

```
albgrd = albsod·(1 − frac_sno_albedo) + albsnd·frac_sno_albedo
albgri = albsoi·(1 − frac_sno_albedo) + albsni·frac_sno_albedo
```

with albsnd/albsni the SNICAR direct/diffuse snow albedos (albsnd_hst/
albsni_hst at the SurfaceRadiation call, F90:877-879). With
frac_sno_albedo = 1 the blend collapses to albgrd = albsnd, albgri = albsni —
the endpoint the dump exercises (every `nourban` row has
frac_sno_albedo = 1).

**Bare-patch shortwave partition (SurfaceRadiation, SurfaceRadiationMod.F90:
688-711, 897-899).** On a bare patch (ptype = 0) the two-stream inputs are
fabd = fabi = 0, ftdd = ftii = 1, ftid = 0 and albd = albgrd, albi = albgri
(SurfaceAlbedoMod.F90:1118-1128), so

```
sabv = 0,  parveg = 0
sabg      = Σ_ib [forc_solad_ib·(1−albgrd_ib) + forc_solai_ib·(1−albgri_ib)]
sabg_soil = Σ_ib [forc_solad_ib·(1−albsod_ib) + forc_solai_ib·(1−albsoi_ib)]
sabg_snow = Σ_ib [forc_solad_ib·(1−albsnd_ib) + forc_solai_ib·(1−albsni_ib)]
fsa = sabv + sabg,  fsr = Σ_ib [albgrd_ib·forc_solad_ib + albgri_ib·forc_solai_ib]
```

(ib = vis, nir). The vegetation (two-stream) branch — fabd/fabi/ftdd/ftii/
ftid and sabv/parveg from the Bonan et al. 2011 two-stream solver — is not
traced by the stage-1 dumps (every clm_surfrad.csv row is a bare patch,
ptype = 0, sabv = parveg = 0) and is not implemented here; see
[Two-stream: what the dump would need](#two-stream-what-the-dump-would-need).

## Dump analysis

### Both kinds of clm_surfalbedo.csv

The 352 rows are (16 nourban + 16 urban) columns × 11 steps. The two kinds
pin different things:

- **nourban rows (176).** frac_sno_albedo = 1 everywhere, albsod = albsoi
  everywhere, and the soil albedo is constant across all rows: the traced
  columns' layer-1 water and color class do not change over the window. The
  soil law pins from these rows; the ground-albedo *blend* appears only at
  its all-snow endpoint (albgrd = albsnd, albgri = albsni).
- **urban rows (176).** frac_sno_albedo = 0 and albgrd ≠ albsod: for urban
  columns the SurfaceAlbedo ground albedo is overridden by the UrbanAlbedo
  sref output (UrbanAlbedoMod.F90:442-463). Each urban row of step s equals
  the clm_urbalb.csv row of step s−1 to round-off on all 160 rows with a
  predecessor (steps 3..12; the 16 step-2 rows have no predecessor because
  doalb(1) used the restart albedo). These rows pin no law of this component
  — the urban albedo law is the UrbanAlbedo component's — and they document
  the SurfaceAlbedo wiring order (the UrbanAlbedo call overwrites albgrd/
  albgri after the soil/snow blend, SurfaceAlbedoMod.F90:1254-1270).

### Wetness input and color class (nourban pins)

The W3 extension adds the law's inputs as columns (h2osoi_vol_1 and
isoicol per row), so the pins take them directly from the dump:

- **h2osoi_vol_1 = 0.1505 m3/m3** — dumped as
  1.50499999999999995E-01 on every nourban row. The pre-extension
  reconstruction cross-checks: clm_soiltemp.csv (step 2, exit, j = 1)
  h2osoi_liq = 0.01 kg/m2, h2osoi_ice = 2.751 kg/m2, dz = 0.02 m,
  evaluated as the Hydrology definition h2osoi_liq/(dz·denh2o) +
  h2osoi_ice/(dz·denice) (HydrologyNoDrainageMod.F90:580, denh2o = 1000,
  denice = 917 kg/m3) gives the same value bit-exactly (the traced columns'
  layer-1 water is bit-identical from soiltemp exit(step) to entry(step+1),
  so this is exactly the value SoilAlbedo reads at doalb(step)).
- **isoicol = 10** — dumped on every nourban row, matching the surfdata
  SOIL_COLOR = 10 of the traced gridcells
  (surfdata_1.9x2.5_hist_2000_78pfts_c251022.nc, mxsoil_color = 20); the
  pre-extension argument (class 10 is the only class consistent with the
  dumped albsod at h2osoi_vol = 0.1505) still holds.

Then inc = max(0.11 − 0.40·0.1505, 0) = 0.0498, and the law pins the dumped
albsod_vis/albsod_nir = 1.89800000000000024E-01 / 3.29800000000000038E-01
to the last bit on all 176 nourban rows (the identical row values across
steps/columns make the law constant; the three soil tests pin steps 2/7/12,
columns 1/8/16).

### The mixing-law endpoint (nourban pins)

The W3 albsnd/albsni columns let the all-snow endpoint of the ground-albedo
blend be pinned with the real SNICAR values: on every nourban row the
dumped albsnd = albgrd and albsni = albgri bit-exactly (all rows carry
frac_sno_albedo = 1, so the F90:816-817 blend collapses to albgrd = albsnd,
albgri = albsni; the equality also confirms there are no further
contributions to albgrd/albgri on the traced bare columns). The two mixing
tests pin steps 2 and 7, columns 1 and 8 — the albsnd/albsni values differ
between steps because the SNICAR snow albedo ages (e.g. albsnd_vis:
9.84056006521105808E-01 at step 2 → 9.84760966905754653E-01 at step 7).

### clm_surfrad.csv (bare patches)

All 192 rows are bare (ptype = 0, sabv = parveg = 0 exactly) and
sabg_snow = sabg on every row (the frac_sno_albedo = 1 degeneracy of the
mixing law: albsnd_hst = albgrd, albsni_hst = albgri). Step s uses the
albedo of doalb(s−1): with the previous step's clm_surfalbedo.csv albedo
columns and the row's forc_solad/forc_solai, the sabg/sabg_soil/fsr law
matches the dump to ≤ 3.2e-16 relative on all 160 rows with a predecessor
(steps 3..12); the step-1/step-2 rows used the restart albedo and are not
pinned. The three SurfaceRadiationBare tests pin steps 3/7/12, patches
1/8/16 (forcing from clm_surfrad.csv, albedos from clm_surfalbedo.csv
steps 2/6/11).

## Tests

All 8 tests use time_span {0,0} (instantaneous law pins), rel 1e-12,
and cite their dump rows in the description:

| Test | Pins | Dump rows |
|---|---|---|
| `clm_groundalbedo_soil_step2_c1` | albsat/albdry table entries (class 10), inc, albsod/albsoi both bands | clm_surfalbedo.csv step 2 col 1 nourban (h2osoi_vol_1, isoicol columns) |
| `clm_groundalbedo_soil_step7_c8` | albsod/albsoi both bands (same law, later step/column) | clm_surfalbedo.csv step 7 col 8 nourban |
| `clm_groundalbedo_soil_step12_c16` | albsod/albsoi both bands (last traced step) | clm_surfalbedo.csv step 12 col 16 nourban |
| `clm_groundalbedo_mix_allsnow_step2_c1` | albgrd/albgri both bands = the dumped albsnd/albsni (all-snow endpoint of the blend) | clm_surfalbedo.csv step 2 col 1 nourban (albsnd/albsni columns) |
| `clm_groundalbedo_mix_allsnow_step7_c8` | same 4 outputs (aged snow albedo, so different values) | clm_surfalbedo.csv step 7 col 8 nourban |
| `clm_surfrad_bare_step3_p1` | sabv, sabg, sabg_soil, sabg_snow, fsa, fsr, parveg | clm_surfrad.csv step 3 patch 1; clm_surfalbedo.csv step 2 col 1 |
| `clm_surfrad_bare_step7_p8` | same 7 outputs | clm_surfrad.csv step 7 patch 8; clm_surfalbedo.csv step 6 col 8 |
| `clm_surfrad_bare_step12_p16` | same 7 outputs | clm_surfrad.csv step 12 patch 16; clm_surfalbedo.csv step 11 col 16 |

The surfrad tests parameterize albsnd/albsni from the albsnd/albsni columns
of the previous step's clm_surfalbedo.csv row (which equal albgrd/albgri
bit-exactly at frac_sno_albedo = 1, and the dump's sabg_snow = sabg equality
corroborates the partition on all 192 rows).

## Two-stream: what the dump would need

The two-stream vegetation branch (SurfaceAlbedoMod.F90:1280-1350,
TwoStream 1351-1510: the Sellers/ Dickinson two-stream with the Bonan et al.
2011 clumping and canopy integration; SurfaceRadiationMod.F90:688-711 for
fabd/fabi/ftdd/ftii/ftid and the sabv/parveg partition) is not traced: every
clm_surfrad.csv row has ptype = 0, sabv = parveg = 0. To pin it, the
stage-1 instrumentation would need to trace vegetated patches
(frac_veg_nosno = 1, elai/esai > 0) with the two-stream inputs — laisun/
laisun_frac, clump_frac, coszen, forc_solad/forc_solai, and the leaf
optics (rhol/rhos/tau/taus by band, from CLM canopy structure and optics
state) — and outputs (albd/albi, fabd/fabi/ftdd/ftii/ftid, sabv/parveg).
The traced case has bare patches only, so the vegetation branch is left to
the stage-2 task that instruments a vegetated patch; the bare branch is
fully pinned here.

## Fork dump extension

To pin the mixing law and the soil-law inputs exactly (rather than through
the reconstruction/degeneracy arguments above), the fork instrumentation
extends `dump_surfalbedo` with the SNICAR snow albedos, the top-layer
volumetric soil water and the soil color class per column — fork clm
`3d7038304` "CLM dumps W3: surfalbedo slice + soil color/volumetric water +
SNICAR snow albedo columns" / superproject `527ab5d`, branch
`instrument/clm-dumps`:

- new columns: `albsnd_vis, albsnd_nir, albsni_vis, albsni_nir`
  (SNICAR direct/diffuse snow albedo, both bands), `h2osoi_vol_1` (the
  SoilAlbedo wetness input) and `isoicol` (soil color class), appended to
  the existing 16 columns;
- call-site wiring in SurfaceAlbedoMod.F90 (the dump call gains
  albsnd/albsni from the SNICAR arrays and h2osoi_vol/isoicol from the
  waterstate/soilstate types; urban columns, which skip SNICAR, write NaN
  in the albsnd/albsni slots — SurfaceAlbedoMod.F90:1134-1139 runs
  SnowAlbedo only for non-urban, non-lake columns).

**Rerun.** The instrumented rerun (in-node `case.submit --no-batch` inside
allocation 12027802, LID 261103-133002) regenerated all 11 dumps; the
shared columns are bit-identical to the previous W3 CSVs (the extension is
dump-only; the model trajectory is unchanged) and the new columns verify
the pins above. The first rerun attempt (LID 261103-131613) aborted in
`dump_surfalbedo`'s urban branch — the module's derived-type NaN parameter
(`nan => shr_infnan_nan`, a `type(shr_infnan_nan_type)` with a LOGICAL
component) was passed directly in a formatted es26.17 write list, which
gfortran rejects at runtime; fixed by materializing a real local through
the module's `assignment(=)` interface (bugs.md TBD entry). The run details
are documented in `clm_instrument.json` (changelog) and CTSM_LANE.md (W2c).

## Coverage gaps

The 8 tests pin exactly what the dumps support; the following branches are
transcribed or documented but not dump-pinned:

- **Mixing-law interior 0 < frac_sno_albedo < 1.** Every nourban row is at
  the all-snow endpoint (frac_sno_albedo = 1) and every urban row at the
  snow-free endpoint is the UrbanAlbedo override, so the interior blend and
  the soil endpoint of albgrd/albgri carry no direct pin (the law itself is
  the simple linear blend of F90:816-817, written as templates). A state
  with partial snow cover (or a nourban snow-free column) would pin it.
- **The soil endpoint of the blend (frac_sno_albedo = 0, nourban).** As
  above; also the only place albgrd = albsod would be visible.
- **albsnd/albsni as SNICAR outputs.** The albsnd/albsni *values* are now
  pinned from the W3 columns (at the frac_sno_albedo = 1 endpoint, and
  corroborated by sabg_snow = sabg on all 192 surfrad rows); the SNICAR
  laws that produce them (SnowSnicarMod) are a separate component.
- **The 8-class albsat/albdry branch (F90:204-208).** The reference run's
  surfdata has mxsoil_color = 20, so only the 20-class tables are
  transcribed.
- **SurfaceRadiation steps 1-2.** Their albedo inputs come from the restart
  (not the traced doalb chain), so they are not pinned.
- **Vegetated two-stream branch.** See
  [Two-stream: what the dump would need](#two-stream-what-the-dump-would-need).
- **coszen ≤ 0 (night) rows.** doalb skips the albedo computation there;
  the dump traces only the daytime window.
