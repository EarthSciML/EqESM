# CICE delta-Eddington shortwave — stage 2

Stage-2 transcription of the CICE `shortwave = 'dEdd'` 3-band
delta-Eddington radiation scheme as the esm 1.3.0 component
`cice_shortwave.esm` (with the shared law library `cice_templates.esm`),
tested against the stage-1 instrumentation dump
`stage1/dumps/cice/cice_dedd.csv` of the D-class reference run
(`stage1/cime/cases/D_cice_inst/`).

Sources (checkout `stage1/EqESM-cesm-cice/`): `icepack_shortwave.F90`
`compute_dEdd_3bd` 1703-3013 and the per-layer solver `solution_dEdd`
3022-3488, plus the SNICAR grain-radius and bulk-aerosol data tables of
`icepack_shortwave_data.F90` (107-263, 331-353).

## Design

`cice_shortwave.esm` is a single model `Shortwave` built entirely from
expression templates imported by reference from `./cice_templates.esm`
(append-only: the 6 fresnel templates of this stage are appended at the
end of the file). The model is organized as the Fortran law, stage by
stage:

| Stage | Model unknowns | Fortran source |
|---|---|---|
| layer geometry | `rnslyr`, `rnilyr`, `dz_s`, `dz_i`, `dz_ssl_s`, `dz_ssl_i`, `dzk`, `frsnw` | F90:2132-2156, 2122-2125 |
| snow IOPs | `Qs`, `ws_snw`, `gs_snw`, `ks_snw`, `ks_l`, `ws_l`, `gs_l` | F90:2396-2424, `icepack_shortwave_data.F90` 116-232; mass extinction ks 2422-2425 |
| chlorophyll / pond / delta-layer factors | `kabs_chl`, `dz_pond`, `fs_dl`, `kabs_i`, `sig_i` | F90:2244-2251, 2550, 2571, 2586-2592 |
| aerosol mixing regions | `mix_region`, `aero_slot_idx`, `aero_scale`, `aero_am` | F90:2451-2460, 2605-2612 |
| aerosol optical depths | `taer`, `waer`, `gaer` (reduce over the 3 species) | F90:2451-2547, 2605-2717 |
| layer optical properties | `tau_base`, `w0_base`, `g_base`, then mixed `tau`, `w0`, `g` | F90:2425-2547, 2562-2717 |
| ocean albedo, fresnel geometry | `albodr`, `mu0`, `mu0nij`, `kfrsnl`, `mu0n` | F90:2772-2775, 3251-3266 |
| delta-Eddington layer solution | `ts`, `ws`, `gs`, `lm`, `ue`, `extins`, `ne`, `rdif_a0`, `tdif_a0`, `trnlay0`, `alp0`, `gam0`, `apg0`, `amg0`, `rdir0`, `tdir0`; gauss quadrature `alp_g`/`gam_g`/`apg_g`/`amg_g`/`trn_g`/`rdr_g`/`tdr_g`, `swt`, `smr`, `smt`, `rdif_g`, `tdif_g` | F90:3294-3351 |
| fresnel (refraction) layer | `R1`-`T2`, `Rf_dir_a`, `Tf_dir_a`, `Rf_dif_a/b`, `Tf_dif_a/b`, `rintfc`, combined `rdir`, `tdir`, `rdif_a`, `rdif_b`, `tdif_a`, `tdif_b`, `trnlay` | F90:3365-3406 |
| downward/upward recurrences | `trndir_f`, `trntdr_f`, `rdndif_f`, `trndif_f`, `active`; guarded `trndir`, `trntdr`, `rdndif`, `trndif`, `m_*`; ocean-up `rupdir`, `rupdif` | F90:3285-3287, 3408-3466 |
| flux fractions and band sums | `dfdir`, `dfdif`, `flux`, `w2`, `wght`, `swdr`, `swdf`, `ksrf`, `klevp`, `avdr`, `avdf`, `aidr`, `aidf`, `fsfc`, `fint`, `fthru`, `fthruvdr`, `fthruvdf`, `fthruidr`, `fthruidf`, `Sabs`, `Iabs`, `fthrul` | F90:2804-2832, 2834-2960 |

### Conventions (model units vs Fortran)

- **Layer axis:** model layer `k` = Fortran layer `k-1`. The sea-ice
  surface scattering layer (SSL) is model layer `N_SLYR+2`, the delta
  layer `N_SLYR+3`, the fresnel layer `kfrsnl = N_SLYR+3` (Fortran
  `nslyr+2`). All layer-branch conditions use `N_SLYR+2`/`N_SLYR+3` (the
  stage-2 Class-A fix: the first transcription used `N_ILYR+1/2`, which
  coincided with the truth only for the ice layers).
- **Downward interface arrays** are top-down: model iface 1 = the top,
  Fortran iface 0. **`rupdir`/`rupdif` are ocean-up**: model iface 1 =
  the ocean base (albodr), iface `klevp = N_SLYR+N_ILYR+3` = the surface
  (the Fortran's `rupdir(klevp)` = albedos). The `dfdir`/`dfdif` flux
  fractions index `rupdir/rupdif` at `klevp+1-k`; `avdr`/`avdf`/`aidr`/
  `aidf` take the surface value at `klevp`.
- **Aerosol mass path** `aero_mp` is species-major: `aero_mp(sp, slot)`
  with slots snow-SSL, snow interior, ice SSL, ice interior — the
  `aero_mp(1:12)` of the Fortran with the species index `na` stepping 4
  (`F90:1476-1485`). The stage-1 dump columns are likewise species-major
  (`aero_mp_{(sp-1)*4+slot}`).
- **`ksrf`** = 2 for snow, `N_SLYR+3` otherwise (Fortran `ksrf = 1` /
  `nslyr+2` + 1 for model units).
- **`mu0n`** is `mu0nij` below the fresnel level and `mu0` above it for
  `srftyp < 2` (F90:3302-3305); the model carries it per (band, layer).

### Arithmetic notes (bit-level)

The Fortran arithmetic was reproduced bit-for-bit; three associations
matter and are documented in the templates:

- `smr/smt` accumulate `mu*(rdr*gwt)` (the gwt multiplies rdr before the
  mu scaling, F90:3347-3348).
- `alpha`/`agamm` denominators are `1 - (e*e)*(mu*mu)` (squared `e` times
  squared `mu`, F90:3321-3322, 3343-3344).
- `exp(-min(10, x))` clamps the extinction/transmission arguments at
  `exp_argmax = 10` (F90:3307-3308, 3339-3340, 3317-3318).

The 8-point Gauss-Legendre angular integration of the diffuse
reflectivity/transmissivity is written as an explicit sum over the 8
abscissas rather than an esm integral: no EarthSciDiscretizations
quadrature rule exists for a bare `mu` axis (an integral reaches
evaluation only when a rewrite rule lowers it — probed on this branch).
Similarly the downward/upward recurrences are explicit unrolled
`faq`/`ifelse` chains over the layer axis (13 layers), not integrals.

The Fortran's `trmin` guard (F90:3285-3287) skips the layer solution
below the activity cutoff. Here the full solution is computed and the
outputs are masked with the same condition (`active = trntdr_f > trmin`
for the downward chain, the `m_*` layer props and the upward recurrence
inputs), which reproduces the Fortran bit-for-bit: the guard condition
is a property of the full solution.

## Transcription bugs caught by the tests

- **tdr_g sign** (`cice_templates.esm` `dEdd_direct_transmissivity`):
  the first transcription computed `apg*T1 + amg*R1*trn - apg*trn + trn`
  as `(apg*trn - trn) + trn = apg*trn`; the Fortran (F90:3346) is
  `- apg*trn + trn`. This corrupted every band-2/3 `tdir0`/`tdif_g`.
- **ksrf**: `fsfc` was identically zero for snow because `ksrf` was
  `N_SLYR+3` for all surface types; Fortran F90:2236-2240 gives `ksrf = 1`
  (model 2) for snow.
- **aerosol table transpose** (`cice_shortwave.esm` `aerosol_3bd`): the
  kaer/waer/gaer data rows were written species-major while the axes are
  declared (band, aero) and the Fortran array is `(nspint_3bd, max_aero)`
  — a transposed data literal. The lookup itself was correct, so every
  layer whose aerosol slot was nonzero got the wrong mixed `tau/w0/g`:
  the snow SSL/interior and ice SSL/DL of the snow group and the ice SSL
  of the bare group (band 3 layer 5). The error was masked by the pin
  curation (all pinned `tau/w0/g`/`m_*` were at the zero-aerosol ice
  interior layers) and surfaced as residual ~1e-8 layer-solution and
  ~5e-8 scalar discrepancies, localized by pinning the intermediate
  `taer/waer/gaer` against the Python reference. The data is now
  band-major: `kaer(ns, na)` with `kaer(1, :) = [11580.61872,
  25798.96479, 196.49772]` (F90:331-337).
- **aero_mp dump layout** (`stage1/analysis/gen_cice_dedd_tests.py`):
  the stage-1 dump writes the aerosol columns species-major (index
  `(sp-1)*4+slot`), not slot-major; the test generator originally read
  them slot-major, which put sp2/sp3 mass paths in the wrong slots.
- **rupdir convention**: the first transcription of `dfdir`/`dfdif` used
  the top-down index into the ocean-up `rupdir`/`rupdif`, giving
  `dfdir(1,1) = 0.9687` instead of `0.0313`; the binding is
  `index(rupdir, b, klevp+1-k)`.

None of these are CICE bugs — all were transcription errors in the esm
files, fixed on the branch. (No `bugs.md` entries.)

## Coverage and gaps

Pinned (189 assertions, all against the instrumented dump):

- layer geometry `dzk` at the snow SSL / snow interior / ice SSL / delta
  layer / ice interior (bare-group snow layers are identically zero and
  unpinned); `frsnw` all 3 snow layers
- snow-table IOPs `Qs`/`ws_snw`/`gs_snw`/`ks_snw` at snow layers {1,3}
  × all bands (snow group only)
- mixed `tau/w0/g` at layer 6 × all bands (snow) / layer 6 × all bands
  (bare; the bare SSL values are identically zero and unpinned)
- the guard-zeroed layer props `m_*` at one active layer per band (band
  1 layer 6, bands 2/3 layer 1 for snow; band 1/2 layer 6 and band 3
  layer 5 for bare — the last active layer before the cutoff), including
  `m_* = 0` pins below the cutoff
- the downward chain `trndir`/`trntdr`/`rdndif`/`trndif` at the top and
  bottom of band 1 plus the first guarded interface of bands 2 and 3
  (expected 0, exercising the trmin guard; the iface-1 `rdndif` boundary
  value 0 is unpinned)
- the upward chain `rupdir`/`rupdif` at the top of the column (iface 14)
  for bands 1 and 3; the ocean-interface (iface 1) values equal the band
  ocean albedo `albodr` by construction and are unpinned
- `dfdir`/`dfdif` straddling the cutoff, the scalar band sums
  `avdr`/`avdf`/`aidr`/`aidf`, `fsfc`/`fint`/`fthru`,
  `fthruvdr`/`fthruvdf` and the band-3 row `fthruidr`/`fthruidf`
- `Sabs` snow layers {1,2,3} (snow only), `Iabs` ice layers {1,2,8},
  `fthrul` levels {1,5,9}

The extractor (`gen_cice_dedd_tests.py`) skips the unfailable pins with
tolerance-aware copy/zero checks (|value − reference| ≤ rel·|value|): the
ocean-interface `rupdir`/`rupdif` (a copy of `albodr`), the iface-1
`rdndif` (= 0) and the bare-group `dzk` layers 1-2 and `tau/w0/g` at
layer 1 (the hs = 0 zero region).

Tolerances: rel 1e-12 except rel 1e-9 on the interface quantities
(`trndir`/`trntdr`/`rdndif`/`trndif`/`rupdir`/`rupdif`/`dfdir`/`dfdif`).
The interface values are accumulated over 13 layers of the upward/downward
recurrences, each step dividing by `1 - rupdif*rdif_b` terms that are
small for the highly scattering snow-covered column, so round-off
amplification of the ~1e-16 per-layer algebra reaches ~1e-9 at the
surface; the layer solutions themselves pin at 1e-12. (The band-3 near
snow-SSL algebra is the worst case: `ue` is small, `extins` is tiny and
`ne = (u+1)²/extins - ...` amplifies relative error by ~1e4.)

Gaps (not pinned, and why):

- **Pond branch** (`srftyp = 2`): transcribed but untested — the
  instrumented D-class run contains no pond rows. `kw`, `ki_p_ssl`,
  `ki_p_int` and the pond `w0_base`/`g_base` zero branches are
  therefore only reviewed by eye.
- **`fthruuvdr`/`fthruuvdf`** (F90:2861-2872, 2924-2933): the
  absorbed-in-surface-visible flux fractions used only by the UV
  shortwave scheme; omitted from the model (the UV branch of the CICE
  shortwave is out of scope here).
- **`pardr`/`pardf`** (F90:2910-2911): the pond area weighted albedo
  adjustments; the pond branch is untested, so the pins stop at
  `aidr`/`aidf`.
- **Bare-group snow IOPs**: `Qs/ws_snw/gs_snw/ks_snw` are not pinned for
  `srftyp = 0` (the snow branch is inactive and the values are
  identically zero); `Sabs` is skipped for bare ice (identically zero).
- **Aerosol species 4-6**: the instrumented run's aerosol fields carry
  only species 1-3 (`n_aero = 3`, the aerosol loop F90:2451-2547 runs
  `na = 1..4*n_aero`); the table holds only those species (the Fortran
  `max_aero = 6` data beyond species 3 is not transcribed).
- **Integral form**: the Gauss angular integration and the layer
  recurrences are explicit sums/chains (see above); the Fortran's
  per-layer `rdir/tdir` products inside the downward chain are exact
  algebraic products of the fresnel-combined layer props, not integrals.

## Testing

`stage1/analysis/gen_cice_dedd_tests.py` regenerates the curated tests
block (the tuple values are extracted from `stage1/dumps/cice/cice_dedd.csv`
into the hand-authored block; the equations, structure and factoring are
hand-authored). The two test groups cover the snow-covered category
(step 1, i=2, j=4, ncat 2, srftyp 1: hs = 0.15364 m, hi = 1.2089 m,
frsnw = 290.36 um, aerosols in all four slots) and the bare-ice category
(same cell, ncat 1, srftyp 0: hs = 0, hi = 0.34181 m). Both run with
`coszen = 0.21184` (mu0 clamped), `fnidr = 0.6458`, `kalg = 0`.

Run:

```
export ESD_ROOT=/projects/illinois/eng/cee/ctessum/ctessum/code/EarthSciDiscretizations
./esm test stage2/cice_shortwave.esm
```

The Python reference used during development (association-level
validation against the dump, `scratch/dedd_reference.py` in the cice-dedd
worktrees) computes the full column bit-for-bit; the final tests pin the
dump directly.
