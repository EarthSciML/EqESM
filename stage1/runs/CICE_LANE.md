# CICE lane: delta-Eddington shortwave (dEdd)

Scope: the CICE sea-ice shortwave radiation physics (`shortwave = 'dEdd'`,
Briegleb & Light 2007 / Holland et al. 2012 three-band scheme as implemented
in `components/cice/src/icepack/columnphysics/icepack_shortwave.F90`) is to
be re-implemented as `stage2/cice_shortwave.esm` in the EqESM model. This
lane delivers (1) the instrumented D-class boundary dump
(`stage1/dumps/cice/cice_dedd.csv`, one row per (cell, category, step) per
spectral band), (2) the identity-law validation script, (3) the stage-2 .esm
component with tests pinned to the dump at rel 1e-12.

Sibling lanes' patterns: CISM (`stage1/runs/CISM_LANE.md`), MOSART, CTSM.

## Ownership / collision guards

- **Ours:** `stage1/runs/CICE_LANE.md`, `stage1/runs/cice_instrument.json`,
  fork `stage1/EqESM-cesm-cice` (branch `instrument/cice-dumps`; cice
  wrapper submodule + `cice/src` submodule + icepack submodule on the same
  branch name), `stage1/dumps/cice/`, `stage1/analysis/*cice*`,
  stage-2 CICE work (`stage2/cice_*.esm`).
- **Coordinate:** `PROGRESS.md` (manager-maintained), the shared
  `stage1/inputdata`, the `ctessum` partition, `bugs.md`.

## 1. The physics: what the dEdd scheme computes

`icepack_step_radiation` (`icepack_shortwave.F90:1298`) calls
`compute_dEdd_3bd` (`:1703`) once per surface type (0 = bare ice, 1 = snow,
2 = pond; pond fraction from `shortwave_dEdd_set_eff`, `:3900`), each call
solving the delta-Eddington two-stream radiative transfer for a 12-layer
column (3 snow layers + snow SSL + 8 ice layers; `klev=12`, `klevp=13`,
`kii=4`) in three spectral bands (visible, near-IR 0.7–1.19 µm, near-IR
1.19–5.0 µm; `nspint_3bd=3`).

Per band the scheme:

1. builds the single-scattering properties per radiation layer — snow
   layers: `frsnw = (fnidr + 0.8·(1−fnidr))·rsnw` grain adjustment,
   `ks = Qs·(ρs/ρi)·3/(4·frsnw·10⁻⁶)`, `τ = 3·Qs·ρs·dz/(4·ρi·frsnw·10⁻⁶)`,
   `w0 = 1 − ks·dz/τ`, `g = gs`; ice SSL/DL/interior and ponded-ice IOPs
   from the tuned `ki_*_mn_3bd` tables, with the `R_ice`/`R_pnd` tuning
   factors (0 in this run → identity); bulk aerosol modification from the
   12 `aero_mp` mass-path entries (n_aero=3: species 1..3 =
   kaer/waer/gaer_3bd(:,1:3); per species na = 1,5,9 the snow SSL, snow
   layers, ice SSL, ice layers use aero_mp(na), aero_mp(na+1)·rnslyr,
   aero_mp(na+2), aero_mp(na+3)·rnilyr), then per layer
   `τ += taer`, `w0 = (w0·τ + waer)/(τ + taer)`,
   `g = (g·w0·τ + gaer)/(w0·τ + waer)`;
2. delta-scales and combines layers down/up (`solution_dEdd`, `:3022`):
   combined layer properties `rdir_l/rdif_a_l/rdif_b_l/tdir_l/tdif_a_l/
   tdif_b_l/trnlay_l`, then the 8-point Gauss quadrature angular integral
   (8 `gauspt/gauswt` points), Fresnel interface at `kfrsnl = nslyr+2`,
   and the downward/upward interface fluxes `trndir/trntdr/trndif/rdndif/
   rupdir/rupdif` (`dfdir/dfdif = downward − upward` per interface);
3. accumulates the band outputs with the near-IR spectral weights
   `wghtns(1)=1`, `wghtns(2)=0.67+0.11·(1−fnidr)`, `wghtns(3)=1−wghtns(2)`:
   `avdr/avdf` = band-1 `rupdir(0)/rupdif(0)`; `aidr/aidf` +=
   `rupdir(0)/rupdif(0)·wghtns`; `fsfc/fint/fthru` +=
   `(flux(0)−flux(ksrf))`, `(flux(ksrf)−flux(13))`, `flux(13)` with `ksrf=1`
   (snow) / 5 (bare, pond); `Sabs` (snow layers only) and `Iabs` (ice
   layers, with the snow-covered top-layer km=4/kp=6 indexing); `fthrul(i)`
   = visible flux at interface 4+i (i = 1..9).

Layer/interface truncation: after the full recurrence, every
`trntdr(k) < trmin = 0.001` value (and the derived downward quantities below
it) is set to 0 — the guard only *skips recomputation*, so the unguarded
recurrence with a post-hoc active mask is bit-exact (the stage-2 model uses
that form, avoiding the EarthSciAST `observed_cycle` rejection of guarded
recurrences).

## 2. The reference configuration (what the D-class run does)

Compset **D** (`2000_DATM%NYF_SLND_SICE_SOCN_DROF%NYF_SGLC_SWAV`) on
**f19_g17** (a%1.9x2.5, oi%null, m%gx1v7), one day, 24 CICE timesteps of
3600 s, from `stage1/runs/cice_instrument.json` (job 12124251, the final
`cice_dedd.csv` layout; rerun history 11864874/12118982/12122398 in the
record). Namelist facts for the dEdd path (all verified from the case's
`ice_in`):

| Option | Value | Effect on the dEdd law |
|---|---|---|
| `shortwave` | `dEdd` | 3-band delta-Eddington; `compute_dEdd_3bd` |
| `tr_aero` | .true., `modal_aero` .false., `n_aero=3` | `aero_mp` 12 entries; species 1..3 = kaer_3bd(:,1:3) |
| `kalg` | 0 | `kabs_chl` = 0 (chlorophyll term off) |
| `R_ice`, `R_pnd` | 0 | IOP tuning identity → `ki_ssl..gi_p_int` = the `*_mn_3bd` tables (carried as inputs) |
| `tr_zaero`, `tr_bgc_N` | .false. | no z-aerosol/bgc shortwave terms |
| `tr_pond_sealvl` | .true. | `shortwave_dEdd_set_eff` pond path (hp-transition block skipped) |
| `frzpnd` | `hlid` | lid-ice ponds |
| `ktherm` | 2 (mushy) | zTin/zTsn real (0 in the dump are the `ktherm /= 2` placeholders) |
| `ncat=5, nilyr=8, nslyr=3, kii=4, kfrsnl=5` | — | layer geometry of the dump |

The dump (PET 0, first 32 tmask cells, first 12 steps) has 5760 rows = 1920
(cell, category, step) triples × 3 bands, 383 columns: state/inputs
(aicen/vicen/vsnon, Tsfcn, pond tracers, zTin/zTsn, fnidr, coszen, the four
incident fluxes, hs/hi/hp/fi/kalg, rsnw/rhosnw/frsnw per snow layer,
dzk_0..12, the 15 tuned IOPs, aero_mp_1..12) and the per-band solution
(tau/w0/g per radiation layer, rdir_l..trnlay_l, trndir..dfdif per
interface) plus the running outputs (final at band 3).

## 3. Validation (identity laws)

`stage1/analysis/check_cice_dedd.py` verifies 27 identities per row × 5760
rows = **155520/155520** at rtol/atol 1e-12:

- `fthrul_1..9` = visible flux `dfdir·swvdr + dfdif·swvdf` at interfaces
  5..13 (from the triple's band-1 row);
- `avdr/avdf` = band-1 `rupdir(0)/rupdif(0)`;
- `aidr/aidf` = running sums of `rupdir(0)/rupdif(0)·wghtns(ns)` over
  bands 2..3;
- `Sabs_1..3` (0 on bare rows), `Iabs_1..8` (snow-covered top layer km=4/
  kp=6), `fsfc/fint/fthru` = running band sums 1..ns(row) of the per-band
  flux differences, `wghtns` applied exactly as the Fortran forms them
  (one rounding per term), so the checks hold to round-off.

Composition notes: srftyp 1 (snow) 1884 triples, srftyp 0 (bare) 36
triples, no srftyp 2 rows (effective pond fraction 0 everywhere — see the
`row_counts` check); `fi` = coverage of the surface type, 1.0 everywhere;
bare rows carry zero-thickness snow layers (`dzk_0..3 = tau_0..3 = 0`);
`fthrul_9 = fthru` only where the near-IR bands transmit nothing to the ice
base (1662/1920 triples; elsewhere IR tails 1e-6..1e-2 W/m²).

## 4. Run record (2026-11-06)

Job **12124251** (final `cice_dedd.csv` layout), fork commits
**04eecbf** (icepack) / **5c879ee** (cice/src) / **8579a7c** (cice wrapper)
/ **a34b33d** (superproject) on branch `instrument/cice-dumps`: recorder
dummy lower bounds fixed per-array (0: for dzk/qs_k..ks_k/tau/w0/g/
rdir_l..trnlay_l/trndir..dfdif, 1: for rsnw/rhosnw/frsnw/ki_ssl..gi_p_int/
Sabs/Iabs/fthrul/aero_mp), fixing the gfortran assumed-shape
lbound-1 renormalization that shifted every layer/interface array in the
earlier dumps (reruns 12118982, 12122398 — see the record). The run is
**solution-neutral**: instrumented vs reference D_f19g17 iceh and restart
files bit-identical modulo title/timestamp attributes. All 14 CSVs pass
the record's checks (row counts, freezing-temp law, EVP stress/momentum
laws, no blank lines) and the dEdd dump passes 155520/155520 identity
checks (above).

## 5. Stage 2 plan

`stage2/cice_shortwave.esm` (model `ShortwaveDEdd`): inputs from the dump;
pointwise `dzk` law, `frsnw` law, `kabs_chl`; snow IOP interpolation
(`Qs/ws/gs` vs `rsnw_tab`, clamped linear) and the aerosol-modified
`tau/w0/g` per (band, srftyp, layer); the single-layer dEdd solution
(8-point Gauss quadrature via `faq` reduce over the 8 `gauspt/gauswt`
constants, Fresnel interface at `kfrsnl`); six coupled recurrences over the
14-interface axis (downward 4-variable pass, upward 2-variable pass) in the
unguarded + post-hoc `active`-mask form (bit-exact vs the Fortran, no
`observed_cycle`); `dfdir/dfdif`; band-accumulation recurrences over the
3-band axis for `avdr/avdf/aidr/aidf/fsfc/fint/fthru/fthru*/Sabs/Iabs/
fthrul`; srftyp aggregation with the fs/fp/fi masks. Tests pin the dump
rows (snow rows incl. band-3 truncation; the 36 bare rows; pond unpinnable
— no srftyp=2 rows, documented). New expression templates (layer solution
pieces, `wghtns`, IOP lookups) appended to `stage2/cice_templates.esm`;
`stage2/README-cice-dedd.md` documents the mapping to
`icepack_shortwave.F90`.

## 6. Risks

- **Pond branch unpinnable** from this run (no srftyp=2 rows) — the pond
  laws are pinned from the source formulas only; a later run with ponding
  (e.g. a summer start) would complete the tests.
- **Recurrence support in EarthSciAST** is new (v0.7.0): `faq`-based
  recurrences need 1-based index sets, an `ifelse(k <= 1, base, recur)`
  guard, and the self-variable excluded from `faq` args; the probe
  `scratch/probe_recur4.esm` documents the working form.
- **gfortran denormal floats** in the dumps lack the 'E' (e.g.
  `5.58...-322`) — parse with the DENORM regex in `check_cice_dedd.py`.
- The shared queue is contended with other lanes' runs; submit early,
  rerun late.

## 7. Ridge + ITD instrumentation (2026-11-09)

Runs 12277990 / 12280008 / 12280215 instrumented the mechanical
redistribution (ridging) and the thickness-category remap for stage 2 of
`stage2/cice_itd.esm`, on the same D_cice_inst case (fork branch
`instrument/cice-dumps`, superproject commit 90b4900; icepack d404867 =
thread-private recorder buffers, 90b4900 = da0_melt/da0_himin init;
cice 425c1db = dump rows, 9e504f1 pointer; wrapper 7cdffa2).

### Dumps (both land in stage1/dumps/cice/, 32 traced cells x 12 steps)

- `cice_ridge.csv` (156 cols): the ridge law
  (`icepack_mechred.F90` ridge_ice + ridge_shift, armed/drained in
  `ice_step_mod.F90` step_dyn_ridge): inputs rdg_conv/rdg_shear/Tf; the
  first-iteration ridge_prep rates closing_net/opning/divu_adv; ridge_itd
  outputs aksum/apartic_0..5/hrmin/hrmax/hrexp/krdg (krdg_partic=1,
  krdg_redist=1, exponential); per-category transfers ardg1n/ardg2n/
  virdgn + totals ardg1/ardg2/virdg; the n=1 redistribution fractions
  aredistn/vredistn; rafting araftn/vraftn; entry/exit state + energies;
  msnow_mlt/esnow_mlt/mpond; niter; hin_max_0..4. niter = 1 on all rows.
- `cice_itd.csv` (136 cols): the linear remap (`icepack_therm_itd.F90`
  linear_itd + shift_ice, armed/drained in step_therm2): pre-thermo and
  pre-ITD entry states; hicen_init/hicen/dhicen; hbnew_0..5; the
  per-category linear fits g0/g1/hL/hR + remap_flag (= 1 on all rows);
  boundary transfers donor/daice/dvice; da0_melt/da0_himin; post- and
  exit states + aice/aice0; the 12 conservation totals; linear_itd_ran.

### Bugs found and fixed in the instrumentation

- **OMP race (cice run 12277990)**: `step_dyn_ridge` and `step_therm2`
  run their iblk loops under `!$OMP PARALLEL DO`; an untraced block's
  `ridge_ice`/`icepack_step_therm2` could see the shared recorder
  `*_armed` flag while a traced cell sat between arm and drain and
  overwrite its row (11/384 ridge rows internally inconsistent, e.g.
  entry snapshot of one cell with the transfers and exit state of
  another). Fixed with `!$OMP THREADPRIVATE` on both recorders' module
  state (icepack d404867); the dEdd recorder is untouched (its dump is
  clean and must stay bit-identical).
- **Uninitialized da0_melt/da0_himin** (run 12280008): the dh0<0 melt
  and hi_min closures assign them only when they fire; uninitialized
  locals leaked denormal garbage into rows where they did not. Fixed by
  initializing both to 0 at the top of linear_itd (icepack 90b4900).

### Validation (final run 12280215)

- `stage1/analysis/check_cice_itd.py`: 55007/55007 identity checks pass
  (71 per row x 768 rows, rtol/atol 1e-12) - every recorded quantity
  recomputed from the others with the exact Fortran formulas: ridge
  participation/aksum/krdg/hrmin/hrexp, ridge_prep rates, per-category
  transfers + totals, aredistn/vredistn from the exponential ridge ITD
  integrals (sum farea = sum fvol = 1), rafting mask, energy/volume/
  snow conservation (msnow_mlt/rhos and esnow_mlt added back); ITD
  hbnew interpolation + top extrapolation, fit_line area/first-moment
  identities + h13/h23 bound shifts, boundary-transfer eta integrals
  with the puny/full-category clamps, post state (da0_melt, transfers,
  hi_min closure at 0.01 m), exit-state aggregate identities, and the
  conservation totals within linear_itd's own thresholds.
- `stage1/analysis/check_cice.py` now also delegates to
  check_cice_itd.py and counts both new CSVs (384 rows each, zero blank
  lines); the pre-existing 15 CSVs re-verified bit-identical in content
  (deterministic trajectory + unchanged dump code, and the checks pass).
- Solution neutrality: `compare_iceh.py` bit-identical vs the D_f19g17
  reference after every rerun (modulo case metadata).

### Notes for stage 2 (`stage2/cice_itd.esm`)

- Ridge rates are the first-iteration ridge_prep values (the loop
  zeroes closing_net/opning on convergence; niter = 1 everywhere in
  this run, so iteration accumulation never engaged).
- `step_dyn_ridge` does not pass `closing` to icepack_step_ridge, so
  ridge_ice takes the ridge_prep branch (rates from rdg_conv/rdg_shear).
- The dumped `hicen` is post-shift_ice (shift_ice recomputes it from
  vicen/aicen); `dhicen` is pre-shift; the fits used the pre-shift
  entry thicknesses.
- `da0_melt`/`da0_himin` fire on all 384 rows in this run (category-1
  melt everywhere; post-remap hicen_1 < 0.01 m everywhere).
- Coverage gaps: no niter > 1 rows (no rate-capped ridging), no
  remap_flag = 0 rows (no failed remap), no srftyp-related ITD variety;
  all rows have ice in every category.
