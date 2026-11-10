# CAM6 RRTMG radiation (cam_radiation.esm)

Stage-2 lane for the RRTMG radiation chain (`radiation_tend`, tphysbc
site 11, `radiation.F90`) as six `.esm` models pinned to the stage-1
instrumentation dumps of run 12255018
(`stage1/dumps/cam/cam_rad_input.csv`, `cam_rad_scal.csv`,
`cam_rad_out.csv`, `cam_rad_out_ifc.csv`, `cam_rad_cldopt_sw.csv`,
`cam_rad_cldopt_lw.csv`, `cam_rad_solver_sw.csv`,
`cam_rad_solver_lw.csv`; 23 radiation steps of the 48-step SCAM run —
every iradsw=iradlw=2 step, 13 day steps 3, 25, …, 47 and 10 night
steps 5–23).  Sources: the instrumented cam6_4_212 fork
(`stage1/EqESM-cesm-cam`, branch `instrument/cam-dumps`, cam submodule
`81922142`): `components/cam/src/physics/rrtmg/radiation.F90`
(`radiation_tend`, the cldfprime assembly), `radsw.F90`
(`rad_rrtmg_sw`; the E_cld_tau → tauc_sw/ssac_sw/asmc_sw/fsfc_sw
conversion :437-475, the McICA call :489-501, the flux/heating
reversal :545-575), `radlw.F90` (`rad_rrtmg_lw`), `rrtmg_state.F90`
(the layer/gas-vmr construction :60-247), the solver internals
`aer_src/mcica_subcol_gen_sw.f90` / `mcica_subcol_gen_lw.f90`, and the
band constants of `radconstants.F90:34-40`.  Boundary/column
documentation and the stage-2 decomposition plan:
`stage1/runs/CAM_LANE.md` (RRTMG dump section), `stage1/runs/
cam_instrument.json`, `stage2/README-cam-rad.md`; dump validation:
`stage1/analysis/validate_cam_rad_dump.py` (ALL RADIATION DUMP CHECKS
PASSED 2026-11-08, run 12255018, incl. the flux-divergence heating
identity verified to 3e-16 W/kg).

## The chain and the six models

```
state (t,q,pmid,pint) + cam_in (ts, albedos) + pbuf (cld, cldfsnow,
      ICLWP, ICIWP, DEI, LAMBDAC, MU, ICSWP, DES) + coszrs/eccf/sfac
  |
  | CamRadState     rrtmg_state_create (rrtmg_state.F90:60-247):
  |                 the RRTMG 33-layer gas state from the 32 CAM layers
  |                 + the extra top layer; pmid, the shared pint_hi/
  |                 T_ifc_hi edges, h2ovmr = (q/(1-q))*amdw,
  |                 cfc22vmr = ccl4vmr = 0
  v
solin, coszrs/eccf/sfac
  |
  | CamRadSolin     solin = sum_b sfac_b*solar_band_irrad_b * eccf * coszrs
  |                 (radsw.F90:395-397; 0 at night)
  v
cloud optics per (band, lev): tau, tau_w, tau_w_g, tau_w_f  (E_cld_tau,
radiation.F90:956-1210, not dumped; the ratios below are)
  |
  | CamRadSolverInputs  ssa = tau_w/tau (tau>0) else 1,
  |                 asm = tau_w_g/tau_w, fsfc = tau_w_f/tau_w
  |                 (tau_w>0) else 0 (radsw.F90:437-475),
  |                 cldfprime = max(cld, cldfsnow) (radiation.F90:1012-1016)
  v
tauc/ssac/asmc/fsfc per (band, lev)  +  cld_sto (RNG draw, not dumped)
  |
  | CamRadMcicaSw   tauc_sto = cld_sto*tau(band(gpt)),
  |                 ssac_sto = cld_sto*ssa + (1-cld_sto),
  |                 asmc_sto = cld_sto*asm, fsfc_sto = cld_sto*fsfc
  |                 (mcica_subcol_gen_sw.f90:439-444)
  | CamRadMcicaLw   tauc_sto_lw = cld_sto*tau_lw(band_lw(gpt_lw))
  |                 (mcica_subcol_gen_lw.f90:408-416)
  v
two-stream solvers (radsw.F90 / radlw.F90 internals, not dumped)
  |
  | CamRadHeating   fns = fds - fus, fnl = ful - fdl at the interfaces;
  |                 qrs = g*D(fns), qrl = -g*D(fnl) (the layer flux
  |                 divergence of radsw.F90:573-574 / radlw.F90:244-245);
  |                 the history scalars fsnt/fsds/fsns/flnt/flns/flut/
  |                 flwds as index picks (radsw.F90:552-562,
  |                 radlw.F90:214-228); sw_net/lw_net as the whole-column
  |                 integrals sum(qrs dp/g), sum(qrl dp/g) (the
  |                 telescoping identities fsnt-fsns, flns-flnt)
  v
qrs, qrl  (W/kg)  ->  heating rates for the state update
```

| Model | Dump pins | Equations |
|---|---|---|
| `CamRadState` | pmid at k=1, 2, 17, 33; pint_hi at k=1, 2, 17, 32; T_ifc_hi at k=17, 32; h2ovmr at k=1, 17, 33; cfc22vmr/ccl4vmr = 0 at k=17 | `pmid(1) = 0.5*pint_hi(1)`, `pmid(k) = 0.5*(pint_lo+pint_hi)` (rrtmg_state.F90:126,140), the shared edges `pint_hi(k) = pint_lo(k+1)` and `T_ifc_hi(k) = T_ifc_lo(k+1)` via the `edge_shift` template (binding n = NK−1, pinned at four of the k=1..32 identities), `h2ovmr = (q/(1-q))*amdw`, `cfc22vmr = ccl4vmr = 0` |
| `CamRadSolin` | solin on 4 day steps | `solin = sum_b sfac_b*solar_band_irrad_b * eccf * coszrs`, the band sum a `faq` reduce over the band axis (solar_band_irrad the 14 radconstants.F90:34-40 values, sum 1368.22 W/m²) |
| `CamRadSolverInputs` | ssa at the tau=0 else branch (ssa=1, esd lev 32 and 16 of the clear step) and ssa/asm/fsfc at the cloudy levels of steps 33/47, bands 1, 8, 14 | the radsw.F90:437-475 ratios with `guarded_ratio` |
| `CamRadMcicaSw` | tauc_sto/ssac_sto/asmc_sto/fsfc_sto at the first cloudy subcolumns (gpts 13, 47, 110 at esd lev 18; 9, 20, 47, 101 at esd lev 13) and the ssac_sto = 1 clear fill | the mcica_subcol_gen_sw.f90:439-444 laws with the `band_gather` + `mcica_subcolumn` templates, `cld_sto` an input (RNG draw) |
| `CamRadMcicaLw` | tauc_sto_lw at the first cloudy subcolumns of bands 4, 8, 14 (gpt_lw 46/89/135 and 42/91/136 at esd lev 13) | `tauc_sto_lw = cld_sto*tau_lw(band_lw(gpt_lw))` with the LW band map (rrtmg_lw_init.f90:308) |
| `CamRadHeating` | qrs/qrl at esd levs 1, 17, 32; fns/fnl at lev_nodes 1, 17, 33; the eight history scalars; sw_net/lw_net | the esd flux divergence `qrs = g*D(fns)`, `qrl = -g*D(fnl)` on `ze = -pint_hi` reversed; scalar index picks; whole-column `integral` of the named observeds `qrs_over_g`/`qrl_over_g` |

## Column and index conventions

- The esd column is **lev 1 = surface … lev NLEV = top**, interfaces
  lev_nodes 1 = surface … NLEV+1 = top, exactly as in
  `cam_vertical_diffusion.esm`; `ze = -pint_hi` reversed (Pa, increasing
  upward), esd interface n <-> CAM interface 34-n.
- The RRTMG state construction and all pointwise dump rows are written
  **CAM top-down**; the extractor reverses the level axes and maps
  assertion coordinates `l -> 33-l` (CAM lev 15 = esd lev 18, CAM lev
  20 = esd lev 13 in the cloudy pins).  The RRTMG layer index `k`
  (1..33, top-down, the extra top layer) is kept **unreversed** inside
  CamRadState — the extra layer is k=1 and the surface layer k=33, and
  the esd assertion coords for that model are the CAM k values.
- Scalars are index picks: `fsnt = fns[NLEV+1]`, `fsns = fns[1]`,
  `flnt = fnl[NLEV+1]`, `flns = fnl[1]`, `fsds = fds[1]` (surface
  downwelling), `flut = ful[NLEV+1]` (TOA upwelling), `flwds = fdl[1]`.
- Metaparameters NLEV=32, NK=33, NBND=14, NGPT=112, NBNDL=16,
  NGPTL=140; index sets lev/lev_nodes/k/band/gpt/band_lw/gpt_lw.

## Templates (cam_templates.esm, appended in this lane)

- `guarded_ratio` — the `tau > 0`-guarded ratios of radsw.F90:437-475
  (also used by cam_vertical_diffusion); bound as `a/(a+b)` with the
  else-branch value.
- `edge_shift` — the shared-edge identity `u_hi(k) = u_lo(k+1)`;
  bound here with n = NK−1 for pint_hi/T_ifc_hi.
- `band_gather` — the k-distribution band map gpt -> band (SW) and
  gpt_lw -> band_lw (LW), the rrtmg init maps.
- `mcica_subcolumn` — the McICA masking of one (gpt, lev) column from
  `cld_sto` and the gathered band optics, `cld_sto` as an input
  (RNG draw).
- The solin band sum uses a `faq` with `"reduce": "+"` and empty
  `output_idx` (the whole-band-axis sum).

## Pinned tests (19 tests, 198 assertions, all rel 1e-12, `time_span {0,0}`)

| model | tests | assertions | character |
|---|---|---|---|
| CamRadState | step3, step41, step47 | 15 each | clear / dense cloud (cld ~ 1.0) / dense cloud, day |
| CamRadSolin | step3, step25, step37, step47 | 1 each | the 4 day steps (coszrs 0.1308, 0.1495, 0.9688, 0.5101; eccf ~ 0.9685; sfac 14 values ~ 0.917-1.027, the same every day step) |
| CamRadSolverInputs | step3, step33, step47 | 4 / 9 / 9 | tau = 0 else-branch (ssa = 1, asm = fsfc = 0); cloudy levels CAM 15 / CAM 20, bands 1, 8, 14 |
| CamRadMcicaSw | step37, step47 | 13 / 17 | thin cloud (cld ~ 0.12) / dense cloud; first cloudy subcolumns of bands 2, 6, 14 / 2, 3, 6, 14 at the cloudy level + the cld_sto = 0 fill branch (ssac_sto = 1) |
| CamRadMcicaLw | step41, step47 | 3 each | first cloudy subcolumns of bands 4, 8, 14 at the cloudy level |
| CamRadHeating | step3, step25, step37, step47, step21 | 19 / 21 / 21 / 19 / 11 | day clear / day clear / thin cloud / dense cloud / night |

- Every model's first test runs on the **step-3 parameter defaults**
  (`parameter_overrides: {}`); the other tests override the dump
  inputs.  Each assertion is guaranteed to fail on a wrong law — the
  extractor copy-checks every tuple against the input rows and drops
  exact input copies (nothing in the tests block is a pass-through of
  an overridden input).
- Step 3 and step 47 have sfc albedo = 1, so `fsns = fns[1] = 0`
  exactly and the fsns assertion is skipped there; the step-21 (night)
  test pins the longwave heating only — the night SW side is the
  all-zero boundary of the radsw early return and is never pinned.
- Units: h2ovmr is declared mol/mol with the dimensionless amdw
  conversion (the units checker treats kg/kg as dimensionless, so the
  q/(1-q) ratio is unit-free); ze in Pa (the CAM pint), dz in Pa, and
  the D operator is the flux divergence dF/dp, so qrs/qrl come out in
  W/kg with g in m/s².

## Known gaps (out of scope of this slice)

- **k-distribution and two-stream solvers**: the RRTMG gpt absorption
  coefficients, the delta-scaled two-stream solution and the
  fds/fus/fdl/ful up/down split are not dumped; CamRadHeating takes
  the interface fluxes (fds/fus/fdl/ful, cam_rad_out_ifc.csv) as
  inputs and reproduces the net-flux layer heating, and the McICA
  models reproduce the subcolumn masking from the band optics.  The
  full SW/LW radiative transfer itself is a later stage.
- **o3vmr top-layer decay**: the P_top-dependent o3vmr extrapolation of
  rrtmg_state.F90 is not pinnable from the dump (only the dumped
  levels appear; the extrapolated layer value is an input to the
  dump's state rows), so CamRadState pins pmid/edges/h2ovmr/CFC laws
  only.
- **cldfprime branch**: `cldfprime = max(cld, cldfsnow)`
  (radiation.F90:1012-1016) is implemented but the dump's cldfsnow = 0
  everywhere, so only the cldfsnow <= cld branch is exercised; the
  cldfsnow > cld branch is untested.
- **EarthSciModels `rrtm_lw` is NOT reused**: the existing
  EarthSciModels RRTM implementation descends from the WRF fork whose
  band structure and constants differ from CAM6's rrtmg_lw (different
  band counts and wavenumber ranges); reusing it would silently
  pin the wrong physics, so this lane is written against the CAM6
  sources directly.
- **asm/fsfc else branch**: the tau = 0 -> ssa = 1 branch is pinned on
  the clear step 3, but the tau_w = 0 -> asm = fsfc = 0 else branch is
  not pinned anywhere (on the cloudy steps the guarded ratio takes the
  main branch; on step 3 only ssa is asserted at the tau = 0 levels).

## Reproducing the test tuples

`stage1/analysis/extract_cam_rad_pins.py` reads the eight
`cam_rad_*.csv` dumps (run 12255018) and emits the per-model test
tuples (defaults = step-3 inputs, assertion coordinates already mapped
to the esd axes, copy-checked) as JSON; `--render` prints the
ready-to-paste `tests` entries.  Scripts paste only numeric tuples;
the model structure, equations and descriptions are hand-authored.

## Testing

    ./esm validate stage2/cam_radiation.esm
    ./esm units stage2/cam_radiation.esm
    ./esm test stage2/cam_radiation.esm     # 198 assertions

Whole-suite stage-2 baseline: `./esm test stage2/` = 3271P/0F/0E
(3073 pre-existing + 198 from this lane).
