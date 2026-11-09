# CAM6 RRTMG radiation — proposed stage-2 `.esm` decomposition

Stage-2 plan for the RRTMG radiation chain (`radiation_tend`, tphysbc
site 11, `radiation.F90`), to be implemented as `.esm` components
pinned by the stage-1 radiation dumps `cam_rad_input.csv`,
`cam_rad_scal.csv`, `cam_rad_out.csv`, `cam_rad_out_ifc.csv`,
`cam_rad_cldopt_sw.csv`, `cam_rad_cldopt_lw.csv`,
`cam_rad_solver_sw.csv`, `cam_rad_solver_lw.csv` (23 radiation steps of
the 48-step SCAM run: 3, 5, …, 47 — every iradsw=iradlw=2 steps; 13
day steps 3, 25, …, 47 and 10 night steps 5–23).  Sources: the
instrumented cam6_4_212 fork (`stage1/EqESM-cesm-cam`, branch
`instrument/cam-dumps`, cam submodule `81922142`):
`components/cam/src/physics/rrtmg/radiation.F90` (`radiation_tend`,
`rrtmg_state_create`, the cldfprime/cloud-optics assembly at
:956-1210), `radsw.F90` (`rad_rrtmg_sw`; the E_cld_tau → tauc_sw/
ssac_sw/asmc_sw/fsfc_sw conversion :437-475, the McICA call :489-501,
the flux/heating reversal :545-575), `radlw.F90` (`rad_rrtmg_lw`),
`rrtmg_state.F90` (the layer/gas-vmr construction :60-247), and the
solver internals `aer_src/mcica_subcol_gen_sw.f90` /
`mcica_subcol_gen_lw.f90`.  Boundary/column documentation:
`stage1/runs/CAM_LANE.md` (RRTMG dump section), `stage1/runs/
cam_instrument.json`; dump validation: `stage1/analysis/
validate_cam_rad_dump.py` (ALL RADIATION DUMP CHECKS PASSED
2026-11-08, run 12255018).

## The chain and its data flow

```
state (t,q,pmid,pint) + cam_in (ts, albedos) + pbuf (cld, cldfsnow,
      ICLWP, ICIWP, DEI, LAMBDAC, MU, ICSWP, DES) + coszrs/eccf/sfac
  |  rrtmg_state_create (rrtmg_state.F90)   -> pmidmb/pintmb/tlay/tlev,
  |     gas vmrs h2ovmr..ccl4vmr (P_top=50 Pa top-layer o3 decay)
  |  cldfprime = max(cld, cldfsnow)         (radiation.F90:1012-1016)
  |  cloud_rad_props (gammadist/mitchell, cloud_rad_props.F90)
  |     -> E_cld_tau, E_cld_tau_w, E_cld_tau_w_g, E_cld_tau_w_f (SW),
  |        c_cld_lw_abs (LW)                [cam_rad_cldopt_*.csv]
  |  rad_rrtmg_sw (radsw.F90):
  |     solin = sum(sfac*solar_band_irrad)*eccf*coszrs     [:395-397]
  |     tauc_sw = E_cld_tau; ssa = tau_w/tau (tau>0 else 1);
  |     asm = tau_w_g/tau_w, fsfc = tau_w_f/tau_w (tau_w>0 else 0)
  |                                                  [:437-475]
  |     mcica_subcol_sw (permuteseed=1) -> cld_sto, tauc_sto,
  |        ssac_sto, asmc_sto, fsfc_sto     [cam_rad_solver_sw.csv]
  |     rrtmg_sw_rad (k-distribution gas optics + two-stream + McICA)
  |        -> swdflx/swuflx, swhr (bottom-up)   [NOT dumped]
  |  rad_rrtmg_lw (radlw.F90):
  |     tauc_lw = c_cld_lw_abs              [cam_rad_cldopt_lw.csv]
  |     mcica_subcol_lw (permuteseed=150) -> cld_sto, tauc_sto
  |                                          [cam_rad_solver_lw.csv]
  |     rrtmg_lw_rad -> dflx/uflx, hr (bottom-up)   [NOT dumped]
  |  reversal to CAM top-down + K/d -> W/kg  (radsw.F90:545-575,
  |     radlw.F90:207-245)                  [cam_rad_out*.csv]
  v
qrs/qrl pbuf (x pdel in radheat_tend, radiation.F90:1427-1433) +
fsnt..solld history fluxes [cam_rad_scal.csv]
```

`cam_rad_input.csv` pins the solver boundary that both radsw and radlw
consume (one row per (step, RRTMG layer k = 1..num_rrtmg_levs = 33);
k = 1 is the extra layer above model top with pint_lo = 0.01 Pa,
cam_lev = max(k−1,1), and the top-layer o3vmr includes the P_top=50 Pa
linear-decay column, rrtmg_state.F90:228-247).  `cam_rad_cldopt_*.csv`
pin the CAM cloud-optics output per band (14 SW / 16 LW), with the
derived ssa/asm/fsfc that the solver actually receives;
`cam_rad_solver_*.csv` pin the McICA output per g-point (112 SW /
140 LW); `cam_rad_out*.csv` pin the solver verdict (heating rates per
level, fluxes per interface, and the history scalars).

## Proposed components (each its own `.esm`, tests from the dumps)

| component | Fortran source | pinned by | internal identities (validator checks) |
|---|---|---|---|
| `cam_rad_state` | rrtmg_state_create (rrtmg_state.F90:60-247) | pmid/pint/T/T_ifc, h2ovmr..ccl4vmr per layer | h2ovmr = (q/(1−q))·amdw, amdw = 1.607793; pmid(k) = ½(pint_lo+pint_hi) k≥2, pmid(1) = ½·pint_hi(1); pint_hi(k) = pint_lo(k+1); T_ifc_hi(k) = T_ifc_lo(k+1); cfc22vmr = ccl4vmr = 0 |
| `cam_rad_cldfprime` | radiation.F90:956-976 | cld, cldfsnow → cldfprime | cldfprime = max(cld, cldfsnow) |
| `cam_rad_cldopt` | cloud_rad_props.F90:123-131, 248-460 (gammadist + mitchell) | iclwp/iciwp/dei/lamc/mu/icswp/des → tau, tau_w, tau_w_g, tau_w_f (SW) and tau_lw (LW) | ssa = tau_w/tau (tau>0) else 1; asm = tau_w_g/tau_w, fsfc = tau_w_f/tau_w (tau_w>0) else 0; 0 ≤ tau_w ≤ tau; tau = 0 when cld = 0 |
| `cam_rad_mcica` | mcica_subcol_gen_sw.f90 / mcica_subcol_gen_lw.f90 (deterministic: permuteseed 1 / 150, single column) | cld, band optics → cld_sto, tauc_sto, ssac_sto, asmc_sto, fsfc_sto | cld_sto ∈ {0,1}; tauc_sto = tau(band(gpt)) when cld_sto=1 else 0; ssac_sto = ssa(band) when cld_sto=1 else 1 (day); asm/fsfc analog; band→gpt mappings from the k-distribution (SW 6,12,8,8,10,10,2,10,8,6,6,8,6,12; LW 10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2) |
| `cam_rad_solver_sw` | rrtmg_sw_rad.f90 (two-stream + McICA) — or EarthSciModels' existing RRTMG SW | cam_rad_out*.csv + cam_rad_scal.csv | fns = fds−fus; qrs = −g·(fns(k+1)−fns(k))/dp (max resid 2.6e-16 W/kg); solin identity; fsnt = fns(1), fsns = fns(33), fsds = fds(33), sols+soll+solsd+solld = fds(33) |
| `cam_rad_solver_lw` | rrtmg_lw_rad.f90 — candidate: reuse of the RRTMG LW in EarthSciModels atmospheric_radiation (EarthSciModels-rrtmlw worktree exists) | cam_rad_out*.csv + cam_rad_scal.csv | fnl = ful−fdl; qrl = +g·(fnl(k+1)−fnl(k))/dp (max resid 4.3e-16); flnt = fnl(1), flns = fnl(33), flut = ful(1), flwds = fdl(33) |
| `cam_rad_driver` | radiation_tend skeleton (radiation.F90) | the full chain: state → qrs/qrl + fluxes | column integrals Σ qrs·dp/g = fsnt−fsns, Σ qrl·dp/g = flns−flnt; fds(1) ≤ solin (extra-layer absorption 2.0..5.3 W/m2); night steps: all-zero SW boundary |

## Unit and scaling conventions that must be reproduced exactly

- All dumped radiation arrays are CAM **top-down**; RRTMG layers
  bottom-up internally (rrtmg_sw_rad.f90 reverses; radsw.F90:566-568
  reverses back).  RRTMG layer k (1..33, top down) holds CAM layer
  max(k−1,1)'s state; CAM interface i has pressure pintmb(i+1).
- Pressures in the dump are **Pa** (r_state%pmidmb/pintmb are mb,
  ×100); T in K; q in kg/kg; gas vmrs in mol/mol; fluxes in W/m2;
  heating rates in **W/kg** (= swhr[K/d]·cpair/86400, radsw.F90:573-574
  with dps = 1/86400).
- The flux-divergence identities: qrs(k) = −g·(fns(k+1)−fns(k))/
  (pintmb(k+2)−pintmb(k+1)), qrl(k) = +g·(fnl(k+1)−fnl(k))/dp with
  g = 9.80616 — the RRTMG heating is exactly the layer flux divergence
  (verified to ~3e-16 W/kg).
- ssa/asm/fsfc use the radsw.F90:437-475 branch structure (tau = 0 →
  ssa = 1, asm = fsfc = 0; tau_w = 0 → asm = fsfc = 0) with the
  `max(x, 1e-80)` guards.
- McICA maps each band's optics to its g-points through the fixed
  k-distribution band→gpt count lists above; the stochastic cld_sto
  profile is generated from the permuted seed sequence (reproducible
  for a single column, permuteseed 1 SW / 150 LW).
- Night steps (coszrs ≤ 0, 10 of 23) have the all-zero SW boundary
  (radsw.F90 early return; the dump flush zero-fills the SW buffers);
  LW always runs.  SW tests must branch on day/night.

## Test plan (mirroring the validator identities)

All tests `time_span {0,0}`, rel 1e-12 (as in the other CAM stage-2
files), citing dump file, step, k/lev, band or gpt:

1. **state layer laws** (pins `cam_rad_state`): h2ovmr = (q/(1−q))·amdw;
   pmid/pint/T_ifc edge sharing and the extra-layer pmid(1) =
   ½·pint_hi(1); cfc22vmr = ccl4vmr = 0 (759 rows).
2. **cloud optics** (pins `cam_rad_cldopt`): ssa/asm/fsfc from
   tau/tau_w/tau_w_g/tau_w_f with the branch structure above
   (10304 SW + 11776 LW rows); tau = 0 where cld = 0.
3. **McICA** (pins `cam_rad_mcica`): cld_sto ∈ {0,1}; tauc_sto =
   cld_sto·tau(band(gpt)); ssac_sto = 1 where cld_sto = 0 (day steps);
   the band→gpt mappings (82432 SW + 103040 LW rows).
4. **solver verdict** (pins `cam_rad_solver_*`): qrs/qrl vs the flux
   divergence; fns = fds−fus, fnl = ful−fdl; the scal↔ifc identities
   (fsnt = fns(1) etc.); solin identity on day steps; column integrals.
5. **night branch**: all-zero SW boundary on the 10 night steps.
6. **k-distribution gap**: the RRTMG gas-optics k-tables are not
   dumped (they are fixed tabulated data, not per-step state); pin the
   solver with the per-g-point boundary and the tabulated
   k-distributions from the upstream RRTMG sources, or reuse the
   EarthSciModels RRTMG SW/LW implementation and pin it with the
   solver verdict tests.

Data volume: 23 steps × (33 + 1 + 32 + 33 + 448 + 512 + 3584 + 4480)
rows per identity; solin spans 172..1277 W/m2 on day steps, coszrs
−0.50..0.97; clouds present on most steps (fsnt drops to ~176 W/m2 on
the cloudy steps 39-47 while fsds stays ~1000 W/m2).

## Coverage gaps / notes

- **The k-distribution gas optics are not dumped** (per-g-point solver
  boundary dumped instead, as agreed in the lane plan): the solver
  input per g-point is pinned, but the gas-absorption optical depths
  inside rrtmg_sw_rad/rrtmg_lw_rad are not directly observable.  The
  band→gpt mappings and the fixed k-tables come from the upstream
  RRTMG sources (aer_src/rrsw_kg*.f90, rrlw_kg*.f90).
- The clear-sky up/down interface fluxes (fusc/fdsc/fsul/fsdl) are not
  dumped; the clear-sky nets fcns/fcnl are pinned at ifc 1 and 33 via
  cam_rad_scal (fsntc/fcns(1), fsnsc/fcns(33), flntc/fcnl(1),
  flnsc/fcnl(33)).
- qrs/qrl are the pre-pdel W/kg values; the pbuf QRS/QRL consumed by
  radheat_tend are qrs·pdel (radiation.F90:1427-1433) — the .esm
  driver must apply the pdel scaling.
- fsn200/fsn200c/fln200/fln200c (fluxes at 200 mb) are dumped but not
  independently verifiable from the other columns (they interpolate
  the interface fluxes to 200 mb, radiation.F90:1271-1272, 1325-1326);
  they pin the interpolation in `cam_rad_driver`.
- fsntoa (net at TOA, edge 34 = 0.01 Pa) differs from fns(1) (CAM
  top) by the extra-layer absorption, 2.0..5.3 W/m2 (the O3 P_top
  column) — do not assert fsntoa = fns(1).
- The SCAM greenhouse gas vmrs (co2/ch4/n2o/cfc) are prescribed by the
  boundary data in this lane (scamMod.F90 "using column value of
  co2vmr from boundary data"), so `cam_rad_state` tests should use the
  dumped vmrs, not CAM-default values.
- cfc22vmr/ccl4vmr are 0 in this configuration (rrtmg_state.F90:222-223)
  — untested by data.
- 2026-11-08 queue note: the case batch script's `--exclusive` was
  removed from `.case.run` for runs 12254141/12255018 (partition
  saturated) and restored afterwards; see CAM_LANE.md.
