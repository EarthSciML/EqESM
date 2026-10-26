# Stage 1 — CESM Subassembly Boundary Catalog

This is the deliverable of PLAN.md §4.2 (boundary catalog): it maps the §3
component inventory to concrete CESM modules, entry routines, and
input/output boundaries, so that Stage-1 instrumentation (§4.3) and Stage-2
`.esm` stubs have an authoritative target list. It follows the EqIce pattern
(`eqice/stage1/boundaries.md`).

Reference implementation: CESM **cesm3.0-alphabranch** (checkout commit
`e59f2512`, grafted), submodule tags per PLAN.md §3 (cmeps1.1.73,
cam6_4_212, ctsm5.4.063, cesm3_cice6_6_3_22, mi_260909, mosart1.1.13,
cismwrap_3_00_002, cdeps1.0.112). All paths below are relative to
`stage1/cesm/components/`.

Legend: **D** = instantaneous/derivative computation (pointwise function of
current state — preferred for tests); **I** = integrated trajectory (advances
state over `dt`, or an accumulation/average over a coupling interval).

## 0. Coupling architecture (CMEPS mediator)

The mediator is `cmeps/mediator/` (NUOPC Mediator). It is a pure data-mediation
component: it maps, merges, accumulates/averages, and (for a few fields)
computes fluxes between the component caps; no component state evolves inside
it. Coupling mode is selected by the `coupling_mode` attribute
(`med.F90:840`): `'cesm'` uses `esmFldsExchange_cesm_mod.F90` (this catalog;
the `ufs*` and `hafs` modes are out of scope).

### 0.1 Run phases (`med.F90` SetServices; `cime_config/runseq/runseq_general.py`)

The generated `nuopc.runseq` nests time loops — glc (outermost, coarsest),
then rof, then ocn, then atm (innermost, shortest) — and interleaves mediator
phases with component runs and connectors (`MED -> COMP :remapMethod=redist`,
`COMP -> MED`). Per coupling step the active phases are:

| Phase (file) | Function | Class |
|---|---|---|
| `med_phases_prep_ocn_accum` (`med_phases_prep_ocn_mod.F90`) | accumulate atm/ice/rof fields toward ocn over the ocn coupling interval | **I** |
| `med_phases_prep_ocn_avg` | time-average the accumulated bundle, map to ocn grid, merge | **D** (average of an integral) |
| `med_phases_prep_ocn_enthalpy_*` | runoff/ice enthalpy corrections (Foxx_hrof*, Foxx_hcond …) | **D** |
| `med_phases_ocnalb_run` (`med_phases_ocnalb_mod.F90`) | **ocean albedo** (Sx_avsdr/avsdf/anidr/anidf), orbital params | **D** |
| `med_phases_aofluxes_run` (`med_phases_aofluxes_mod.F90`) | **atm/ocn surface fluxes** (see §0.4) | **D** |
| `med_fraction_set` (`med_fraction_mod.F90`) | update ocean/ice/land fractions (FBFrac: lfrac, ofrac, ifrac) | **D** |
| `med_phases_prep_{lnd,ice,rof,atm,glc,wav}[_accum\|_avg]` | map+merge fields from the relevant sources to each component's grid | **D** |
| `med_phases_post_{atm,ocn,ice,lnd,rof,glc,wav}` | map each component's export to other components' grids; `post_lnd` also calls `med_phases_prep_rof_accum` | **D** |
| `med_phases_cdeps_run` (`med_phases_cdeps_mod.F90`) | inline CDEPS data-component handling | **D** |
| `med_phases_diag_*`, `med_phases_history_write`, `med_phases_restart_write`, `med_phases_profile` | diagnostics/I-O | — |

Run-sequence ordering (atm-inner loop, `runseq_general.py:99–225`):
`aofluxes_run` → `prep_ocn_accum` → `ocnalb_run` → (ocn loop) `prep_ocn_avg`
→ `prep_lnd` → `prep_ice` → `prep_wav` → `prep_rof` → `prep_ocn_avg` →
`ICE` → `LND` → `ROF` → `WAV` → `OCN` → `post_ocn` → `aofluxes_run` →
`LND -> MED` → `post_lnd` (+`prep_rof_accum`) → `ICE -> MED` → `post_ice` →
`ROF -> MED` → `post_rof` → `OCN -> MED` → `post_ocn` → `GLC -> MED` →
`prep_glc` → `GLC` → `post_glc` → `ATM -> MED` → `post_atm` → `prep_atm` →
`ATM`.

### 0.2 Coupling intervals

Per-component coupling periods are set in `cmeps/cime_config/buildnml:180`
from the CIME XML variables: `NCPL_BASE_PERIOD` (`hour`/`day`/`year`/`decade`)
and `<COMP>_NCPL` for each component class; the resulting
`<comp>_cpl_dt = basedt/NCPL` is written into `drv_in`. The mediator clock
inherits the driver clock time step (`med.F90:2315 SetRunClock`), i.e. the
shortest coupling interval; the run sequence generator enforces
`glc_cpl_time >= rof,ocn,atm_cpl_time`. Typical CESM values: `ATM_NCPL=48`
(30 min), `OCN_NCPL=1` (daily), `ROF_NCPL=8` (3-hourly), `GLC_NCPL=1`.

### 0.3 Field bundles and the exchange table (`esmFldsExchange_cesm_mod.F90`)

Fields are registered per (source, destination) pair with
`addfld_from`/`addfld_to` + mapping (`addmap_from`: `mapbilnr`, `mapconsf`,
`mapconsd`, `mappatch`, `mapfcopy`, `mapnstod`; `idmap` when two components
share a mesh) and merging (`addmrg_to`: merge types `copy`, `merge`, `sum`,
fraction-weighted via `mrg_fracname_*`). The `med -> comp` bundles
(`fldListTo(comp)`, assembled in `med_phases_prep_*` by
`med_merge_auto_*` in `med_merge_mod.F90`) and `comp -> med` bundles
(`fldListFr(comp)`) are:

| Bundle | # | Field families (names as in code) |
|---|---|---|
| atm → med | 36 | state `Sa_z,u,v,tbot,ptem,shum,pbot,dens,pslv,o3,topo,u10m,v10m`; radiation `Faxa_swdn,swvdr,swvdf,swndr,swndf,swnet,lwdn`; precip `Faxa_rainc,rainl,snowc,snowl`; aerosols/chem `Faxa_bcph(3),ocph(3),dstwet(4),dstdry(4),ndep(2)`, `Sa_co2prog,co2diag,lightning` |
| med → atm | 44 | fractions `Sl_lfrac,Si_ifrac,So_ofrac`; surface state `Sx_tref,qref,t,So_t,So_ssq,So_re,So_ustar,Sx_u10,So_ugustOut,So_u10withGust,Sl_fv,ram1,snowh(lnd+ice),ddvel,soilw,fztop`; albedo `Sx_avsdr,avsdf,anidr,anidf`; fluxes `Faxx_taux,tauy,lat,sen,lwup,evap`, `Fall_flxdst(4),fco2_lnd,voc,fire`, `Faoo_fco2_ocn` |
| lnd → med | 35 | fluxes `Fall_taux,tauy,lat,sen,lwup,evap,swnet,flxdst(4),fco2_lnd,voc,fire`; state `Sl_t,tref,qref,avsdr,anidr,avsdf,anidf,snowh,u10,fv,ram1,lfrin,ddvel,soilw,fztop`; runoff `Flrl_rofsur,rofsub,rofgwl,rofi,irrig,rofsur_nonh2o`; glc-elev `Sl_tsrf_elev,topo_elev`, `Flgl_qice_elev` |
| med → lnd | 36 | atm state/rad/precip as in atm→med minus swnet; `Flrr_flood,volr,volrmch`, `Sr_tdepth,tdepth_max` (rof); `Sg_icemask,icemask_coupled_fluxes,ice_covered_elev,topo_elev` (glc); `Flgg_hflx_elev` |
| ice → med | 39 | state `Si_imask,ifrac,ifrac_n,t,vice,vsno,tref,qref,snowh,u10,avsdr,avsdf,anidr,anidf,thick,floediam,bpress`; atm/ice fluxes `Faii_taux,tauy,lat,sen,lwup,evap,swnet`; ice/ocn fluxes `Fioi_melth,swpen(+vdr,vdf,idr,idf,ifrac_n),meltw,salt,taux,tauy,bcpho,bcphi,flxdst` |
| med → ice | 29 | atm state `Sa_z,u,v,shum,tbot,pbot,ptem,dens`; rad/precip `Faxa_swvdr,swvdf,swndr,swndf,lwdn,rain,snow,bcph(3),dstwet(4),dstdry(4)`; ocn state `So_dhdx,dhdy,t,s,u,v,roce`, `Fioo_q,frazil`; waves `Sw_elevation_spectrum` |
| ocn → med | 18 | state `So_omask,t,s,u,v,dhdx,dhdy,bldepth,s_depth,t_depth,roce`; fluxes `Fioo_q,frazil`, `Faoo_fco2_ocn,fbrf_ocn,fdms_ocn,fn2o_ocn,fnh3_ocn` |
| med → ocn | 80 | atm state/rad/precip (subset); merged fluxes `Foxx_taux,tauy,sen,lat,lwup,lwnet,evap,swnet(+vdr,vdf,idr,idf,afracr),rofl,rofi,hrain,hsnow,hevap,hcond,hrofl,hrofi(+glc)`; `Forr_rofl,rofi(+glc),rofl_nonh2o`; ice `Fioi_*`; fracs `Sf_afrac,afracr,Si_ifrac,ifrac_n,bpress`; waves `Sw_*` (~30 Stokes/spectral fields); `Sa_pslv,co2`, `So_duu10n` |
| rof → med | 10 | `Forr_rofl,rofi,rofl_glc,rofi_glc,rofl_nonh2o`, `Flrr_flood,volr,volrmch`, `Sr_tdepth,tdepth_max` |
| med → rof | 8 | `Flrl_rofsur,rofsub,rofgwl,rofi,irrig`, `Fgrg_rofl,rofi`, (+scalar) |
| glc → med | 8 | `Sg_area,ice_covered,topo,icemask,icemask_coupled_fluxes`, `Flgg_hflx`, `Fgrg_rofl,rofi` |
| med → glc | 4 | `Sl_tsrf`, `Flgl_qice`, `So_t_depth,s_depth` (ocn2glc coupling) |
| wav ↔ med | 14 / 53 | `Fwxx_taux,tauy`, `Sa_*`, `Si_*`, `So_*` → wav; `Sw_*` spectra/statistics → med (ww3; out of scope) |

Field-name prefixes encode the pairs: `S` = state, `F` = flux;
`a`=atm, `o`=ocn, `i`=ice, `l`=lnd, `r`=rof, `g`=glc, `x`=exchange/average
surface; e.g. `Faxa_*` atm→x, `Fall_*` lnd→atm, `Forr_*` rof→ocn, `Sx_*`
mediator-computed surface state.

### 0.4 Flux computations inside the mediator

- **Atm/ocn fluxes** — `med_phases_aofluxes_mod.F90: med_aofluxes_update`
  maps input to the `aoflux_grid` (`'ogrid'` default; `'agrid'`/`'xgrid'`
  options) and calls `flux_atmocn_driver` (`cmeps/cesm/flux_atmocn/`), which
  dispatches to `flux_atmOcn_Large` (default Large & Pond),
  `flux_atmOcn_COARE`, or `flux_atmOcn_UA`, plus
  `flux_atmocn_Diurnal` (diurnal-cycle correction) and `shr_flux_mod`
  constants. Outputs: `Faxx_taux,tauy,sen,lat,lwup,evap`, diagnostics
  `Sx_tref,qref,So_ustar,re,ssq,duu10n,ugustOut,u10res`. All **D**.
- **Ocean albedo** — `med_phases_ocnalb_mod.F90: med_phases_ocnalb_run`
  (with orbital updates): `Sx_avsdr/avsdf/anidr/anidf`. **D**.
- **Sea-ice fluxes** — NOT in the mediator: CICE computes `Faii_*`
  (atm/ice) and `Fioi_*` (ice/ocn) internally (§1, CICE).
- **Enthalpy/heat-content corrections** — `med_phases_prep_ocn_enthalpy_*`
  (runoff `Foxx_hrof*`, `Foxx_hcond`, `Foxx_hrain/hsnow/hevap`) and
  `med_phases_prep_atm_enthalpy_*`. **D**.
- **Fraction updates** — `med_fraction_mod.F90: med_fraction_set`:
  `So_ofrac`, `Si_ifrac`, `Sl_lfrac` used for fraction-weighted merging. **D**.
- **Runoff mediation** — `med_phases_post_rof` creates the rof→lnd/ocn/ice
  bundles (`Forr_*`, `Flrr_*`, `Foxx_rofl/rofi`), removes negative runoff
  (`med_phases_post_rof_remove_negative_runoff`); `med_phases_prep_rof`
  accumulates/merges lnd (+glc) runoff to rof. **D** (accumulation **I**).

## 1. Component import/export boundaries (NUOPC caps)

All caps follow the same pattern: `advertise_fields` / `realize_fields` /
`import_fields` / `export_fields` in `<comp>_import_export.F90`, driven by
`<comp>_comp_nuopc.F90`. The cap's mesh comes from the `mesh_<comp>`
attribute (§3). Classification: state fields are instantaneous (**D**);
flux fields are coupling-interval averages/accumulations (**I**); the
component's own time advance is **I**.

### 1.1 CAM — atmosphere (`cam/src/cpl/nuopc/`; `atm_comp_nuopc.F90`, `atm_import_export.F90`)

Boundary types: `cam_in_t`/`cam_out_t` (`cam/src/control/camsrfexch.F90`);
`import_fields(cam_in)`, `export_fields(cam_out)`; physics state in the CAM
physics buffer. Run: `cam_run1..4` (**I**).

| Direction | Field groups |
|---|---|
| Export (atm→med) | lowest-level state `Sa_z,u,v,tbot,ptem,shum,pbot,dens,pslv,topo,o3` (+`co2prog/co2diag`, `lightning`); radiation `Faxa_swdn,swvdr,swvdf,swndr,swndf,lwdn,swnet(diag)`; precip `Faxa_rainc,rainl,snowc,snowl`; aerosol/deposition `Faxa_bcph(3),ocph(3),dstwet(4),dstdry(4),ndep(2)` |
| Import (med→atm) | surface state/albedo `Sx_tref,qref,t,u10`, `Sx_avsdr,avsdf,anidr,anidf`, `So_t,ssq,re,ustar,ugustOut,u10withGust`, `Sl_lfrac,ifrac(ofrac via Si/So),fv,ram1,snowh(lnd/ice),ddvel,soilw,fztop`; surface fluxes `Faxx_taux,tauy,lat,sen,lwup,evap`, `Fall_flxdst(4),fco2_lnd,voc,fire`, `Faoo_fco2_ocn` |

### 1.2 CTSM — land (`clm/src/cpl/nuopc/`; `lnd_comp_nuopc.F90`, `lnd_import_export.F90`)

Run: `clm_drv(doalb, nextsw_cday, declin, ...)` (`lnd_comp_nuopc.F90:913`),
**I**. Imports feed CLM's atm-forcing type; exports come from CLM surface
flux/runoff buffers.

| Direction | Field groups |
|---|---|
| Export (lnd→med) | surface state `Sl_lfrin,t,tref,qref,avsdr,anidr,avsdf,anidf,snowh,u10,fv,ram1,soilw,fztop`; fluxes `Fall_taux,tauy,lat,sen,lwup,evap,swnet`, `Fall_flxdst(4),fco2_lnd,voc,fire`, `Sl_ddvel`; runoff `Flrl_rofsur,rofsub,rofgwl,rofi,irrig,rofsur_nonh2o`; glc-elevation `Sl_tsrf_elev,topo_elev,Flgl_qice_elev` (nec+1) |
| Import (med→lnd) | atm state `Sa_z,topo,u,v,ptem,pbot,tbot,shum,o3`, `Sa_co2prog/co2diag`; rad/precip `Faxa_lwdn,rainc,rainl,snowc,snowl,swndr,swvdr,swndf,swvdf`, `Faxa_bcph(3),ocph(3),dstwet(4),dstdry(4),ndep(2)`; rof feedback `Flrr_flood,volr,volrmch`, `Sr_tdepth,tdepth_max`; glc `Sg_icemask,icemask_coupled_fluxes,ice_covered_elev,topo_elev` |

### 1.3 CICE — sea ice (`cice/src/cicecore/drivers/nuopc/cmeps/`; `ice_comp_nuopc.F90`, `ice_import_export.F90`)

Run: `CICE_Run()` (`ice_comp_nuopc.F90:1269`), **I** (EVP dynamics +
Icepack column thermodynamics). CICE computes its own atm/ice and ice/ocn
fluxes (unlike MOM, which receives them from the mediator).

| Direction | Field groups |
|---|---|
| Import (med→ice) | atm state `Sa_z,u,v,shum,tbot,pbot,ptem,dens`; rad/precip `Faxa_swvdr,swvdf,swndr,swndf,lwdn,rain,snow`, `Faxa_bcph(3),dstwet(4),dstdry(4)`; ocn state `So_dhdx,dhdy,t,s,u,v`, `Fioo_q`; waves `Sw_elevation_spectrum` |
| Export (ice→med) | state `Si_imask,ifrac,ifrac_n(ncat),t,vice,vsno,tref,qref,snowh,u10,avsdr,avsdf,anidr,anidf`, `Si_thick,floediam` (wav); atm/ice fluxes `Faii_taux,tauy,lat,sen,lwup,evap,swnet`; ice/ocn fluxes `Fioi_melth,meltw,salt,taux,tauy,swpen(+vdr,vdf,idr,idf,ifrac_n)`, `Fioi_bcpho,bcphi,flxdst` |

### 1.4 MOM — ocean (`mom/MOM6/config_src/drivers/nuopc_cap/`; `mom_cap.F90`)

Run: `ModelAdvance` → `update_ocean_model` (internal subcycling), **I**.
Imports fill `Ice_ocean_boundary` (FMS coupler type); exports from
`ocean_public`. MARBL BGC and wave fields are optional.

| Direction | Field groups |
|---|---|
| Import (med→ocn) | merged fluxes `Foxx_taux,tauy,sen,lat,lwup,lwnet,evap,swnet(+vdr,vdf,idr,idf,afracr)`, `Faxa_rain,snow`, `Sa_pslv`; runoff `Foxx_rofl,rofi(+glc)`, `Forr_rofl_glc,rofi_glc`; heat-content `Foxx_hrain,hsnow,hevap,hcond,hrofl,hrofi(+glc)`; ice `Fioi_salt,meltw,melth,swpen_ifrac_n,flxdst,bcphi,bcpho`, `Si_ifrac,ifrac_n`, `Sf_afracr`; `So_duu10n`; MARBL `Faxa_ndep,dstwet(4),dstdry(4),bcph(3)`, `Sa_co2prog/co2diag`; waves `Sw_lamult`, `Sw_pstokes_x/y` |
| Export (ocn→med) | state `So_omask,t,s,u,v,dhdx,dhdy,bldepth,s_depth,t_depth`; flux `Fioo_q`; MARBL `Faoo_fco2_ocn` (and `fbrf_ocn,fdms_ocn,fn2o_ocn,fnh3_ocn` in mediator table) |

### 1.5 MOSART — river routing (`mosart/src/cpl/nuopc/`; `rof_comp_nuopc.F90`, `rof_import_export.F90`)

Run: `mosart_run(begr,endr,ntracers,rstwr,nlend,rdate)` (`rof_comp_nuopc.F90:714`),
**I** (see §2). Imports converted kg/m²/s → m³/s in `import_fields`
(`state_getimport` with area correction).

| Direction | Field groups |
|---|---|
| Import (med→rof) | `Flrl_rofsur` (surface), `Flrl_rofsub` (subsurface), `Flrl_rofgwl` (groundwater/glacier-wetland-lake), `Flrl_rofi` (frozen → routed to outlet), `Flrl_irrig` (demand), `Fgrg_rofl,rofi` (glc runoff, if `rof_from_glc`) |
| Export (rof→med) | `Forr_rofl,rofi` (to ocn), `Forr_rofl_glc,rofi_glc`, `Flrr_flood` (back to lnd), `Flrr_volr,volrmch` (channel storage for lnd), `Sr_tdepth,tdepth_max` (stream depths, if `flds_r2l_stream_channel_depths`) |

### 1.6 CISM — ice sheets (`cism/drivers/cpl/nuopc/`; `glc_comp_nuopc.F90`, `glc_import_export.F90`)

Run: `glc_run` → SIA/SSA dycore (`source_cism/libglide`/`libglissade`), **I**.
One instance per ice sheet (`mesh_glc` colon-delimited list).

| Direction | Field groups |
|---|---|
| Import (med→glc) | `Sl_tsrf` (surface temperature), `Flgl_qice` (SMB); optional `So_t_depth,s_depth` (ocean depth profiles, ocn2glc) |
| Export (glc→med) | `Sg_area,ice_covered,topo,icemask,icemask_coupled_fluxes`, `Flgg_hflx` (basal heat to lnd), `Figg_rofi` (to ice), `Fgrg_rofl,rofi` (to ocn, via rof) |

### 1.7 CDEPS data components (`cdeps/`)

Shared machinery in `dshr/` (`dshr_fldlist_*`, `dshr_dfield_mod`,
`dshr_mod`); each cap is a thin NUOPC wrapper around per-mode datamodes.
All are data-driven; classify **D** (data interpolation at coupling time;
no state) except docn SOM mode (slab-ocean SST integration, **I**).

| Cap | Datamodes (file) | Export fields (typical) | Import fields |
|---|---|---|---|
| datm (`datm/atm_comp_nuopc.F90`) | `datm_datamode_{clmncep,jra,core2,cplhist,era5,gefs,simple}_mod.F90` | `Sa_z,u,v,u10m,v10m,ptem,dens,pslv,tbot,pbot,shum`, `Faxa_rainc,rainl,snowc,snowl,swndr,swvdr,swndf,swvdf,swnet,lwdn,swdn` | optional `Faxx_*` (cplhist) |
| dlnd (`dlnd/lnd_comp_nuopc.F90`) | `dlnd_datamode_rof_forcing_mod.F90`, `dlnd_datamode_glc_forcing_mod.F90` | `Flrl_rofsur_nonh2o`, or `Sl_tsrf_elev,topo_elev,Flgl_qice_elev` (glc forcing) | — |
| docn (`docn/ocn_comp_nuopc.F90`) | `docn_datamode_{sstdata,multilev,multilev_sstdata,cplhist,multilev_cplhist,aquaplanet,som}_mod.F90` | `So_omask,t,u,v,bldepth,s` (+`dhdx,dhdy,Fioo_q` som; `So_fswpen` som) | som: `Foxx_swnet,lwup,sen,lat`, `Faxa_lwdn,snow`, `Fioi_melth`, `Foxx_rofi` … |
| dice (`dice/ice_comp_nuopc.F90`) | `dice_datamode_{ssmi,cplhist}_mod.F90`, `dice_flux_atmice_mod.F90` | `Si_ifrac,imask,t,vice,vsno`, `Faii_taux,tauy,lat,sen,lwup,evap`, `Si_avsdr,avsdf,anidr,anidf` | atm/ocn state + radiation |

## 2. MOSART internals (detailed — first implementation target)

All routing logic lives in `mosart/src/riverroute/`. The cap
(`src/cpl/nuopc/rof_comp_nuopc.F90`) calls
`mosart_read_namelist` / `mosart_init1` / `mosart_init2` / `mosart_run`
(`mosart_driver.F90`). Namelist `mosart_inparm` (`mosart_driver.F90:91`):
`frivinp` (network file), `coupling_period`, `delt_mosart` (default 3600 s),
`bypass_routing_option` (default `'direct_in_place'`; also
`'direct_to_outlet'`, `'none'`), `qgwl_runoff_option` (default `'threshold'`;
also `'all'`, `'negative'`), `mosart_tracers` (default `'LIQ:ICE'`),
`mosart_euler_calc`, `ice_runoff`, `budget_frq`. Control defaults
(`mosart_tctl_type.F90`): `RoutingMethod=1` (kinematic wave / variable
storage), `DLevelH2R=5`, `DLevelR=3`.

### 2.1 `mosart_run` orchestration (`mosart_driver.F90:375`) — **I**

Per coupling step: (1) load `ctl%qsur/qsub/qgwl` from imports (m³/s);
(2) **irrigation** — remove `qirrig` demand from main-channel storage
`wr`, overage returns to `qsub`, actual rate `qirrig_actual`; (3) **flood**
— if `volr(liq) > fthresh`, send `Flrr_flood` back to lnd (subtract from
forcing); (4) **direct transfer** — glc runoff, all frozen runoff, and
(per `qgwl_runoff_option`/`bypass_routing_option`) qgwl and negative runoff
are routed straight to the outlet via the `rh_direct` SMM, or kept
in-place; water outside the basin mask is always direct; (5) **subcycling**
— `nsub = coupling_period/delt_mosart`, `delt = delt_coupling/nsub`;
(6) per substep `Euler`; (7) average `flow,erout_prev,eroutup_avg,erlat_avg`
over subcycles; (8) update `volr = wt + wr + wh*area`,
`dvolrdt`, and split `runofflnd` (mask=1) vs `runoffocn` (mask≥2, outlet),
`runofftot = direct + runoffocn + direct_glc`; (9) budget check, history,
restart. Inputs→outputs: imports (`ctl%qsur,qsub,qgwl,qirrig,qglc_*`) and
state (`wt,wr,wh`) → exports (`ctl%runofftot,flood,volr,...`).

### 2.2 `Euler` (`mosart_physics.F90:51`) — **I** (state advance over `delt`)

Per tracer `nt` with `euler_calc(nt)` (frozen water is direct-only):
1. **hillslope** — `hillslopeRouting` per cell: overland outflow
   `ehout = -CREHT_nosqrt(...)` (Manning overland flow from storage `yh`),
   `dwh = qsur + ehout`, `wh += dwh*delt`; inflow to subnetwork
   `etin = (-ehout + qsub)*area*frac`. **I** (storage update; the law is **D**).
2. **subnetwork** (tributary channels) — per cell, `numDT_t(nr)` substeps of
   `subnetworkRouting`: velocity `vt = CRVRMAN_nosqrt(tslpsqrt, nt, rt)`
   (Manning), outflow `etout = -vt*mt`, `dwt = etin + etout`,
   `wt += dwt*localdelt`; dry-cell clamp keeps `wt ≥ TINYVALUE`. **I**.
3. **upstream interaction** — `erout` of each cell is mapped to the
   downstream cell's `eroutUp` by the `rh_eroutUp` SMM route handle. **D**
   (sparse-matrix gather).
4. **main channel** — per cell, `numDT_r(nr)` substeps of
   `mainchannelRouting` → `Routing_KW`: inflow `erin = -eroutUp`,
   velocity `vr = CRVRMAN_nosqrt(rslpsqrt, nr, rr)`, outflow
   `erout = -vr*mr`, storage tendency
   `dwr = erlateral + erin + erout + qgwl*area*frac`, `wr += dwr*localdelt`;
   `flow -= erout`; negative-storage check (`negchan < -1e-10` warns). **I**.
5. `erout_prev` (for budget), `eroutup_avg`, `erlat_avg` accumulated; all
   divided by `DLevelH2R` at the end.

`DLevelH2R` outer levels; per-cell CFL substepping from
`set_subtimesteps` (`mosart_tspatialunit_type.F90:624`):
`phi_r = areaTotal2*sqrt(rslp)/(rlen*rwidth)`,
`phi_t = area*sqrt(tslp)/(tlen*twidth)` → `numDT_r/t` (log10 scaling for
`phi ≥ 10`), computed once at init. **D** (diagnostic).

### 2.3 Pointwise flow laws and geometry functions (`mosart_physics.F90:406-574`) — all **D**

| Function | Law |
|---|---|
| `CRVRMAN_nosqrt(sqrtslp, n, rr)` | Manning channel velocity `v = (rr²)^(1/3)·√s/n` (rr = hydraulic radius, n = roughness) |
| `CREHT_nosqrt(sqrthslp, nh, Gxr, yh)` | hillslope overland flow `eht = (yh²)^(1/3)·√s·Gxr/nh` (drainage density Gxr) |
| `GRMR(wr, rlen)` | channel cross-section area `m = wr/rlen` |
| `GRHR(mr, rwidth, rwidth0, rdepth)` | channel depth from storage (bankfull/trapezoid geometry, `rwid0/rdep` from network file) |
| `GRPR(hr, rwidth, rwidth0, rdepth)` | wetted perimeter |
| `GRRR(mr, pr)` | hydraulic radius `r = m/p` |
| `GRHT(mt, twid)`, `GRPT(yt, twid)` | tributary depth / wetted perimeter |

`updateState_hillslope/subnetwork/mainchannel` (`:352-405`) recompute
`yh, mt,yt,pt,rt / mr,yr,pr,rr` from storage — **D**.

### 2.4 Network / grid handling (`mosart_tspatialunit_type.F90`)

`Tspatialunit_type::Init` (`:92`) reads the routing network file `frivinp`
(NetCDF, PIO) with variables: `frac` (land fraction), `fdir` (flow
direction → `mask`: <0 ocean, 0 outlet, >0 land), `ID`, `dnID`
(downstream-cell IDs — the network topology), `area`, `areaTotal`,
`hslp,nh` (hillslope), `tslp,twid,tlen,nt` (tributary), `rslp,rwid,rwid0,
rdep,rlen,nr` (main channel), `gxr` (drainage density). Routing is
**ID/dnID-network based**, not lat-lon; the lat-lon grid enters only via
the ESMF mesh (`mesh_rof` attribute, `rof_comp_nuopc.F90:526`) and the
cap's `nx_global/ny_global` scalars. `set_routehandles` (`:430`) builds the
two ESMF SMM route handles from `gindex`/`outletg`: `rh_direct`
(gridcell→basin outlet, for direct transfer) and `rh_eroutUp`
(cell→downstream cell). `set_areatotal2` (`:533`) accumulates upstream
catchment area along `dnID`.

### 2.5 Lakes, reservoirs, groundwater — TBD (verify at instrumentation time)

This MOSART version (mosart1.1.13) has **no separate lake/reservoir or
groundwater module**: lakes/glaciers/wetlands enter only as the `qgwl`
runoff stream from CTSM (`Flrl_rofgwl`), with the `qgwl_runoff_option`
direct-handling (`'all'`/`'negative'`/`'threshold'` — the threshold option
protects a minimum channel storage `river_depth_minimum=1e-4 m`). Channel
storage regulation exists only via the flood term (§2.1). If a reference
compset needs lake retention or groundwater dynamics, confirm whether they
live in CTSM (they do not appear in MOSART) and mark `mosart_lakes.esm` /
`mosart_groundwater.esm` as TBD in the §3 inventory.

### 2.6 MOSART subassembly → `.esm` mapping (PLAN.md §3 row)

| Subassembly | Files / entry | D/I | Proposed `.esm` |
|---|---|---|---|
| Cap boundary | `src/cpl/nuopc/rof_import_export.F90` (import_fields/export_fields, kg/m²s↔m³/s) | D (unit conversion) | `mosart_coupling.esm` |
| Direct transfer | `mosart_driver.F90:mosart_run` steps 2–4 | I (over coupling period) | `mosart_coupling.esm` |
| Irrigation/flood | `mosart_run` steps 2–3 | D (rate rules) | `mosart_routing.esm` |
| Routing (hillslope+subnetwork+channel) | `mosart_physics.F90:Euler`, `*Routing` | I (storage advance) | `mosart_routing.esm` |
| Flow laws / geometry | `CRVRMAN_nosqrt, CREHT_nosqrt, GR*` | **D — first test targets** | `mosart_routing.esm` (expression templates) |
| Network topology | `mosart_tspatialunit_type.F90:Init, set_routehandles, set_subtimesteps` | D (init-time) | data source + discretization note |
| Lakes/reservoirs/groundwater | — (not present; §2.5) | — | TBD |
| State/budget | `mosart_tstatusflux_type.F90` (TRunoff), `mosart_budget_type.F90`, `mosart_restfile.F90` | I (state) | harness |

## 3. Grids and precision notes

### 3.1 Grids

Each cap realizes its fields on an ESMF mesh given by the `mesh_<comp>`
mediator attribute (`cmeps/cime_config/namelist_definition_drv.xml:4079+`),
set from the CIME XML variables `ATM_DOMAIN_MESH`, `LND_DOMAIN_MESH`,
`ICE_DOMAIN_MESH`, `OCN_DOMAIN_MESH`, `ROF_DOMAIN_MESH`, `GLC_DOMAIN_MESH`,
`WAV_DOMAIN_MESH`. The mediator builds pairwise regridding between every
connected pair (`med_map_mod.F90`): `mapbilnr` (state), `mapconsf/mapconsd`
(fluxes), `mappatch`, `mapfcopy`, `mapnstod`, plus `idmap` when two
components share a mesh (`esmFldsExchange_cesm_mod.F90:196`); uv vectors
mapped with `mappatch_uv3d`/`mapuv_with_cart3d`. Grid conventions:

| Component | Grid | Notes |
|---|---|---|
| CAM | SE cubed sphere (e.g. `ne30np4`) or FV lat-lon (`f19`) | mesh from `mesh_atm`; single-column mode from `PTS_LON/LAT` |
| CTSM | same grid as CAM (idmap) | `mesh_lnd` |
| CICE | same grid as CAM in coupled CESM (`ICE_GRID`) | `mesh_ice`; `ICE_NX/NY` scalars |
| MOM | tripolar (`tx0.66v1`, `gx1v7`, …) | `mesh_ocn`; `So_omask` ocean mask |
| MOSART | lat-lon routing grid (`r05`, `r0125`) | `mesh_rof`; network topology from `frivinp` (§2.4) |
| CISM | on the lnd grid, one mesh per ice sheet (`mesh_glc` list) | elevation classes (nec) in glc↔lnd fields |
| WW3 | global tripole (`mesh_wav`) | out of scope |

The mediator's own flux grids: `aoflux_grid` attribute — `'ogrid'`
(default, fluxes on the ocean grid), `'agrid'`, `'xgrid'` (exchange grid).
Fractions (`med_fraction_mod.F90`) are computed on the `med` (atm) grid and
mapped with the fraction maps (`mrg_fracname_*`).

### 3.2 Precision

All components run real64 in the default CESM configuration; the mediator
fields are `ESMF_TYPEKIND_R8` throughout:

- **CAM/CTSM/MOSART/CISM/CDEPS:** hard-coded `R8`/`r8` kinds (real64).
- **CICE:** `dbl_kind = selected_real_kind(13)` (real64)
  (`cice/src/cicecore/shared/ice_kinds_mod.F90` → `icepack_kinds.F90:30`);
  `NO_R16` CPP controls `r16_kind` only. CICE's standalone build
  (`cice.setup`) historically offers single-precision; the cime config in
  this checkout does not expose a `RUNTIME_PRECISION` option — **TBD
  (verify at instrumentation time)** whether the single-precision code path
  can be selected in this CESM configuration.
- **MOM6/FMS:** MOM6 source uses the default real kind; the FMS library is
  built with `-fdefault-real-8 -fdefault-double-8`
  (`libraries/FMS/Makefile.cesm:91`), fixing real64. (PLAN.md §4.2 note:
  if any component were run single-precision, reference values would have a
  real32 floor — EqWeFiC finding. Not applicable to the default config.)

Per PLAN.md §4.2, record the exact run configuration (grid, timestep,
coupling interval, namelist, precision) with every instrumented run.

## 4. Instrumentation notes

- **Boundary hooks:** each cap's `import_fields`/`export_fields` (the NUOPC
  importState/exportState field bundles) and the mediator's
  `FBImp/FBExp/FBMed` field bundles (`med_internalstate_mod.F90`) are the
  natural dump points for the §1 tables. Dump per coupling step, with the
  run-sequence phase name and timestamp.
- **Preferred D boundaries for first instrumentation (MOSART):**
  `CRVRMAN_nosqrt`, `CREHT_nosqrt`, `GRMR/GRHR/GRPR/GRRR/GRHT/GRPT`,
  `updateState_*`, `Routing_KW`'s outflow/storage-tendency computation,
  irrigation/flood rate rules, `med_phases_ocnalb_run`,
  `med_phases_aofluxes_run` (per-point flux kernels), `med_fraction_set`.
- **I boundaries needing (t, dt) + state dumps:** `mosart_run` (per
  coupling step, with `wt/wr/wh` before/after and the subcycle-averaged
  flows), `Euler` (per substep), each component's `*_run` entry
  (`cam_run*`, `clm_drv`, `CICE_Run`, `update_ocean_model`, `glc_run`),
  and the mediator accumulation phases (`prep_ocn_accum`, `prep_rof_accum`).
- **Budget check:** `mosart_budget_type.F90` already computes a water
  budget per coupling step — reuse it as an instrumented conservation test.
