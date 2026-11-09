#!/usr/bin/env python3
"""Compute the derived (non-dump-column) SoilFluxes outputs for the stage-2
clm_soilfluxes.esm tests.

The stage-1 CLM dumps trace only the columns listed in the clm_soilfluxes.csv
header (t_grnd, eflx_soil_grnd, eflx_lh_grnd, eflx_sh_grnd, eflx_lwrad_net,
eflx_sh_veg, eflx_sh_tot, eflx_lh_tot, qflx_evap_soi, cgrnds, cgrndl,
t_skin). The remaining outputs of SoilFluxesMod.F90 (SoilFluxes 38-518) --
the previous-time-step ground-temperature blend t_grnd0, the temperature
increment tinc, the longwave-weighted T^4 lw_grnd, the corrected raw fluxes
qflx_ev_{snow,soil,h2osfc}_raw / qflx_evap_soi_raw / eflx_sh_grnd_raw, the
evaporation limit, the partition fluxes and the black-body emission terms --
are NOT dump columns; this script computes them by evaluating the transcribed
equations at the pinned dump inputs (the same f64 arithmetic as the esm
replay, cross-checked on every dump column they touch to max rel err 7.9e-15,
see stage2/README-clm-soilfluxes.md), and emits the esm assertion tuples for
the rows the .esm tests pin.

Inputs are read from the stage-1 dumps exactly as the .esm parameter_overrides
do (see the input-source table in stage2/README-clm-soilfluxes.md):
  t_grnd, cgrnds, cgrndl              clm_soilfluxes.csv
  eflx_sh_grnd_pre, qflx_evap_soi_pre,
  qflx_ev_snow_pre, qflx_ev_soil_pre   clm_bareground.csv
  tssbef_snow, tssbef_soil (entry), h2osoi_liq/ice_top (exit),
  frac_sno_fluxes, snl                clm_soiltemp.csv (j = snl+1 / j = 1)
  sabg_soil, sabg_snow                clm_surfrad.csv
  emg, htvp                           clm_run_imports.csv (preflux stage)
  forc_lwrad                          cdeps/datm_gswp3_exports.csv Faxa_lwdn
                                      (datm step s-1, n = 1; columns 1-5 only)
and the run constants dtime = 1800 s (get_step_size_real), sb = SHR_CONST_
STEBOL = 5.67e-8 (shr_const_mod.F90:23), hvap = SHR_CONST_LATVAP = 2.501e6,
tfrz = SHR_CONST_TKFRZ = 273.15, dlrad = ulrad = 0 (BareGroundFluxesMod.F90:
313-314), frac_veg_nosno = 0 and qflx_evap_veg = qflx_tran_veg = 0 (bare
traced patches), frac_h2osfc = 0 and qflx_ev_h2osfc_pre = 0 (not dumped).

Usage: python3 clm_soilfluxes_derived.py [--all] [step patch ...]
Default: the four rows pinned by the .esm tests -- (1,1) (12,1) (7,8) (6,16)
in the order of the tests array; --all prints every step/patch row.
"""
import csv
import math
import os
import sys

# The dumps are gitignored and live in the main checkout ($M/stage1/dumps),
# not in the worktree; fall back to the relative location when run there.
DUMPS = os.environ.get('DUMPS',
                       '/projects/illinois/eng/cee/ctessum/ctessum/code/EqESM/'
                       'stage1/dumps')
if not os.path.isdir(DUMPS):
    DUMPS = os.path.join(os.path.dirname(__file__), '..', 'dumps')

DTIME = 1800.0      # get_step_size_real, SoilFluxesMod.F90:165
SB = 5.67e-8        # SHR_CONST_STEBOL, shr_const_mod.F90:23
HVAP = 2.501e6      # SHR_CONST_LATVAP
TFREEZE = 273.15    # SHR_CONST_TKFRZ


def read_rows(path, key):
    """Read a dump csv into {key_tuple: row_dict}.

    The dump column names live in the '# columns = ...' comment line (the
    data rows start immediately after it), so that line is parsed as the
    header instead of the first data row.
    """
    rows = {}
    with open(path) as f:
        fields = None
        data = []
        for line in f:
            if line.startswith('# columns'):
                fields = line.split('=', 1)[1].strip().split(',')
            elif not line.startswith('#'):
                data.append(line)
        reader = csv.DictReader(data, fieldnames=fields)
        for r in reader:
            def cast(v):
                try:
                    return int(v)
                except ValueError:
                    return v
            rows[tuple(cast(r[k]) for k in key)] = r
    return rows


def main():
    dumps = os.environ.get('DUMPS', DUMPS)
    soilfluxes = read_rows(os.path.join(dumps, 'clm', 'clm_soilfluxes.csv'),
                           ('step', 'p'))
    bareground = read_rows(os.path.join(dumps, 'clm', 'clm_bareground.csv'),
                           ('step', 'p'))
    soiltemp = read_rows(os.path.join(dumps, 'clm', 'clm_soiltemp.csv'),
                         ('step', 'stage', 'c', 'j'))
    surfrad = read_rows(os.path.join(dumps, 'clm', 'clm_surfrad.csv'),
                        ('step', 'p'))
    runimports = read_rows(os.path.join(dumps, 'clm', 'clm_run_imports.csv'),
                           ('step', 'stage', 'p'))
    datm = read_rows(os.path.join(dumps, 'cdeps', 'datm_gswp3_exports.csv'),
                     ('step', 'n'))

    if '--all' in sys.argv:
        wanted = [(int(r['step']), int(r['p']))
                  for r in soilfluxes.values()]
    else:
        pairs = sys.argv[1:]
        wanted = [(int(pairs[i]), int(pairs[i + 1]))
                  for i in range(0, len(pairs), 2)] if pairs else \
            [(1, 1), (12, 1), (7, 8), (6, 16)]

    for step, p in wanted:
        sf = soilfluxes[(step, p)]
        bg = bareground[(step, p)]
        ri = runimports[(step, 'preflux', p)]
        sr = surfrad[(step, p)]
        c = int(sf['c'])
        snl = int(float(soiltemp[(step, 'entry', c, -1)]['snl']))
        jtop = snl + 1
        st_e = soiltemp[(step, 'entry', c, jtop)]
        st_x = soiltemp[(step, 'exit', c, jtop)]
        soil1 = soiltemp[(step, 'entry', c, 1)]

        t_grnd = float(sf['t_grnd'])
        cgrnds = float(sf['cgrnds'])
        cgrndl = float(sf['cgrndl'])
        tssbef_snow = float(st_e['t_soisno'])
        tssbef_soil = float(soil1['t_soisno'])
        frac_sno_fluxes = float(st_e['frac_sno_fluxes'])
        frac_h2osfc = 0.0            # not dumped; 0 on all traced rows
        h2osoi_liq = float(st_x['h2osoi_liq'])
        h2osoi_ice = float(st_x['h2osoi_ice'])
        eflx_sh_grnd_pre = float(bg['eflx_sh_grnd'])
        qflx_evap_soi_pre = float(bg['qflx_evap_soi'])
        qflx_ev_snow_pre = float(bg['qflx_ev_snow'])
        qflx_ev_soil_pre = float(bg['qflx_ev_soil'])
        sabg_soil = float(sr['sabg_soil'])
        sabg_snow = float(sr['sabg_snow'])
        emg = float(ri['emg'])
        htvp = float(ri['htvp'])
        # DATM Faxa_lwdn at datm step s-1, point n = 1 (columns 1-5 only);
        # columns 6-16 have no forc_lwrad (the eflx_soil_grnd / eflx_lwrad_*
        # pins are left out there).
        forc_lwrad = (float(datm[(step - 1, 1)]['Faxa_lwdn'])
                      if c <= 5 else None)
        dlrad = 0.0                  # BareGroundFluxesMod.F90:313
        ulrad = 0.0                  # BareGroundFluxesMod.F90:314
        frac_veg_nosno = 0.0         # bare traced patches
        qflx_evap_veg = 0.0
        qflx_tran_veg = 0.0
        eflx_sh_veg = float(sf['eflx_sh_veg'])
        eflx_sh_stem = 0.0           # bare traced patches
        t_h2osfc_bef = 0.0           # 0-weight (frac_h2osfc = 0)

        # --- the transcribed equations (SoilFluxesMod.F90) ---
        if snl < 0:
            t_grnd0 = (frac_sno_fluxes * tssbef_snow
                       + (1 - frac_sno_fluxes - frac_h2osfc) * tssbef_soil
                       + frac_h2osfc * t_h2osfc_bef)                 # :176-178
        else:
            t_grnd0 = ((1 - frac_h2osfc) * tssbef_soil
                       + frac_h2osfc * t_h2osfc_bef)                 # :180
        tinc = t_grnd - t_grnd0                                      # :183
        lw_grnd = (frac_sno_fluxes * tssbef_snow**4
                   + (1 - frac_sno_fluxes - frac_h2osfc) * tssbef_soil**4
                   + frac_h2osfc * t_h2osfc_bef**4)                  # :344-346

        eflx_sh_grnd_raw = eflx_sh_grnd_pre + tinc * cgrnds          # :192
        qflx_evap_soi_raw = qflx_evap_soi_pre + tinc * cgrndl        # :193
        qflx_ev_snow_raw = qflx_ev_snow_pre + tinc * cgrndl          # :202
        qflx_ev_soil_raw = qflx_ev_soil_pre + tinc * cgrndl          # :203
        qflx_ev_h2osfc = 0.0 + tinc * cgrndl                         # :204

        evaporation_limit = ((h2osoi_liq + h2osoi_ice)
                             / (frac_sno_fluxes * DTIME))            # :287
        qflx_ev_snow_excess = (max(qflx_ev_snow_raw - evaporation_limit, 0.0)
                               if snl < 0 else 0.0)                  # :288-296
        qflx_ev_snow = qflx_ev_snow_raw - qflx_ev_snow_excess        # :291
        qflx_ev_soil = qflx_ev_soil_raw                              # :203
        qflx_evap_soi = qflx_evap_soi_raw - frac_sno_fluxes * qflx_ev_snow_excess
        eflx_sh_grnd = (eflx_sh_grnd_raw
                        + frac_sno_fluxes * qflx_ev_snow_excess * htvp)

        h2osoi_top = h2osoi_liq + h2osoi_ice                         # :95-96
        if qflx_ev_snow >= 0.0:                                      # :238-243
            if h2osoi_top > 0.0:
                qflx_liqevap = max(qflx_ev_snow * h2osoi_liq / h2osoi_top, 0.0)
            else:
                qflx_liqevap = 0.0
            qflx_solidevap = qflx_ev_snow - qflx_liqevap
            qflx_soliddew = 0.0
            qflx_liqdew = 0.0
        else:                                                        # :244-248
            qflx_liqevap = 0.0
            qflx_solidevap = 0.0
            if t_grnd < TFREEZE:
                qflx_soliddew = abs(qflx_ev_snow)
                qflx_liqdew = 0.0
            else:
                qflx_soliddew = 0.0
                qflx_liqdew = abs(qflx_ev_snow)

        lwrad_emission = emg * SB * lw_grnd                          # :350,471
        lwrad_increment = 4.0 * emg * SB * t_grnd0**3 * tinc          # :350,472
        # columns 6-16 have no forc_lwrad; their forc_lwrad-dependent
        # outputs are unpinned (a documented gap, README-clm-soilfluxes.md).
        if forc_lwrad is not None:
            eflx_soil_grnd = ((1 - frac_sno_fluxes) * sabg_soil
                              + frac_sno_fluxes * sabg_snow + dlrad
                              + (1 - frac_veg_nosno) * emg * forc_lwrad
                              - lwrad_emission - lwrad_increment
                              - (eflx_sh_grnd + qflx_evap_soi * htvp))  # :348-351
            eflx_lwrad_out = (ulrad + (1 - frac_veg_nosno) * (1 - emg)
                              * forc_lwrad + (1 - frac_veg_nosno)
                              * lwrad_emission + lwrad_increment)   # :469-472
            eflx_lwrad_net = eflx_lwrad_out - forc_lwrad            # :485
        else:
            eflx_soil_grnd = eflx_lwrad_out = eflx_lwrad_net = None

        eflx_sh_tot = eflx_sh_veg + eflx_sh_grnd + eflx_sh_stem      # :375-376
        qflx_evap_tot = qflx_evap_veg + qflx_evap_soi                # :377
        eflx_lh_tot = HVAP * qflx_evap_veg + htvp * qflx_evap_soi    # :379
        qflx_evap_can = qflx_evap_veg - qflx_tran_veg                # :390
        eflx_lh_vege = qflx_evap_can * HVAP                          # :391
        eflx_lh_vegt = qflx_tran_veg * HVAP                          # :392
        eflx_lh_grnd = qflx_evap_soi * htvp                          # :393

        if frac_veg_nosno < 0.5:                                     # :483
            t_skin = math.sqrt(math.sqrt(lw_grnd))
        else:
            t_skin = None  # t_skin_veg, set in CanopyFluxes (untested)

        print(f'# step {step}, patch {p} (column {c}): snl = {snl}, '
              f'tinc = {tinc:.6e}, forc_lwrad = {forc_lwrad}')
        for name, value in [
                ('t_grnd0', t_grnd0), ('tinc', tinc), ('lw_grnd', lw_grnd),
                ('qflx_ev_snow_raw', qflx_ev_snow_raw),
                ('qflx_ev_soil_raw', qflx_ev_soil_raw),
                ('qflx_evap_soi_raw', qflx_evap_soi_raw),
                ('eflx_sh_grnd_raw', eflx_sh_grnd_raw),
                ('evaporation_limit', evaporation_limit),
                ('qflx_ev_snow_excess', qflx_ev_snow_excess),
                ('qflx_ev_snow', qflx_ev_snow), ('qflx_ev_soil', qflx_ev_soil),
                ('qflx_evap_soi', qflx_evap_soi),
                ('eflx_sh_grnd', eflx_sh_grnd),
                ('qflx_liqevap_from_top_layer', qflx_liqevap),
                ('qflx_solidevap_from_top_layer', qflx_solidevap),
                ('qflx_soliddew_to_top_layer', qflx_soliddew),
                ('qflx_liqdew_to_top_layer', qflx_liqdew),
                ('lwrad_emission', lwrad_emission),
                ('lwrad_increment', lwrad_increment),
                ('eflx_soil_grnd', eflx_soil_grnd),
                ('eflx_sh_tot', eflx_sh_tot), ('qflx_evap_tot', qflx_evap_tot),
                ('eflx_lh_tot', eflx_lh_tot), ('qflx_evap_can', qflx_evap_can),
                ('eflx_lh_vege', eflx_lh_vege), ('eflx_lh_vegt', eflx_lh_vegt),
                ('eflx_lh_grnd', eflx_lh_grnd),
                ('eflx_lwrad_out', eflx_lwrad_out),
                ('eflx_lwrad_net', eflx_lwrad_net), ('t_skin', t_skin)]:
            print(f'  {name:34s} {value!r}')


if __name__ == '__main__':
    main()
