#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2018-2024 Institute of Computer Science of the Czech Academy of
# Sciences, Prague, Czech Republic. Authors: Pavel Krc, Martin Bures, Jaroslav
# Resler.
#
# This file is part of PALM-METEO.
#
# PALM-METEO is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM-METEO is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM-METEO. If not, see <https://www.gnu.org/licenses/>.

import glob
import numpy as np
import scipy.ndimage as ndimage
import netCDF4
import cftime

from palmmeteo.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from palmmeteo.logging import die, warn, log, verbose, log_output
from palmmeteo.config import cfg, ConfigError, parse_duration
from palmmeteo.runtime import rt
from palmmeteo.utils import ensure_dimension, ax_, rad, NotWholeTimestep, utcdefault
from palmmeteo.library import PalmPhysics, TriRegridder, verify_palm_hinterp, AssimCycle, HorizonSelection
from palmmeteo.vinterp import get_vinterp

barom_pres = PalmPhysics.barom_lapse0_pres
barom_gp = PalmPhysics.barom_lapse0_gp
g = PalmPhysics.g

def log_dstat_on(desc, delta):
    """Calculate and log delta statistics if enabled."""
    log_output('{0} ({1:8g} ~ {2:8g}): bias = {3:8g}, MAE = {4:8g}\n'.format(
        desc, delta.min(), delta.max(), delta.mean(), np.abs(delta).mean()))

def log_dstat_off(desc, delta):
    """Do nothing (log disabled)"""
    pass

def get_3dvar(f, varname):
    # ICON NetCDF files have generated varnames and dimnames from sequences in
    # random order (grib messages), therefore we need to find the correct var
    # among matching names. We want one with z dimension 'height*' and not
    # 'plev*'
    v = f.variables[varname]
    num = 1
    while not v.dimensions[1].startswith('height'):
        num += 1
        v = f.variables[f'{varname}_{num}']
    return v

def assign_time(coll, key, value):
    """
    Assigns value = coll[key] only if that was previously None,
    otherwise raises an error.
    """
    if coll[key] is not None:
        raise RuntimeError(f'Time index {key} was already loaded!')
    coll[key] = value

class IconPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def check_config(self, *args, **kwargs):
        if cfg.icon.assim_cycles.cycles_used == 'all':
            raise ConfigError('ICON must use either cycle interval or a '
                    'single cycle, not all cycles', cfg.icon.assim_cycles,
                    'cycles_used')
        if cfg.radiation.enabled and not rt.paths.icon.static_data:
            raise ConfigError('For radiation with ICON, the static data file '
                    'with surface emissivity must be specified', cfg.path,
                    'icon_static_data')
        if (rt.timestep_rad is not None
                and rt.timestep_rad != rt.simulation.timestep):
            raise ConfigError('For radiation with ICON, the radiation timestep '
                    'can not be configured to a different value than the main '
                    'timestep', cfg.radiation, 'timestep')
        if cfg.output.geostrophic_wind:
            raise ConfigError('The ICON plugin does not support geostrophic wind',
                              cfg.output, 'geostrophic_wind')

    def import_data(self, fout, *args, **kwargs):
        log('Importing ICON data...')

        # Prepare time indices
        cycles = AssimCycle(cfg.icon.assim_cycles)
        hor0 = parse_duration(cfg.icon.assim_cycles, 'earliest_horizon')
        hselect = HorizonSelection(cycles, hor0)
        cycle_first_wanted = hselect.dt_from_idx(0)[0]
        cycle_last_wanted = hselect.dt_from_idx(rt.nt)[0]
        rt.times = [None] * rt.nt

        verbose('ICON horizon range for instantaneous values: {}...{}',
                hselect.horiz_first, hselect.horiz_last)

        if cfg.radiation.enabled:
            # Radiation uses same timestep as IBC in ICON
            rt.timestep_rad = rt.simulation.timestep
            rt.nt_rad = rt.tindex(rt.simulation.end_time_rad) + 2 #=nt+1 unless nested

            # Because we disaggregate the time intervals between timesteps, the
            # values represent the centres of the intervals. We start with half
            # a timestep before the simulation and end with half a timestep
            # after (PALM deals with this correctly).
            rt.times_rad_sec = np.arange(-0.5, rt.nt_rad-1) * rt.timestep_rad.total_seconds()

            # Prepare aggregated values
            aggr_start = [None] * (rt.nt_rad)
            aggr_end   = [None] * (rt.nt_rad)
            hselect_rad_left = HorizonSelection(cycles, hor0, -1, rt.nt_rad-1)
            hselect_rad_right = HorizonSelection(cycles, hor0+rt.timestep_rad, 0, rt.nt_rad)

            verbose('ICON horizon range for aggregated values: {}...{}',
                    hselect_rad_left.horiz_first, hselect_rad_right.horiz_last)

        # HHL is only present at time 0 of the run, so we need to cache it
        hhl_d = {}
        run_d = {}

        # Process input files
        verbose('Parsing ICON files from {}', rt.paths.icon.file_mask)
        first = True
        for fn in sorted(glob.glob(rt.paths.icon.file_mask)):
            verbose('Parsing ICON file {}', fn)
            with netCDF4.Dataset(fn) as fin:
                # Decode time and locate timestep
                tvar = fin.variables['time']
                tcycle = utcdefault(cftime.num2date(0, tvar.units, tvar.calendar,
                                                    only_use_cftime_datetimes=False,
                                                    only_use_python_datetimes=True))
                if not cycles.is_selected(tcycle):
                    verbose('Cycle {} not among selected ICON cycles - skipping', tcycle)
                    continue

                assert tvar.shape == (1,)
                t = utcdefault(cftime.num2date(tvar[0], tvar.units, tvar.calendar,
                                               only_use_cftime_datetimes=False,
                                               only_use_python_datetimes=True))
                try:
                    it = rt.tindex(t)
                except NotWholeTimestep:
                    verbose('Time {} is not within timestep intervals - skipping', t)
                    continue
                thoriz = t - tcycle

                # Find indices, decide on using
                use = False
                it_main = hselect.get_idx(thoriz, it)
                if it_main is not False:
                    verbose('File is used for meteorology.')
                    use = True
                if (cycle_first_wanted <= tcycle <= cycle_last_wanted) and not thoriz:
                    verbose('File is used for HHL.')
                    use = True

                if cfg.radiation.enabled:
                    it_rad_left = hselect_rad_left.get_idx(thoriz, it)
                    if it_rad_left is not False:
                        verbose('File is used for left-side deaggregation of radiation.')
                        use = True
                    it_rad_right = hselect_rad_right.get_idx(thoriz, it)
                    if it_rad_right is not False:
                        verbose('File is used for right-side deaggregation of radiation.')
                        use = True

                if not use:
                    verbose('File is not used for anything.')
                    continue

                if first:
                    # coordinate projection
                    verbose('Preparing regridder')

                    clat = fin.variables['CLAT'][0]
                    clon = fin.variables['CLON'][0]
                    rt.regrid_icon = TriRegridder(clat, clon,
                                                  rt.palm_grid_lat, rt.palm_grid_lon,
                                                  cfg.icon.point_selection_buffer)

                    if cfg.hinterp.validate:
                        verbose('Validating horizontal inteprolation.')
                        verify_palm_hinterp(rt.regrid_icon,
                                            rt.regrid_icon.loader(clat)[...],
                                            rt.regrid_icon.loader(clon)[...])

                    # dimensions
                    ensure_dimension(fout, 'time', rt.nt)

                    ntime, nz, ncells_full = get_3dvar(fin, 'P').shape
                    ensure_dimension(fout, 'ncells', rt.regrid_icon.npt) #selection of points
                    ensure_dimension(fout, 'z_meteo', nz)
                    ensure_dimension(fout, 'zw_meteo', nz+1)

                    # Soil layers
                    zsoil_tso_dim = fin.variables['T_SO'].dimensions[1]
                    zsoil_tso = fin.variables[zsoil_tso_dim][:]
                    zsoil_wso_dim = fin.variables['W_SO'].dimensions[1]
                    zsoil_wso = fin.variables[zsoil_wso_dim][:]
                    zsoil_dwso = (zsoil_wso[1:] - zsoil_wso[:-1])[:,ax_]

                    assert zsoil_tso[-2] <= zsoil_wso[-1]

                    rt.nz_soil = len(zsoil_tso) - 1
                    ensure_dimension(fout, 'zsoil', rt.nz_soil)
                    rt.z_soil_levels = zsoil_tso[:-1].data.tolist()
                    verbose('Z soil levels: {}', rt.z_soil_levels)
                    wsoil_coords = np.searchsorted(zsoil_wso, zsoil_tso[:-1], 'right') - 1

                    verbose('Initializing import file variables.')
                    for varname in cfg.icon.vars_3d:
                        fout.createVariable(varname, 'f4', ['time', 'z_meteo', 'ncells'])
                    for varname in cfg.icon.vars_3dw + ['HHL']:
                        fout.createVariable(varname, 'f4', ['time', 'zw_meteo', 'ncells'])
                    for varname in cfg.icon.vars_2d:
                        fout.createVariable(varname, 'f4', ['time', 'ncells'])
                    for varname in cfg.icon.vars_soil:
                        fout.createVariable(varname, 'f4', ['time', 'zsoil', 'ncells'])

                    if cfg.radiation.enabled:
                        verbose('Building list of indices for radiation smoothing.')

                        deg_range = cfg.icon.radiation_smoothing_distance / (PalmPhysics.radius*rad)
                        rad_mask = np.hypot((clon-rt.cent_lon)*rt.regrid_icon.lon_coef,
                                             clat-rt.cent_lat) <= deg_range

                        verbose('Using {} points for radiation', rad_mask.sum())

                        log('Converting net longwave radiation to downward using surface '
                            'emissivity from a static file. CAUTION - this is only valid '
                            'without snow cover!')
                        with netCDF4.Dataset(rt.paths.icon.static_data) as fst:
                            vemis = fst.variables['EMIS_RAD']
                            assert vemis.shape == clat.shape
                            emis = vemis[:][rad_mask]
                            emis_sigma = emis * PalmPhysics.sigma_sb
                            emis_r = 1. / emis
                            del emis
                        if zsoil_tso[0] != 0.0:
                            die('Conversion of net longwave radiation '
                                    'requries surface soil temperature.')

                    first = False
                else:
                    # Verify soil layers
                    zsoil_tso_dim = fin.variables['T_SO'].dimensions[1]
                    zsl2 = fin.variables[zsoil_tso_dim][:-1]
                    assert(np.abs(zsl2-rt.z_soil_levels).max() < 0.01)

                # Process aggregate values
                if cfg.radiation.enabled:
                    # Verified: ASWDIR_S + ASWDIFD_S - ASWDIFU_S (upward) == ASOB_S (net),
                    # therefore ASWDIR_S must be on horizontal (not normal) plane
                    swdir = fin.variables['ASWDIR_S'][0][rad_mask].mean()
                    swdif = fin.variables['ASWDIFD_S'][0][rad_mask].mean()

                    # LWnet cannot be averaged yet, we need to convert to
                    # LWdown after temporal de-aggregation
                    lwnet = fin.variables['ATHB_S'][0][rad_mask]
                    tsurf = fin.variables['T_SO'][0,0,:][rad_mask] # zero-height soil layer
                    h = thoriz.total_seconds()
                    if it_rad_left is not False:
                        assign_time(aggr_start, it_rad_left, (h, swdir, swdif, lwnet, tsurf))
                    if it_rad_right is not False:
                        assign_time(aggr_end, it_rad_right, (h, swdir, swdif, lwnet, tsurf))

                # Save HHL
                if not thoriz: #zero horizon
                    hhl = rt.regrid_icon.loader(fin.variables['HHL'])[0,...]
                    hsurf = rt.regrid_icon.loader(fin.variables['HSURF'])[0,...]

                    # Layers in ICON NetCDF are top->bottom. Make sure that
                    # they are not mixed up in convertor, and that terrain
                    # matches lowest boundary
                    hhld = hhl[1:] - hhl[:-1]
                    assert hhld.max() < 0.
                    assert np.abs(hsurf - hhl[-1]).max() < 0.1 #10 cm
                    hhl_d[tcycle] = hhl

                if it_main is False:
                    continue

                assign_time(rt.times, it, t)
                verbose('Importing time {}, timestep {}', t, it)

                run_d[it] = tcycle

                # ICON netcdf dimension names are just generated, we cannot rely
                # on the exact names.
                for varname in cfg.icon.vars_3d:
                    v_icon = rt.regrid_icon.loader(get_3dvar(fin, varname))
                    fout.variables[varname][it] = v_icon[0,...][::-1]

                for varname in cfg.icon.vars_3dw:
                    v_icon = rt.regrid_icon.loader(get_3dvar(fin, varname))
                    fout.variables[varname][it] = v_icon[0,...][::-1]

                for varname in cfg.icon.vars_2d:
                    v_icon = rt.regrid_icon.loader(fin.variables[varname])
                    fout.variables[varname][it] = v_icon[0,...]

                # soil layers
                fout.variables['T_SO'][it] = rt.regrid_icon.loader(fin.variables['T_SO'])[0,:-1,...]

                wso = rt.regrid_icon.loader(fin.variables['W_SO'])[0,...]
                wso_gradient = (wso[1:] - wso[:-1]) / zsoil_dwso # kg/m2 agg -> kg/m3
                wso_gradient = wso_gradient / 999.0 # -> m3/m3
                qsoil = wso_gradient[wsoil_coords, :]
                fout.variables['QSOIL'][it] = qsoil

        verbose('All provided ICON files imported.')

        # write cached HHL
        v_out = fout.variables['HHL']
        for it in range(rt.nt):
            try:
                dtrun = run_d[it]
            except KeyError:
                c, h, dt = hselect.dt_from_idx(it)
                die('Missing timestep #{} ({}, cycle {}, horizon {}) in loaded data!',
                    it, dt, c, h)
            try:
                v_out[it] = hhl_d[dtrun][::-1]
            except KeyError:
                die('Base ICON run {} not loaded!', dtrun)

        # De-aggregate values
        if cfg.radiation.enabled:
            verbose('De-aggregating SW+LW radiation.')
            if aggr_start[0] is None and cfg.icon.allow_skip_first_disaggr:
                log('Extrapolating ICON disaggregated radiation values for '
                    'the first timestep!')
                aggr_start[0] = aggr_start[1]
                aggr_end[0] = aggr_end[1]
                rt.times_rad_sec[0] = 0.
            if aggr_end[-1] is None and cfg.icon.allow_skip_last_disaggr:
                log('Extrapolating ICON disaggregated radiation values for '
                    'the last timestep!')
                aggr_start[-1] = aggr_start[-2]
                aggr_end[-1] = aggr_end[-2]
                rt.times_rad_sec[-1] = rt.timestep_rad.total_seconds()
            swdown = []
            swdif = []
            lwdown = []
            for i in range(rt.nt_rad):
                if aggr_start[i] is None:
                    c, h, dt = hselect_rad_left.dt_from_idx(i)
                    die('Missing time {}, cycle {}, horizon {} for left side of radiation deaggregation!',
                        dt, c, h)
                h0, sr0, sf0, l0, t0 = aggr_start[i]

                if aggr_end[i] is None:
                    c, h, dt = hselect_rad_right.dt_from_idx(i)
                    die('Missing time {}, cycle {}, horizon {} for right side of radiation deaggregation!',
                        dt, c, h)
                h1, sr1, sf1, l1, t1 = aggr_end[i]

                # We need to de-aggregate the temporal average. Since:
                # mean_end = mean_start * horiz_start / horiz_end
                #            + mean_cur * (horiz_end - horiz_start) / horiz_end
                # mean_cur = (mean_end * horiz_end - mean_start * horiz_start)
                #            / (horiz_end - horiz_start)
                rdh = 1.0 / (h1 - h0)
                deag_swdir = (sr1 * h1 - sr0 * h0) * rdh
                deag_swdif = (sf1 * h1 - sf0 * h0) * rdh
                deag_lwnet = (l1 * h1 - l0 * h0) * rdh

                swdown.append(deag_swdir+deag_swdif)
                swdif.append(deag_swdif)

                # The surface temperature is at inst. times, we need interval
                # centres so we take consecutive means
                mean_tsurf = (t0 + t1) * 0.5

                # Now we need to convert the net radiation to downward radiation
                # using surface temperature and emissivity:
                # Lnet = Ldown - Lup = Ldown*emis - emis*sigma*T^4
                # Ldown = (Lnet + emis*sigma*T^4) / emis
                lwdown.append(((deag_lwnet + mean_tsurf**4 * emis_sigma) * emis_r).mean())

            rt.rad_swdown = swdown
            rt.rad_swdiff = swdif
            rt.rad_lwdown = lwdown
            rt.has_rad_diffuse = True

        log('ICON import finished.')

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing horizontal interpolation')

        verbose('Preparing output file')
        with netCDF4.Dataset(rt.paths.intermediate.import_data) as fin:
            # Create dimensions
            for d in ['time', 'z_meteo', 'zw_meteo', 'zsoil']:
                fout.createDimension(d, len(fin.dimensions[d]))
            fout.createDimension('x', rt.nx)
            fout.createDimension('y', rt.ny)

            # Create variables
            vars = (cfg.icon.vars_2d + cfg.icon.vars_3d + cfg.icon.vars_3dw
                    + ['HHL'] + cfg.icon.vars_soil)
            for varname in vars:
                v_icon = fin.variables[varname]
                if v_icon.dimensions[-1] != 'ncells':
                    raise RuntimeError('Unexpected dimensions for '
                            'variable {}: {}!'.format(varname,
                                v_icon.dimensions))
                fout.createVariable(varname, 'f4', v_icon.dimensions[:-1]
                        + ('y', 'x'))

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # regular vars
                for varname in vars:
                    v_icon = fin.variables[varname]
                    v_out = fout.variables[varname]
                    v_out[it] = rt.regrid_icon.regrid(v_icon[it])

    def interpolate_vert(self, fout, *args, **kwargs):
        verbose_dstat = log_dstat_on if cfg.verbosity >= 2 else log_dstat_off

        log('Performing vertical interpolation')

        with netCDF4.Dataset(rt.paths.intermediate.hinterp) as fin:
            verbose('Preparing output file')
            for dimname in ['time', 'y', 'x']:
                ensure_dimension(fout, dimname, len(fin.dimensions[dimname]))
            ensure_dimension(fout, 'z', rt.nz)
            ensure_dimension(fout, 'zw', rt.nz-1)
            ensure_dimension(fout, 'zsoil', rt.nz_soil)

            fout.createVariable('init_atmosphere_qv', 'f4', ('time', 'z', 'y', 'x'))
            fout.createVariable('init_atmosphere_pt', 'f4', ('time', 'z', 'y', 'x'))
            fout.createVariable('init_atmosphere_u', 'f4', ('time', 'z', 'y', 'x'))
            fout.createVariable('init_atmosphere_v', 'f4', ('time', 'z', 'y', 'x'))
            fout.createVariable('init_atmosphere_w', 'f4', ('time', 'zw', 'y', 'x'))
            fout.createVariable('palm_hydrostatic_pressure', 'f4', ('time', 'z'))
            fout.createVariable('palm_hydrostatic_pressure_stag', 'f4', ('time', 'zw'))
            fout.createVariable('init_soil_t', 'f4', ('time', 'zsoil', 'y', 'x'))
            fout.createVariable('init_soil_m', 'f4', ('time', 'zsoil', 'y', 'x'))
            fout.createVariable('zsoil', 'f4', ('zsoil',))
            fout.createVariable('z', 'f4', ('z',))
            fout.createVariable('zw', 'f4', ('zw',))

            fout.variables['z'][:] = rt.z_levels
            fout.variables['zw'][:] = rt.z_levels_stag
            fout.variables['zsoil'][:] = rt.z_soil_levels #depths of centers of soil layers

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                hhl = fin.variables['HHL'][it,:,:,:]
                gp_w = hhl * g
                iconterr = hhl[0]

                if cfg.vinterp.terrain_smoothing:
                    verbose('Smoothing PALM terrain for the purpose of '
                            'dynamic driver with sigma={0} grid '
                            'points.', cfg.vinterp.terrain_smoothing)
                    target_terrain = ndimage.gaussian_filter(rt.terrain,
                            sigma=cfg.vinterp.terrain_smoothing, order=0)
                else:
                    target_terrain = rt.terrain

                verbose('Morphing ICON terrain ({0} ~ {1}) to PALM terrain ({2} ~ {3})',
                    iconterr.min(), iconterr.max(), target_terrain.min(), target_terrain.max())
                verbose_dstat('Terrain shift [m]', iconterr - target_terrain[:,:])

                t_u = fin.variables['T'][it,:,:,:]
                tair_surf = t_u[0,:,:]

                p = fin.variables['P'][it,:,:,:]
                p_surf = fin.variables['PS'][it,:,:]

                # Save 1-D hydrostatic pressure
                print('lev0shift', rt.origin_z - iconterr) #TODO DEBUG
                p_lev0 = PalmPhysics.barom_ptn_pres(p_surf, rt.origin_z - iconterr, tair_surf).mean()
                tsurf_ref = tair_surf.mean()
                fout.variables['palm_hydrostatic_pressure'][it,:,] = PalmPhysics.barom_ptn_pres(p_lev0, rt.z_levels, tsurf_ref)
                fout.variables['palm_hydrostatic_pressure_stag'][it,:,] = PalmPhysics.barom_ptn_pres(p_lev0, rt.z_levels_stag, tsurf_ref)

                gp_new_surf = target_terrain * g

                # Calculate transition pressure level using horizontal
                # domain-wide pressure average
                z_trans = rt.origin_z + rt.z_levels_stag[rt.canopy_top] + cfg.vinterp.transition_level
                gp_trans = z_trans * g
                p_trans = barom_pres(p_surf, gp_trans, gp_w[0,:,:], tair_surf).mean()
                verbose('Vertical stretching transition level: {} m ASL = {} Pa', z_trans, p_trans)

                # Convert the geopotentials to pressure naively using barometric equation
                p_orig_w = barom_pres(p_surf, gp_w, gp_w[0,:,:], tair_surf)

                # Mass (half) levels should be calculated from full
                # levels by halving pressure, not geopotential
                p_orig_u = (p_orig_w[:-1] + p_orig_w[1:]) * 0.5

                # Calculate terrain pressure shift ratio
                p_surf_new = barom_pres(p_surf, gp_new_surf, gp_w[0,:,:], tair_surf)
                terrain_ratio = (p_surf_new - p_trans) / (p_surf - p_trans)

                # TODO: this may be optimized by finding highest stretched level and
                # caclulating only below that, or by using numexpr
                p_str_u = (p_orig_u[:,:,:] - p_trans) * terrain_ratio + p_trans
                p_str_w = (p_orig_w[:,:,:] - p_trans) * terrain_ratio + p_trans
                del terrain_ratio

                # Stretch levels to match terrain and keep everthing above transition level
                p_new_u = np.where(p_orig_u > p_trans, p_str_u, p_orig_u)
                p_new_w = np.where(p_orig_w > p_trans, p_str_w, p_orig_w)

                # Calculate new geopotentials
                gp_new_u = barom_gp(gp_w[0,:,:], p_new_u, p_surf, tair_surf)
                gp_new_w = barom_gp(gp_w[0,:,:], p_new_w, p_surf, tair_surf)
                # Verified: gp differences in levels above p_trans
                # (~0.03) are only due to float32 precision

                # Calculate new heights
                z_w = gp_new_w * (1./g) - rt.origin_z
                z_u = gp_new_u * (1./g) - rt.origin_z

                # Report
                gpdelta = gp_new_w - gp_w
                for k in range(gp_w.shape[0]):
                    verbose_dstat('GP shift level {:3d}'.format(k), gpdelta[k])

                # Standard heights
                vinterp, vinterp_wind = get_vinterp(rt.z_levels, z_u, True, True)

                var = fin.variables['QV'][it]
                fout.variables['init_atmosphere_qv'][it,:,:,:], = vinterp(var)

                var = fin.variables['T'][it]*PalmPhysics.exner_inv(p)
                fout.variables['init_atmosphere_pt'][it,:,:,:], = vinterp(var)

                var = fin.variables['U'][it]
                fout.variables['init_atmosphere_u'][it,:,:,:], = vinterp_wind(var)

                var = fin.variables['V'][it]
                fout.variables['init_atmosphere_v'][it,:,:,:],  = vinterp_wind(var)

                del vinterp, vinterp_wind

                # Z staggered
                vinterp_wind, = get_vinterp(rt.z_levels_stag, z_w, False, True)

                var = fin.variables['W'][it] #z staggered!
                fout.variables['init_atmosphere_w'][it,:,:,:], = vinterp_wind(var)

                del vinterp_wind

                # Other vars w/o vinterp

                var = fin.variables['T_SO'][it] #soil temperature
                fout.variables['init_soil_t'][it,:,:,:] = var

                var = fin.variables['QSOIL'][it] #soil moisture
                fout.variables['init_soil_m'][it,:,:,:] = var
