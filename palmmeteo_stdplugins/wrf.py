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

import os
import glob
from datetime import datetime, timezone
import numpy as np
import scipy.ndimage as ndimage
import netCDF4
from pyproj import transform
from metpy.interpolate import interpolate_1d

from palmmeteo.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from palmmeteo.logging import die, warn, log, verbose, log_output
from palmmeteo.config import cfg, ConfigError
from palmmeteo.runtime import rt
from palmmeteo.utils import ensure_dimension
from .wrf_utils import WRFCoordTransform, BilinearRegridder, calc_ph_hybrid, \
    calc_ph_sigma, wrf_t, palm_wrf_gw, WrfPhysics
from palmmeteo.library import PalmPhysics, verify_palm_hinterp

barom_pres = PalmPhysics.barom_lapse0_pres
barom_gp = PalmPhysics.barom_lapse0_gp
g = PalmPhysics.g


def lpad(var):
    """Pad variable in first dimension by repeating lowest layer twice"""
    return np.r_[var[0:1], var]

def log_dstat_on(desc, delta):
    """Calculate and log delta statistics if enabled."""
    log_output('{0} ({1:8g} ~ {2:8g}): bias = {3:8g}, MAE = {4:8g}\n'.format(
        desc, delta.min(), delta.max(), delta.mean(), np.abs(delta).mean()))

def log_dstat_off(desc, delta):
    """Do nothing (log disabled)"""
    pass

class WRFPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def check_config(self, *args, **kwargs):
        if cfg.wrf.vertical_stretching in ['hybrid', 'sigma']:
            warn('The configuration setting wrf.vertical_stretching="{}" is '
                 'deprecated, the "universal" vertical stretching provides '
                 'smoother upper level values. For the "hybrid" or "sigma" '
                 'settings, the user must also manually verify that the '
                 'selected method matches the vertical coordinate system '
                 'used in the provided WRFOUT files.',
                 cfg.wrf.vertical_stretching)

    def import_data(self, fout, *args, **kwargs):
        log('Importing WRF data...')

        # Process input files
        verbose('Parsing WRF files from {}', rt.paths.wrf.file_mask)
        rt.times = [None] * rt.nt
        first = True
        for fn in glob.glob(rt.paths.wrf.file_mask):
            verbose('Parsing WRF file {}', fn)
            with netCDF4.Dataset(fn) as fin:
                # Decode time and locate timestep
                ts = fin.variables['Times'][:].tobytes().decode('utf-8')
                t = datetime.strptime(ts, '%Y-%m-%d_%H:%M:%S')
                t = t.replace(tzinfo=timezone.utc)
                try:
                    it = rt.tindex(t)
                except ValueError:
                    verbose('Time {} is not within timestep intervals - skipping', t)
                    continue
                if not (0 <= it < rt.nt):
                    verbose('Time {} is out of range - skipping', t)
                    continue
                if rt.times[it] is not None:
                    die('Time {} has been already loaded!', t)
                rt.times[it] = t
                verbose('Importing time {}, timestep {}', t, it)

                if first:
                    # coordinate projection
                    verbose('Loading projection and preparing regridder')
                    rt.trans_wrf = WRFCoordTransform(fin)
                    if cfg.hinterp.validate:
                        rt.trans_wrf.verify(fin)

                    palm_in_wrf_y, palm_in_wrf_x = rt.trans_wrf.latlon_to_ji(
                                                    rt.palm_grid_lat, rt.palm_grid_lon)
                    rt.regrid_wrf = BilinearRegridder(palm_in_wrf_x, palm_in_wrf_y, preloaded=True)
                    rt.regrid_wrf_u = BilinearRegridder(palm_in_wrf_x+.5, palm_in_wrf_y, preloaded=True)
                    rt.regrid_wrf_v = BilinearRegridder(palm_in_wrf_x, palm_in_wrf_y+.5, preloaded=True)
                    del palm_in_wrf_y, palm_in_wrf_x
                    if cfg.hinterp.validate:
                        verbose('Validating horizontal inteprolation.')
                        verify_palm_hinterp(rt.regrid_wrf,
                                            rt.regrid_wrf.loader(fin.variables['XLAT' ])[0],
                                            rt.regrid_wrf.loader(fin.variables['XLONG'])[0])

                    # dimensions
                    ensure_dimension(fout, 'time', rt.nt)
                    ensure_dimension(fout, 'z', rt.nz) #final z-coord for UG, VG
                    for orig_dim, new_dim in cfg.wrf.dimensions:
                        if new_dim == 'time':
                            pass
                        elif new_dim == 'x_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_wrf.xlen)
                        elif new_dim == 'xu_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_wrf_u.xlen)
                        elif new_dim == 'y_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_wrf.ylen)
                        elif new_dim == 'yv_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_wrf_v.ylen)
                        else:
                            ensure_dimension(fout, new_dim,
                                    len(fin.dimensions[orig_dim]))

                # 1D vars are copied as-is
                for varname in cfg.wrf.vars_1d:
                    v_wrf = fin.variables[varname]
                    v_out = (fout.createVariable(varname, 'f4', 
                                [cfg.wrf.dimensions[d] for d in v_wrf.dimensions])
                            if first else fout.variables[varname])
                    v_out[it] = v_wrf[0]

                # for hinterp vars, only the requested region is copied
                for varname in cfg.wrf.hinterp_vars:
                    v_wrf = fin.variables[varname]
                    v_out = (fout.createVariable(varname, 'f4', 
                                [cfg.wrf.dimensions[d] for d in v_wrf.dimensions])
                            if first else fout.variables[varname])
                    v_out[it] = v_wrf[0,...,rt.regrid_wrf.ys,rt.regrid_wrf.xs]

                # U and V (staggered coords)
                v_wrf = fin.variables['U']
                v_out = (fout.createVariable('U', 'f4', 
                            [cfg.wrf.dimensions[d] for d in v_wrf.dimensions])
                        if first else fout.variables['U'])
                v_out[it] = v_wrf[0,...,rt.regrid_wrf_u.ys,rt.regrid_wrf_u.xs]
                v_wrf = fin.variables['V']
                v_out = (fout.createVariable('V', 'f4', 
                            [cfg.wrf.dimensions[d] for d in v_wrf.dimensions])
                        if first else fout.variables['V'])
                v_out[it] = v_wrf[0,...,rt.regrid_wrf_v.ys,rt.regrid_wrf_v.xs]

                # calculated SPECHUM
                if first:
                    shvars = sorted(set(cfg.wrf.spechum_vars
                        ).intersection(fin.variables.keys()))
                    verbose('Hydro variables in wrf files: {}', ', '.join(shvars))
                v_out = (fout.createVariable('SPECHUM', 'f4',
                                ('time', 'z_meteo', 'y_meteo', 'x_meteo'))
                            if first else fout.variables['SPECHUM'])
                vdata = fin.variables[shvars[0]][0,...,rt.regrid_wrf.ys,rt.regrid_wrf.xs]
                for vname in shvars[1:]:
                    vdata += fin.variables[vname][0,...,rt.regrid_wrf.ys,rt.regrid_wrf.xs]
                v_out[it] = vdata
                del vdata

                # calculated geostrophic wind
                ug, vg = palm_wrf_gw(fin, rt.cent_lon, rt.cent_lat, rt.z_levels, 0)
                v_out = (fout.createVariable('UG', 'f4', ('time', 'z')) if first
                    else fout.variables['UG'])
                v_out[it] = ug
                v_out = (fout.createVariable('VG', 'f4', ('time', 'z')) if first
                    else fout.variables['VG'])
                v_out[it] = vg

                # soil layers
                if first:
                    if 'ZS' in fin.variables.keys():
                        rt.z_soil_levels = fin.variables['ZS'][0].data.tolist()
                    else:
                        rt.z_soil_levels = []
                    rt.nz_soil = len(rt.z_soil_levels)
                    verbose('Z soil levels: {}', rt.z_soil_levels)

                first = False

        if not all(rt.times):
            die('Some times are missing: {}', rt.times)
        log('All WRF files imported.')

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing horizontal interpolation')

        with netCDF4.Dataset(rt.paths.intermediate.import_data) as fin:
            verbose('Preparing output file')
            # Create dimensions
            for d in ['time', 'z_meteo', 'zw_meteo', 'z', 'zsoil_meteo']:
                ensure_dimension(fout, d, len(fin.dimensions[d]))
            ensure_dimension(fout, 'x', rt.nx)
            ensure_dimension(fout, 'y', rt.ny)

            # Create variables
            for varname in cfg.wrf.hinterp_vars + ['SPECHUM']:
                v_wrf = fin.variables[varname]
                if v_wrf.dimensions[-2:] != ('y_meteo', 'x_meteo'):
                    raise RuntimeError('Unexpected dimensions for '
                            'variable {}: {}!'.format(varname,
                                v_wrf.dimensions))
                fout.createVariable(varname, 'f4', v_wrf.dimensions[:-2]
                        + ('y', 'x'))
            fout.createVariable('U', 'f4', ('time', 'z_meteo', 'y', 'x'))
            fout.createVariable('V', 'f4', ('time', 'z_meteo', 'y', 'x'))
            for varname in cfg.wrf.vars_1d + ['UG', 'VG']:
                v_wrf = fin.variables[varname]
                fout.createVariable(varname, 'f4', v_wrf.dimensions)

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # regular vars
                for varname in cfg.wrf.hinterp_vars + ['SPECHUM']:
                    v_wrf = fin.variables[varname]
                    v_out = fout.variables[varname]
                    v_out[it] = rt.regrid_wrf.regrid(v_wrf[it])

                # U and V have special treatment (unstaggering)
                fout.variables['U'][it] = rt.regrid_wrf_u.regrid(
                        fin.variables['U'][it])
                fout.variables['V'][it] = rt.regrid_wrf_v.regrid(
                        fin.variables['V'][it])

                # direct copy
                for varname in cfg.wrf.vars_1d + ['UG', 'VG']:
                    fout.variables[varname][it] = fin.variables[varname][it]

    def interpolate_vert(self, fout, *args, **kwargs):
        verbose_dstat = log_dstat_on if cfg.verbosity >= 2 else log_dstat_off

        log('Performing vertical interpolation')

        with netCDF4.Dataset(rt.paths.intermediate.hinterp) as fin:
            verbose('Preparing output file')
            for dimname in ['time', 'y', 'x', 'zsoil_meteo']:
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
            fout.createVariable('ls_forcing_ug', 'f4', ('time', 'z'))
            fout.createVariable('ls_forcing_vg', 'f4', ('time', 'z'))
            fout.createVariable('zsoil', 'f4', ('zsoil',))
            fout.createVariable('z', 'f4', ('z',))
            fout.createVariable('zw', 'f4', ('zw',))

            fout.variables['z'][:] = rt.z_levels
            fout.variables['zw'][:] = rt.z_levels_stag
            fout.variables['zsoil'][:] = rt.z_soil_levels #depths of centers of soil layers

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # Use hybrid ETA levels in WRF and stretch them so that the WRF terrain
                # matches either PALM terrain or flat terrain at requested height
                gp_w = fin.variables['PH'][it,:,:,:] + fin.variables['PHB'][it,:,:,:]
                wrfterr = gp_w[0]*(1./g) #verified: equals HGT

                if cfg.vinterp.terrain_smoothing:
                    verbose('Smoothing PALM terrain for the purpose of '
                            'dynamic driver with sigma={0} grid '
                            'points.', cfg.vinterp.terrain_smoothing)
                    target_terrain = ndimage.gaussian_filter(rt.terrain,
                            sigma=cfg.vinterp.terrain_smoothing, order=0)
                else:
                    target_terrain = rt.terrain

                verbose('Morphing WRF terrain ({0} ~ {1}) to PALM terrain ({2} ~ {3})',
                    wrfterr.min(), wrfterr.max(), target_terrain.min(), target_terrain.max())
                verbose_dstat('Terrain shift [m]', wrfterr - target_terrain[:,:])

                # Load real temperature
                t_u = wrf_t(fin, it)
                tair_surf = t_u[0,:,:]

                # Load original dry air column pressure
                mu = fin.variables['MUB'][it,:,:] + fin.variables['MU'][it,:,:]
                p_top = fin.variables['P_TOP'][it]
                p_surf = mu + p_top

                # Save 1-D hydrostatic pressure
                print('lev0shift', rt.origin_z - wrfterr) #TODO DEBUG
                p_lev0 = PalmPhysics.barom_ptn_pres(p_surf, rt.origin_z - wrfterr, tair_surf).mean()
                tsurf_ref = tair_surf.mean()
                fout.variables['palm_hydrostatic_pressure'][it,:,] = PalmPhysics.barom_ptn_pres(p_lev0, rt.z_levels, tsurf_ref)
                fout.variables['palm_hydrostatic_pressure_stag'][it,:,] = PalmPhysics.barom_ptn_pres(p_lev0, rt.z_levels_stag, tsurf_ref)

                gp_new_surf = target_terrain * g

                if cfg.wrf.vertical_stretching == 'universal':
                    # Calculate transition pressure level using horizontal
                    # domain-wide pressure average
                    gp_trans = (rt.origin_z + cfg.wrf.transition_level) * g
                    p_trans = barom_pres(p_surf, gp_trans, gp_w[0,:,:], tair_surf).mean()
                    verbose('Vertical stretching transition level: {} Pa', p_trans)

                    # Convert the geopotentials to pressure naively using barometric equation
                    p_orig_w = barom_pres(p_surf, gp_w, gp_w[0,:,:], tair_surf)

                    # Mass (half) levels should be calculated from full
                    # levels by halving pressure, not geopotential, because
                    # ZNU = (ZNW[:-1]+ZNW[1:])/2 (verified)
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
                else:
                    # Sigma or hybrid
                    # Shift column pressure so that it matches PALM terrain
                    mu2 = barom_pres(p_surf, gp_new_surf, gp_w[0,:,:], tair_surf) - p_top

                    # Calculate original and shifted 3D dry air pressure
                    if cfg.wrf.vertical_stretching == 'hybrid':
                        p_orig_w, p_orig_u = calc_ph_hybrid(fin, it, mu)
                        p_new_w, p_new_u = calc_ph_hybrid(fin, it, mu2)
                    else:
                        p_orig_w, p_orig_u = calc_ph_sigma(fin, it, mu)
                        p_new_w, p_new_u = calc_ph_sigma(fin, it, mu2)

                    t_w = np.concatenate((t_u, t_u[-1:,:,:]), axis=0) # repeat highest layer

                    # Shift 3D geopotential according to delta dry air pressure
                    gp_new_w = barom_gp(gp_w, p_new_w, p_orig_w, t_w)
                    # For half-levs, originate from gp full levs rather than less accurate gp halving
                    gp_new_u = barom_gp(gp_w[:-1,:,:], p_new_u, p_orig_w[:-1,:,:], t_u)

                # Calculate new heights
                z_w = gp_new_w * (1./g) - rt.origin_z
                z_u = gp_new_u * (1./g) - rt.origin_z

                # Report
                gpdelta = gp_new_w - gp_w
                for k in range(gp_w.shape[0]):
                    verbose_dstat('GP shift level {:3d}'.format(k), gpdelta[k])

                # Because we require levels below the lowest level from WRF, we will always
                # add one layer at zero level with repeated values from the lowest level.
                # WRF-python had some special treatment for theta in this case.
                height = np.zeros((z_u.shape[0]+1,) + z_u.shape[1:], dtype=z_u.dtype)
                height[0,:,:] = -999. #always below terrain
                height[1:,:,:] = z_u
                heightw = np.zeros((z_w.shape[0]+1,) + z_w.shape[1:], dtype=z_w.dtype)
                heightw[0,:,:] = -999. #always below terrain
                heightw[1:,:,:] = z_w

                var = lpad(fin.variables['SPECHUM'][it])
                fout.variables['init_atmosphere_qv'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['T'][it] + WrfPhysics.base_temp) #from perturbation pt to standard
                fout.variables['init_atmosphere_pt'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['U'][it])
                fout.variables['init_atmosphere_u'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['V'][it])
                fout.variables['init_atmosphere_v'][it,:,:,:]  = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['W'][it]) #z staggered!
                fout.variables['init_atmosphere_w'][it,:,:,:] = interpolate_1d(rt.z_levels_stag, heightw, var)

                var = fin.variables['TSLB'][it] #soil temperature
                fout.variables['init_soil_t'][it,:,:,:] = var

                var = fin.variables['SMOIS'][it] #soil moisture
                fout.variables['init_soil_m'][it,:,:,:] = var

                var = fin.variables['UG'][it]
                fout.variables['ls_forcing_ug'][it,:] = var

                var = fin.variables['VG'][it]
                fout.variables['ls_forcing_vg'][it,:] = var


class WRFRadPlugin(ImportPluginMixin):
    def import_data(self, *args, **kwargs):
        log('Importing WRF radiation data...')
        verbose('Parsing WRF radiation files from {}', rt.paths.wrf.rad_file_mask)

        rad_data = []
        first = True
        for fn in glob.glob(rt.paths.wrf.rad_file_mask):
            verbose('Parsing WRF radiation file {}', fn)
            with netCDF4.Dataset(fn) as fin:
                # Decode time
                ts = fin.variables['Times'][:].tobytes().decode('utf-8')
                t = datetime.strptime(ts, '%Y-%m-%d_%H:%M:%S')
                t = t.replace(tzinfo=timezone.utc)
                if not (rt.simulation.start_time <= t <= rt.simulation.end_time_rad):
                    verbose('Time {} is out of range - skipping', t)
                    continue

                verbose('Importing radiation for time {}', t)
                if not rad_data:
                    verbose('Building list of indices for radiation smoothig.')

                    # Find mask using PALM projection
                    lons = fin.variables['XLONG'][0]
                    lats = fin.variables['XLAT'][0]
                    xs, ys = transform(rt.lonlatproj, rt.inproj, lons, lats)
                    #TODO: improve - change to circle
                    mask = (np.abs(xs-rt.xcent) <= cfg.wrf.radiation_smoothing_distance
                            ) & (np.abs(ys-rt.ycent) <= cfg.wrf.radiation_smoothing_distance)
                    del lons, lats, xs, ys

                    # Detect bounding box of the mask, prepare slices for
                    # faster loading
                    xmask = np.logical_or.reduce(mask, axis=0)
                    ymask = np.logical_or.reduce(mask, axis=1)
                    xfrom = np.argmax(xmask)
                    yfrom = np.argmax(ymask)
                    xto = len(xmask) - np.argmax(xmask[::-1])
                    yto = len(ymask) - np.argmax(ymask[::-1])
                    assert not any(xmask[:xfrom]) #TODO comment out
                    assert not any(xmask[xto:])
                    assert all(xmask[xfrom:xto])
                    assert not any(ymask[:yfrom])
                    assert not any(ymask[yto:])
                    assert all(ymask[yfrom:yto])
                    mask = ~mask[yfrom:yto,xfrom:xto]

                # load radiation
                if first:
                    rad_vars = [cfg.wrf.rad_vars.sw_tot_h, cfg.wrf.rad_vars.lw_tot_h]
                    if cfg.wrf.rad_vars.sw_dif_h in fin.variables:
                        verbose('WRF file does contain {} variable, adding diffuse component.',
                                cfg.wrf.rad_vars.sw_dif_h)
                        rt.has_rad_diffuse = True
                        rad_vars.append(cfg.wrf.rad_vars.sw_dif_h)
                    else:
                        verbose('WRF file does not contain {} variable, diffuse compoment will be estimated in PALM.',
                                cfg.wrf.rad_vars.sw_dif_h)
                        rt.has_rad_diffuse = False

                entry = [t]
                for varname in rad_vars:
                    arr = fin.variables[varname][0,yfrom:yto,xfrom:xto]
                    arr.mask &= mask
                    entry.append(arr.mean())
                rad_data.append(entry)
            first = False

        verbose('Processing loaded radiation values')
        rad_data.sort()
        rad_data_uz = list(zip(*rad_data)) #unzip/transpose
        rad_times = rad_data_uz[0]

        # Determine timestep and check consistency
        rt.times_rad = list(rad_times)
        rt.nt_rad = len(rt.times_rad)
        if rt.times_rad[0] != rt.simulation.start_time:
            die('Radiation must start with start time ({}), but they start with '
                    '{}!', rt.simulation.start_time, rt.times_rad[0])
        if rt.times_rad[-1] != rt.simulation.end_time_rad:
            die('Radiation must start with end time ({}), but they end with '
                    '{}!', rt.simulation.end_time_rad, rt.times_rad[-1])
        rt.timestep_rad = rt.times_rad[1] - rt.times_rad[0]
        for i in range(1, rt.nt_rad-1):
            step = rt.times_rad[i+1] - rt.times_rad[i]
            if step != rt.timestep_rad:
                die('Time delta between steps {} and {} ({}) is different from '
                        'radiation timestep ({})!', i, i+1, step, rt.timestep_rad)
        rt.times_rad_sec = np.arange(rt.nt_rad) * rt.timestep_rad.total_seconds()
        verbose('Using detected radiation timestep {} with {} times.',
                rt.timestep_rad, rt.nt_rad)

        # Store loaded data
        # TODO: move to netCDF (opened once among plugins)
        rt.rad_swdown = list(rad_data_uz[1])
        rt.rad_lwdown = list(rad_data_uz[2])
        if rt.has_rad_diffuse:
            rt.rad_swdiff = list(rad_data_uz[3])
