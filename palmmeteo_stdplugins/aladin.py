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

# trial of aladin plugin for palm_meteo

import os
import glob
from datetime import datetime, timezone, timedelta
import sys
import numpy as np
import cfgrib as cg
import pyproj
import netCDF4
import math
import metpy
metpy_version_master = int(metpy.__version__.split('.', 1)[0])
import metpy.calc as mpcalc
from metpy.interpolate import log_interpolate_1d
from metpy.units import units
from pyproj import transform
from metpy.interpolate import interpolate_1d
from palmmeteo.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from palmmeteo.logging import die, warn, log, verbose, log_output
from palmmeteo.config import cfg, ConfigError
from palmmeteo.runtime import rt
from palmmeteo.utils import ensure_dimension
import scipy.ndimage as ndimage
from .wrf_utils import BilinearRegridder, calc_ph_hybrid, \
    calc_ph_sigma, WrfPhysics
from palmmeteo.library import PalmPhysics

barom_pres = PalmPhysics.barom_lapse0_pres
barom_gp = PalmPhysics.barom_lapse0_gp
g = PalmPhysics.g

ax_ = np.newaxis

# Settings for geostrophic wind FIXME: change them to Aladin
gw_gfs_margin_deg = 5.  # smoothing area in degrees lat/lon
gw_wrf_margin_km = 10.  # smoothing area in km
# gw_alpha = .143 #GW vertical interpolation by power law
gw_alpha = 1.  # ignore wind power law, interpolate linearly

def variable2parameter():
    filter_indicator = {}
    filter_indicator['pressure'] = {'indicatorOfParameter' : 1, 'shortname' : 'pres',
                                    'typeOfLevel' : 'hybrid', 'netCDFname' : 'P'}
    filter_indicator['geopotential'] = {'indicatorOfParameter' : 6, 'shortname' : 'z',
                                        'typeOfLevel' : 'hybrid', 'netCDFname' : 'PH'}
    filter_indicator['terrain_height'] = {'indicatorOfParameter': 6, 'shortname': 'z',
                                          'typeOfLevel': 'heightAboveGround', 'netCDFname': 'HGT'}
    filter_indicator['temperature'] = {'indicatorOfParameter' : 11, 'shortname' : 't',
                                       'typeOfLevel' : 'hybrid', 'netCDFname' : 'T'}
    filter_indicator['u'] = {'indicatorOfParameter' : 33, 'shortname' : 'u',
                             'typeOfLevel' : 'hybrid', 'netCDFname' : 'U'}
    filter_indicator['v'] = {'indicatorOfParameter' : 34, 'shortname' : 'v',
                             'typeOfLevel' : 'hybrid', 'netCDFname' : 'V'}
    filter_indicator['w'] = {'indicatorOfParameter' : 40, 'shortname' : 'unknown',
                             'typeOfLevel' : 'hybrid', 'netCDFname' : 'W'}
    filter_indicator['specific humidity'] = {'indicatorOfParameter' : 51, 'shortname' : 'q',
                                             'typeOfLevel' : 'hybrid', 'netCDFname' : 'Q'}
    filter_indicator['soil moisture content depthBelowLand'] = {'indicatorOfParameter' : 86, 'shortname' : 'sm',
                                                 'typeOfLevel' : 'depthBelowLand', 'netCDFname' : 'SMOIS'}
    filter_indicator['soil moisture content heightAboveGround'] = {'indicatorOfParameter': 86, 'shortname': 'sm',
                                                 'typeOfLevel': 'heightAboveGround', 'netCDFname': 'SMOIS1'}
    #filter_indicator['surface_forcing_surface_pressure'] = {'indicatorOfParameter': 1, 'shortname': 'pres',
    #                         'typeOfLevel': 'heightAboveGround', 'netCDFname': 'PSFC'}
    filter_indicator['init_soil_t'] = {'indicatorOfParameter': 11, 'shortname': 't',
                             'typeOfLevel': 'heightAboveGround', 'netCDFname': 'TSLB'}

    # filter_indicator['thermal downward flux'] = {'indicatorOfParameter' : 153, 'shortname' : 'unknown',
    #                                              'typeOfLevel' : , 'netCDFname' : }
    # filter_indicator['solar downward flux'] = {'indicatorOfParameter' : 154, 'shortname' : 'unknown',
    #                                            'typeOfLevel' : , 'netCDFname' : }
    # filter_indicator['solar direct horizontal'] = {'indicatorOfParameter' : 159, 'shortname' : 'unknown',
    #                                                'typeOfLevel' : , 'netCDFname' : }
    # filter_indicator['total cloud cover'] = {'indicatorOfParameter' : 171, 'shortname' : 'unknown',
    #                                          'typeOfLevel' : , 'netCDFname' : }

    typeOfLevel = ['hybrid', 'depthBelowLand', 'heightAboveGround']
    ## depthBelowLand only :: sm
    ## heightAboveGround   :: {u10, q, v10, t2m, r} these have more levels
    ##                     :: sf (snow fall), t (air temperature), z (geopotential),
    ##                           unknown (GRIB_paramId 0), p3112, p3111, pres, p3064, lsm, p3067, sm
    ## hybrid              :: q, t, pres, r, u, v, z, unknown
    # Ondra email: P.S. V seznamu byl uveden i 176:~:Surface geopotential, spectral [m^2/s^2] #missing.
    # Ten podle mne neni treba. V gribu je geopotencial (parametr 6) v hladine 0, coz odpovida vysce terenu.
    return filter_indicator

def findnearest(xlon, xlat, point):
    import scipy.spatial
    from math import floor
    if len(xlon.shape) == 3:
        xlon = xlon[0]
        xlat = xlat[0]

    mytree = scipy.spatial.KDTree(list(zip(xlon.ravel(order='F'),
                                           xlat.ravel(order='F'))))
    dist, index = mytree.query(point)
    ncols = xlon.shape[0]

    y = index % ncols
    x = int(floor(index / ncols))

    # print(point, x, y)
    return (x, y)

# class AladinPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
class AladinPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def check_config(self, *args, **kwargs):
        if cfg.aladin.vertical_stretching in ['hybrid', 'sigma']:
            warn('The configuration setting aladin.vertical_stretching="{}" is '
                 'deprecated, the "universal" vertical stretching provides '
                 'smoother upper level values. For the "hybrid" or "sigma" '
                 'settings, the user must also manually verify that the '
                 'selected method matches the vertical coordinate system '
                 'used in the provided WRFOUT files.',
                 cfg.aladin.vertical_stretching)


    def import_data(self, fout, *args, **kwargs):
        log('Importing Aladin data...')
        # Process input files
        if cfg.aladin.load_gribs:
            ierr = transform_from_grib(rt.paths.aladin.grib_file_mask, rt.paths.aladin.output, cfg, **kwargs)

        log('Tranform from grid, done ...')

        verbose('Parsing Aladin files from {}', rt.paths.aladin.output)

        rt.times = [None] * rt.nt
        first = True
        # for fn in glob.glob(rt.paths.aladin.output):
        fn = rt.paths.aladin.output

        with netCDF4.Dataset(fn, "r", format='NETCDF4') as fin:
            verbose('Parsing Aladin file {}', fn)
            for ta, ts in enumerate(fin.variables['Times']):
                # Decode time and locate timestep
                # ts = fin.variables['Times'][it]
                t = datetime.utcfromtimestamp(ts)
                t = t.replace(tzinfo=timezone.utc)
                # check if time is in selected times

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
                verbose('Importing time {}, timestep {}, [fin timestep: {}]', t, it, ta)

                if first:
                    # coordinate projection
                    verbose('Loading projection and preparing regridder')
                    trans_alad = AladinCoordTransform(fin)
                    trans_alad.verify(fin)
                    palm_in_alad_y, palm_in_alad_x = trans_alad.latlon_to_ji(
                        rt.palm_grid_lat, rt.palm_grid_lon)
                    rt.regrid_alad   = BilinearRegridder(palm_in_alad_x, palm_in_alad_y, preloaded=True)
                    rt.regrid_alad_u = BilinearRegridder(palm_in_alad_x + .5, palm_in_alad_y, preloaded=True)
                    rt.regrid_alad_v = BilinearRegridder(palm_in_alad_x, palm_in_alad_y + .5, preloaded=True)
                    verbose('Checking lat, lon of slices')
                    verbose('LONGITUDE\n\t{}'.format(fin.variables['XLONG'][rt.regrid_alad.ys, rt.regrid_alad.xs]))
                    verbose('LATITUDE\n\t{}'.format(fin.variables['XLAT'][rt.regrid_alad.ys, rt.regrid_alad.xs]))
                    del palm_in_alad_y, palm_in_alad_x

                    # dimensions
                    # ensure_dimension(fout, 'time', rt.nt)
                    ensure_dimension(fout, 'z', rt.nz)
                    for orig_dim, new_dim in cfg.aladin.dimensions:
                        if new_dim == 'time':
                            ensure_dimension(fout, new_dim,rt.nt)
                        elif new_dim == 'x_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_alad.xlen)
                        elif new_dim == 'xu_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_alad_u.xlen)
                        elif new_dim == 'y_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_alad.ylen)
                        elif new_dim == 'yv_meteo':
                            ensure_dimension(fout, new_dim, rt.regrid_alad_v.ylen)
                        else:
                            ensure_dimension(fout, new_dim,
                                    len(fin.dimensions[orig_dim]))

                # 1D vars are copied as-is
                for varname in cfg.aladin.vars_1d:
                    v_alad = fin.variables[varname]
                    v_out = (fout.createVariable(varname, 'f4',
                            [cfg.aladin.dimensions[d] for d in v_alad.dimensions])
                             if first else fout.variables[varname])
                    v_out[it] = v_alad[ta]

                # for hinterp vars, only the requested region is copied
                for varname in cfg.aladin.hinterp_vars:
                    v_alad = fin.variables[varname]
                    v_out = (fout.createVariable(varname, 'f4',
                            [cfg.aladin.dimensions[d] for d in v_alad.dimensions])
                             if first else fout.variables[varname])
                    v_out[it,...] = v_alad[ta,..., rt.regrid_alad.ys, rt.regrid_alad.xs]

                # U and V (staggered coords)
                v_alad = fin.variables['U']
                v_out = (fout.createVariable('U', 'f4',
                            [cfg.aladin.dimensions[d] for d in v_alad.dimensions])
                         if first else fout.variables['U'])
                v_out[it] = v_alad[ta, ..., rt.regrid_alad_u.ys, rt.regrid_alad_u.xs]

                v_alad = fin.variables['V']
                v_out = (fout.createVariable('V', 'f4',
                            [cfg.aladin.dimensions[d] for d in v_alad.dimensions])
                         if first else fout.variables['V'])
                v_out[it] = v_alad[ta, ..., rt.regrid_alad_v.ys, rt.regrid_alad_v.xs]

                # calculated SPECHUM
                if first:
                    shvars = sorted(set(cfg.aladin.spechum_vars
                                        ).intersection(fin.variables.keys()))
                    verbose('Hydro variables in alad files: {}', ', '.join(shvars))
                v_out = (fout.createVariable('SPECHUM', 'f4',
                                ('time', 'z_meteo', 'y_meteo', 'x_meteo'))
                         if first else fout.variables['SPECHUM'])
                vdata = fin.variables[shvars[0]][ta, ..., rt.regrid_alad.ys, rt.regrid_alad.xs]
                # for vname in shvars[1:]:
                #     vdata += fin.variables[vname][0, ..., rt.regrid_alad.ys, rt.regrid_alad.xs]

                v_out[it] = vdata
                del vdata

                # calculated geostrophic wind
                ug, vg = palm_alad_gw(fin, rt.cent_lon, rt.cent_lat, rt.z_levels, ta)
                v_out = (fout.createVariable('UG', 'f4', ('time', 'z')) if first
                         else fout.variables['UG'])
                v_out[it] = ug
                v_out = (fout.createVariable('VG', 'f4', ('time', 'z')) if first
                         else fout.variables['VG'])
                v_out[it] = vg

                # soil layers
                if first:
                    if 'soil_layers_stag' in fin.variables.keys():
                        rt.z_soil_levels = np.asarray(fin.variables['soil_layers_stag'])
                        rt.nz_soil = np.asarray(fin.variables['soil_layers_stag']).size
                    else:
                        rt.z_soil_levels = []
                        rt.nz_soil = 0
                    verbose('Z soil levels: {}', rt.z_soil_levels)

                first = False

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing horizontal interpolation')

        verbose('Preparing output file')
        with netCDF4.Dataset(rt.paths.intermediate.imported) as fin:
            # Create dimensions
            for d in ['time', 'z_meteo', 'zw_meteo', 'z', 'zsoil_meteo']:
                ensure_dimension(fout, d, len(fin.dimensions[d]))
            ensure_dimension(fout, 'x', rt.nx)
            ensure_dimension(fout, 'y', rt.ny)

            # Create variables
            for varname in cfg.aladin.hinterp_vars + ['SPECHUM']:
                v_wrf = fin.variables[varname]
                if v_wrf.dimensions[-2:] != ('y_meteo', 'x_meteo'):
                    raise RuntimeError('Unexpected dimensions for '
                            'variable {}: {}!'.format(varname,
                                v_wrf.dimensions))
                fout.createVariable(varname, 'f4', v_wrf.dimensions[:-2]
                        + ('y', 'x'))
            fout.createVariable('U', 'f4', ('time', 'z_meteo', 'y', 'x'))
            fout.createVariable('V', 'f4', ('time', 'z_meteo', 'y', 'x'))
            for varname in cfg.aladin.vars_1d + ['UG', 'VG']:
                v_wrf = fin.variables[varname]
                fout.createVariable(varname, 'f4', v_wrf.dimensions)

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # regular vars
                for varname in cfg.aladin.hinterp_vars + ['SPECHUM']:
                    v_wrf = fin.variables[varname]
                    v_out = fout.variables[varname]
                    v_out[it] = rt.regrid_alad.regrid(v_wrf[it])

                # U and V have special treatment (unstaggering)
                fout.variables['U'][it] = rt.regrid_alad_u.regrid(
                        fin.variables['U'][it])
                fout.variables['V'][it] = rt.regrid_alad_v.regrid(
                        fin.variables['V'][it])

                # direct copy
                for varname in cfg.aladin.vars_1d + ['UG', 'VG']:
                    fout.variables[varname][it] = fin.variables[varname][it]

    def interpolate_vert(self, fout, *args, **kwargs):
        verbose_dstat = log_dstat_on if cfg.verbosity >= 2 else log_dstat_off

        log('Performing vertical interpolation')

        verbose('Preparing output file')
        with netCDF4.Dataset(rt.paths.intermediate.hinterp) as fin:
            for dimname in ['time', 'y', 'x', 'zsoil_meteo']:
                ensure_dimension(fout, dimname, len(fin.dimensions[dimname]))
            ensure_dimension(fout, 'z', rt.nz)
            ensure_dimension(fout, 'zw', rt.nz - 1)
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

            # DEBUG:
            fout.createDimension('z_meteo', len(fin.dimensions['z_meteo']))
            fout.createDimension('zw_meteo', len(fin.dimensions['zw_meteo']))
            fout.createVariable('height', 'f4', ('time', 'z_meteo', 'y', 'x'))
            fout.createVariable('heightw', 'f4', ('time', 'zw_meteo', 'y', 'x'))

            fout.variables['z'][:] = rt.z_levels
            fout.variables['zw'][:] = rt.z_levels_stag
            fout.variables['zsoil'][:] = rt.z_soil_levels #depths of centers of soil layers

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # Use hybrid ETA levels in WRF and stretch them so that the WRF terrain
                # matches either PALM terrain or flat terrain at requested height
                gp_w = fin.variables['PH'][it,:,:,:] # FIXME, does not have PHB + fin.variables['PHB'][it,:,:,:]
                terr = gp_w[0]*(1./g)

                if cfg.vinterp.terrain_smoothing:
                    verbose('Smoothing PALM terrain for the purpose of '
                            'dynamic driver with sigma={0} grid '
                            'points.', cfg.vinterp.terrain_smoothing)
                    target_terrain = ndimage.gaussian_filter(rt.terrain,
                            sigma=cfg.vinterp.terrain_smoothing, order=0)
                else:
                    target_terrain = rt.terrain

                verbose('Morphing ALADIN terrain ({0} ~ {1}) to PALM terrain ({2} ~ {3})',
                    terr.min(), terr.max(), target_terrain.min(), target_terrain.max())
                verbose_dstat('Terrain shift [m]', terr - target_terrain[:,:])

                # Load real temperature
                t_u = aladin_t(fin, it)
                tair_surf = t_u[0, :, :]

                # Load original dry air column pressure
                mu = fin.variables['SPECHUM'][it, :, :] # + fin.variables['MU'][it, :, :]
                # p_top = np.mean(fin.variables['P'][it, -1, :, :])
                # print('p top : {}'.format(p_top))
                # p_surf = mu + p_top
                # p_top = fin.variables['P'][it, -1, :, :]
                p_surf = fin.variables['P'][it, 0, :, :]

                # Save 1-D hydrostatic pressure
                print('lev0shift', rt.origin_z - terr) #TODO DEBUG
                p_lev0 = PalmPhysics.barom_ptn_pres(p_surf, rt.origin_z - terr, tair_surf).mean()
                tsurf_ref = tair_surf.mean()
                fout.variables['palm_hydrostatic_pressure'][it,:,] = PalmPhysics.barom_ptn_pres(p_lev0, rt.z_levels, tsurf_ref)
                fout.variables['palm_hydrostatic_pressure_stag'][it,:,] = PalmPhysics.barom_ptn_pres(p_lev0, rt.z_levels_stag, tsurf_ref)

                gp_new_surf = target_terrain * g

                if cfg.aladin.vertical_stretching == 'universal':
                    # Calculate transition pressure level using horizontal
                    # domain-wide pressure average
                    gp_trans = (rt.origin_z + cfg.aladin.transition_level) * g
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
                    p_str_w = (p_orig_w[:,:,:] - p_trans) * terrain_ratio + p_trans  # is {z + } 1 higher that p_str_u
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
                    # TODO: test it and complete also this option
                    # Sigma or hybrid
                    # Shift column pressure so that it matches PALM terrain
                    # mu2 = barom_pres(p_surf, gp_new_surf, gp_w[0,:,:], tair_surf) - p_top

                    # Calculate original and shifted 3D dry air pressure
                    # if cfg.aladin.vertical_stretching == 'hybrid':
                    #     p_orig_w, p_orig_u = calc_ph_hybrid(fin, it, mu)
                    #     p_new_w, p_new_u = calc_ph_hybrid(fin, it, mu2)
                    # else:
                    #     p_orig_w, p_orig_u = calc_ph_sigma(fin, it, mu)
                    #     p_new_w, p_new_u = calc_ph_sigma(fin, it, mu2)

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

                #DEBUG:
                fout.variables['height'][it] = height[:,:,:]
                fout.variables['heightw'][it] = heightw[1:,:,:]

                var = lpad(fin.variables['SPECHUM'][it])
                fout.variables['init_atmosphere_qv'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['T'][it]) #from perturbation pt to standard
                fout.variables['init_atmosphere_pt'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['U'][it])
                fout.variables['init_atmosphere_u'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['V'][it])
                fout.variables['init_atmosphere_v'][it,:,:,:]  = interpolate_1d(rt.z_levels, height, var)

                var = lpad(fin.variables['W'][it]) #z staggered!
                fout.variables['init_atmosphere_w'][it,:,:,:] = interpolate_1d(rt.z_levels_stag, heightw, var)

                #var = fin.variables['PSFC'][it]
                #fout.variables['surface_forcing_surface_pressure'][it,:,:] = p_surf_new

                var = fin.variables['TSLB'][it] #soil temperature
                fout.variables['init_soil_t'][it,:,:,:] = var

                var = fin.variables['SMOIS'][it] #soil moisture
                fout.variables['init_soil_m'][it,:,:,:] = var

                var = fin.variables['UG'][it]
                fout.variables['ls_forcing_ug'][it,:] = var

                var = fin.variables['VG'][it]
                fout.variables['ls_forcing_vg'][it,:] = var


class AladinCoordTransform(object):
    'Coordinate transformer'
    def __init__(self,ncf):
        verbose('Coord transform')
        attr = lambda a: getattr(ncf, a)

        # Define grids

        # see http://www.pkrc.net/wrf-lambert.html
        latlon_sphere = pyproj.Proj(proj='latlong',
                                    ellps='WGS84',
                                    towgs84='0,0,0', no_defs=True)

        lambert_grid = pyproj.Proj(proj='lcc',
                                   lat_1=attr('lat_1'),
                                   lat_2=attr('lat_2'),
                                   lat_0=attr('lat_0'),
                                   lon_0=attr('lon_0'),
                                   a=attr('a_sphere'), b=attr('b_sphere'),
                                   towgs84='0,0,0', no_defs=True)

        # resoltion in m
        self.dx = dx = attr('DX')
        self.dy = dy = attr('DY')

        # number of mass grid points
        self.nx = nx = attr('nx')
        self.ny = ny = attr('ny')

        # distance between centers of mass grid points at edges
        extent_x = (nx - 1) * dx
        extent_y = (ny - 1) * dy

        # grid first grid point
        # ? left lower corner?
        i0_x, j0_y = pyproj.transform(latlon_sphere, lambert_grid,
                attr('lon_first'),
                attr('lat_first'))

        # grid center in lambert
        #center_x, center_y = pyproj.transform(latlon_wgs84, lambert_grid,
        #    attr('CEN_LON'), attr('CEN_LAT'))

        # grid origin coordinates in lambert
        #i0_x = center_x - extent_x*.5
        #j0_y = center_y - extent_y*.5

        # Define fast transformation methods

        def latlon_to_ji(lat, lon):
            x, y = pyproj.transform(latlon_sphere, lambert_grid,
                    lon, lat)
            return (y-j0_y)/dy, (x-i0_x)/dx
        self.latlon_to_ji = latlon_to_ji

        def ji_to_latlon(j, i):
            lon, lat = pyproj.transform(lambert_grid, latlon_sphere,
                i*dx+i0_x, j*dy+j0_y)
            return lat, lon
        self.ji_to_latlon = ji_to_latlon

    def verify(self, ncf):
        verbose('verifing Coords')
        lat = np.asarray(ncf.variables['XLAT'])
        lon = np.asarray(ncf.variables['XLONG'])
        j, i = np.mgrid[0:self.ny, 0:self.nx]

        jj, ii = self.latlon_to_ji(lat, lon)
        d = np.hypot(jj-j, ii-i)
        verbose('error for ll->ji: max {0} m, avg {1} m.'.format(d.max(), d.mean()))

        llat, llon = self.ji_to_latlon(j, i)
        d = np.hypot(llat - lat, llon - lon)
        verbose('error for ji->ll: max {0} deg, avg {1} deg.'.format(d.max(), d.mean()))

        # lat = ncf.variables['XLAT_U'][0]
        # lon = ncf.variables['XLONG_U'][0]
        # j, i = np.mgrid[0:self.ny, 0:self.nx+1]
        # jj, ii = self.latlon_to_ji(lat, lon)
        # ii = ii + .5
        # d = np.hypot(jj-j, ii-i)
        # print('error for U-staggered ll->ji: max {0} m, avg {1} m.'.format(d.max(), d.mean()))

class AladinRadPlugin(ImportPluginMixin):
    def import_data(self, *args, **kwargs):
        log('Importing Aladin radiation data...')
        verbose('Parsing Aladin radiation files from {}', rt.paths.aladin.grib_file_mask)

        rad_data = []

        # TODO: insert specification into yaml config file if possible
        grib_filter = {}
        grib_filter['solar downward flux'] = {
            'indicatorOfParameter': 154, 'shortname': 'unknown',
            'typeOfLevel': 'heightAboveGround', 'netCDFname': 'Nan'}
        grib_filter['solar direct horizontal'] = {
            'indicatorOfParameter': 159, 'shortname': 'unknown',
            'typeOfLevel': 'heightAboveGround', 'netCDFname': 'Nan1'}
        # grib_filter['Accumulated solar upward flux at the surface [J/m^2]'] = \
        #     {'indicatorOfParameter': 152, 'shortname': 'unknown',
        #      'typeOfLevel': 'heightAboveGround', 'netCDFname': 'Nan1'}
        # grib_filter['Accumulated thermal downward flux at the surface [J/m^2]'] =    \
        #     {'indicatorOfParameter': 153, 'shortname': 'unknown',
        #      'typeOfLevel': 'heightAboveGround', 'netCDFname': 'Nan2'}
        grib_filter['Accumulated thermal downward flux at the surface [J/m^2]'] = {
            'indicatorOfParameter': 153, 'shortname': 'unknown',
            'typeOfLevel': 'heightAboveGround', 'netCDFname': 'Nan1'}
        # grib_filter['total cloud cover'] = {'indicatorOfParameter' : 171, 'shortname' : 'unknown',
        #                                          'typeOfLevel' : 'heightAboveGround', 'netCDFname' : 'Nan2'}
        origin_time = False
        for fn in glob.glob(rt.paths.aladin.grib_file_mask, recursive=True):
            verbose('Parsing Aladin radiation file {}', fn)
            if 'soil_depth' in fn:
                verbose('soil depth file, skip')
                continue
            first = True
            for grib_var in grib_filter.keys():
                # loop over grib vars
                filter_keys = {}
                filter_keys['typeOfLevel'] = grib_filter[grib_var]['typeOfLevel']
                filter_keys['indicatorOfParameter'] = grib_filter[grib_var]['indicatorOfParameter']
                shortname = grib_filter[grib_var]['shortname']
                # try:
                ds = cg.open_dataset(fn,
                             backend_kwargs=dict(filter_by_keys
                             =filter_keys))
                # except:
                #     print('error in opening data set of grib file {}'.format(fn))


                if first:
                    first = False
                    dt64 = ds['valid_time']
                    # attr = ds[shortname].attrs
                    ts = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
                    t = datetime.utcfromtimestamp(ts)
                    t = t.replace(tzinfo=timezone.utc)
                    entry = [t]
                    time_ok = True
                    if not (rt.simulation.start_time <= t <= rt.simulation.end_time_rad):
                        if (rt.simulation.start_time - timedelta(hours=1) == t ):
                            verbose('Loading time 1h before origin_time: {}', t)
                            time_ok = False
                            origin_time = True
                        else:
                            verbose('Time {} is out of range - skipping', t)
                            time_ok = False
                            origin_time = False
                            break

                    verbose('Importing radiation for time {}', t)

                if not rad_data:
                    verbose('Building list of indices for radiation smoothig.')
                    # Find mask using PALM projection

                    lons = np.asarray(np.mod(ds['longitude']+180,360)-180)
                    lats = np.asarray(ds['latitude'])
                    xs, ys = transform(rt.lonlatproj, rt.inproj, lons, lats)
                    # TODO: improve - change to circle
                    mask = (np.abs(xs - rt.xcent) <= cfg.aladin.radiation_smoothing_distance
                            ) & (np.abs(ys - rt.ycent) <= cfg.aladin.radiation_smoothing_distance)
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
                    # mask = ~mask[yfrom:yto,xfrom:xto]


                    # load radiation

                arr = np.asarray(ds[shortname])[mask]
                entry.append(arr.mean())

            if origin_time:
                ot_time, rad_1_or, rad_2_or, rad_3_or = entry
                verbose('Origin rad 1: {}, rad 2: {}, rad 3: {}', rad_1_or, rad_2_or, rad_3_or)
                origin_time = False
            if time_ok:
                rad_data.append(entry)

        verbose('Processing loaded radiation values')
        rad_data.sort()
        rad_times, rad_1, rad_2, rad_3 = zip(*rad_data) #unzip
        rad_times, rad_1, rad_2, rad_3 = list(rad_times), list(rad_1), list(rad_2), list(rad_3)
        # rad_times, rad_swdown, rad_lwdown, rad_swdiff = zip(*rad_data)  # unzip
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
        rad_1_c, rad_2_c, rad_3_c = [], [], []
        for i in range(rt.nt_rad):
            verbose('Time: {}, rad 1: {}, rad 2: {}, rad 3: {}', rt.times_rad[i], rad_1[i], rad_2[i], rad_3[i])
            if i == 0:
                rad_1_c.append((rad_1[i] - rad_1_or) / rt.timestep_rad.total_seconds())
                rad_2_c.append((rad_2[i] - rad_2_or) / rt.timestep_rad.total_seconds())
                rad_3_c.append((rad_3[i] - rad_3_or) / rt.timestep_rad.total_seconds())
            else:
                rad_1_c.append((rad_1[i] - rad_1[i-1]) / rt.timestep_rad.total_seconds())
                rad_2_c.append((rad_2[i] - rad_2[i-1]) / rt.timestep_rad.total_seconds())
                rad_3_c.append((rad_3[i] - rad_3[i-1]) / rt.timestep_rad.total_seconds())

        rt.times_rad_sec = np.arange(rt.nt_rad) * rt.timestep_rad.total_seconds()
        verbose('Using detected radiation timestep {} with {} times.',
                rt.timestep_rad, rt.nt_rad)
        # Improve overlaps
        for i in range(rt.nt_rad):
            if rt.times_rad[i].hour in [7, 13, 19, 1]:
                if i > 0:
                    rad_1_c[i] = (rad_1_c[i-1] + rad_1_c[i+1]) / 2.0
                    rad_2_c[i] = (rad_2_c[i-1] + rad_2_c[i+1]) / 2.0
                    rad_3_c[i] = (rad_3_c[i-1] + rad_3_c[i+1]) / 2.0
                else:
                    rad_1_c[i] = rad_1_c[i + 1]
                    rad_2_c[i] = rad_2_c[i + 1]
                    rad_3_c[i] = rad_3_c[i + 1]
            verbose('Time: {}, rad 1: {}, rad 2: {}, rad 3: {}', rt.times_rad[i], rad_1_c[i], rad_2_c[i], rad_3_c[i])

        # Store loaded data
        # TODO: move to netCDF (opened once among plugins)
        rt.rad_swdown = rad_1_c
        rt.rad_lwdown = rad_3_c
        rt.rad_swdiff = rad_2_c

        rt.has_rad_diffuse = True

        # rt.rad_swdown = list(rad_swdown)
        # rt.rad_lwdown = list(rad_lwdown)
        # rt.rad_swdiff = list(rad_swdiff)

class BilinearRegridder(object):
    '''Bilinear regridder for multidimensional data.

    By standard, the last two dimensions are always Y,X in that order.
    '''
    def __init__(self, projected_x, projected_y, preloaded=False):
        projected_x = np.asanyarray(projected_x)
        projected_y = np.asanyarray(projected_y)
        self.shape = projected_x.shape
        self.rank = len(self.shape)
        assert self.shape == projected_y.shape

        y0 = np.floor(projected_y)
        yd = projected_y - y0
        ydd = 1. - yd
        self.y0 = y0.astype('i8')

        x0 = np.floor(projected_x)
        xd = projected_x - x0
        xdd = 1. - xd
        self.x0 = x0.astype('i8')

        if preloaded:
            # Prepare slices for preloading from NetCDF files (for cases where
            # the range of loaded Y, X coordinates is much less than total
            # size. The regrid method then expects preloaded data.

            ybase = self.y0.min()
            ytop = self.y0.max()+2
            assert 0 <= ybase <= ytop, (0, ybase, ytop)
            self.ys = slice(ybase, ytop)
            self.y0 -= ybase
            self.ylen = ytop - ybase

            xbase = self.x0.min()
            xtop = self.x0.max()+2
            assert 0 <= xbase <= xtop, (0, xbase, xtop)
            self.xs = slice(xbase, xtop)
            self.x0 -= xbase
            self.xlen = xtop - xbase

        self.y1 = self.y0 + 1
        self.x1 = self.x0 + 1

        self.weights = np.array([
            ydd * xdd, #wy0x0
            ydd * xd , #wy0x1
            yd  * xdd, #wy1x0
            yd  * xd , #wy1x1
            ])

    def regrid(self, data):
        # data may contain additional dimensions (before Y,X)
        dshape = data.shape[:-2]
        drank = len(dshape)

        # Prepare array for selected data
        sel_shape = (4,) + dshape + self.shape
        selection = np.empty(sel_shape, dtype=data.dtype)

        selection[0, ...] = data[..., self.y0, self.x0]
        selection[1, ...] = data[..., self.y0, self.x1]
        selection[2, ...] = data[..., self.y1, self.x0]
        selection[3, ...] = data[..., self.y1, self.x1]

        # Slice weights to match the extra dimensions
        wslice = ((slice(None),) +      #weights
            (ax_,) * drank +            #data minus Y,X
            (slice(None),) * self.rank) #regridded shape

        w = selection * self.weights[wslice]
        return w.sum(axis=0)

def transform_from_grib(filename,fileout, cfg, **kwargs):
    #log('Importing WRF data...')
    # data are stored in xy slices with a few data def as
    # Heighaboveground, and other as hybrid levels

    filter_indicator = variable2parameter()

    it = 0
    first = True
    with netCDF4.Dataset(fileout, 'w', format='NETCDF4') as fout:
        # for file in glob.glob(filename):
        for file in glob.glob(filename, recursive=True):
            if os.path.isdir(file):
                verbose('{} si dir, skip'.format(file))
                continue
            elif os.path.isfile(file):
                verbose('{} is file, take it'.format(file))
                if 'soil_depth' in file:
                    verbose('soil depth file, skip')
                    continue
            log('loading grib file: {}'.format(file))

            filter_keys = {}
            key0 = list(filter_indicator.keys())[0]
            filter_keys['typeOfLevel'] = filter_indicator[key0]['typeOfLevel']
            filter_keys['indicatorOfParameter'] = filter_indicator[key0]['indicatorOfParameter']
            shortname0 = filter_indicator[key0]['shortname']
            # print('\t ',shortname0, filter_keys)
            # get attributes
            try:
                ds = cg.open_dataset(file,
                    backend_kwargs=dict(filter_by_keys
                    =filter_keys))
            except:
                die('error in opening data set of grib file')
                exit(1)

            if first:
                # cut original domain into smaller / around the center
                attr = ds[shortname0].attrs
                from pyproj import Proj, transform, CRS, Transformer

                transformer_2_latlon = Transformer.from_crs(CRS.from_user_input(cfg.srid_palm),
                                                            CRS.from_user_input(cfg.srid_wgs84))

                lat_orig = np.asarray(ds['latitude'])
                lon_orig = np.asarray(np.mod(ds['longitude'][:, :] + 180, 360) - 180)

                # From static driver
                # xcent = rt.origin_x + rt.nx * rt.dx / 2.0
                # ycent = rt.origin_y + rt.ny * rt.dy / 2.0
                # lat_cent, lon_cent = transformer_2_latlon.transform(xcent, ycent)
                #
                # print('lat lon, cen', lat_cent, lon_cent)
                # j_cent, i_cent = findnearest(lon_orig, lat_orig, [lon_cent, lat_cent])
                # print(lon_orig[j_cent, i_cent], lat_orig[j_cent, i_cent])
                # print('centers', j_cent, i_cent)
                # print(lon_orig[i_cent, j_cent], lat_orig[i_cent, j_cent])
                # i_cent, j_cent = findnearest(lon_orig, lat_orig, [lon_cent, lat_cent])
                # print(lon_orig[j_cent, i_cent], lat_orig[j_cent, i_cent])
                # print(lon_orig[i_cent, j_cent], lat_orig[i_cent, j_cent])
                # print('////')


                # xleft   = xcent - 4.0 * rt.nx * rt.dx
                # xright  = xcent + 4.0 * rt.nx * rt.dx
                # ybottom = ycent - 4.0 * rt.ny * rt.dy
                # ytop    = ycent + 4.0 * rt.ny * rt.dy
                #
                # domain_x = rt.nx * rt.dx * 5.0 # 10 times bigger (10.0 / 2.0)
                # domain_nx = max([20, int(round(domain_x / attr['GRIB_DxInMetres']))])
                # print('domain nx multiplicat', domain_nx)
                #
                # domain_y = rt.ny * rt.dy * 5.0 # 10 times bigger (10.0 / 2.0)
                # domain_ny = max([20, int(round(domain_y / attr['GRIB_DyInMetres']))])
                # print('domain ny multiplicat', domain_ny)
                #
                # # south west border
                # lat_sw, lon_sw = transformer_2_latlon.transform(xleft, ybottom)
                # # north east border
                # lat_ne, lon_ne = transformer_2_latlon.transform(xright, ytop)
                #
                # print('lat lon, sw', lat_sw, lon_sw)
                #
                # print('lat lon, ne', lat_ne, lon_ne)
                #
                # lat_mask = np.logical_and(lat_sw < lat_orig, lat_orig < lat_ne)
                # lon_mask = np.logical_and(lon_sw < lon_orig, lon_orig < lon_ne)
                #
                #
                # # j_start = int(np.argwhere(lat_mask[:,i_cent])[0]) - 40
                # # j_end   = int(np.argwhere(lat_mask[:,i_cent])[1]) + 40
                # # i_start = int(np.argwhere(lon_mask[j_cent,:])[0]) - 40
                # # i_end   = int(np.argwhere(lon_mask[j_cent,:])[1]) + 40
                #
                # j_start = j_cent - domain_ny
                # j_end   = j_cent + domain_ny
                # i_start = i_cent - domain_nx
                # i_end   = i_cent + domain_nx

                # print(j_start, j_end, i_start, i_end)
                # j_start = 0
                # j_end = lat_mask.shape[0]
                # i_start = 0
                # i_end = lat_mask.shape[1]

                j_start = 0
                j_end = lat_orig.shape[0]
                i_start = 0
                i_end = lat_orig.shape[1]

                nx = i_end - i_start
                ny = j_end - j_start

                # print('Cutted lat, lon')
                # print(lat_orig[j_start:j_end, i_start:i_end],
                #       lon_orig[j_start:j_end, i_start:i_end])
                #
                # print('Cutted lat, lon .. switched')
                # print(lat_orig[i_start:i_end, j_start:j_end],
                #       lon_orig[i_start:i_end, j_start:j_end])



                keys_levels_hy = np.asarray(ds.hybrid)
                # print(attr)
                fout.ny = ny
                fout.nx = nx
                fout.DY = dy = attr['GRIB_DyInMetres']
                fout.DX = dx = attr['GRIB_DxInMetres']
                fout.lat_1 = lat_1 = attr['GRIB_Latin2InDegrees']
                fout.lat_2 = lat_2 = attr['GRIB_Latin1InDegrees']
                fout.lon_0 = lon_0 = attr['GRIB_LoVInDegrees']
                fout.lat_0 = lat_0 = attr['GRIB_LaDInDegrees']
                fout.a_sphere = a_sphere = 6371229
                fout.b_sphere = b_sphere = 6371229
                fout.lat_first = lat_firts = attr['GRIB_latitudeOfFirstGridPointInDegrees']
                fout.lon_first = lon_first = attr['GRIB_longitudeOfFirstGridPointInDegrees']

                nz = np.size(keys_levels_hy)-1

                # create netCDF dimension
                fout.createDimension('x', nx)
                fout.createDimension('xs', nx)
                fout.createDimension('y', ny)
                fout.createDimension('ys', ny)
                fout.createDimension('z', nz)
                fout.createDimension('zs', nz)
                fout.createDimension('soil_layers_stag', 2)
                fout.createDimension('time', None)  # FIXME: unlimited?
                verbose('\t dims created')

                # create netCDF variables
                # lat, lon
                fout.createVariable('Times','f8','time')
                fout.createVariable('XLAT','f4',('y','x'))
                fout.createVariable('XLONG','f4',('y', 'x'))
                fout.createVariable('soil_layers_stag', 'f4', ('soil_layers_stag'))


                fout.variables['XLAT'][:,:] = lat = np.asarray(ds['latitude'])
                fout.variables['XLONG'][:,:] = lon = np.asarray(np.mod(ds['longitude']+180,360)-180)
                #
                # fout.variables['XLAT'][:, :]  = lat_orig[j_start:j_end, i_start:i_end]
                # fout.variables['XLONG'][:, :] = lon_orig[j_start:j_end, i_start:i_end]
                #
                # fout.lat_first = lat_firts = lat_orig[j_start, i_start]
                # fout.lon_first = lon_first = lon_orig[j_start, i_start]
                #
                # print(lat_orig[j_start, i_start], lon_orig[j_start, i_start])
                # # print(lat_orig[i_start, j_start], lon_orig[i_start, j_start])
                # print('end')
                # print(lat_orig[j_end, i_end], lon_orig[j_end, i_end])
                # # print(lat_orig[i_end, j_end], lon_orig[i_end, j_end])
                # # exit(1)

                try:
                    plon = cfg.aladin.soil_avg_point[0], cfg.aladin.soil_avg_point[1]
                    plat = cfg.aladin.soil_avg_point[0], cfg.aladin.soil_avg_point[1]
                    # FIXME: is it correct?
                    i, j = findnearest(lon_orig, lat_orig, (plon, plat))
                    ds_soil = np.asarray(cg.open_dataset(rt.paths.aladin.soil_depth)['z'])
                    ds_soil_avg = ds_soil[i,j]
                except:
                    ds_soil_avg = cfg.aladin.soil_depth_default

                fout.variables['soil_layers_stag'][:] = [0.005, ds_soil_avg/2.0]
                verbose('Soil depths are in this case {} m'.format([0.005, ds_soil_avg/2.0]))

                #  creating variables
                for key in filter_indicator.keys():
                    var_name = filter_indicator[key]['netCDFname']
                    # print('\t creating: ',var_name, ' from ', key)
                    if filter_indicator[key]['typeOfLevel'] == 'hybrid':
                        # create 3D variable
                        if var_name == 'U':
                            fout.createVariable('U', 'f4', ('time', 'z', 'y', 'xs'))
                        elif var_name == 'V':
                            fout.createVariable('V', 'f4', ('time', 'z', 'ys', 'x'))
                        elif var_name == 'W':
                            fout.createVariable('W', 'f4', ('time', 'zs', 'y', 'x'))
                        else:
                            fout.createVariable(var_name, 'f4', ('time', 'z', 'y', 'x'))
                    elif filter_indicator[key]['typeOfLevel'] == 'depthBelowLand':
                        # create 2D variable
                        fout.createVariable(var_name, 'f4', ('time', 'soil_layers_stag', 'y', 'x'))
                    elif filter_indicator[key]['typeOfLevel'] == 'heightAboveGround':
                        # TODO: create 2D variable
                        fout.createVariable(var_name, 'f4', ('time', 'y', 'x'))
                verbose('\t vars created')

                first = False

            dt64 = ds['valid_time']
            dt = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
            fout.variables['Times'][it] = dt
            ds.close()

            # print(dt)
            # t = datetime.utcfromtimestamp(dt)
            # print(t)
            # t = t.replace(tzinfo=timezone.utc)



            ### load the data
            for key in filter_indicator.keys():
                verbose('\t loading: {}, {}', key, filter_indicator[key])
                if filter_indicator[key]['typeOfLevel'] == 'hybrid':
                    filter_keys = {}
                    filter_keys['typeOfLevel'] = 'hybrid'
                    shortname = filter_indicator[key]['shortname']
                    netcdfname = filter_indicator[key]['netCDFname']
                    filter_keys['shortName'] = shortname
                    # print('\t loading: ', key, shortname, netcdfname)
                    fout.variables[netcdfname][it,...] = cg.open_dataset(file,
                                         backend_kwargs=dict(filter_by_keys
                                                             =filter_keys))[shortname][::-1,j_start:j_end,i_start:i_end][1:, ...]

                elif filter_indicator[key]['typeOfLevel'] == 'depthBelowLand':
                    filter_keys = {}
                    filter_keys['typeOfLevel'] = 'depthBelowLand'
                    shortname = filter_indicator[key]['shortname']
                    netcdfname = filter_indicator[key]['netCDFname']
                    filter_keys['shortName'] = shortname
                    # print('\t loading: ', key, shortname, netcdfname)
                    fout.variables[netcdfname][it, ...] = cg.open_dataset(file,
                                                    backend_kwargs=dict(filter_by_keys
                                                    =filter_keys))[shortname][j_start:j_end,i_start:i_end]
                elif filter_indicator[key]['typeOfLevel'] == 'heightAboveGround':
                    filter_keys = {}
                    filter_keys['typeOfLevel'] = 'heightAboveGround'
                    shortname = filter_indicator[key]['shortname']
                    netcdfname = filter_indicator[key]['netCDFname']
                    filter_keys['shortName'] = shortname
                    # print('\t loading: ', key, shortname, netcdfname)
                    if netcdfname == 'HGT':
                        # calculate terrain height from surface geopotential
                        fout.variables[netcdfname][it, ...] = cg.open_dataset(file,
                            backend_kwargs=dict(filter_by_keys
                            =filter_keys))[shortname][j_start:j_end,i_start:i_end] / 9.81
                    elif netcdfname == 'PSFC' or netcdfname == 'TSLB':
                        filter_keys['typeOfLevel'] = 'hybrid'
                        fout.variables[netcdfname][it, ...] = cg.open_dataset(file,
                            backend_kwargs=dict(filter_by_keys
                            =filter_keys))[shortname][-1, j_start:j_end,i_start:i_end]
                    elif netcdfname == 'SMOIS1':
                        fout.variables[netcdfname][it, ...] = cg.open_dataset(file,
                            backend_kwargs=dict(filter_by_keys
                            =filter_keys))[shortname][j_start:j_end,i_start:i_end]
                    else:
                        fout.variables[netcdfname][it, ...] = cg.open_dataset(file,
                            backend_kwargs=dict(filter_by_keys
                            =filter_keys))[shortname][-1, j_start:j_end,i_start:i_end]

            # transform temperature to potential temperature
            rdivcp = 0.286
            pressure_ref = 1.0e5
            fout.variables['T'][it, ...] = fout.variables['T'][it, ...]*\
                                           (pressure_ref/fout.variables['P'][it, ...])**rdivcp
            verbose('done with iteration: {}, time: {:.1f}'.format(it, dt.data))
            it += 1

    verbose('loading data done')
    ierr = 1
    return ierr

def log_dstat_on(desc, delta):
    """Calculate and log delta statistics if enabled."""
    log_output('{0} ({1:8g} ~ {2:8g}): bias = {3:8g}, MAE = {4:8g}\n'.format(
        desc, delta.min(), delta.max(), delta.mean(), np.abs(delta).mean()))

def log_dstat_off(desc, delta):
    """Do nothing (log disabled)"""
    pass

def lpad(var):
    """Pad variable in first dimension by repeating lowest layer twice"""
    return np.r_[var[0:1], var]

def palm_alad_gw(f, lon, lat, levels, tidx=0):
    '''Calculate geostrophic wind from WRF using metpy'''

    hgts, ug, vg = calcgw_wrf(f, lat, lon, levels, tidx)

    # extrapolate at the bottom
    hgts = np.r_[np.array([0.]), hgts]
    ug = np.r_[ug[0], ug]
    vg = np.r_[vg[0], vg]

    return minterp(levels, hgts, ug, vg)

def minterp(interp_heights, data_heights, u, v):
    '''Interpolate wind using power law for agl levels'''

    pdata = data_heights ** gw_alpha
    pinterp = interp_heights ** gw_alpha
    hindex = np.searchsorted(data_heights, interp_heights, side='right')
    lindex = hindex - 1
    assert lindex[0] >= 0
    assert hindex[-1] < len(data_heights)
    lbound = pdata[lindex]
    hcoef = (pinterp - lbound) / (pdata[hindex] - lbound)
    #print(data_heights)
    #print(lindex)
    #print(hcoef)
    lcoef = 1. - hcoef
    iu = u[lindex] * lcoef + u[hindex] * hcoef
    iv = v[lindex] * lcoef + v[hindex] * hcoef
    return iu, iv

def get_wrf_dims(f, lat, lon, xlat, xlong):
    '''A crude method, yet satisfactory for approximate WRF surroundings'''

    sqdist = (xlat - lat)**2 + (xlong - lon)**2
    coords = np.unravel_index(sqdist.argmin(), sqdist.shape)

    xmargin = int(math.ceil(gw_wrf_margin_km * 1000 / f.DX)) #py2 ceil produces float
    ymargin = int(math.ceil(gw_wrf_margin_km * 1000 / f.DY))
    y0, y1 = coords[0] - ymargin, coords[0] + ymargin
    x0, x1 = coords[1] - xmargin, coords[1] + xmargin
    assert 0 <= y0 < y1 < sqdist.shape[0], "Point {0} + surroundings not inside domain".format(coords[0])
    assert 0 <= x0 < x1 < sqdist.shape[1], "Point {0} + surroundings not inside domain".format(coords[1])

    return coords, (slice(y0, y1+1), slice(x0, x1+1)), (ymargin, xmargin)

def calcgw_wrf(f, lat, lon, levels, tidx=0):
    # MFDataset removes the time dimension from XLAT, XLONG
    xlat = f.variables['XLAT']
    xlslice = (0,) * (len(xlat.shape)-2) + (slice(None), slice(None))
    xlat = xlat[xlslice]
    xlong = f.variables['XLONG'][xlslice]

    (iy, ix), area, (iby, ibx) = get_wrf_dims(f, lat, lon, xlat, xlong)
    areat = (tidx,) + area
    areatz = (tidx, slice(None)) + area
    #print('wrf coords', lat, lon, xlat[iy,ix], xlong[iy,ix])
    #print(xlat[area][iby,ibx], xlong[area][iby,ibx], areat)

    # # load area
    # hgt = (f.variables['PH'][areatz] + f.variables['PHB'][areatz]) / 9.81
    # hgtu = (hgt[:-1] + hgt[1:]) * .5
    # pres = f.variables['P'][areatz] + f.variables['PB'][areatz]
    # terrain = f.variables['HGT'][areat]
    #
    # # find suitable pressure levels
    # yminpres, xminpres = np.unravel_index(pres[0].argmin(), pres[0].shape)
    # pres1 = pres[0, yminpres, xminpres] - 1.
    #
    # aglpt = hgtu[:,iby,ibx] - terrain[iby,ibx]
    # pres0 = pres[np.searchsorted(aglpt, levels[-1]), iby, ibx]
    # plevels = np.arange(pres1, min(pres0, pres1)-1, -1000.)
    #
    # # interpolate wrf into pressure levels
    # phgt = log_interpolate_1d(plevels, pres, hgtu, axis=0)

    # load area
    hgt = (f.variables['PH'][areatz]) / 9.81  #FIXME: + f.variables['PHB'][areatz]) / 9.81
    hgtu = (hgt[:-1] + hgt[1:]) * .5
    pres = f.variables['P'][areatz] # FIXME: + f.variables['PB'][areatz]
    terrain = f.variables['HGT'][areat]

    # find suitable pressure levels
    yminpres, xminpres = np.unravel_index(pres[0].argmin(), pres[0].shape)
    pres1 = pres[0, yminpres, xminpres] - 1.

    aglpt = hgtu[:,iby,ibx] - terrain[iby,ibx]
    pres0 = pres[np.searchsorted(aglpt, levels[-1]), iby, ibx]
    plevels = np.arange(pres1, min(pres0, pres1)-1, -1000.)

    # interpolate wrf into pressure levels
    pres = pres[:-1,:,:]
    phgt = log_interpolate_1d(plevels, pres, hgtu, axis=0)

    # lat_lon_grid_deltas doesn't work under py2, but for WRF grid it is still
    # not very accurate, better use direct values.
    #dx, dy = mpcalc.lat_lon_grid_deltas(xlong[area], xlat[area])
    dx = f.DX * units.m
    dy = f.DY * units.m

    if metpy_version_master >= 1:
        mylat = np.deg2rad(xlat[area])
        my_geostrophic_wind = lambda sh: mpcalc.geostrophic_wind(sh, dx=dx,
                dy=dy, latitude=mylat)
    else:
        coriol = mpcalc.coriolis_parameter(np.deg2rad(xlat[area])).to('1/s')
        my_geostrophic_wind = lambda sh: mpcalc.geostrophic_wind(sh, coriol,
                dx, dy)

    # Smooth height data. Sigma=1.5 for gfs 0.5deg
    res_km = f.DX / 1000.

    ug = np.zeros(plevels.shape, 'f8')
    vg = np.zeros(plevels.shape, 'f8')
    for i in range(len(plevels)):
        sh = ndimage.gaussian_filter(phgt[i,:,:], sigma=1.5*50/res_km, order=0)
        ugl, vgl = my_geostrophic_wind(sh * units.m)
        ug[i] = ugl[iby, ibx].magnitude
        vg[i] = vgl[iby, ibx].magnitude

    return phgt[:,iby,ibx], ug, vg

#TODO suspicious - check!!!
def aladin_t(f, it):
    p = f.variables['P'][it,:,:,:] # from wrf_utils+ f.variables['PB'][it,:,:,:]
    return (f.variables['T'][it,:,:,:] + WrfPhysics.base_temp) * WrfPhysics.exner(p)


# for dump purpouse, will be deleted
def description_key():
    describtion = {}
    # level 0
    describtion['p3067'] = 'mixed layer depth'
    describtion['paramId_0'] = 'original GRIB paramId_0'
    describtion['p311'] = 'Net short-wave radiation flux'
    describtion['z'] = 'Geopotential'
    describtion['lsm'] = 'land sea mark'
    describtion['pres'] = 'pressure'
    describtion['sm'] = 'soil moisture'
    describtion['sf'] = 'snow fall water equivalent'
    describtion['p3064'] = 'snow fall rate water equivalent'
    describtion['t'] = 'temperature'
    # level 2
    describtion['r'] = 'relative humidity'
    describtion['q'] = 'specific humidity'
    describtion['2t'] = '2 meter temperature'   ##also t2m
    # level 10
    describtion['10u'] = '10 meter U wind component'
    describtion['10v'] = '10 meter V wind component'
    # level hybrid
    describtion['z'] = 'geopotential'
    describtion['r'] = 'relative humidity'
    describtion['q'] = 'specific humidity'
    describtion['paramId_0'] = 'original Grib paramId: 0'
    describtion['pres'] = 'pressure'
    describtion['t'] = 'temperature'
    describtion['u'] = 'U wind component'
    describtion['v'] = 'V wind component'

    filter_indicator = {}
    filter_indicator['pressure'] = {'indicatorOfParameter' : 1, 'shortname' : 'pres'}
    filter_indicator['geopotential'] = {'indicatorOfParameter' : 6, 'shortname' : 'z'}
    filter_indicator['temperature'] = {'indicatorOfParameter' : 11, 'shortname' : 't'}
    filter_indicator['u'] = {'indicatorOfParameter' : 33, 'shortname' : 'u'}
    filter_indicator['v'] = {'indicatorOfParameter' : 34, 'shortname' : 'v'}
    filter_indicator['w'] = {'indicatorOfParameter' : 40, 'shortname' : 'unknown'}
    filter_indicator['specific humidity'] = {'indicatorOfParameter' : 51, 'shortname' : 'q'}
    filter_indicator['soil moisture content'] = {'indicatorOfParameter' : 86, 'shortname' : 'sm'}
    filter_indicator['thermal downward flux'] = {'indicatorOfParameter' : 153, 'shortname' : 'unknown'}
    filter_indicator['solar downward flux'] = {'indicatorOfParameter' : 154, 'shortname' : 'unknown'}
    filter_indicator['solar direct horizontal'] = {'indicatorOfParameter' : 159, 'shortname' : 'unknown'}
    filter_indicator['total cloud cover'] = {'indicatorOfParameter' : 171, 'shortname' : 'unknown'}

    typeOfLevel = {'depthBelowLand', 'heightAboveGround', 'hybrid'}

    path = '/home/bures/palm/palm_meteo/data/Aladin_trial/alad_grib4camx_19890101_01.grb'
    filter_keys = {}
    # filter_keys['typeOfLevel'] = 'hybrid'
    for key in filter_indicator.keys():
        filter_keys['indicatorOfParameter'] = filter_indicator[key]['indicatorOfParameter']
        var_name = filter_indicator[key]['shortname']
        ds = cg.open_dataset(path,
                             backend_kwargs=dict(filter_by_keys
                                                 =filter_keys))
        verbose(var_name, ds[var_name].GRIB_DxInMetres)
