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

import numpy as np
import netCDF4
from pyproj import Proj, transform

from .plugins import SetupPluginMixin, WritePluginMixin
from .logging import die, warn, log, verbose
from .config import cfg
from .runtime import rt
from .utils import find_free_fname, tstep, td0
from .library import PalmPhysics

ax_ = np.newaxis

class SetupPlugin(SetupPluginMixin):
    def setup_model(self, *args, **kwargs):
        log('Setting up model domain...')

        # absolute terrain needed for vertical interpolation of wrf data
        rt.terrain = rt.terrain_rel + rt.origin_z

        # print domain parameters and check ist existence in caso of setup from config
        verbose('Domain parameters:')
        verbose('nx={}, ny={}, nz={}', rt.nx, rt.ny, rt.nz)
        verbose('dx={}, dy={}, dz={}', rt.dx, rt.dy, rt.dz)
        verbose('origin_x={}, origin_y={}', rt.origin_x, rt.origin_y)
        verbose('Base of domain is in level origin_z={}', rt.origin_z)

        # centre of the domain (needed for ug,vg calculation)
        rt.xcent = rt.origin_x + rt.nx * rt.dx / 2.0
        rt.ycent = rt.origin_y + rt.ny * rt.dy / 2.0
        # WGS84 projection for transformation to lat-lon
        rt.inproj = Proj('+init='+cfg.domain.proj_palm)
        rt.lonlatproj = Proj('+init='+cfg.domain.proj_wgs84)
        rt.cent_lon, rt.cent_lat = transform(rt.inproj, rt.lonlatproj,
                rt.xcent, rt.ycent)
        verbose('xcent={}, ycent={}', rt.xcent, rt.ycent)
        verbose('cent_lon={}, cent_lat={}', rt.cent_lon, rt.cent_lat)
        # prepare target grid
        irange = rt.origin_x + rt.dx * (np.arange(rt.nx, dtype='f8') + .5)
        jrange = rt.origin_y + rt.dy * (np.arange(rt.ny, dtype='f8') + .5)
        rt.palm_grid_y, rt.palm_grid_x = np.meshgrid(jrange, irange, indexing='ij')
        rt.palm_grid_lon, rt.palm_grid_lat = transform(rt.inproj, rt.lonlatproj,
                rt.palm_grid_x, rt.palm_grid_y)

        ######################################
        # build structure of vertical layers
        # remark:
        # PALM input requires nz=ztop in PALM
        # but the output file in PALM has max z higher than z in PARIN.
        # The highest levels in PALM are wrongly initialized !!!
        #####################################
        if rt.stretching:
            if cfg.domain.dz_stretch_level < 0:
                raise ConfigError('Stretch level has to be specified for '
                    'stretching', cfg.domain, 'dz_stretch_level')
            if cfg.domain.dz_max < rt.dz:
                raise ConfigError('dz_max has to be higher or equal than than '
                        'dz (={})'.format(rt.dz), cfg.domain, 'dz_max')
        # fill out z_levels
        rt.z_levels = np.zeros(rt.nz, dtype=float)
        rt.z_levels_stag = np.zeros(rt.nz-1, dtype=float)
        dzs = rt.dz
        rt.z_levels[0] = dzs/2.0
        for i in range(rt.nz-1):
            rt.z_levels[i+1] = rt.z_levels[i] + dzs
            rt.z_levels_stag[i] = (rt.z_levels[i+1]+rt.z_levels[i])/2.0
            if rt.stretching and rt.z_levels[i+1] + dzs >= cfg.domain.dz_stretch_level:
                dzs = min(dzs * cfg.domain.dz_stretch_factor, cfg.domain.dz_max)
        rt.ztop = rt.z_levels[-1] + dzs / 2.
        verbose('z: {}', rt.z_levels)
        verbose('zw: {}', rt.z_levels_stag)

        # configure times
        rt.simulation.end_time_rad = rt.simulation.start_time + rt.simulation.length
        rt.tindex = lambda dt: tstep(dt-rt.simulation.start_time, rt.simulation.timestep)
        if rt.nested_domain:
            log('Nested domain - preparing only initialization (1 timestep).')
            rt.nt = 1
            rt.simulation.duration = td0
            rt.simulation.end_time = rt.simulation.start_time
        else:
            rt.simulation.end_time = rt.simulation.end_time_rad
            rt.nt = rt.tindex(rt.simulation.end_time) + 1

        rt.times_sec = np.arange(rt.nt) * rt.simulation.timestep.total_seconds()
        verbose('PALM meteo data extent {} - {} ({} timesteps).',
                rt.simulation.start_time, rt.simulation.end_time, rt.nt)


class WritePlugin(WritePluginMixin):
    def write_data(self, *args, **kwargs):
        log('Writing data to dynamic driver')

        fn_out = find_free_fname(rt.paths.dynamic_driver, cfg.output.overwrite)
        dtdefault = cfg.output.default_precision
        filldefault = cfg.output.fill_value


        log('Preparing dynamic driver file {}.', fn_out)
        with netCDF4.Dataset(fn_out, 'w', format='NETCDF4') as fout:

            # shorthands
            fov = fout.variables
            def mkvar(vname, dims, lod=None, units=None, dtype=dtdefault,
                    fill=filldefault, attrs_from=None):
                v = fout.createVariable(vname, dtype, dims, fill_value=fill)
                if attrs_from is None:
                    attrs = {}
                else:
                    attrs = {a: attrs_from.getncattr(a)
                                    for a in attrs_from.ncattrs()}
                if lod is not None:
                    attrs['lod'] = lod
                if units is not None:
                    attrs['units'] = units
                if attrs:
                    v.setncatts(attrs)
                return v

            # Create dimensions
            fout.createDimension('time',  rt.nt     )
            fout.createDimension('z',     rt.nz     )
            fout.createDimension('zw',    rt.nz-1   )
            fout.createDimension('zsoil', rt.nz_soil)
            fout.createDimension('x',     rt.nx     )
            fout.createDimension('xu',    rt.nx-1   )
            fout.createDimension('y',     rt.ny     )
            fout.createDimension('yv',    rt.ny-1   )

            # Create and write dimension variables
            mkvar('time',  ('time',) )[:] = rt.times_sec
            mkvar('z',     ('z',)    )[:] = rt.z_levels[:]
            mkvar('zw',    ('zw',)   )[:] = rt.z_levels_stag[:]
            mkvar('zsoil', ('zsoil',))[:] = rt.z_soil_levels[:]
            mkvar('y',     ('y',)    )[:] = rt.dy/2 + rt.dy*np.arange(rt.ny,dtype=dtdefault)
            mkvar('x',     ('x',)    )[:] = rt.dx/2 + rt.dx*np.arange(rt.nx,dtype=dtdefault)
            mkvar('yv',    ('yv',)   )[:] = rt.dy*np.arange(1,rt.ny,dtype=dtdefault)
            mkvar('xu',    ('xu',)   )[:] = rt.dx*np.arange(1,rt.nx,dtype=dtdefault)

            # Create init variables
            mkvar('init_atmosphere_pt', ('z', 'y', 'x'),     2)
            mkvar('init_atmosphere_qv', ('z', 'y', 'x'),     2)
            mkvar('init_atmosphere_u',  ('z', 'y', 'xu'),    2)
            mkvar('init_atmosphere_v',  ('z', 'yv', 'x'),    2)
            mkvar('init_atmosphere_w',  ('zw', 'y', 'x'),    2)
            mkvar('init_soil_t',        ('zsoil', 'y', 'x'), 2)
            mkvar('init_soil_m',        ('zsoil', 'y', 'x'), 2)

            # Create forcing variables
            if not rt.nested_domain:
                # surface pressure (scalar)
                mkvar('surface_forcing_surface_pressure', ('time',))

                # boundary - vertical slices from left, right, south, north, top
                mkvar('ls_forcing_left_pt',  ('time', 'z',  'y' ), 2)
                mkvar('ls_forcing_right_pt', ('time', 'z',  'y' ), 2)
                mkvar('ls_forcing_south_pt', ('time', 'z',  'x' ), 2)
                mkvar('ls_forcing_north_pt', ('time', 'z',  'x' ), 2)
                mkvar('ls_forcing_top_pt',   ('time', 'y',  'x' ), 2)

                mkvar('ls_forcing_left_qv',  ('time', 'z',  'y' ), 2)
                mkvar('ls_forcing_right_qv', ('time', 'z',  'y' ), 2)
                mkvar('ls_forcing_south_qv', ('time', 'z',  'x' ), 2)
                mkvar('ls_forcing_north_qv', ('time', 'z',  'x' ), 2)
                mkvar('ls_forcing_top_qv',   ('time', 'y',  'x' ), 2)

                mkvar('ls_forcing_left_u',   ('time', 'z',  'y' ), 2)
                mkvar('ls_forcing_right_u',  ('time', 'z',  'y' ), 2)
                mkvar('ls_forcing_south_u',  ('time', 'z',  'xu'), 2)
                mkvar('ls_forcing_north_u',  ('time', 'z',  'xu'), 2)
                mkvar('ls_forcing_top_u',    ('time', 'y',  'xu'), 2)

                mkvar('ls_forcing_left_v',   ('time', 'z',  'yv'), 2)
                mkvar('ls_forcing_right_v',  ('time', 'z',  'yv'), 2)
                mkvar('ls_forcing_south_v',  ('time', 'z',  'x' ), 2)
                mkvar('ls_forcing_north_v',  ('time', 'z',  'x' ), 2)
                mkvar('ls_forcing_top_v',    ('time', 'yv', 'x' ), 2)

                mkvar('ls_forcing_left_w',   ('time', 'zw', 'y' ), 2)
                mkvar('ls_forcing_right_w',  ('time', 'zw', 'y' ), 2)
                mkvar('ls_forcing_south_w',  ('time', 'zw', 'x' ), 2)
                mkvar('ls_forcing_north_w',  ('time', 'zw', 'x' ), 2)
                mkvar('ls_forcing_top_w',    ('time', 'y',  'x' ), 2)

            # prepare influx/outflux area sizes
            zstag_all = np.r_[0., rt.z_levels_stag, rt.ztop]
            zwidths = zstag_all[1:] - zstag_all[:-1]
            verbose('Z widths: {}', zwidths)
            areas_xb = (zwidths * rt.dy)[:,ax_]
            areas_yb = (zwidths * rt.dx)[:,ax_]
            areas_zb = rt.dx * rt.dy
            area_boundaries = (areas_xb.sum()*rt.ny*2
                    + areas_yb.sum()*rt.nx*2
                    + areas_zb*rt.nx*rt.ny)

            log('Writing values for initialization variables')
            with netCDF4.Dataset(rt.paths.vinterp) as fin:
                fiv = fin.variables

                # geostrophic wind (1D)
                if 'ls_forcing_ug' in fiv:
                    mkvar('ls_forcing_ug', ('time', 'z'))
                    mkvar('ls_forcing_vg', ('time', 'z'))

                # write values for initialization variables
                fov['init_atmosphere_pt'][:,:,:] = fiv['init_atmosphere_pt'][0, :, :, :]
                fov['init_atmosphere_qv'][:,:,:] = fiv['init_atmosphere_qv'][0, :, :, :]
                fov['init_atmosphere_u'][:,:,:] = fiv['init_atmosphere_u'][0, :, :, 1:] #TODO fix stag
                fov['init_atmosphere_v'][:,:,:] = fiv['init_atmosphere_v'][0, :, 1:, :] #TODO fix stag
                fov['init_atmosphere_w'][:,:,:] = fiv['init_atmosphere_w'][0, :, :, :]
                fov['init_soil_t'][:,:,:] = fiv['init_soil_t'][0,:,:,:]
                fov['init_soil_m'][:,:,:] = (fiv['init_soil_m'][0,:,:,:]
                        * rt.soil_moisture_adjust[ax_,:,:])

                # Write values for time dependent values
                if not rt.nested_domain:
                    for it in range(rt.nt):
                        verbose('Processing timestep {}', it)

                        # surface pressure: in PALM, surface pressure is actually level 0 (zb) pressure!!!
                        fov['surface_forcing_surface_pressure'][it] = fiv['palm_hydrostatic_pressure_stag'][it,0]

                        # boundary conditions
                        fov['ls_forcing_left_pt' ][it,:,:] = fiv['init_atmosphere_pt'][it, :, :, 0]
                        fov['ls_forcing_right_pt'][it,:,:] = fiv['init_atmosphere_pt'][it, :, :, rt.nx-1]
                        fov['ls_forcing_south_pt'][it,:,:] = fiv['init_atmosphere_pt'][it, :, 0, :]
                        fov['ls_forcing_north_pt'][it,:,:] = fiv['init_atmosphere_pt'][it, :, rt.ny-1, :]
                        fov['ls_forcing_top_pt'  ][it,:,:] = fiv['init_atmosphere_pt'][it, rt.nz-1, :, :]

                        fov['ls_forcing_left_qv' ][it,:,:] = fiv['init_atmosphere_qv'][it, :, :, 0]
                        fov['ls_forcing_right_qv'][it,:,:] = fiv['init_atmosphere_qv'][it, :, :, rt.nx-1]
                        fov['ls_forcing_south_qv'][it,:,:] = fiv['init_atmosphere_qv'][it, :, 0, :]
                        fov['ls_forcing_north_qv'][it,:,:] = fiv['init_atmosphere_qv'][it, :, rt.ny-1, :]
                        fov['ls_forcing_top_qv'  ][it,:,:] = fiv['init_atmosphere_qv'][it, rt.nz-1, :, :]

                        # Perform mass balancing for U, V, W
                        uxleft = fiv['init_atmosphere_u'][it, :, :, 0]
                        uxright = fiv['init_atmosphere_u'][it, :, :, rt.nx-1]
                        vysouth = fiv['init_atmosphere_v'][it, :, 0, :]
                        vynorth = fiv['init_atmosphere_v'][it, :, rt.ny-1, :]
                        wztop = fiv['init_atmosphere_w'][it, rt.nz-2, :, :]#nzw=nz-1
                        mass_disbalance = ((uxleft * areas_xb).sum()
                            - (uxright * areas_xb).sum()
                            + (vysouth * areas_yb).sum()
                            - (vynorth * areas_yb).sum()
                            - (wztop * areas_zb).sum())
                        mass_corr_v = mass_disbalance / area_boundaries
                        log('Mass disbalance: {0:8g} m3/s (avg = {1:8g} m/s)',
                            mass_disbalance, mass_corr_v)
                        uxleft -= mass_corr_v
                        uxright += mass_corr_v
                        vysouth -= mass_corr_v
                        vynorth += mass_corr_v
                        wztop += mass_corr_v

                        # Verify mass balance
                        if cfg.output.check_mass_balance and cfg.verbosity >= 1:
                            mass_disbalance = ((uxleft * areas_xb).sum()
                                - (uxright * areas_xb).sum()
                                + (vysouth * areas_yb).sum()
                                - (vynorth * areas_yb).sum()
                                - (wztop * areas_zb).sum())
                            mass_corr_v = mass_disbalance / area_boundaries
                            log('Mass balanced:   {0:8g} m3/s (avg = {1:8g} m/s)',
                                mass_disbalance, mass_corr_v)

                        # Write U, V, W
                        fov['ls_forcing_left_u' ][it,:,:] = uxleft
                        fov['ls_forcing_right_u'][it,:,:] = uxright
                        fov['ls_forcing_south_u'][it,:,:] = fiv['init_atmosphere_u'][it, :, 0, 1:]
                        fov['ls_forcing_north_u'][it,:,:] = fiv['init_atmosphere_u'][it, :, rt.ny-1, 1:]
                        fov['ls_forcing_top_u'  ][it,:,:] = fiv['init_atmosphere_u'][it, rt.nz-1, :, 1:]

                        fov['ls_forcing_left_v' ][it,:,:] = fiv['init_atmosphere_v'][it, :, 1:, 0]
                        fov['ls_forcing_right_v'][it,:,:] = fiv['init_atmosphere_v'][it, :, 1:, rt.nx-1]
                        fov['ls_forcing_south_v'][it,:,:] = vysouth
                        fov['ls_forcing_north_v'][it,:,:] = vynorth
                        fov['ls_forcing_top_v'  ][it,:,:] = fiv['init_atmosphere_v'][it, rt.nz-1, 1:, :]

                        fov['ls_forcing_left_w' ][it,:,:] = fiv['init_atmosphere_w'][it, :, :, 0]
                        fov['ls_forcing_right_w'][it,:,:] = fiv['init_atmosphere_w'][it, :, :, rt.nx-1]
                        fov['ls_forcing_south_w'][it,:,:] = fiv['init_atmosphere_w'][it, :, 0, :]
                        fov['ls_forcing_north_w'][it,:,:] = fiv['init_atmosphere_w'][it, :, rt.ny-1, :]
                        fov['ls_forcing_top_w'  ][it,:,:] = wztop

                        # geostrophic wind (1D)
                        if 'ls_forcing_ug' in fiv:
                            fov['ls_forcing_ug'][it] = fiv['ls_forcing_ug'][it]
                            fov['ls_forcing_vg'][it] = fiv['ls_forcing_vg'][it]

                # Write chemical boundary conds
                if cfg.chem_species:
                    log('Writing values for chemistry variables')

                    convert_to_ppmv = set()
                    for vn in cfg.chem_species:
                        vin = fiv[vn]

                        unit = vin.units
                        if (unit == cfg.chem_units.targets.kgm3
                                and not getattr(vin, 'non_gasphase', False)):
                            if not hasattr(vin, 'molar_mass'):
                                die('Variable {} needs to be converted to ppmv but it '
                                        'is missing molar mass!', vn)
                            convert_to_ppmv.add(vn)
                            unit = cfg.chem_units.targets.ppmv

                        mkvar('init_atmosphere_'+vn, ('z',), 1, unit, attrs_from=vin)
                        if not rt.nested_domain:
                            mkvar('ls_forcing_left_'+vn,  ('time','z','y'), 2, unit, attrs_from=vin)
                            mkvar('ls_forcing_right_'+vn, ('time','z','y'), 2, unit, attrs_from=vin)
                            mkvar('ls_forcing_south_'+vn, ('time','z','x'), 2, unit, attrs_from=vin)
                            mkvar('ls_forcing_north_'+vn, ('time','z','x'), 2, unit, attrs_from=vin)
                            mkvar('ls_forcing_top_'+vn,   ('time','y','x'), 2, unit, attrs_from=vin)

                    if convert_to_ppmv:
                        log('Chemical quantities converted from kg/m3 to ppmv: {}', convert_to_ppmv)

                    for it in range(1 if rt.nested_domain else rt.nt):
                        verbose('Processing timestep {}', it)

                        if convert_to_ppmv:
                            # calculate molar volume V/n = R*T/p
                            pres = fiv['palm_hydrostatic_pressure'][it,:][:,ax_,ax_]
                            t = fiv['init_atmosphere_pt'][it,:,:,:] * PalmPhysics.exner(pres)
                            mol_vol = t * PalmPhysics.R / pres # m3/mol

                        for vn in cfg.chem_species:
                            # Load timestep
                            v = fiv[vn][it,:,:,:]

                            if vn in convert_to_ppmv:
                                v *= mol_vol * (1e9 / fiv[vn].molar_mass) # kg/g*1e9 = 1e6

                            # PALM doesn't support 3D LOD=2 init for chem yet, we have
                            # to average the field
                            fov['init_atmosphere_'+vn][:] = v.mean(axis=(1,2))

                            if not rt.nested_domain:
                                fov['ls_forcing_left_' +vn][it] = v[:,:,0]
                                fov['ls_forcing_right_'+vn][it] = v[:,:,-1]
                                fov['ls_forcing_south_'+vn][it] = v[:,0,:]
                                fov['ls_forcing_north_'+vn][it] = v[:,-1,:]
                                fov['ls_forcing_top_'  +vn][it] = v[-1,:,:]

            if cfg.radiation:
                # Separate time dimension for radiation
                fout.createDimension('time_rad', rt.nt_rad)
                var = mkvar('time_rad', ('time_rad',))
                var[:] = rt.times_rad_sec

                # radiation variables
                var = mkvar('rad_sw_in', ('time_rad',), 1, 'W/m2')
                var[:] = rt.rad_swdown

                var = mkvar('rad_lw_in', ('time_rad',), 1, 'W/m2')
                var[:] = rt.rad_lwdown

                if rt.has_rad_diffuse:
                    var = mkvar('rad_sw_in_dif', ('time_rad',), 1, 'W/m2')
                    var[:] = rt.rad_swdiff

        log('Dynamic driver written successfully.')
