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
import re
import glob
from datetime import datetime, timedelta, timezone
import numpy as np
import netCDF4
from metpy.interpolate import interpolate_1d

from palmmeteo.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from palmmeteo.logging import die, warn, log, verbose
from palmmeteo.config import cfg
from palmmeteo.runtime import rt
from palmmeteo.utils import ensure_dimension
from palmmeteo.library import QuantityCalculator, LatLonRegularGrid, verify_palm_hinterp
from .wrf_utils import BilinearRegridder, WrfPhysics
import pyproj

ax_ = np.newaxis
re_num = re.compile(r'[0-9\.]+')

class CAMSPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing CAMS data...')

        filled = [False] * rt.nt
        timesteps = [QuantityCalculator.new_timestep() for t in rt.times]
        zcoord = [None] * rt.nt

        # Process input files
        verbose('Parsing CAMS file {}', rt.paths.cams.file)
        with netCDF4.Dataset(rt.paths.cams.file, 'r') as fin:
            # Decode time and locate timesteps
            origin_time = fin.variables['time'].long_name.split(' ')[-1]
            origin_time = datetime.strptime(origin_time, '%Y%m%d')
            tflag = fin.variables['time'][:].data
            times = [origin_time + timedelta(hours=float(h)) for h in tflag]

            dts = []
            for i in range(len(times)):
                dt = times[i]
                dt = dt.replace(tzinfo=timezone.utc)
                ireq = rt.tindex(dt)
                if not 0 <= ireq < rt.nt:
                    continue
                dts.append((ireq, i))

            if not dts:
                verbose('Skipping CAMS file time: '.format(dt))

            # locate layer heights
            height = fin.variables['level'][:]

            # coordinate projection
            verbose('Loading projection and preparing regridder')
            lats = fin.variables['latitude'][:]
            lons = fin.variables['longitude'][:]
            transform = LatLonRegularGrid(lats, lons)
            palm_in_cams_y, palm_in_cams_x = transform.latlon_to_ji(rt.palm_grid_lat, rt.palm_grid_lon)
            rt.regrid_cams = BilinearRegridder(palm_in_cams_x, palm_in_cams_y, preloaded=True)
            del palm_in_cams_y, palm_in_cams_x

            if cfg.hinterp.validate:
                verbose('Validating horizontal inteprolation.')
                llats, llons = np.broadcast_arrays(lats[:,ax_], lons[ax_,:])
                verify_palm_hinterp(rt.regrid_cams,
                                    rt.regrid_cams.loader(llats)[...],
                                    rt.regrid_cams.loader(llons)[...])
                del llats, llons

            convertor = QuantityCalculator(cfg.chem_species,
                                      cfg.cams.output_var_defs, cfg.cams.preprocessors,
                                      rt.regrid_cams)

            # dimensions
            ensure_dimension(fout, 'time', rt.nt)
            ensure_dimension(fout, 'z_chem', height.size)
            ensure_dimension(fout, 'y_chem', rt.regrid_cams.ylen)
            ensure_dimension(fout, 'x_chem', rt.regrid_cams.xlen)
            chem_dims = ('time', 'z_chem', 'y_chem', 'x_chem')

            vz_out = fout.createVariable('height_chem', 'f4', chem_dims)
            vz_out.units = 'm'

            for itout, itf in dts:
                verbose('Importing timestep {} -> {}', itf, itout)
                verbose('\tProcessing CAMS time {0}.'.format(dts[itout]))
                zcoord[itout] = np.tile(height[:, ax_, ax_], (1, rt.regrid_cams.ylen, rt.regrid_cams.xlen))

                filled[itout] = convertor.load_timestep_vars(fin, itf,
                                                             timesteps[itout])

                vz_out[itout,:,:,:] = np.tile(height[:, ax_, ax_], (1, rt.regrid_cams.ylen, rt.regrid_cams.xlen))

            if not all(filled):
                die('Could not find all CAMS variables for all times.\n'
                    'Missing variables in times:\n{}',
                    '\n'.join('{}: {}'.format(dt, ', '.join(vn
                                                            for vn in (convertor.loaded_vars - tsdata)))
                              for dt, fil, tsdata in zip(rt.times, filled, timesteps)
                              if not fil))

            for i, tsdata in enumerate(timesteps):
                # Save heights
                vz_out[i, :, :, :] = zcoord[i]

                # Save computed variables
                convertor.validate_timestep(tsdata)
                for sn, v, unit, attrs in convertor.calc_timestep_species(tsdata):
                    v_out = (fout.variables[sn] if i
                             else fout.createVariable(sn, 'f4', chem_dims))
                    v_out.units = unit
                    if attrs:
                        v_out.setncatts(attrs)
                    v_out[i, :, :, :] = v

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing CAMS horizontal interpolation')
        hvars = ['height_chem'] + cfg.chem_species

        with netCDF4.Dataset(rt.paths.intermediate.import_data) as fin:
            verbose('Preparing output file')
            # Create dimensions
            for d in ['time', 'z_chem']:
                ensure_dimension(fout, d, len(fin.dimensions[d]))
            ensure_dimension(fout, 'x', rt.nx)
            ensure_dimension(fout, 'y', rt.ny)

            # Create variables
            for varname in hvars:
                v_in = fin.variables[varname]
                if v_in.dimensions[-2:] != ('y_chem', 'x_chem'):
                    raise RuntimeError('Unexpected dimensions for '
                            'variable {}: {}!'.format(varname,
                                v_in.dimensions))
                v_out = fout.createVariable(varname, 'f4', v_in.dimensions[:-2]
                        + ('y', 'x'))
                v_out.setncatts({a: v_in.getncattr(a) for a in v_in.ncattrs()})
            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # regular vars
                for varname in hvars:
                    v_in = fin.variables[varname]
                    v_out = fout.variables[varname]
                    v_out[it] = rt.regrid_cams.regrid(v_in[it])

    def interpolate_vert(self, fout, *args, **kwargs):
        log('Performing CAMS vertical interpolation')
        terrain_rel = rt.terrain_rel[ax_,:,:]

        with netCDF4.Dataset(rt.paths.intermediate.hinterp) as fin:
            agl_chem = fin.variables['height_chem']
            chem_heights = np.zeros((agl_chem.shape[1]+1,) + agl_chem.shape[2:], dtype=agl_chem.dtype)
            chem_heights[0,:,:] = -999.

            verbose('Preparing output file')
            for dimname in ['time', 'y', 'x']:
                ensure_dimension(fout, dimname, len(fin.dimensions[dimname]))
            ensure_dimension(fout, 'z', rt.nz)

            for vn in cfg.chem_species:
                v_in = fin.variables[vn]
                var = fout.createVariable(vn, 'f4', ('time', 'z', 'y', 'x'))
                var.setncatts({a: v_in.getncattr(a) for a in v_in.ncattrs()})

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # Calc CAMS layer heights
                chem_heights[1:,:,:] = agl_chem[it] + terrain_rel

                # Load all variables for the timestep
                vardata = []
                for vn in cfg.chem_species:
                    data = fin.variables[vn][it]
                    data = np.r_[data[0:1], data]
                    vardata.append(data)

                # Perform vertical interpolation on all currently loaded vars at once
                vinterp = interpolate_1d(rt.z_levels, chem_heights, *vardata,
                        return_list_always=True)
                del vardata

                for vn, vd in zip(cfg.chem_species, vinterp):
                    v = fout.variables[vn]
                    v[it] = vd

