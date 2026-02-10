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

import re
import glob
import datetime
import numpy as np
import netCDF4

from palmmeteo.plugins import (ImportPluginMixin, HInterpPluginMixin,
                               VInterpPluginMixin)
from palmmeteo.logging import die, warn, log, verbose
from palmmeteo.config import cfg
from palmmeteo.runtime import rt
from palmmeteo.utils import ensure_dimension
from palmmeteo.library import (QuantityCalculator, LatLonRegularGrid,
                               verify_palm_hinterp, NCDates, InputGatherer,
                               HorizonSelection)
from palmmeteo.vinterp import get_vinterp
from .wrf_utils import BilinearRegridder

ax_ = np.newaxis
re_num = re.compile(r'[0-9\.]+')
re_lev = re.compile(r'^.* at ([\d\.]+) meters above the surface.*$')
re_time_nonst = re.compile(r'^.* time from (\d{8})$')

class CAMSPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing CAMS data...')
        hselect = HorizonSelection.from_cfg(cfg.cams.assim_cycles)
        lvars = QuantityCalculator.get_loaded_vars(cfg.chem_species,
                                                   cfg.cams.output_var_defs)
        ig = InputGatherer(fout, sorted(lvars), rt.nt,
                           ('time', 'z_chem', 'y_chem', 'x_chem'),
                           copy_attrs=['units'])

        first = True
        for fn in glob.glob(rt.paths.cams.file_mask):
            verbose('Parsing CAMS file {}', fn)
            with netCDF4.Dataset(fn, 'r') as fin:
                # Decode times
                vtime = fin.variables['time']
                if vtime.units == 'hours':
                    # CF non-compliant units
                    m = re_time_nonst.match(vtime.long_name)
                    if not m:
                        die('Cannot decode CF non-compliant time units at {}', fn)
                    origin_dt = datetime.datetime.strptime(m.group(1), '%Y%m%d')
                    ncdates = NCDates.from_origin(origin_dt, vtime.units, vtime[:])
                else:
                    ncdates = NCDates(vtime)

                if not ncdates.match_hselect(hselect):
                    verbose('No matching times, skipping file')
                    continue

                # Decode levels
                try:
                    vlev = fin.variables['level']
                except KeyError:
                    levs_in = None
                else:
                    levs_in = vlev[:]

                for vn in ig.varnames:
                    try:
                        var = fin.variables[vn]
                    except KeyError:
                        continue
                    verbose('Loading variable {}', vn)

                    if levs_in is None:
                        m = re_lev.match(var.source)
                        if not m:
                            die('Cannot autodetect level for variable {}', vn)
                        lev = float(m.group(1))

                    if first:
                        # We only setup projection etc. after we know for sure
                        # that we will be loading something from the file

                        verbose('Loading projection and preparing regridder')
                        lats = fin.variables[cfg.cams.varnames.latitude][:]
                        lons = fin.variables[cfg.cams.varnames.longitude][:]
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

                        first = False

                    for it_file, it_sel in ncdates.matching:
                        if levs_in is None:
                            ig.add_single_lev(vn, it_sel, lev, var[it_file], var)
                        else:
                            for il, lev in enumerate(levs_in):
                                ig.add_single_lev(vn, it_sel, lev, var[it_file,il], var)

        ig.finalize()

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing CAMS conversion and horizontal interpolation')
        qc = QuantityCalculator(cfg.chem_species, cfg.cams.output_var_defs,
                                cfg.cams.preprocessors, rt.regrid_cams)
        outdims = ['time', 'z_chem', 'y', 'x']

        with netCDF4.Dataset(rt.paths.intermediate.import_data) as fin:
            levs = fin.variables['z_chem'][:]
            lev_order = np.argsort(levs)

            verbose('Preparing output file')
            # Create dimensions
            for d in ['time', 'z_chem']:
                ensure_dimension(fout, d, len(fin.dimensions[d]))
            ensure_dimension(fout, 'x', rt.nx)
            ensure_dimension(fout, 'y', rt.ny)
            fout.createVariable('z_chem', levs.dtype, ('z_chem',))[:] = levs[lev_order]

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # Load, validate and convert here instead of on import so that
                # we can reorder levels that were loaded sequentially on import
                ts = QuantityCalculator.new_timestep()
                qc.load_timestep_vars(fin, it, ts)
                qc.validate_timestep(ts)

                for sn, v, unit, attrs in qc.calc_timestep_species(ts):
                    if it:
                        v_out = fout.variables[sn]
                    else:
                        v_out = fout.createVariable(sn, v.dtype, outdims)
                        v_out.units = unit
                        if attrs:
                            v_out.setncatts(attrs)

                    v_out[it,:,:,:] = rt.regrid_cams.regrid(v[lev_order,:,:])

    def interpolate_vert(self, fout, *args, **kwargs):
        log('Performing CAMS vertical interpolation')

        with netCDF4.Dataset(rt.paths.intermediate.hinterp) as fin:
            chem_heights = rt.terrain[ax_,:,:] + fin.variables['z_chem'][:][:,ax_,ax_]

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

                # Load all variables for the timestep
                vardata = [fin.variables[vn][it] for vn in cfg.chem_species]

                # Perform vertical interpolation on all currently loaded vars at once
                vinterpolator, = get_vinterp(rt.z_levels_msl, chem_heights, True, False)
                vinterp = vinterpolator(*vardata)
                del vardata, vinterpolator

                for vn, vd in zip(cfg.chem_species, vinterp):
                    v = fout.variables[vn]
                    v[it] = vd

