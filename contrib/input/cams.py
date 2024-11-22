import os
import re
import glob
from datetime import datetime, timedelta, timezone
import numpy as np
import netCDF4
from metpy.interpolate import interpolate_1d

from core.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose, log_output
from core.config import cfg
from core.runtime import rt
from core.utils import ensure_dimension
from .wrf_utils import BilinearRegridder, radius
from .camx import CAMxConvertor
import pyproj

ax_ = np.newaxis
re_num = re.compile(r'[0-9\.]+')

rdivcp = 0.286
pressure_ref = 1.0e5
R = 8.314
ppm_conv = 1e6


class CAMSPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing CAMS data...')

        filled = [False] * rt.nt
        timesteps = [CAMxConvertor.new_timestep() for t in rt.times]
        zcoord = [None] * rt.nt

        # Process input files
        verbose('Parsing CAMS file {}', cfg.paths.cams_output)
        with netCDF4.Dataset(cfg.paths.cams_output, 'r') as fin:
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
            rt.trans_cams = CAMSCoordTransform(fin, rt.palm_epsg)
            palm_in_cams_y, palm_in_cams_x = rt.trans_cams.latlon_to_ji(rt.palm_grid_lat, rt.palm_grid_lon)
            rt.regrid_cams = BilinearRegridder(palm_in_cams_x, palm_in_cams_y, preloaded=True)
            del palm_in_cams_y, palm_in_cams_x

            convertor = CAMxConvertor(cfg.chem_species,
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
                die('Could not find all CAMx variables for all times.\n'
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
                for sn, v, unit in convertor.calc_timestep_species(tsdata):
                    v_out = (fout.variables[sn] if i
                             else fout.createVariable(sn, 'f4', chem_dims))
                    v_out.units = unit
                    v_out[i, :, :, :] = v

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing CAMS horizontal interpolation')
        hvars = ['height_chem'] + cfg.chem_species

        with netCDF4.Dataset(rt.paths.imported) as fin:
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
                v_out.units = v_in.units
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

        with netCDF4.Dataset(rt.paths.hinterp) as fin:
            agl_chem = fin.variables['height_chem']
            chem_heights = np.zeros((agl_chem.shape[1]+1,) + agl_chem.shape[2:], dtype=agl_chem.dtype)
            chem_heights[0,:,:] = -999.

            verbose('Preparing output file')
            for dimname in ['time', 'y', 'x']:
                ensure_dimension(fout, dimname, len(fin.dimensions[dimname]))
            ensure_dimension(fout, 'z', rt.nz)

            for vn in cfg.chem_species:
                var = fout.createVariable(vn, 'f4', ('time', 'z', 'y', 'x'))
                var.units = fin.variables[vn].units

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

class CAMSCoordTransform(object):
    'Coordinate transformer for CAMx files running from WRF'

    def __init__(self, ncf, palm_epsg):
        attr = lambda a: getattr(ncf, a)

        # Define grids

        latlon_sphere = pyproj.Proj(proj='latlong',
            a=radius, b=radius,
            towgs84='0,0,0', no_defs=True)

        lambert_grid = pyproj.Proj(init='EPSG:{}'.format(palm_epsg))

        # number of mass grid points
        self.nx = nx = ncf.variables['longitude'][:].size
        self.ny = ny = ncf.variables['latitude'][:].size

        # Define fast transformation methods
        lat = ncf.variables['latitude'][:]
        lon = ncf.variables['longitude'][:]

        # create lat, lon grid
        LON, LAT = np.meshgrid(lon, lat)
        x, y = pyproj.transform(latlon_sphere, lambert_grid,
                                LON, LAT)
        dx, dy = x[0,1]-x[0,0], y[1,0]-y[0,0]
        i0_x, j0_y = x[0,0], y[0,0]

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
        lat = ncf.variables['latitude'][:]
        lon = ncf.variables['longitude'][:]
        j, i = np.mgrid[0:self.ny, 0:self.nx]

        jj, ii = self.latlon_to_ji(lat, lon)
        d = np.hypot(jj-j, ii-i)
        print('error for ll->ji: max {0} m, avg {1} m.'.format(d.max(), d.mean()))

        llat, llon = self.ji_to_latlon(j, i)
        d = np.hypot(llat - lat, llon - lon)
        print('error for ji->ll: max {0} deg, avg {1} deg.'.format(d.max(), d.mean()))
