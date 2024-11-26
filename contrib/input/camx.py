import os
import re
import glob
from datetime import datetime
import numpy as np
import netCDF4
from metpy.interpolate import interpolate_1d

from core.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose, log_output
from core.config import cfg
from core.runtime import rt
from core.utils import ensure_dimension
from core.library import QuantityCalculator
from .wrf_utils import CAMxCoordTransform, BilinearRegridder

ax_ = np.newaxis
re_num = re.compile(r'[0-9\.]+')


class CAMxPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing CAMx data...')

        filled = [False] * rt.nt
        timesteps = [QuantityCalculator.new_timestep() for t in rt.times]
        zcoord = [None] * rt.nt

        # Process input files
        camxglob = os.path.join(rt.paths.base,
                cfg.paths.camx_output.format(**rt.paths.expand),
                cfg.paths.camx_file_mask)
        verbose('Parsing CAMx files from {}', camxglob)
        first = True
        for fn in glob.glob(camxglob):
            verbose('Parsing CAMx file {}', fn)
            with netCDF4.Dataset(fn) as fin:
                # Decode time and locate timesteps
                tflag = fin.variables['TFLAG'][:]
                assert len(tflag.shape) == 3 and tflag.shape[2] == 2
                xdate = tflag[:,0,0]
                xtime = tflag[:,0,1]

                # Verify that dates are equal for each variable
                assert (tflag[:,:,0] == xdate[:,ax_]).all()
                assert (tflag[:,:,1] == xtime[:,ax_]).all()

                dts = []
                for i in range(len(xdate)):
                    dt = datetime.strptime(
                            '{0:07d} {1:06d} +0000'.format(xdate[i], xtime[i]),
                            '%Y%j %H%M%S %z')
                    ireq = rt.tindex(dt)
                    if not 0 <= ireq < rt.nt:
                        continue
                    dts.append((ireq, i))

                if not dts:
                    verbose('Skipping CAMx file {0}: no requested times.'.format(fn))
                    continue

                verbose('Processing CAMx file {0}.'.format(fn))

                # locate layer heights
                try:
                    vz = fin.variables['z']
                except KeyError:
                    verbose('Loading heights from separate file')
                    vz = None
                    with open(fn+'.heights') as fh:
                        fix_hgt = np.array(list(map(float, re_num.findall(fh.read())))) * 1000. #orig in km
                        fix_hgt = fix_hgt[:,ax_,ax_] #convert to 3D for broadcasting
                else:
                    verbose('Loading heights from variable z')

                if first:
                    # coordinate projection
                    verbose('Loading projection and preparing regridder')
                    rt.trans_camx = CAMxCoordTransform(fin)
                    palm_in_camx_y, palm_in_camx_x = rt.trans_camx.latlon_to_ji(
                                                        rt.palm_grid_lat, rt.palm_grid_lon)
                    rt.regrid_camx = BilinearRegridder(palm_in_camx_x, palm_in_camx_y, preloaded=True)
                    del palm_in_camx_y, palm_in_camx_x

                    convertor = QuantityCalculator(cfg.chem_species,
                            cfg.camx.output_var_defs, cfg.camx.preprocessors,
                            rt.regrid_camx)

                    # dimensions
                    ensure_dimension(fout, 'time', rt.nt)
                    ensure_dimension(fout, 'z_chem', fix_hgt.shape[0]
                            if vz is None else vz.shape[1])
                    ensure_dimension(fout, 'y_chem', rt.regrid_camx.ylen)
                    ensure_dimension(fout, 'x_chem', rt.regrid_camx.xlen)
                    chem_dims = ('time', 'z_chem', 'y_chem', 'x_chem')

                for itout, itf in dts:
                    verbose('Importing timestep {} -> {}', itf, itout)

                    # TODO: check for dicrepancy among input files
                    zcoord[itout] = fix_hgt if vz is None else rt.regrid_camx.loader(vz)[itf, :]

                    filled[itout] = convertor.load_timestep_vars(fin, itf,
                            timesteps[itout])
            first = False

        if first:
            die('No CAMx files found under {}.', camxglob)

        if not all(filled):
            die('Could not find all CAMx variables for all times.\n'
                    'Missing variables in times:\n{}',
                    '\n'.join('{}: {}'.format(dt, ', '.join(vn
                            for vn in(convertor.loaded_vars-tsdata)))
                        for dt, fil, tsdata in zip(rt.times, filled, timesteps)
                        if not fil))

        vz_out = fout.createVariable('height_chem', 'f4', chem_dims)
        vz_out.units = 'm'

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
        log('Performing CAMx horizontal interpolation')
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
                v_out.setncatts({a: v_in.getncattr(a) for a in v_in.ncattrs()})
            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # regular vars
                for varname in hvars:
                    v_in = fin.variables[varname]
                    v_out = fout.variables[varname]
                    v_out[it] = rt.regrid_camx.regrid(v_in[it])

    def interpolate_vert(self, fout, *args, **kwargs):
        log('Performing CAMx vertical interpolation')
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
                v_in = fin.variables[vn]
                var = fout.createVariable(vn, 'f4', ('time', 'z', 'y', 'x'))
                var.setncatts({a: v_in.getncattr(a) for a in v_in.ncattrs()})

            for it in range(rt.nt):
                verbose('Processing timestep {}', it)

                # Calc CAMx layer heights
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

