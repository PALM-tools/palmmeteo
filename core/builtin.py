import os
import re
import numpy as np
import netCDF4
from pyproj import Proj, transform

from .plugins import SetupPluginMixin, WritePluginMixin
from .logging import die, warn, log, verbose
from .config import cfg
from .runtime import rt

fext_re = re.compile(r'\.(\d{3})$')

def find_free_fname(fpath):
    if not os.path.exists(fpath):
        return fpath

    if cfg.output.overwrite:
        log('Existing file {} will be overwritten.', fpath)
        return fpath

    # Try to find free fpath.###
    path, base = os.path.split(fpath)
    nbase = len(base)
    maxnum = -1
    for fn in os.listdir(path):
        if not fn.startswith(base):
            continue
        m = fext_re.match(fn[nbase:])
        if not m:
            continue
        maxnum = max(maxnum, int(m.group(1)))
    if maxnum >= 999:
        raise RuntimeError('Cannot find free filename starting with ' + fpath)

    newpath = '{}.{:03d}'.format(fpath, maxnum+1)
    log('Filename {} exists, using {}.', fpath, newpath)
    return newpath


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
        inproj = Proj('+init='+cfg.domain.proj_palm)
        lonlatproj = Proj('+init='+cfg.domain.proj_wgs84)
        rt.cent_lon, rt.cent_lat = transform(inproj, lonlatproj, rt.xcent, rt.ycent)
        verbose('xcent={}, ycent={}', rt.xcent, rt.ycent)
        verbose('cent_lon={}, cent_lat={}', rt.cent_lon, rt.cent_lat)
        # prepare target grid
        irange = rt.origin_x + rt.dx * (np.arange(rt.nx, dtype='f8') + .5)
        jrange = rt.origin_y + rt.dy * (np.arange(rt.ny, dtype='f8') + .5)
        rt.palm_grid_y, rt.palm_grid_x = np.meshgrid(jrange, irange, indexing='ij')
        rt.palm_grid_lon, rt.palm_grid_lat = transform(inproj, lonlatproj,
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
                die('domain:dz_stretch_level has to be set for stretching.', cfg.domain.dz_max, rt.dz)
            if cfg.domain.dz_max < rt.dz:
                die('domain:dz_max (={}) has to be >= than dz (={}).', cfg.domain.dz_max, rt.dz)
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


class WritePlugin(WritePluginMixin):
    def write_data(self, *args, **kwargs):
        log('Writing data to dynamic driver')

        fn_out = find_free_fname(rt.paths.dynamic_driver)
        log('Preparing dynamic driver file {}.', fn_out)
        with netCDF4.Dataset(fn_out, 'w', format='NETCDF4') as fout:
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
            fout.createVariable('time',  'f4', ('time',) )[:] = rt.times_sec
            fout.createVariable('z',     'f4', ('z',)    )[:] = rt.z_levels[:]
            fout.createVariable('zw',    'f4', ('zw',)   )[:] = rt.z_levels_stag[:]
            fout.createVariable('zsoil', 'f4', ('zsoil',))[:] = rt.z_soil_levels[:]
            fout.createVariable('y',     'f4', ('y',)    )[:] = rt.dy/2 + rt.dy*np.arange(rt.ny)
            fout.createVariable('x',     'f4', ('x',)    )[:] = rt.dx/2 + rt.dx*np.arange(rt.nx)

            # Create init variables
            fout.createVariable('init_atmosphere_pt', 'f4', ('z', 'y', 'x')).lod = 2
            fout.createVariable('init_atmosphere_qv', 'f4', ('z', 'y', 'x')).lod = 2
            fout.createVariable('init_atmosphere_u', 'f4', ('z', 'y', 'xu')).lod = 2
            fout.createVariable('init_atmosphere_v', 'f4', ('z', 'yv', 'x')).lod = 2
            fout.createVariable('init_atmosphere_w', 'f4', ('zw', 'y', 'x')).lod = 2
            fout.createVariable('init_soil_t', 'f4', ('zsoil', 'y', 'x')).lod = 2
            fout.createVariable('init_soil_m', 'f4', ('zsoil', 'y', 'x')).lod = 2

            # Create forcing variables
            if not rt.nested_domain:
                # surface pressure (scalar)
                fout.createVariable('surface_forcing_surface_pressure', 'f4', ('time',))

                # boundary - vertical slices from left, right, south, north, top
                fout.createVariable('ls_forcing_left_pt',  'f4', ('time', 'z',  'y' ))
                fout.createVariable('ls_forcing_right_pt', 'f4', ('time', 'z',  'y' ))
                fout.createVariable('ls_forcing_south_pt', 'f4', ('time', 'z',  'x' ))
                fout.createVariable('ls_forcing_north_pt', 'f4', ('time', 'z',  'x' ))
                fout.createVariable('ls_forcing_top_pt',   'f4', ('time', 'y',  'x' ))

                fout.createVariable('ls_forcing_left_qv',  'f4', ('time', 'z',  'y' ))
                fout.createVariable('ls_forcing_right_qv', 'f4', ('time', 'z',  'y' ))
                fout.createVariable('ls_forcing_south_qv', 'f4', ('time', 'z',  'x' ))
                fout.createVariable('ls_forcing_north_qv', 'f4', ('time', 'z',  'x' ))
                fout.createVariable('ls_forcing_top_qv',   'f4', ('time', 'y',  'x' ))

                fout.createVariable('ls_forcing_left_u',   'f4', ('time', 'z',  'y' ))
                fout.createVariable('ls_forcing_right_u',  'f4', ('time', 'z',  'y' ))
                fout.createVariable('ls_forcing_south_u',  'f4', ('time', 'z',  'xu'))
                fout.createVariable('ls_forcing_north_u',  'f4', ('time', 'z',  'xu'))
                fout.createVariable('ls_forcing_top_u',    'f4', ('time', 'y',  'xu'))

                fout.createVariable('ls_forcing_left_v',   'f4', ('time', 'z',  'yv'))
                fout.createVariable('ls_forcing_right_v',  'f4', ('time', 'z',  'yv'))
                fout.createVariable('ls_forcing_south_v',  'f4', ('time', 'z',  'x' ))
                fout.createVariable('ls_forcing_north_v',  'f4', ('time', 'z',  'x' ))
                fout.createVariable('ls_forcing_top_v',    'f4', ('time', 'yv', 'x' ))

                fout.createVariable('ls_forcing_left_w',   'f4', ('time', 'zw', 'y' ))
                fout.createVariable('ls_forcing_right_w',  'f4', ('time', 'zw', 'y' ))
                fout.createVariable('ls_forcing_south_w',  'f4', ('time', 'zw', 'x' ))
                fout.createVariable('ls_forcing_north_w',  'f4', ('time', 'zw', 'x' ))
                fout.createVariable('ls_forcing_top_w',    'f4', ('time', 'y',  'x' ))

                # geostrophic wind (1D)
                fout.createVariable('ls_forcing_ug', 'f4', ('time', 'z'))
                fout.createVariable('ls_forcing_vg', 'f4', ('time', 'z'))

            # prepare influx/outflux area sizes
            zstag_all = np.r_[0., rt.z_levels_stag, rt.ztop]
            zwidths = zstag_all[1:] - zstag_all[:-1]
            verbose('Z widths: {}', zwidths)
            areas_xb = (zwidths * rt.dy)[:,np.newaxis]
            areas_yb = (zwidths * rt.dx)[:,np.newaxis]
            areas_zb = rt.dx * rt.dy
            area_boundaries = (areas_xb.sum()*rt.ny*2
                    + areas_yb.sum()*rt.nx*2
                    + areas_zb*rt.nx*rt.ny)

            log('Writing values for initialization variables')
            with netCDF4.Dataset(rt.paths.vinterp) as fin:
                # write values for initialization variables
                fout.variables['init_atmosphere_pt'][:,:,:] = fin.variables['init_atmosphere_pt'][0, :, :, :]
                fout.variables['init_atmosphere_qv'][:,:,:] = fin.variables['init_atmosphere_qv'][0, :, :, :]
                fout.variables['init_atmosphere_u'][:,:,:] = fin.variables['init_atmosphere_u'][0, :, :, 1:]
                fout.variables['init_atmosphere_v'][:,:,:] = fin.variables['init_atmosphere_v'][0, :, 1:, :]
                fout.variables['init_atmosphere_w'][:,:,:] = fin.variables['init_atmosphere_w'][0, :, :, :]
                fout.variables['init_soil_t'][:,:,:] = fin.variables['init_soil_t'][0,:,:,:]
                fout.variables['init_soil_m'][:,:,:] = (fin.variables['init_soil_m'][0,:,:,:]
                        * rt.soil_moisture_adjust[np.newaxis,:,:])

                # Write values for time dependent values
                if not rt.nested_domain:
                    for it in range(rt.nt):
                        verbose('Processing timestep {}', it)

                        # surface pressure
                        fout.variables['surface_forcing_surface_pressure'][it] = (
                                fin.variables['surface_forcing_surface_pressure'][it,:,:].mean())

                        # boundary conditions
                        fout.variables['ls_forcing_left_pt' ][it,:,:] = fin.variables['init_atmosphere_pt'][0, :, :, 0]
                        fout.variables['ls_forcing_right_pt'][it,:,:] = fin.variables['init_atmosphere_pt'][0, :, :, rt.nx-1]
                        fout.variables['ls_forcing_south_pt'][it,:,:] = fin.variables['init_atmosphere_pt'][0, :, 0, :]
                        fout.variables['ls_forcing_north_pt'][it,:,:] = fin.variables['init_atmosphere_pt'][0, :, rt.ny-1, :]
                        fout.variables['ls_forcing_top_pt'  ][it,:,:] = fin.variables['init_atmosphere_pt'][0, rt.nz-1, :, :]

                        fout.variables['ls_forcing_left_qv' ][it,:,:] = fin.variables['init_atmosphere_qv'][0, :, :, 0]
                        fout.variables['ls_forcing_right_qv'][it,:,:] = fin.variables['init_atmosphere_qv'][0, :, :, rt.nx-1]
                        fout.variables['ls_forcing_south_qv'][it,:,:] = fin.variables['init_atmosphere_qv'][0, :, 0, :]
                        fout.variables['ls_forcing_north_qv'][it,:,:] = fin.variables['init_atmosphere_qv'][0, :, rt.ny-1, :]
                        fout.variables['ls_forcing_top_qv'  ][it,:,:] = fin.variables['init_atmosphere_qv'][0, rt.nz-1, :, :]

                        # Perform mass balancing for U, V, W
                        uxleft = fin.variables['init_atmosphere_u'][0, :, :, 0]
                        uxright = fin.variables['init_atmosphere_u'][0, :, :, rt.nx-1]
                        vysouth = fin.variables['init_atmosphere_v'][0, :, 0, :]
                        vynorth = fin.variables['init_atmosphere_v'][0, :, rt.ny-1, :]
                        wztop = fin.variables['init_atmosphere_w'][0, rt.nz-2, :, :]#nzw=nz-1
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
                        fout.variables['ls_forcing_left_u' ][it,:,:] = uxleft
                        fout.variables['ls_forcing_right_u'][it,:,:] = uxright
                        fout.variables['ls_forcing_south_u'][it,:,:] = fin.variables['init_atmosphere_u'][0, :, 0, 1:]
                        fout.variables['ls_forcing_north_u'][it,:,:] = fin.variables['init_atmosphere_u'][0, :, rt.ny-1, 1:]
                        fout.variables['ls_forcing_top_u'  ][it,:,:] = fin.variables['init_atmosphere_u'][0, rt.nz-1, :, 1:]

                        fout.variables['ls_forcing_left_v' ][it,:,:] = fin.variables['init_atmosphere_v'][0, :, 1:, 0]
                        fout.variables['ls_forcing_right_v'][it,:,:] = fin.variables['init_atmosphere_v'][0, :, 1:, rt.nx-1]
                        fout.variables['ls_forcing_south_v'][it,:,:] = vysouth
                        fout.variables['ls_forcing_north_v'][it,:,:] = vynorth
                        fout.variables['ls_forcing_top_v'  ][it,:,:] = fin.variables['init_atmosphere_v'][0, rt.nz-1, 1:, :]

                        fout.variables['ls_forcing_left_w' ][it,:,:] = fin.variables['init_atmosphere_w'][0, :, :, 0]
                        fout.variables['ls_forcing_right_w'][it,:,:] = fin.variables['init_atmosphere_w'][0, :, :, rt.nx-1]
                        fout.variables['ls_forcing_south_w'][it,:,:] = fin.variables['init_atmosphere_w'][0, :, 0, :]
                        fout.variables['ls_forcing_north_w'][it,:,:] = fin.variables['init_atmosphere_w'][0, :, rt.ny-1, :]
                        fout.variables['ls_forcing_top_w'  ][it,:,:] = wztop

                        # copy geostrophic wind
                        fout.variables['ls_forcing_ug'][it, :] = fin.variables['ls_forcing_ug'][it, :]
                        fout.variables['ls_forcing_vg'][it, :] = fin.variables['ls_forcing_vg'][it, :]

                ### # Write chemical boundary conds
                ### if camx_interp_fname:
                ###     f_camx = netCDF4.Dataset(camx_interp_fname)
                ###     try:
                ###         for vname, vval in f_camx.variables.items():
                ###             # PALM doesn't support 3D LOD=2 init for chem yet, we have to average the field
                ###             var = fout.createVariable('init_atmosphere_'+vname,
                ###                     'f4', ('z',), )
                ###             var.units = vval.units
                ###             var.lod = 1
                ###             var[:] = vval[0,:,:,:].mean(axis=(1,2))

                ###             var = fout.createVariable('ls_forcing_left_'+vname,
                ###                     'f4', ('time','z','y'), )
                ###             var.units = vval.units
                ###             var[:] = vval[:,:,:,0]

                ###             var = fout.createVariable('ls_forcing_right_'+vname,
                ###                     'f4', ('time','z','y'), )
                ###             var.units = vval.units
                ###             var[:] = vval[:,:,:,-1]

                ###             var = fout.createVariable('ls_forcing_south_'+vname,
                ###                     'f4', ('time','z','x'), )
                ###             var.units = vval.units
                ###             var[:] = vval[:,:,0,:]

                ###             var = fout.createVariable('ls_forcing_north_'+vname,
                ###                     'f4', ('time','z','x'), )
                ###             var.units = vval.units
                ###             var[:] = vval[:,:,-1,:]

                ###             var = fout.createVariable('ls_forcing_top_'+vname,
                ###                     'f4', ('time','y','x'), )
                ###             var.units = vval.units
                ###             var[:] = vval[:,-1,:,:]
                ###     finally:
                ###         f_camx.close()

            ### if len(rad_times_proc) > 0:
            ###     # process radiation inputs
            ###     # radiation time dimension and variable
            ###     fout.createDimension('time_rad', len(rad_times_proc))
            ###     _val_times =  fout.createVariable('time_rad',"f4", ("time_rad"))
            ###     _val_times[:] = rad_times_proc[:]
            ###     # radiation variables
            ###     var = fout.createVariable('rad_sw_in', "f4", ("time_rad"), )
            ###     var.setncattr('lod', 1)
            ###     var.units = 'W/m2'
            ###     var[:] = rad_values_proc[0][:]
            ###     var = fout.createVariable('rad_lw_in', "f4", ("time_rad"), )
            ###     var.setncattr('lod', 1)
            ###     var.units = 'W/m2'
            ###     var[:] = rad_values_proc[1][:]
            ###     var = fout.createVariable('rad_sw_in_dif', "f4", ("time_rad"), )
            ###     var.setncattr('lod', 1)
            ###     var.units = 'W/m2'
            ###     var[:] = rad_values_proc[2][:]

        log('Dynamic driver written successfully.')
