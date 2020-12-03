import os
import glob
from datetime import datetime, timedelta, timezone
import netCDF4

from core.plugins import Plugin, ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose
from core.config import cfg
from core.runtime import rt
from .palm_wrf_utils import WRFCoordTransform, BilinearRegridder

available_meteo_vars = {
    'tas':  {'desc': 'temperature at surface', 'units': 'K'},
    'ta':   {'desc': '3D temperature', 'units': 'K'},
    'qas':  {'desc': 'specific humidity at the surface', 'units': 'kg/kg'},
    'qa':   {'desc': '3D specific humidity', 'units': '1'},
    'rsds': {'desc': 'surface incident SW radiation for BVOC', 'units': 'W/m2'},
    'pa':   {'desc': '3D pressure', 'units': 'Pa'},
    'zf':   {'desc': 'layer interface heights', 'units': 'm'},
    'ua':   {'desc': 'U-wind', 'units': 'm/s'},
    'va':   {'desc': 'V-wind', 'units': 'm/s'}
}

required_variables = set()


class RequiresMeteoPluginMixin(Plugin):
    """
    Set a list of required meteorological variables in plugin metainformation:

    class Requires:
        meteo_vars = [ ... ]

    Global available_meteo_vars holds names of all variables known to the
    processor.
    """

    def __init__(self, *args, **kwargs):
        if hasattr(self, 'Requires') and hasattr(self.Requires, 'meteo_vars'):
            my_required_vars = set()
            for v in self.Requires.meteo_vars:
                if v not in available_meteo_vars:
                    raise ValueError(
                        'Unknown meteorological variable required by plugin {}.'
                        .format(self))
                else:
                    my_required_vars.add(v)

            required_variables.update(my_required_vars)
        else:
            raise AttributeError(
                'Missing Requires.meteo_vars in plugin {} derived from '
                'RequiresMeteoPluginMixin.'.format(self))

td0 = timedelta(hours=0)
td1h = timedelta(hours=1)
def hrdiff(td):
    d, m = divmod(td, td1h)
    if m != td0:
        raise ValueError('Not a whole hour!')
    return d

class WRFPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, *args, **kwargs):
        log('Importing WRF data...')

        # Set up WRF times
        if rt.nested_domain:
            log('Nested domain - processing only initialization (1 timestep).')
            rt.nt = 1
        else:
            rt.nt = cfg.simulation.length_hours+1

        rt.tindex = lambda dt: hrdiff(dt-rt.start_time)
        rt.end_time = rt.start_time + timedelta(hours=cfg.simulation.length_hours)
        rt.end_time_rad = rt.end_time
        verbose('PALM simulation extent {} - {} ({} timesteps).', rt.start_time,
                rt.end_time, rt.nt)

        with netCDF4.Dataset(rt.paths.imported, 'w', format='NETCDF4') as fout:
            # Process input files
            wrfglob = os.path.join(rt.paths.base,
                    cfg.paths.wrf_output.format(**rt.paths.expand),
                    cfg.paths.wrf_file_mask)
            verbose('Parsing WRF files from {}', wrfglob)
            rt.times = [None] * rt.nt
            first = True
            for fn in glob.glob(wrfglob):
                verbose('Parsing WRF file {}', fn)
                with netCDF4.Dataset(fn) as fin:
                    # Decode time and locate timestep
                    ts = fin.variables['Times'][:].tobytes().decode("utf-8")
                    t = datetime.strptime(ts, '%Y-%m-%d_%H:%M:%S')
                    t = t.replace(tzinfo=timezone.utc)
                    try:
                        it = rt.tindex(t)
                    except ValueError:
                        verbose('Time {} is not within hourly intervals - skipping', t)
                        continue
                    if not (0 <= it < rt.nt):
                        verbose('Time {} is out of range - skipping', t)
                        continue
                    if rt.times[it] is not None:
                        die('Time {} has been already loaded!', t)
                    rt.times[it] = t
                    verbose('Importing time {}, timestep {}', t, it)

                    # coordinate projection
                    rt.trans_wrf = WRFCoordTransform(fin)

                    # dimensions
                    if first:
                        fout.createDimension('Time', rt.nt)
                        for d in '''west_east west_east_stag south_north
                                south_north_stag bottom_top bottom_top_stag
                                soil_layers_stag'''.split():
                            fout.createDimension(d, len(fin.dimensions[d]))

                    # copied vars
                    for varname in cfg.wrf.copy_vars + ['U', 'V', 'P_TOP']:
                        v_wrf = fin.variables[varname]
                        v_out = (fout.createVariable(varname, 'f4', v_wrf.dimensions)
                                if first else fout.variables[varname])
                        v_out[it] = v_wrf[0]

                    # calculated SPECHUM
                    if first:
                        shvars = sorted(set(cfg.wrf.spechum_vars
                            ).intersection(fin.variables.keys()))
                        verbose('Hydro variables in wrf files: {}', ', '.join(shvars))
                    v_out = (fout.createVariable('SPECHUM', 'f4', ('Time',
                        'bottom_top', 'south_north', 'west_east')) if first
                        else fout.variables['SPECHUM'])
                    vdata = fin.variables[shvars[0]][0]
                    for vname in shvars[1:]:
                        vdata += fin.variables[vname][0]
                    v_out[it] = vdata
                    del vdata

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

    def interpolate_horiz(self, *args, **kwargs):
        log('Performing horizontal interpolation')

        verbose('Preparing regridder indices and weights')
        palm_in_wrf_y, palm_in_wrf_x = rt.trans_wrf.latlon_to_ji(
                                        rt.palm_grid_lat, rt.palm_grid_lon)
        regridder = BilinearRegridder(palm_in_wrf_x, palm_in_wrf_y, preloaded=True)
        regridder_u = BilinearRegridder(palm_in_wrf_x+.5, palm_in_wrf_y, preloaded=True)
        regridder_v = BilinearRegridder(palm_in_wrf_x, palm_in_wrf_y+.5, preloaded=True)

        verbose('Preparing output file')
        with netCDF4.Dataset(rt.paths.imported) as fin:
            with netCDF4.Dataset(rt.paths.hinterp, 'w', format='NETCDF4') as fout:
                for d in '''Time  bottom_top bottom_top_stag
                        soil_layers_stag'''.split():
                    fout.createDimension(d, len(fin.dimensions[d]))
                fout.createDimension('west_east', rt.nx)
                fout.createDimension('south_north', rt.ny)

                for varname in cfg.wrf.copy_vars + ['SPECHUM', 'P_TOP']:
                    v_wrf = fin.variables[varname]
                    v_out = fout.createVariable(varname, 'f4', v_wrf.dimensions)
                v_out = fout.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
                v_out = fout.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))

                for it in range(rt.nt):
                    verbose('Processing timestep {}', it)

                    # regular vars
                    for varname in cfg.wrf.copy_vars + ['SPECHUM']:
                        v_wrf = fin.variables[varname]
                        v_out = fout.variables[varname]
                        v_out[it] = regridder.regrid(v_wrf[it,...,regridder.ys,regridder.xs])

                    # U and V have special treatment (unstaggering)
                    fout.variables['U'][it] = regridder_u.regrid(
                            fin.variables['U'][it,...,regridder_u.ys,regridder_u.xs])
                    fout.variables['V'][it] = regridder_v.regrid(
                            fin.variables['V'][it,...,regridder_v.ys,regridder_v.xs])

                    # direct copy
                    fout.variables['P_TOP'][it] = fin.variables['P_TOP'][it]


def palm_wrf_vertical_interp(
        origin_z, terrain, wrf_hybrid_levs, vinterp_terrain_smoothing):

    def interpolate_vert(self, *args, **kwargs):
        log('Performing vertical interpolation')

        verbose('Preparing output file')
        with netCDF4.Dataset(rt.paths.imported) as fin:
            with netCDF4.Dataset(rt.paths.hinterp, 'w', format='NETCDF4') as fout:
                for dimname in ['Time', 'west_east', 'south_north', 'soil_layers_stag']:
                    fout.createDimension(dimname, len(fin.dimensions[dimname]))
                fout.createDimension('z', rt.nz)
                fout.createDimension('zw', rt.nz-1)
                fout.createDimension('zsoil', rt.nz_soil)




                for varname in cfg.wrf.copy_vars + ['SPECHUM']:
                    v_wrf = fin.variables[varname]
                    v_out = fout.createVariable(varname, 'f4', v_wrf.dimensions)
                v_out = fout.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
                v_out = fout.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))

                for it in range(rt.nt):
                    verbose('Processing timestep {}', it)

                    # Use hybrid ETA levels in WRF and stretch them so that the WRF terrain
                    # matches either PALM terrain or flat terrain at requested height
                    gpf = fin.variables['PH'][it,:,:,:] + fin.variables['PHB'][it,:,:,:]
                    wrfterr = gpf[0]*(1./g)

                    if vinterp_terrain_smoothing is None:
                        target_terrain = terrain
                    else:
                        print('Smoothing PALM terrain for the purpose of dynamic driver with sigma={0} grid points.'.format(
                            vinterp_terrain_smoothing))
                        target_terrain = ndimage.gaussian_filter(terrain, sigma=vinterp_terrain_smoothing, order=0)
                    print('Morphing WRF terrain ({0} ~ {1}) to PALM terrain ({2} ~ {3})'.format(
                        wrfterr.min(), wrfterr.max(), target_terrain.min(), target_terrain.max()))
                    print_dstat('terrain shift', wrfterr - target_terrain[:,:])

                    # Load original dry air column pressure
                    mu = fin.variables['MUB'][it,:,:] + fin.variables['MU'][it,:,:]
                    pht = fin.variables['P_TOP'][it]

                    # Shift column pressure so that it matches PALM terrain
                    t = wrf_t(fin)
                    mu2 = barom_pres(mu+pht, target_terrain*g, gpf[0,:,:], t[0,:,:])-pht

                    # Calculate original and shifted 3D dry air pressure
                    if wrf_hybrid_levs:
                        phf, phh = calc_ph_hybrid(nc_wrf, mu)
                        phf2, phh2 = calc_ph_hybrid(nc_wrf, mu2)
                    else:
                        phf, phh = calc_ph_sigma(nc_wrf, mu)
                        phf2, phh2 = calc_ph_sigma(nc_wrf, mu2)

                    # Shift 3D geopotential according to delta dry air pressure
                    tf = np.concatenate((t, t[-1:,:,:]), axis=0) # repeat highest layer
                    gpf2 = barom_gp(gpf, phf2, phf, tf)
                    # For half-levs, originate from gp full levs rather than less accurate gp halving
                    gph2 = barom_gp(gpf[:-1,:,:], phh2, phf[:-1,:,:], t)
                    zf = gpf2 * (1./g) - origin_z
                    zh = gph2 * (1./g) - origin_z

                    # Report
                    gpdelta = gpf2 - gpf
                    print('GP deltas by level:')
                    for k in range(gpf.shape[0]):
                        print_dstat(k, gpdelta[k])

                    # Because we require levels below the lowest level from WRF, we will always
                    # add one layer at zero level with repeated values from the lowest level.
                    # WRF-python had some special treatment for theta in this case.
                    height = np.zeros((zh.shape[0]+1,) + zh.shape[1:], dtype=zh.dtype)
                    height[0,:,:] = -999. #always below terrain
                    height[1:,:,:] = zh
                    heightw = np.zeros((zf.shape[0]+1,) + zf.shape[1:], dtype=zf.dtype)
                    heightw[0,:,:] = -999. #always below terrain
                    heightw[1:,:,:] = zf

                    # ======================== SPECIFIC HUMIDITY ==============================
                    qv_raw = nc_infile.variables['SPECHUM'][it]
                    qv_raw = np.r_[qv_raw[0:1], qv_raw]
                    
                    # Vertical interpolation to grid height levels (specified in km!)
                    # Levels start at 50m (below that the interpolation looks very sketchy)
                    init_atmosphere_qv = interpolate_1d(z_levels, height, qv_raw)
                    vdata = nc_outfile.createVariable('init_atmosphere_qv', "f4", ("Time", "z","south_north","west_east"))
                    vdata[it,:,:,:] = init_atmosphere_qv
                    
                    # ======================= POTENTIAL TEMPERATURE ==========================
                    pt_raw = nc_infile.variables['T'][it] + 300.   # from perturbation pt to standard
                    pt_raw = np.r_[pt_raw[0:1], pt_raw]
                    
                    #plt.figure(); plt.contourf(pt[0]) ; plt.colorbar() ; plt.show()
                    vdata = nc_outfile.createVariable('init_atmosphere_pt', "f4", ("Time", "z","south_north","west_east"))
                    vdata[it,:,:,:] = init_atmosphere_pt
                    
                    # ======================= Wind ==========================================
                    u_raw = nc_infile.variables['U'][it]
                    u_raw = np.r_[u_raw[0:1], u_raw]
                    init_atmosphere_u = interpolate_1d(z_levels, height, u_raw)
                    
                    vdata = nc_outfile.createVariable('init_atmosphere_u', "f4", ("Time", "z","south_north","west_east"))
                    vdata[it,:,:,:] = init_atmosphere_u

                    v_raw = nc_infile.variables['V'][it]
                    v_raw = np.r_[v_raw[0:1], v_raw]
                    init_atmosphere_v = interpolate_1d(z_levels, height, v_raw)
                    
                    vdata = nc_outfile.createVariable('init_atmosphere_v', "f4", ("Time", "z","south_north","west_east"))
                    #vdata.coordinates = "XLONG_V XLAT_V XTIME"
                    vdata[it,:,:,:] = init_atmosphere_v
                    
                    w_raw = nc_infile.variables['W'][it]
                    w_raw = np.r_[w_raw[0:1], w_raw]
                    init_atmosphere_w = interpolate_1d(z_levels_stag, heightw, w_raw)
                    
                    vdata = nc_outfile.createVariable('init_atmosphere_w', "f4", ("Time", "zw","south_north","west_east"))
                    #vdata.coordinates = "XLONG XLAT XTIME"
                    vdata[it,:,:,:] = init_atmosphere_w

                    # ===================== SURFACE PRESSURE ==================================
                    surface_forcing_surface_pressure = nc_infile.variables['PSFC']
                    vdata = nc_outfile.createVariable('surface_forcing_surface_pressure', "f4", ("Time", "south_north","west_east"))
                    vdata[it,:,:] = surface_forcing_surface_pressure[it,:,:]
                    
                    # ======================== SOIL VARIABLES (without vertical interpolation) =============
                    # soil temperature
                    init_soil_t = nc_infile.variables['TSLB']
                    vdata = nc_outfile.createVariable('init_soil_t', "f4", ("Time", "zsoil","south_north","west_east"))
                    vdata[it,:,:,:] = init_soil_t[it,:,:,:]

                    # soil moisture
                    init_soil_m = nc_infile.variables['SMOIS']
                    vdata = nc_outfile.createVariable('init_soil_m', "f4", ("Time","zsoil","south_north","west_east"))
                    vdata[it,:,:,:] = init_soil_m[it,:,:,:]

                    # zsoil
                    zsoil = nc_wrf.variables['ZS']    #ZS:description = "DEPTHS OF CENTERS OF SOIL LAYERS" ;
                    vdata = nc_outfile.createVariable('zsoil', "f4", ("zsoil"))
                    vdata[:] = zsoil[it,:]

                    # coordinates z, zw
                    vdata = nc_outfile.createVariable('z', "f4", ("z"))
                    vdata[:] = list(z_levels)

                    vdata = nc_outfile.createVariable('zw', "f4", ("zw"))
                    vdata[:] = list (z_levels_stag) 

                    # zsoil is taken from wrf - not need to define it

                    nc_infile.close()
                    nc_wrf.close()
                    nc_outfile.close()


    class Provides:
        meteo_vars = ['tas', 'pa']


class EmisPlugin(RequiresMeteoPluginMixin, ImportPluginMixin):
    def import_data(self, *args, **kwargs):
        print('Import emission data')

    class Requires:
        meteo_vars = ['tas', 'pa']
