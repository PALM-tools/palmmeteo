from core.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose, log_output
from core.config import cfg, ConfigError
from core.runtime import rt

class CAMxPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing CAMx data...')

        # Process input files
        camxglob = os.path.join(rt.paths.base,
                cfg.paths.camx_output.format(**rt.paths.expand),
                cfg.paths.camx_file_mask)
        verbose('Parsing CAMx files from {}', camxglob)
        rt.times = [None] * rt.nt
        first = True
        for fn in glob.glob(camxglob):
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
                    palm_in_wrf_y, palm_in_wrf_x = rt.trans_wrf.latlon_to_ji(
                                                    rt.palm_grid_lat, rt.palm_grid_lon)
                    rt.regrid_wrf = BilinearRegridder(palm_in_wrf_x, palm_in_wrf_y, preloaded=True)
                    rt.regrid_wrf_u = BilinearRegridder(palm_in_wrf_x+.5, palm_in_wrf_y, preloaded=True)
                    rt.regrid_wrf_v = BilinearRegridder(palm_in_wrf_x, palm_in_wrf_y+.5, preloaded=True)
                    del palm_in_wrf_y, palm_in_wrf_x

def tflag(data, req_dts):
    assert len(data.shape) == 3 and data.shape[2] == 2
    xdate = data[:,0,0]
    xtime = data[:,0,1]

    # Verify that dates are equal for each variable
    assert (data[:,:,0] == xdate[:,_na]).all()
    assert (data[:,:,1] == xtime[:,_na]).all()

    dts = []
    for i in range(len(xdate)):
        dt = datetime.datetime.strptime(
            '{0:07d} {1:06d}'.format(xdate[i], xtime[i]), '%Y%j %H%M%S')
        try:
            ireq = req_dts[dt]
        except KeyError:
            continue
        dts.append((ireq, i))
    return dts

def load_conversion(spc, f, sl, hlp, camx_vars, camx_units, formula, **kwargs):
    loaded_vars = []
    for varname, unit in zip(camx_vars, camx_units):
        try:
            v = f.variables[varname]
        except KeyError:
            print('Skipping {0}, because input variable {1}({2}) is missing.'.format(
                spc, varname, unit))
            return None
        if getattr(v, 'units', None) != unit:
            print('Skipping {0}, because input variable {1} has wrong unit ({2} <> {3}).'.format(
                spc, varname, getattr(v, 'units', None), unit))
            return None

        print('Loading variable {0}.'.format(varname))
        loaded_vars.append(v[sl])
    try:
        val = formula(*(loaded_vars + [hlp]))
    except HelperNotFound as e:
        print('Skipping {0} - missing helper variable {1}.'.format(spc, e.name))
        return None

    return val

def process_tstep(f, itf, regridder, lay_height, fout, itout, z_levels,
        vars_remaining, filled, conversions, helpers):

    # Load helper vars for this timestep
    hlp = Helpers()
    for helper_name, helper in helpers:
        data = load_conversion(helper_name, f,
                (itf,slice(None),regridder.ys,regridder.xs), hlp, **helper)
        if data is not None:
            setattr(hlp, helper_name, data)

    # Load all usable vars for this timestep, regrid horizontally
    varmeta = []
    vardata = []
    for spc in list(vars_remaining):
        conv = conversions[spc]
        data = load_conversion(spc, f, (itf,slice(None),regridder.ys,regridder.xs),
                hlp, **conv)
        if data is None:
            continue

        data = regridder.regrid(data)
        vardata.append(np.r_[data[0:1], data]) #add peg below
        varmeta.append((spc, conv['output_unit']))

    # Perform vertical interpolation on all currently loaded vars at once
    print('Interpolating vertically...')
    vinterp = interpolate_1d(z_levels, lay_height, *vardata)
    if len(vardata) == 1:
        # return_list_always=True argument is only in later versions of MetPy
        vinterp = [vinterp]
    del vardata
    for (vn, vu), vd in zip(varmeta, vinterp):
        v = fout.variables[vn]
        v[itout] = vd
        v.units = vu
        filled[vn][itout] = True

def process_files(camx_file_list, camx_interp_fname, palm_grid_lat,
        palm_grid_lon, terrain_rel, z_levels, times, species_names,
        conversions, helpers):

    terrain_shift = terrain_rel[_na,:,:]
    lowest_layer = np.zeros(((1,) + palm_grid_lat.shape), dtype='f4')
    lowest_layer[:] = -999.
    tindex = dict((dt, i) for i, dt in enumerate(times))
    filled = {}

    with netCDF4.Dataset(camx_interp_fname, 'w', format='NETCDF4') as fout:
        fout.createDimension('time', len(times))
        fout.createDimension('z', len(z_levels))
        fout.createDimension('y', palm_grid_lat.shape[0])
        fout.createDimension('x', palm_grid_lat.shape[1])
        for vn in species_names:
            fout.createVariable(vn, 'f4', ('time', 'z', 'y', 'x'))
            filled[vn] = [False] * len(times)

        for fname in sorted(camx_file_list):
            with netCDF4.Dataset(fname) as f:
                dts = tflag(f.variables['TFLAG'][:], tindex)
                if dts:
                    print('Processing CAMx file {0}.'.format(fname))

                    # preprare projection
                    trans = palm_wrf_utils.CAMxCoordTransform(f)
                    palm_in_camx_y, palm_in_camx_x = trans.latlon_to_ji(
                                                    palm_grid_lat, palm_grid_lon)
                    regridder = palm_wrf_utils.BilinearRegridder(
                            palm_in_camx_x, palm_in_camx_y, preloaded=True)

                    # locate layer heights
                    try:
                        vz = f.variables['z']
                    except KeyError:
                        print('Loading heights from separate file')
                        with open(fname+'.heights') as fh:
                            fix_hgt = np.array(list(map(float, re_num.findall(fh.read())))) * 1000. #orig in km
                            fix_hgt = fix_hgt[:,_na,_na]
                    else:
                        print('Loading heights from variable z')
                        fix_hgt = None

                    for itout, itf in dts:
                        print('Timestep {0}'.format(itout))
                        vars_remaining = [vn for vn, vf in filled.items() if not vf[itout]]

                        lay_height = fix_hgt if fix_hgt is not None else (
                                regridder.regrid(vz[itf,:,regridder.ys,regridder.xs]))
                        lay_height = np.r_[lowest_layer, lay_height + terrain_shift]
                                                #add 1 pegging layer always below
                        process_tstep(f, itf, regridder, lay_height, fout, itout, z_levels,
                                vars_remaining, filled, conversions, helpers)
                else:
                    print('Skipping CAMx file {0} - no required times.'.format(fname))

    if not all(all(vf) for vf in filled.values()):
        sys.exit('CAMx data not complete - missing some variables/timesteps: {0}'
                .format(filled))
