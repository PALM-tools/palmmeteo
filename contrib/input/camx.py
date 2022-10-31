from core.plugins import ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose, log_output
from core.config import cfg, ConfigError
from core.runtime import rt

class CAMxPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing CAMx data...')

        filled = [False] * rt.nt
        timesteps = [convertor.new_timestep() for t in rt.times]
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
                # Decode time and locate timestep
                dts = tflag(f.variables['TFLAG'][:], rt.tindex)
                if not dts:
                    print('Skipping CAMx file {0}: no requested times.'.format(fname))
                    continue

                print('Processing CAMx file {0}.'.format(fname))

                if first:
                    # coordinate projection
                    verbose('Loading projection and preparing regridder')
                    rt.trans_camx = CAMxCoordTransform(fin)
                    palm_in_camx_y, palm_in_camx_x = rt.trans_camx.latlon_to_ji(
                                                        rt.palm_grid_lat, rt.palm_grid_lon)
                    rt.regrid_camx = BilinearRegridder(palm_in_camx_x, palm_in_camx_y, preloaded=True)
                    del palm_in_camx_y, palm_in_camx_x

                    convertor = CAMxConvertor(cfg.chem_species,
                            cfg.camx.output_var_defs, cfg.camx.preprocessors,
                            rt.regrid_camx)

                # locate layer heights
                try:
                    vz = fin.variables['z']
                except KeyError:
                    print('Loading heights from separate file')
                    with open(fname+'.heights') as fh:
                        fix_hgt = np.array(list(map(float, re_num.findall(fh.read())))) * 1000. #orig in km
                        zcoord[itout] = fix_hgt[:,_na,_na]
                else:
                    print('Loading heights from variable z')
                    # TODO: check for dicrepancy among input files
                    zcoord[itout] = rt.regrid_camx.regrid(vz[0, :, rt.regrid_camx.ys, rt.regrid_camx.xs])

                if first:
                    # dimensions
                    ensure_dimension(fout, 'time', rt.nt)
                    ensure_dimension(fout, 'z_chem', zcoord[itout].shape[0])
                    ensure_dimension(fout, 'y_chem', rt.regrid_camx.ylen)
                    ensure_dimension(fout, 'x_chem', rt.regrid_camx.xlen)
                    chem_dims = ('time', 'z_chem', 'y_chem', 'x_chem')

                for itout, itf in dts:
                    times[itout] = True
                    verbose('Importing timestep {} -> {}', itf, itout)

                    filled[itout] = convertor.load_timestep_vars(fin, itf,
                            timesteps[itout])
            first = False

        if not all(filled):
            die('Could not find all CAMx variables for all times: {}', filled)

        for i, tsdata in enumerate(timesteps):
            # Save heights
            v_out = (fout.createVariable('height_chem', 'f4', chem_dims) if i
                    else fout.variables['height_chem'])
            v_out.units = 'm'
            v_out[i, :, :, :] = zcoord[i]

            # Save computed variables
            convertor.validate_timestep(tsdata)
            for sn, v, unit in convertor.calc_timestep_species(tsdata):
                v_out = (fout.createVariable(spn, 'f4', chem_dims) if i
                        else fout.variables[spn])
                v_out.units = unit
                v_out[i, :, :, :] = v

    def interpolate_horiz(self, fout, *args, **kwargs):
        log('Performing horizontal interpolation')

        with netCDF4.Dataset(rt.paths.imported) as fin:
            verbose('Preparing output file')
            # Create dimensions
            for d in ['time', 'z_chem']:
                ensure_dimension(fout, d, len(fin.dimensions[d]))
            ensure_dimension(fout, 'x', rt.nx)
            ensure_dimension(fout, 'y', rt.ny)

            # Create variables
            for varname in cfg.chem_species:
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
                for varname in cfg.chem_species:
                    v_in = fin.variables[varname]
                    v_out = fout.variables[varname]
                    v_out[it] = rt.regrid_camx.regrid(v_in[it])

    def interpolate_vert(self, *args, **kwargs):
                    lay_height = np.r_[lowest_layer, lay_height + terrain_shift]
        verbose_dstat = log_dstat_on if cfg.verbosity >= 2 else log_dstat_off

        log('Performing vertical interpolation')

        verbose('Preparing output file')
        with netCDF4.Dataset(rt.paths.hinterp) as fin:
            with netCDF4.Dataset(rt.paths.vinterp, 'w', format='NETCDF4') as fout:
                for dimname in ['time', 'y', 'x', 'zsoil_meteo']:
                    fout.createDimension(dimname, len(fin.dimensions[dimname]))
                fout.createDimension('z', rt.nz)
                fout.createDimension('zw', rt.nz-1)
                fout.createDimension('zsoil', rt.nz_soil)

                fout.createVariable('init_atmosphere_qv', 'f4', ('time', 'z', 'y', 'x'))
                fout.createVariable('init_atmosphere_pt', 'f4', ('time', 'z', 'y', 'x'))
                fout.createVariable('init_atmosphere_u', 'f4', ('time', 'z', 'y', 'x'))
                fout.createVariable('init_atmosphere_v', 'f4', ('time', 'z', 'y', 'x'))
                fout.createVariable('init_atmosphere_w', 'f4', ('time', 'zw', 'y', 'x'))
                fout.createVariable('surface_forcing_surface_pressure', 'f4', ('time', 'y', 'x'))
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

                    var = lpad(fin.variables['T'][it] + wrf_base_temp) #from perturbation pt to standard
                    fout.variables['init_atmosphere_pt'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                    var = lpad(fin.variables['U'][it])
                    fout.variables['init_atmosphere_u'][it,:,:,:] = interpolate_1d(rt.z_levels, height, var)

                    var = lpad(fin.variables['V'][it])
                    fout.variables['init_atmosphere_v'][it,:,:,:]  = interpolate_1d(rt.z_levels, height, var)

                    var = lpad(fin.variables['W'][it]) #z staggered!
                    fout.variables['init_atmosphere_w'][it,:,:,:] = interpolate_1d(rt.z_levels_stag, heightw, var)

                    var = fin.variables['PSFC'][it]
                    fout.variables['surface_forcing_surface_pressure'][it,:,:] = var

                    var = fin.variables['TSLB'][it] #soil temperature
                    fout.variables['init_soil_t'][it,:,:,:] = var

                    var = fin.variables['SMOIS'][it] #soil moisture
                    fout.variables['init_soil_m'][it,:,:,:] = var

                    var = fin.variables['UG'][it]
                    fout.variables['ls_forcing_ug'][it,:] = var

                    var = fin.variables['VG'][it]
                    fout.variables['ls_forcing_vg'][it,:] = var


class CAMxUnitsInfo:
    pass

class CAMxConvertor:
    def __init__(self, species, var_defs, preprocessors, regridder):
        self.regridder = regridder
        self.loaded_vars = set()
        self.preprocessors = {}
        self.validations = {}
        self.vars = {}

        for spn, sp in species:
            try:
                vdef = var_defs[sp]
            except KeyError:
                die('Requested CAMx variable {} not found in configured '
                        'variable definitions.', sp)

            self.loaded_vars.add(sp.loaded_vars)

            for pp in vdef.preprocessors:
                if pp not in self.preprocessors:
                    try:
                        prs = preprocessors[pp]
                    except KeyError:
                        die('Requested CAMx preprocessor {} not found in '
                                'configured variable definitions.', pp)
                    try:
                        self.preprocessors[pp] = compile(prs,
                                '<camx_preprocessor_{}>'.format(pp), 'exec')
                    except SyntaxError:
                        die('Syntax error in definition of the CAMx '
                                'preprocessor {}: "{}".', pp, prs)

            for val in vdef.validations:
                if val not in self.validations:
                    try:
                        self.validations[val] = compile(val,
                                '<camx_validation>', 'eval')
                    except SyntaxError:
                        die('Syntax error in definition of the CAMx '
                                'validation: "{}".', val)

            try:
                fml = compile(sp.formula, '<camx_formula_{}>'.format(spn),
                        'eval')
            except SyntaxError:
                die('Syntax error in definition of the CAMx '
                        'variable {} formula: "{}".', spn, sp.formula)
            self.species[spn] = (fml, sp.unit)

    @staticmethod
    def new_timestep():
        return {'_units': CAMxUnitsInfo()}

    def load_timestep_vars(self, f, tindex, tsdata):
        complete = True

        for vn in self.loaded_vars:
            if vn in tsdata:
                if vn in f.variables:
                    die('Error: duplicate CAMx variable {}.', vn)
            else:
                try:
                    var = f.variables[vn]
                except KeyError:
                    complete = False
                    continue

                tsdata[vn] = var[tindex, ..., self.regridder.ys,
                        self.regridder.xs]

        return complete

    def validate_timestep(self, tsdata):
        for vs, val in self.validations.items():
            if not eval(val, tsdata):
                die('CAMx validation {} failed!', vs)

    def calc_timestep_species(self, tsdata):
        for pp in self.preprocessors.values():
            exec(pp, tsdata)

        for sn, (fml, unit) in self.species:
            v = eval(fml, tsdata)
            yield sn, v, unit


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
