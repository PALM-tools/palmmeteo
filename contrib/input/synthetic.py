from datetime import datetime, timezone, timedelta
import numpy as np
import scipy.ndimage as ndimage
from metpy.interpolate import interpolate_1d

from core.plugins import ImportPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose, log_output
from core.config import cfg, ConfigError
from core.runtime import rt
from core.utils import ensure_dimension
from core.library import barom_pres, barom_gp, utcdefault


rd_cp = 2./7. #from WRF v4 technote (R_d / c_p)
g = 9.81 #m/s2
_na = np.newaxis

def expand(var):
    """expand horizontal dimensions and pad at top and bottom"""
    v = np.zeros((var.shape[0]+2, rt.ny, rt.nx), dtype=var.dtype)
    v[0,:,:] = var[0]
    v[1:-1,:,:] = var[:,_na,_na]
    v[-1,:,:] = var[-1,_na,_na]
    return v

class SyntheticPlugin(ImportPluginMixin, VInterpPluginMixin):
    def import_data(self, fout, *args, **kwargs):
        log('Importing synthetic profiles')
        rt.times = rt.simulation.start_time + np.arange(rt.nt) * rt.simulation.timestep

        rt.synth_profiles = {vn: ProfileInterpolator(vn, vc)
                for vn, vc in cfg.synthetic.prof_vars}

        rt.z_soil_levels = rt.synth_profiles['soil_t'].heights
        rt.nz_soil = len(rt.z_soil_levels)

    def interpolate_vert(self, fout, *args, **kwargs):
        log('Performing vertical interpolation of synthetic profiles')

        verbose('Preparing output file')
        ensure_dimension(fout, 'time', rt.nt)
        ensure_dimension(fout, 'x', rt.nx)
        ensure_dimension(fout, 'y', rt.ny)
        ensure_dimension(fout, 'z', rt.nz)
        ensure_dimension(fout, 'zw', rt.nz-1)
        ensure_dimension(fout, 'zsoil', rt.nz_soil)

        fout.createVariable('init_atmosphere_qv', 'f4', ('time', 'z', 'y', 'x'))
        fout.createVariable('init_atmosphere_pt', 'f4', ('time', 'z', 'y', 'x'))
        fout.createVariable('init_atmosphere_u', 'f4', ('time', 'z', 'y', 'x'))
        fout.createVariable('init_atmosphere_v', 'f4', ('time', 'z', 'y', 'x'))
        fout.createVariable('init_atmosphere_w', 'f4', ('time', 'zw', 'y', 'x'))
        fout.createVariable('surface_forcing_surface_pressure', 'f4', ('time', 'y', 'x'))
        fout.createVariable('init_soil_t', 'f4', ('time', 'zsoil', 'y', 'x'))
        fout.createVariable('init_soil_m', 'f4', ('time', 'zsoil', 'y', 'x'))
        fout.createVariable('zsoil', 'f4', ('zsoil',))
        fout.createVariable('z', 'f4', ('z',))
        fout.createVariable('zw', 'f4', ('zw',))

        fout.variables['z'][:] = rt.z_levels
        fout.variables['zw'][:] = rt.z_levels_stag
        fout.variables['zsoil'][:] = rt.z_soil_levels

        if cfg.vinterp.terrain_smoothing:
            verbose('Smoothing PALM terrain for the purpose of '
                    'dynamic driver with sigma={0} grid '
                    'points.', cfg.vinterp.terrain_smoothing)
            target_terrain = ndimage.gaussian_filter(rt.terrain,
                    sigma=cfg.vinterp.terrain_smoothing, order=0)
        else:
            target_terrain = rt.terrain

        gp0 = rt.origin_z*g
        gp_trans = (rt.origin_z + cfg.synthetic.transition_level) * g
        gp_new_surf = target_terrain * g

        for it in range(rt.nt):
            verbose('Processing timestep {}', it)

            prof = rt.synth_profiles['pt']
            pt = prof.interp_next_time(rt.times[it])

            if cfg.synthetic.p_surf is None:
                p_surf = barom_pres(cfg.synthetic.p_sl, gp0, 0., pt[0])
            else:
                p_surf = cfg.synthetic.p_surf

            t0 = pt[0] * np.power(0.00001*p_surf, rd_cp)

            # Calculate transition pressure level using horizontal
            # domain-wide pressure average
            p_trans = barom_pres(p_surf, gp_trans, gp0, t0)
            verbose('Vertical stretching transition level: {} Pa', p_trans)

            # Calculate terrain pressure shift ratio
            p_surf_new = barom_pres(p_surf, gp_new_surf, gp0, t0)
            terrain_ratio = (p_surf_new - p_trans) / (p_surf - p_trans)

            z = stretch_heigts(prof.heights, gp_new_surf, gp0, t0, p_surf, p_trans)
            fout.variables['init_atmosphere_pt'][it,:,:,:] = interpolate_1d(rt.z_levels, z, expand(pt))
            del pt

            prof = rt.synth_profiles['qv']
            var = expand(prof.interp_next_time(rt.times[it]))
            z = stretch_heigts(prof.heights, terrain_ratio, gp0, t0, p_surf, p_trans)
            fout.variables['init_atmosphere_qv'][it,:,:,:] = interpolate_1d(rt.z_levels, z, var)

            prof = rt.synth_profiles['u']
            var = expand(prof.interp_next_time(rt.times[it]))
            z = stretch_heigts(prof.heights, terrain_ratio, gp0, t0, p_surf, p_trans)
            fout.variables['init_atmosphere_u'][it,:,:,:] = interpolate_1d(rt.z_levels, z, var)

            prof = rt.synth_profiles['v']
            var = expand(prof.interp_next_time(rt.times[it]))
            z = stretch_heigts(prof.heights, terrain_ratio, gp0, t0, p_surf, p_trans)
            fout.variables['init_atmosphere_v'][it,:,:,:] = interpolate_1d(rt.z_levels, z, var)

            prof = rt.synth_profiles['w']
            var = expand(prof.interp_next_time(rt.times[it]))
            z = stretch_heigts(prof.heights, terrain_ratio, gp0, t0, p_surf, p_trans)
            fout.variables['init_atmosphere_w'][it,:,:,:] = interpolate_1d(rt.z_levels_stag, z, var)

            fout.variables['surface_forcing_surface_pressure'][it,:,:] = p_surf_new

            fout.variables['init_soil_t'][it,:,:,:] = (
                    rt.synth_profiles['soil_t'].interp_next_time(rt.times[it])[:,_na,_na])

            fout.variables['init_soil_m'][it,:,:,:] = (
                    rt.synth_profiles['soil_m'].interp_next_time(rt.times[it])[:,_na,_na])

def stretch_heigts(heights, terrain_ratio, gp0, t0, p_surf, p_trans):
    # Convert the geopotentials to pressure naively using barometric equation
    gp_prof = heights*g + gp0
    p_orig = barom_pres(p_surf, gp_prof, gp0, t0)[:,_na,_na]

    # TODO: this may be optimized by finding highest stretched level and
    # caclulating only below that, or by using numexpr
    p_str = (p_orig - p_trans) * terrain_ratio + p_trans

    # Stretch levels to match terrain and keep everthing above transition level
    p_new = np.where(p_orig > p_trans, p_str, p_orig)

    # Calculate new geopotentials
    gp_new = barom_gp(gp0, p_new, p_surf, t0)

    # Calculate new heights
    z = (gp_new-gp0) * (1./g)

    # Pad at top and bottom
    height = np.zeros((z.shape[0]+2,) + z.shape[1:], dtype=z.dtype)
    height[0,:,:] = -999. #always below terrain
    height[1:-1,:,:] = z
    height[-1,:,:] = 100000.
    return height

class ProfileInterpolator:
    def __init__(self, varname, vc):
        if not vc.profiles:
            die('Variable {} has no vertical profiles configured', varname)
        prof = np.array(vc.profiles)
        nprof, nz = prof.shape

        self.heights = np.array(vc.heights)
        if self.heights.shape != (nz,):
            die('Bad length/shape for heights for variable {}', varname)

        if vc.timeseries:
            prof = prof[vc.timeseries]

        times = vc.times
        if times:
            times = list(map(utcdefault, times))
            if times[0] < rt.times[0]:
                # Extrapolate at start
                times[0:0] = rt.times[0:1]
                prof = np.r_[prof[0:1], prof]
        else:
            if nprof == 1:
                # No times = constant profile
                times = [rt.times[0]]
            else:
                die('Times is only allowed to be unspecified if there is '
                    'exactly one profile, variable {} has {}.', varname, nprof)

        if times[-1] < rt.times[-1] or len(times) == 1:
            # Extrapolate at end
            times.append(rt.times[-1] + timedelta(hours=1)) # avoid div-by-0
            prof = np.r_[prof, prof[-1:]]

        self.cur_id = 0
        self.cur_dt = times[0]
        self.next_dt = times[1]
        self.delta = prof[1] - prof[0]

        self.times = times
        self.prof = prof

    def interp_next_time(self, dt):
        if dt > self.next_dt:
            while True:
                self.cur_dt = self.next_dt
                self.cur_id += 1
                self.next_dt = self.times[self.cur_id+1]
                if dt <= self.times[self.next_dt]:
                    break
            self.delta = self.prof[self.cur_id+1] - self.prof[self.cur_id]

        ratio = (dt-self.cur_dt) / (self.next_dt-self.cur_dt)
        return self.prof[self.cur_id] + self.delta*ratio
