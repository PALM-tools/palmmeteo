import numpy as np
from pyproj import Proj, transform

from .plugins import SetupPluginMixin
from .logging import die, warn, log, verbose
from .config import cfg
from .runtime import rt


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
        palm_grid_y, palm_grid_x = np.meshgrid(jrange, irange, indexing='ij')
        palm_grid_lon, palm_grid_lat = transform(inproj, lonlatproj, palm_grid_x, palm_grid_y)

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
        rt.z_levels = np.zeros(rt.nz,dtype=float)
        rt.z_levels_stag = np.zeros(rt.nz-1,dtype=float)
        dzs = rt.dz
        rt.z_levels[0] = dzs/2.0
        for i in range(rt.nz-1):
            rt.z_levels[i+1] = rt.z_levels[i] + dzs
            rt.z_levels_stag[i] = (rt.z_levels[i+1]+rt.z_levels[i])/2.0
            if rt.stretching and rt.z_levels[i+1] + dzs >= cfg.domain.dz_stretch_level:
                dzs = min(dzs * cfg.domain.dz_stretch_factor, cfg.domain.dz_max)
        ztop = rt.z_levels[-1] + dzs / 2.
        verbose('z: {}', rt.z_levels)
        verbose('zw: {}', rt.z_levels_stag)

