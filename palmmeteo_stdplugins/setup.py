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

import numpy as np
from pyproj import Proj, transform

from palmmeteo.plugins import SetupPluginMixin
from palmmeteo.logging import die, warn, log, verbose
from palmmeteo.config import cfg, ConfigError
from palmmeteo.runtime import rt
from palmmeteo.utils import tstep, td0, DTIndexer

ax_ = np.newaxis

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
        rt.inproj = Proj('+init='+cfg.domain.proj_palm)
        rt.lonlatproj = Proj('+init='+cfg.domain.proj_wgs84)
        rt.cent_lon, rt.cent_lat = transform(rt.inproj, rt.lonlatproj,
                rt.xcent, rt.ycent)
        verbose('xcent={}, ycent={}', rt.xcent, rt.ycent)
        verbose('cent_lon={}, cent_lat={}', rt.cent_lon, rt.cent_lat)
        # prepare target grid
        irange = rt.origin_x + rt.dx * (np.arange(rt.nx, dtype='f8') + .5)
        jrange = rt.origin_y + rt.dy * (np.arange(rt.ny, dtype='f8') + .5)
        rt.palm_grid_y, rt.palm_grid_x = np.meshgrid(jrange, irange, indexing='ij')
        rt.palm_grid_lon, rt.palm_grid_lat = transform(rt.inproj, rt.lonlatproj,
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
                raise ConfigError('Stretch level has to be specified for '
                    'stretching', cfg.domain, 'dz_stretch_level')
            if cfg.domain.dz_max < rt.dz:
                raise ConfigError('dz_max has to be higher or equal than than '
                        'dz (={})'.format(rt.dz), cfg.domain, 'dz_max')
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

        # configure times
        rt.simulation.end_time_rad = rt.simulation.start_time + rt.simulation.length
        rt.tindex = DTIndexer(rt.simulation.start_time, rt.simulation.timestep)
        if rt.nested_domain:
            log('Nested domain - preparing only initialization (1 timestep).')
            rt.nt = 1
            rt.simulation.duration = td0
            rt.simulation.end_time = rt.simulation.start_time
        else:
            rt.simulation.end_time = rt.simulation.end_time_rad
            rt.nt = rt.tindex(rt.simulation.end_time) + 1

        rt.times_sec = np.arange(rt.nt) * rt.simulation.timestep.total_seconds()
        verbose('PALM meteo data extent {} - {} ({} timesteps).',
                rt.simulation.start_time, rt.simulation.end_time, rt.nt)

