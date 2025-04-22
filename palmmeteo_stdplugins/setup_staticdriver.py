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

import os
from datetime import datetime, timezone
import numpy as np
import netCDF4

from palmmeteo.config import cfg, ConfigError
from palmmeteo.runtime import rt
from palmmeteo.plugins import SetupPluginMixin
from palmmeteo.logging import die, warn, log, verbose

na_ = np.newaxis

class StaticDriverPlugin(SetupPluginMixin):
    """Default setup plugin for loading domain info from static driver file."""

    def setup_model(self, *args, **kwargs):
        log('Loading domain info from static driver file {}...', rt.paths.palm_input.static_driver)
        try:
            ncs = netCDF4.Dataset(rt.paths.palm_input.static_driver, 'r')
        except Exception as err:
            die("Error opening static driver file {}: {}", rt.paths.palm_input.static_driver, err)

        # get horizontal structure of the domain
        rt.nx = ncs.dimensions['x'].size
        rt.ny = ncs.dimensions['y'].size
        rt.dx = ncs.variables['x'][:][1] - ncs.variables['x'][:][0]
        rt.dy = ncs.variables['y'][:][1] - ncs.variables['y'][:][0]
        rt.origin_x = ncs.getncattr('origin_x')
        rt.origin_y = ncs.getncattr('origin_y')
        rt.origin_z = ncs.getncattr('origin_z')

        # start_time may be provided in configuration or read from static driver
        if cfg.simulation.origin_time:
            rt.simulation.start_time = cfg.simulation.origin_time
        else:
            dt = ncs.origin_time
            dts = dt.split()
            if len(dts) == 3:
                # extra timezone string
                if len(dts[2]) == 3:
                    # need to add zeros for minutes, otherwise datetime refuses
                    # to parse
                    dt += '00'
                rt.simulation.start_time = datetime.strptime(dt, '%Y-%m-%d %H:%M:%S %z')
            else:
                rt.simulation.start_time = datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')
        if rt.simulation.start_time.tzinfo is None:
            rt.simulation.start_time = rt.simulation.start_time.replace(
                    tzinfo=timezone.utc)

        # create vertical structure of the domain
        rt.dz = cfg.domain.dz
        if not rt.dz:
            log('dz not set: using dx value ({}) as dz.', rt.dx)
            rt.dz = rt.dx
        rt.nz = cfg.domain.nz
        if not rt.nz:
            raise ConfigError('nz > 0 needs to be specified', cfg.domain, 'nz')

        # read terrain height (relative to origin_z) and
        # calculate and check the height of the surface canopy layer
        if 'zt' in ncs.variables.keys():
            rt.terrain_rel = ncs.variables['zt'][:]
        else:
            rt.terrain_rel = np.zeros([rt.ny,rt.nx])

        # Check terrain
        terrain_min = rt.terrain_rel.min()
        if rt.nested_domain and terrain_min != 0:
            warn('The lowest point of the terrain variable zt in the parent '
                'domain is {} (relative to origin_z={}). Please check '
                'that the lowest point of ALL domains equals zero, otherwise PALM '
                'shifts the terrain to ensure that, which can lead to vertical '
                'mismatching with the dynamic driver.', rt.origin_z, terrain_min)

        # Calculate terrain height in integral grid points, which is also the
        # k-coordinate of the lowest air-cell.
        # NOTE: PALM assigns terrain to those grid cells whose center lies on
        # or below terrain (assuming that it is not shifted due to the lowest
        # point not being 0).
        rt.th = np.floor(rt.terrain_rel / rt.dz + 0.5).astype('i8')
        rt.terrain_mask = np.arange(rt.nz)[:,na_,na_] < rt.th[na_,:,:]

        # building height
        if 'buildings_3d' in ncs.variables.keys():
            #print(np.argmax(a != 0, axis=0)) #### FIXME: what is this?
            bh3 = ncs.variables['buildings_3d'][:]
            # minimum index of nonzeo value along inverted z
            rt.bh = np.argmax(bh3[::-1], axis=0)
            # inversion back and masking grids with no buildings
            rt.bh = bh3.shape[0] - rt.bh
            rt.bh[np.max(bh3, axis=0) == 0] = 0
        elif 'buildings_2d' in ncs.variables.keys():
            rt.bh = ncs.variables['buildings_2d'][:]
            rt.bh[rt.bh.mask] = 0
            rt.bh = np.ceil(rt.bh / rt.dz)
        else:
            rt.bh = np.zeros([rt.ny,rt.nx])

        # plant canopy height
        if 'lad' in ncs.variables.keys():
            lad3 = ncs.variables['lad'][:]
            # replace non-zero values with 1
            lad3[lad3 != 0] = 1
            # minimum index of nonzeo value along inverted z
            rt.lad = np.argmax(lad3[::-1], axis=0)
            # inversion back and masking grids with no buildings
            rt.lad = lad3.shape[0] - rt.lad
            rt.lad[np.max(lad3, axis=0) == 0] = 0
        else:
            rt.lad = np.zeros([rt.ny,rt.nx])

        # calculate maximum of surface canopy layer
        nscl = max(np.amax(rt.th+rt.bh),np.amax(rt.th+rt.lad))

        # check nz with ncl
        if rt.nz < nscl + cfg.domain.nscl_free:
            die('nz has to be higher than {}.\nnz={}, dz={}, number of '
                    'scl={}, nscl_free={}', nscl + cfg.domain.nscl_free, rt.nz,
                    rt.dz,  nscl, cfg.domain.nscl_free)
        if (rt.stretching and cfg.domain.dz_stretch_level
                < (nscl + cfg.domain.nscl_free) * rt.dz):
            die('stretching has to start in level above '
                    '{}.\ndz_stretch_level={}, nscl={}, nscl_free={}, dz={}',
                    (nscl + cfg.domain.nscl_free) * rt.dz,
                    cfg.domain.dz_stretch_level, nscl, cfg.domain.nscl_free,
                    rt.dz)
        if 'soil_moisture_adjust' in ncs.variables.keys():
            rt.soil_moisture_adjust = ncs.variables['soil_moisture_adjust'][:]
        else:
            rt.soil_moisture_adjust = np.ones(shape=(rt.ny, rt.nx), dtype=float)

        # geospatial information from static driver
        rt.palm_epsg = int(ncs.variables['crs'].epsg_code.split(':')[-1])

        # close static driver nc file
        ncs.close()
