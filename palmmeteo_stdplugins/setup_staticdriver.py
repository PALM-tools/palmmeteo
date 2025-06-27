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

from datetime import datetime, timezone
import numpy as np
import netCDF4

from palmmeteo.config import cfg, ConfigError
from palmmeteo.runtime import rt
from palmmeteo.plugins import SetupPluginMixin
from palmmeteo.logging import die, warn, log, verbose
from palmmeteo.utils import ax_, where_range

class Building:
    def __init__(self, bid, bids):
        self.id = bid
        building_mask = (bids == bid)
        self.x0, self.x1 = where_range(np.any(building_mask, axis=0))
        self.y0, self.y1 = where_range(np.any(building_mask, axis=1))
        self.slices = (slice(self.y0, self.y1), slice(self.x0, self.x1))
        self.mask = building_mask[self.slices]
        verbose('Building ID {}: {} points between [{}:{},{}:{}].',
                bid, self.mask.sum(),
                self.y0, self.y1, self.x0, self.x1)

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
        rt.th = np.floor(rt.terrain_rel / rt.dz + 0.5).astype('i8') #terrain top
        rt.terrain_mask = np.arange(rt.nz)[:,ax_,ax_] < rt.th[ax_,:,:]

        # Detect individual buildings
        if 'building_id' in ncs.variables:
            building_ids = []
            bids = ncs.variables['building_id'][:]
            building_ids = [Building(bid, bids) for bid
                            in np.ma.unique(bids).compressed()]
            del bids

        # Load buildings
        rt.building_mask = rt.terrain_mask.copy()
        btop = -1 #global max building top
        if 'buildings_3d' in ncs.variables:
            b3ds = ncs.variables['buildings_3d']
            nz_b3d = b3ds.shape[0]
            for bld in building_ids:
                b3d = b3ds[(slice(None),)+bld.slices]
                thmax = rt.th[bld.slices][bld.mask].max() #building terrain top
                bcol = np.any(b3d,axis=(1,2))
                bmax = len(bcol) - np.argmax(bcol[::-1]) + thmax
                verbose('Building ID {} top: {}', bld.id, bmax)
                btop = max(btop, bmax)

                # Put building mask on top of terrain
                rt.building_mask[(slice(0,thmax),)+bld.slices] = 1 #terrain below
                b3d_height = min(rt.nz-thmax, nz_b3d)
                rt.building_mask[(slice(thmax,thmax+b3d_height),)+bld.slices] = b3d[:b3d_height,:,:]
        elif 'buildings_2d' in ncs.variables:
            b2ds = ncs.variables['buildings_2d'][:]
            for bld in building_ids:
                b2d = b2ds[bld.slices]
                bh = np.floor(b2d / rt.dz + 0.5).astype('i8')
                thmax = rt.th[bld.slices][bld.mask].max() #building terrain top
                bmax = bh.max() + thmax
                verbose('Building ID {} top: {}', bld.id, bmax)
                btop = max(btop, bmax)

                # Put building mask on top of terrain
                rt.building_mask[(slice(None),)+bld.slices] = (np.arange(rt.nz)[:,ax_,ax_] <
                                                               (bh+thmax)[ax_,:,:])
        else:
            rt.building_mask = rt.terrain_mask
            btop = rt.th.max()

        # plant canopy height
        if 'lad' in ncs.variables:
            lads = ncs.variables['lad'][:]
            lad_mask = lads > 0
            # minimum index of nonzeo value along inverted z
            rt.lad = lads.shape[0] - np.argmax(lad_mask[::-1], axis=0)
            rt.lad[~np.any(lad_mask, axis=0)] = 0
        else:
            rt.lad = np.zeros([rt.ny,rt.nx])

        # calculate maximum of surface canopy layer
        nscl = max(btop, (rt.th+rt.lad).max())

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
