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

from palmmeteo.plugins import WritePluginMixin
from palmmeteo.logging import die, warn, log, verbose
from palmmeteo.config import cfg
from palmmeteo.runtime import rt

class WindDampPlugin(WritePluginMixin):
    """
    A plugin which provides damping of wind near vertical walls in the initial
    conditions. This helps to avoid instabilities in pressure solver in the
    first timestep.
    """
    def write_data(self, fout, *args, **kwargs):
        log('Processing wind damping near walls')

        defint = np.dtype(int) #system's default integer
        maxval = np.iinfo(defint).max-1 #almost maximum (can be increased)
        ddist = cfg.postproc.wind_damping_dist

        distances = np.empty((rt.nz, rt.ny, rt.nx), dtype=defint)
        distances[:] = maxval
        distances[rt.building_mask] = -1 #inside buildings

        # Iteratively find distance to nearest building
        for i in range(ddist):
            verbose('Finding nearest walls, step {}', i)

            # one step towards:
            distances[:,1:,:]  = np.minimum(distances[:,1:,:],  distances[:,:-1,:]+1) #north
            distances[:,:-1,:] = np.minimum(distances[:,:-1,:], distances[:,1:,:] +1) #south
            distances[:,:,1:]  = np.minimum(distances[:,:,1:],  distances[:,:,:-1]+1) #east
            distances[:,:,:-1] = np.minimum(distances[:,:,:-1], distances[:,:,1:] +1) #west
        np.clip(distances, 0, ddist, out=distances) #remaining distances capped

        verbose('Applying wind damping')

        # Damping formula: directly adjacent to building (distance 0) has
        # factor 0, then linearly towards 1 (which is at points further than
        # max distance).
        dampfact = distances.astype(cfg.output.default_precision) / ddist

        # Apply to wind components using staggered coordinates
        u = fout.variables['init_atmosphere_u']
        u[:] = u[:] * ((dampfact[:,:,:-1]+dampfact[:,:,1:]) * .5)
        v = fout.variables['init_atmosphere_v']
        v[:] = v[:] * ((dampfact[:,:-1,:]+dampfact[:,1:,:]) * .5)
        w = fout.variables['init_atmosphere_w']
        w[:] = w[:] * ((dampfact[:-1,:,:]+dampfact[1:,:,:]) * .5)

        verbose('Wind damping finished.')
