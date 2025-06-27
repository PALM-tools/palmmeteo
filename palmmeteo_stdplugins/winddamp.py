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
    conditions. This helps to avoid instabilities in the pressure solver in the
    first timestep.
    """
    def __init__(self):
        match cfg.winddamp.stagger_method:
            case 'minimum':
                self.stagger = np.minimum
            case 'average':
                self.stagger = lambda a, b: (a+b)*.5
            case other:
                die('Unknown winddamp:stagger_method "{}".', other)

    def check_config(self, *args, **kwargs):
        if cfg.winddamp.num_zeroed > cfg.winddamp.damping_dist:
            die('Configuration winddamp:num_zeroed = {} must not be higher '
                'than winddamp:damping_dist = {}.',
                cfg.winddamp.num_zeroed, cfg.winddamp.damping_dist)

    def write_data(self, fout, *args, **kwargs):
        """
        Decrease the written initial condition wind using a calculated factor.
        Exactly num_zeroed cells outside the wall will have a factor of 0, then
        linearly inceasing and all cells beyond damping_dist will have a factor
        of 1.
        """
        log('Processing wind damping near walls')

        ddist = cfg.winddamp.damping_dist
        nzero = cfg.winddamp.num_zeroed

        defint = np.dtype(int) #system's default integer
        maxval = np.iinfo(defint).max-1 #almost maximum (can be increased)

        # Prepare values for damping formula. 
        inside = -nzero #value inside buildings, will increase with
                        #each step towards outside
        maxdist = ddist - nzero + 1 #value for cells beyond damping_dist
        factor = 1. / maxdist

        distances = np.empty((rt.nz, rt.ny, rt.nx), dtype=defint)
        distances[:] = maxval
        distances[rt.building_mask] = inside

        # Iteratively find distance to nearest building
        for i in range(ddist):
            verbose('Finding nearest walls, step {}', i)

            # one step towards:
            distances[:,1:,:]  = np.minimum(distances[:,1:,:],  distances[:,:-1,:]+1) #north
            distances[:,:-1,:] = np.minimum(distances[:,:-1,:], distances[:,1:,:] +1) #south
            distances[:,:,1:]  = np.minimum(distances[:,:,1:],  distances[:,:,:-1]+1) #east
            distances[:,:,:-1] = np.minimum(distances[:,:,:-1], distances[:,:,1:] +1) #west
        np.clip(distances, 0, maxdist, out=distances) #remaining distances capped

        verbose('Applying wind damping')
        dampfact = distances.astype(cfg.output.default_precision) * factor

        # Apply to wind components using staggered coordinates
        u = fout.variables['init_atmosphere_u']
        u[:] = u[:] * self.stagger(dampfact[:,:,:-1], dampfact[:,:,1:])
        v = fout.variables['init_atmosphere_v']
        v[:] = v[:] * self.stagger(dampfact[:,:-1,:], dampfact[:,1:,:])
        w = fout.variables['init_atmosphere_w']
        w[:] = w[:] * self.stagger(dampfact[:-1,:,:], dampfact[1:,:,:])

        verbose('Wind damping finished.')
