#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2018-2026 Institute of Computer Science of the Czech Academy of
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
import numpy as np
import netCDF4
from matplotlib import pyplot as plt
from matplotlib import cm

from palmmeteo.utils import parse_pos, nearest_gridpt
from palmmeteo.plugins import WritePluginMixin
from palmmeteo.logging import die, warn, log, verbose
from palmmeteo.config import cfg
from palmmeteo.runtime import rt

deg_ticks = np.linspace(0, 360, 9)
deg_labels = [f'{t:.0f}°' for t in deg_ticks]

class PlotPlugin(WritePluginMixin):
    """
    A plugin for plotting time series of vertical profiles of various
    meteorological quantities from the inputs to the dynamic driver.
    """

    def check_config(self, *args, **kwargs):
        rt.plot_positions = []
        for pos_cfg in cfg.plot.positions:
            p = [pos_cfg.get('name')]
            if not p[0]:
                die('All plot positions must have a name')

            for pos, ngrid, resol in [(pos_cfg.get('xpos'), rt.nx, rt.dx),
                                      (pos_cfg.get('ypos'), rt.ny, rt.dy)]:
                try:
                    vpos, is_deg = parse_pos(pos, ngrid, resol)
                except ValueError:
                    die('The expression "{}" is not a valid position within domain!', pos)
                if is_deg:
                    die('Specification in degrees is not yet supported for plotting.')
                p.append(nearest_gridpt(vpos, ngrid))

            p.append(pos_cfg.get('vars', []))
            rt.plot_positions.append(p)

    def write_data(self, fout, *args, **kwargs):
        """
        Loads data from the vertical interpolation step (where there are still
        full 3D fields for all timesteps, not just the boundaries) and plots
        them as configured.
        """
        log('Plotting time-series of vertical profiles.')
        dpath = rt.paths.plot.output_dir
        if not os.path.isdir(dpath):
            os.makedirs(dpath)

        with netCDF4.Dataset(rt.paths.intermediate.vinterp) as fin:
            fiv = fin.variables

            for pos_name, xpos, ypos, plot_vars in rt.plot_positions:
                for vn in plot_vars:
                    verbose('Loading profiles of {} at position {} (i={}, j={}).',
                            vn, pos_name, xpos, ypos)
                    if vn=='wind':
                        u = fiv['init_atmosphere_u'][:,:,ypos,xpos]
                        v = fiv['init_atmosphere_v'][:,:,ypos,xpos]
                        wspd = np.hypot(u, v).T
                        wdir = ((np.degrees(np.arctan2(u, v)) + 180.) % 360).T
                        del u, v
                    elif vn=='temp':
                        val = fiv['init_atmosphere_pt'][:,:,ypos,xpos] - 273.15
                        title = f'Potential temperature at {pos_name} (°C)'
                    else:
                        val = fiv['init_atmosphere_'+vn][:,:,ypos,xpos]
                        title = f'{vn} at {pos_name}'

                    verbose('Plotting.')

                    if vn=='wind':
                        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
                        m = ax1.pcolormesh(rt.times, rt.z_levels, wspd)
                        plt.colorbar(m, label='Wind speed (m/s)')
                        ax1.set_ylabel('Height (m above origin_z)')
                        m = ax2.pcolormesh(rt.times, rt.z_levels, wdir,
                                #alpha=np.minimum(1., wspd),
                                cmap=cm.twilight)
                        cb = plt.colorbar(m, label='Wind direction')
                        cb.set_ticks(ticks=deg_ticks, labels=deg_labels)
                        ax2.tick_params('x', rotation=45)
                        ax2.set_ylabel('Height (m above origin_z)')
                        fig.suptitle(f'Wind at {pos_name}')
                    else:
                        fig, ax = plt.subplots(figsize=(8, 6))
                        m = ax.pcolormesh(rt.times, rt.z_levels, val.T)
                        plt.colorbar(m)
                        ax.tick_params('x', rotation=45)
                        ax.set_ylabel('Height (m above origin_z)')
                        ax.set_title(title)

                    fig.savefig(os.path.join(dpath,
                        f'plot_{pos_name}_{vn}.{cfg.plot.format}'),
                        dpi=cfg.plot.dpi)
                    plt.close(fig)

        verbose('Plotting finished.')
