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

import netCDF4

from . import plugins as plg
from .logging import die, warn, log, verbose, configure_log
from .config import load_config, cfg
from .runtime import rt, basic_init
from .utils import assert_dir


def build_exec_queue(event, from_plugins):
    # logika vytvareni fronty muze byt slozitejsi nez jen prosty seznam (strom, mozna paralelizace....)
    queue = []
    for plugin in from_plugins:
        if isinstance(plugin, getattr(plg, plg.event_hooks[event]['class'])):
            queue.append(plugin)

    return queue


def execute_event(event):
    queue = build_exec_queue(event, rt.plugins)

    kwargs = {}
    common_files = []
    try:
        # Prepare common files or other common processing for specific events
        if event == 'import':
            assert_dir(rt.paths.intermediate.imported)
            f = netCDF4.Dataset(rt.paths.intermediate.imported, 'w', format='NETCDF4')
            common_files.append(f)
            kwargs['fout'] = f
        elif event == 'hinterp':
            assert_dir(rt.paths.intermediate.hinterp)
            f = netCDF4.Dataset(rt.paths.intermediate.hinterp, 'w', format='NETCDF4')
            common_files.append(f)
            kwargs['fout'] = f
        elif event == 'vinterp':
            assert_dir(rt.paths.intermediate.vinterp)
            f = netCDF4.Dataset(rt.paths.intermediate.vinterp, 'w', format='NETCDF4')
            common_files.append(f)
            kwargs['fout'] = f

        # Execute each plugin in a queue
        for plugin in queue:
            getattr(plugin, plg.event_hooks[event]['method'])(**kwargs)
    finally:
        for f in common_files:
            try:
                f.close()
            except:
                warn('Error closing file {}!', f)


def run(argv):
    # Set initial verbosity from commandline, so that we can log the
    # configuration progress appropriately.
    cfg._settings['verbosity'] = (argv.verbosity_arg
                                  if argv.verbosity_arg is not None else 1)
    configure_log(cfg)

    # Load all configfiles and apply commandline config
    load_config(argv)

    # Configure logging according to final config
    configure_log(cfg)

    # Runtime data
    basic_init(rt)

    # Load plugins as configured
    rt.plugins = [plg.plugin_factory(p, cfg=cfg, rt=rt)
                      for p in cfg.plugins]

    # Execute all stages in the workflow
    for event in cfg.workflow:
        execute_event(event)
    log('Finished all stages in the workflow.')
