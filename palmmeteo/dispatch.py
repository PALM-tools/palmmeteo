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

import sys
from argparse import ArgumentParser
import netCDF4

from . import __doc__, __version__
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


def execute_event(event, from_plugins):
    queue = build_exec_queue(event, from_plugins)

    kwargs = {}
    common_files = []
    try:
        # Prepare common files or other common processing for specific events
        try:
            fn_out = getattr(rt.paths.intermediate, event)
        except AttributeError: pass
        else:
            # Output filename is defined for this stage
            assert_dir(fn_out)
            f = netCDF4.Dataset(fn_out, 'w', format='NETCDF4')
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

    # Save snapshot if applicable
    try:
        fn_snapshot = getattr(rt.paths.snapshot, event)
    except AttributeError: pass
    else:
        assert_dir(fn_snapshot)
        rt._save(fn_snapshot)

def run(argv):
    # Set initial verbosity from commandline, so that we can log the
    # configuration progress appropriately.
    cfg._settings['verbosity'] = (argv.verbosity_arg
                                  if argv.verbosity_arg is not None else 1)
    configure_log(cfg)

    # Load all configfiles and apply commandline config
    workflow = load_config(argv)

    # Configure logging according to final config
    configure_log(cfg)

    # Runtime data
    basic_init(rt)

    # Load plugins as configured
    plugins = [plg.plugin_factory(p, cfg=cfg, rt=rt)
                      for p in cfg.plugins]

    if workflow.snapshot_from:
        try:
            fn_snapshot = getattr(rt.paths.snapshot, workflow.snapshot_from)
        except AttributeError: pass
        else:
            # Workflow not from start - load snapshot
            rt._load(fn_snapshot)

    # Execute all stages in the workflow
    for event in workflow:
        execute_event(event, plugins)

    log('Finished all stages in the workflow.')

def main():
    argp = ArgumentParser(prog='pmeteo', description=__doc__)
    argp.add_argument('-c', '--config', nargs='+', help='configuration file (s)', required=True)
    argp.add_argument('-w', '--workflow-from', help='start workflow at STAGE', metavar='STAGE')
    argp.add_argument('-W', '--workflow-to', help='stop workflow at STAGE', metavar='STAGE')
    argp.add_argument('--version', action='version', version=f'PALM-meteo version {__version__}')
    verbosity = argp.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_const',
            dest='verbosity_arg', const=2, help='increase verbosity')
    verbosity.add_argument('-s', '--silent', action='store_const',
            dest='verbosity_arg', const=0, help='print only errors')

    if len(sys.argv) <= 1:
        argp.print_help()
        sys.exit(2)

    argv = argp.parse_args()
    run(argv)
