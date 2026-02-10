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
import os
import pickle

from . import signature
from .logging import die, warn, log, verbose
from .config import cfg, parse_duration, ConfigObj

zstd = None

def myopen(fpath, *args, **kwargs):
    global zstd

    if fpath.endswith('.zst'):
        if zstd is None:
            try:
                from compression import zstd
            except ImportError:
                import pyzstd as zstd
        return zstd.open(fpath, *args, **kwargs)
    else:
        return open(fpath, *args, **kwargs)

class RuntimeObj:
    """An object for holding runtime-related values.

    May be nested.
    """

    def _get_child(self, child_name):
        try:
            return self.__dict__[child_name]
        except KeyError:
            newchild = self.__class__()
            self.__dict__[child_name] = newchild
            return newchild

    def _save(self, fpath):
        log('Saving snapshot to {}', fpath)
        with myopen(fpath, 'wb') as f:
            p = pickle.Pickler(f, protocol=cfg.intermediate_files.pickle_protocol,
                    fix_imports=False)
            p.dump(signature)
            p.dump(self.__dict__)
            if rt.debug.snapshots:
                for k, v in self.__dict__.items():
                    if hasattr(v, 'nbytes'):
                        verbose('{}:\t{} bytes (np).', k, v.nbytes)
                    else:
                        try:
                            l = len(v)
                        except:
                            verbose('{}:\t{} bytes (sys).', k, sys.getsizeof(v))
                        else:
                            verbose('{}:\tlength {}, {} bytes (sys).', k, l, sys.getsizeof(v))
        verbose('Snapshot saved.')

    def _load(self, fpath):
        log('Loading snapshot from {}', fpath)
        with myopen(fpath, 'rb') as f:
            p = pickle.Unpickler(f, fix_imports=False)
            sig_loaded = p.load()
            loaded = p.load()
        if sig_loaded == signature:
            verbose('Loaded snapshot version: {}', sig_loaded)
        else:
            warn('Loaded snapshot version "{}" does not match current '
                 'version "{}", errors may follow!', sig_loaded, signature)
        assert(isinstance(loaded, dict))
        self.__dict__.update(loaded)

def basic_init(rt):
    """Performs initializaiton of basic values from config."""

    # Times
    simulation = rt._get_child('simulation')
    simulation.timestep = parse_duration(cfg.simulation, 'timestep')
    simulation.length = parse_duration(cfg.simulation, 'length')
    if cfg.radiation.timestep == 'auto':
        rt.timestep_rad = None
    else:
        rt.timestep_rad = parse_duration(cfg.radiation, 'timestep')
    simulation.spinup_rad = parse_duration(cfg.radiation, 'spinup_length')

    # Paths
    rt.path_strings = {}
    for key, func in cfg.path_strings:
        try:
            code = compile(func, '<path_string_{}>'.format(key), 'eval')
        except SyntaxError as e:
            die('Syntax error in definition of the path '
                'string {} "{}": {}.', key, func, e)
        try:
            val = eval(code, cfg._settings)
        except Exception as e:
            die('Error while evaluating the path string {} "{}": {}.', key, func, e)
            raise
        rt.path_strings[key] = val

    paths = rt._get_child('paths')
    paths.base = cfg.paths.base.format(**rt.path_strings)
    for sect_name, sect_cfg in cfg.paths:
        if isinstance(sect_cfg, ConfigObj):
            path_sect = paths._get_child(sect_name)
            for key, val in sect_cfg:
                if isinstance(val, str):
                    setattr(path_sect, key,
                        os.path.join(paths.base, val.format(**rt.path_strings)))

    # Domain
    rt.nested_domain = (cfg.dnum > 1)
    rt.stretching = (cfg.domain.dz_stretch_factor != 1.0)

    # Debugging
    debug = rt._get_child('debug')
    isall = cfg.debug.all
    for dname, don in cfg.debug:
        setattr(debug, dname, isall or don)

# Global object
rt = RuntimeObj()
