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
from .config import cfg, parse_duration

class RuntimeObj(object):
    """An object for holding runtime-related values.

    May be nested.
    """
    pass

def basic_init(rt):
    """Performs initializaiton of basic values from config."""

    # Times
    rt.simulation = RuntimeObj()
    rt.simulation.timestep = parse_duration(cfg.simulation, 'timestep')
    rt.simulation.length = parse_duration(cfg.simulation, 'length')

    # Paths
    rt.paths = RuntimeObj()
    rt.paths.expand = dict(
            case=cfg.case,
            scenario='.'+cfg.scenario if cfg.scenario else '',
            domain='' if cfg.dnum==1 else '_N{:02d}'.format(cfg.dnum))

    rt.paths.base = cfg.paths.base.format(**rt.paths.expand)
    rt.paths.palm_input = os.path.join(rt.paths.base,
            cfg.paths.palm_input.format(**rt.paths.expand))
    rt.paths.dynamic_driver = os.path.join(rt.paths.palm_input,
            cfg.paths.dynamic_driver.format(**rt.paths.expand))
    rt.paths.intermediate = os.path.join(rt.paths.base,
            cfg.paths.intermediate.format(**rt.paths.expand))
    rt.paths.imported = os.path.join(rt.paths.intermediate,
            cfg.paths.imported.format(**rt.paths.expand))
    rt.paths.hinterp = os.path.join(rt.paths.intermediate,
            cfg.paths.hinterp.format(**rt.paths.expand))
    rt.paths.vinterp = os.path.join(rt.paths.intermediate,
            cfg.paths.vinterp.format(**rt.paths.expand))

    # Domain
    rt.nested_domain = (cfg.dnum > 1)
    rt.stretching = (cfg.domain.dz_stretch_factor != 1.0)

rt = RuntimeObj()
