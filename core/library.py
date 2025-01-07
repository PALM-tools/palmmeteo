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

"""Library functions for plugins"""

import datetime
import numpy as np

rd = 287. #dry air gas constant (J/kg/K), value used directly in WRF

def barom_pres(p0, gp, gp0, t0):
    """Converts pressure based on geopotential using barometric equation"""

    barom = 1. / (rd * t0)
    return p0 * np.exp((gp0-gp)*barom)

def barom_gp(gp0, p, p0, t0):
    """Converts geopotential based on pressure using barometric equation"""

    baromi = rd * t0
    return gp0 - np.log(p/p0) * baromi

utc = datetime.timezone.utc
utcdefault = lambda dt: dt.replace(tzinfo=utc) if dt.tzinfo is None else dt
