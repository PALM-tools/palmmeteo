#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
