#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Library functions for plugins"""

import datetime
import re
import numpy as np
from core.logging import die, warn, log, verbose, log_output
from core.config import cfg

class PalmPhysics:
    """Physics calculations with defined constants

    This class contains various physics calculations, implemented as class
    methods, which use defined physical constants. It can be subclassed with
    different constants to suit processing of data related to models with
    constants defined otherwise.
    """

    # Constants directly from PALM code:
    g = 9.81                # gravitational acceleration (m s-2)
    c_p = 1005.             # heat capacity of dry air (J kg-1 K-1)
    r_d = 287.              # sp. gas const. dry air (J kg-1 K-1)
                            # (identical in PALM and WRF)
    p0  = 1e5               # standard pressure reference state
    sigma_sb = 5.67037e-08  # Stefan-Boltzmann constant

    d_p0 = 1e-5             # precomputed 1 / p0
    g_d_cp  = g   / c_p     # precomputed g / c_p
    cp_d_rd = c_p / r_d     # precomputed c_p / r_d
    rd_d_cp = r_d / c_p     # precomputed r_d / c_p

    # Other generic constants
    radius = 6371.          # mean radius of earth
    R = 8.314               # ideal gas constant (m3⋅Pa⋅K−1⋅mol−1)

    @classmethod
    def barom_lapse0_pres(cls, p0, gp, gp0, t0):
        """Converts pressure based on geopotential using barometric equation"""

        barom = 1. / (cls.r_d * t0)
        return p0 * np.exp((gp0-gp)*barom)

    @classmethod
    def barom_lapse0_gp(cls, gp0, p, p0, t0):
        """Converts geopotential based on pressure using barometric equation"""

        baromi = cls.r_d * t0
        return gp0 - np.log(p/p0) * baromi

    @classmethod
    def barom_ptn_pres(cls, p0, z, t0):
        """Compute the barometric formula for 1-D array arguments.

        The calculation is based on the assumption of a polytropic atmosphere
        and neutral stratification, where the temperature lapse rate is g/cp.
        """

        return  p0 * (1. - z*(cls.g_d_cp/t0))**cls.cp_d_rd

    @classmethod
    def exner(cls, p):
        """Exner function"""

        return (p*cls.d_p0)**cls.rd_d_cp

    @classmethod
    def exner_inv(cls, p):
        """Reciprocal of the Exner function"""

        return (cls.p0/p)**cls.rd_d_cp


utc = datetime.timezone.utc
utcdefault = lambda dt: dt.replace(tzinfo=utc) if dt.tzinfo is None else dt
