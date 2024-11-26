#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Library functions for plugins"""

import datetime
import re
import numpy as np
from .logging import die, warn, log, verbose, log_output
from .config import cfg
from .runtime import rt
from .utils import ax_, rad, SliceExtender
from scipy.spatial import Delaunay

utc = datetime.timezone.utc
utcdefault = lambda dt: dt.replace(tzinfo=utc) if dt.tzinfo is None else dt

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


class UnitConverter:
    loaded = None

    def __init__(self):
       self.re_ppmv = re.compile(cfg.chem_units.regexes.ppmv)
       self.re_ppbv = re.compile(cfg.chem_units.regexes.ppbv)
       self.re_ugm3 = re.compile(cfg.chem_units.regexes.ugm3)
       self.re_gm3  = re.compile(cfg.chem_units.regexes.gm3)
       self.re_kgm3 = re.compile(cfg.chem_units.regexes.kgm3)

    def convert_auto(self, name, value, unit):
        # volumetric fractional
        if self.re_ppmv.match(unit):
            verbose('Unit {} for variable {} understood as ppmv', unit, name)
            return value, cfg.chem_units.targets.ppmv
        if self.re_ppbv.match(unit):
            verbose('Converting {} from {} (understood as ppbv) to ppmv', name, unit)
            return value*1e-3, cfg.chem_units.targets.ppmv

        # mass per volume
        if self.re_ugm3.match(unit):
            verbose('Converting {} from {} (understood as ug/m3) to kg/m3', name, unit)
            return value*1e-9, cfg.chem_units.targets.kgm3
        if self.re_gm3.match(unit):
            verbose('Converting {} from {} (understood as g/m3) to kg/m3', name, unit)
            return value*1e-3, cfg.chem_units.targets.kgm3
        if self.re_kgm3.match(unit):
            verbose('Unit {} for variable {} understood as kg/m3', unit, name)
            return value, cfg.chem_units.targets.kgm3

        # default
        warn('Unknown unit {} for variable {} - keeping.', unit, name)
        return value, unit

    @classmethod
    def convert(cls, name, value, unit):
        if cls.loaded is None:
            cls.loaded = cls()
        return cls.loaded.convert_auto(name, value, unit)

class InputUnitsInfo:
    pass

class LoadedQuantity:
    __slots__ = ['name', 'formula', 'code', 'unit', 'attrs']

class QuantityCalculator:
    def __init__(self, quantities, var_defs, preprocessors, regridder):
        self.regridder = regridder
        self.loaded_vars = set()
        self.preprocessors = {}
        self.validations = {}
        self.quantities = []

        for qname in quantities:
            try:
                vdef = var_defs[qname]
            except KeyError:
                die('Requested quantity {} not found in configured '
                        'quantity definitions.', qname)

            q = LoadedQuantity()
            q.name = qname

            self.loaded_vars.update(vdef.loaded_vars)
            if len(vdef.loaded_vars) == 1 and 'formula' not in vdef:
                # When we have exactly 1 loaded variable and we do not define
                # an explicit formula, we assume that the formula is identity
                # for that variable and the unit is taken from the input
                # variable unless specified otherwise.
                q.formula = vdef.loaded_vars[0]
                q.unit = getattr(vdef, 'unit', None)
                                #None = taken later from loaded variable
            else:
                q.formula = vdef.formula
                q.unit = vdef.unit

            for pp in getattr(vdef, 'preprocessors', []):
                if pp not in self.preprocessors:
                    try:
                        prs = preprocessors[pp]
                    except KeyError:
                        die('Requested input preprocessor {} not found in '
                                'configured variable definitions.', pp)
                    try:
                        self.preprocessors[pp] = compile(prs,
                                '<quantity_preprocessor_{}>'.format(pp), 'exec')
                    except SyntaxError:
                        die('Syntax error in definition of the input '
                                'preprocessor {}: "{}".', pp, prs)

            for val in getattr(vdef, 'validations', []):
                if val not in self.validations:
                    try:
                        self.validations[val] = compile(val,
                                '<quantity_validation>', 'eval')
                    except SyntaxError:
                        die('Syntax error in definition of the input '
                                'validation: "{}".', val)

            q.attrs = {}
            if 'molar_mass' in vdef:
                q.attrs['molar_mass'] = vdef.molar_mass
            for flag in getattr(vdef, 'flags', []):
                q.attrs[flag] = np.int8(1)

            try:
                q.code = compile(q.formula, '<quantity_formula_{}>'.format(qname), 'eval')
            except SyntaxError:
                die('Syntax error in definition of the quantity '
                        '{} formula: "{}".', qname, q.formula)
            self.quantities.append(q)

    @staticmethod
    def new_timestep():
        return {'_units': InputUnitsInfo()}

    def load_timestep_vars(self, f, tindex, tsdata):
        complete = True
        units = tsdata['_units']

        for vn in self.loaded_vars:
            if vn in tsdata:
                if vn in f.variables:
                    die('Error: duplicate input variable {}.', vn)
            else:
                try:
                    var = f.variables[vn]
                except KeyError:
                    complete = False
                    continue

                tsdata[vn] = self.regridder.loader(var)[tindex,...]
                setattr(units, vn, var.units)

        return complete

    def validate_timestep(self, tsdata):
        for vs, val in self.validations.items():
            if not eval(val, tsdata):
                die('Input validation {} failed!', vs)

    def calc_timestep_species(self, tsdata):
        for pp in self.preprocessors.values():
            exec(pp, tsdata)

        for q in self.quantities:
            value = eval(q.code, tsdata)
            unit = q.unit

            # Assign default unit with directly loaded variables
            if unit is None:
                unit = getattr(tsdata['_units'], q.formula)

            # Check for necessary unit conversion
            value, unit = UnitConverter.convert(q.name, value, unit)

            yield q.name, value, unit, q.attrs


def barycentric(tri, pt, isimp):
    """Calculate barycentric coordinates of a multi-dimensional point set
    within a triangulation.

    :param pt:      selection of points (multi-dimensional)
    :param isimp:   selection of simplices (same dims as pt)
    """
    sel_transform = tri.transform[isimp,:,:] #transform(simp, bary, cart) -> (pt, bary, cart)

    # based on help(Delaunay), changing np.dot to selection among dims, using
    # selected simplices
    fact2 = (pt - sel_transform[...,2,:])[...,ax_,:]
    bary0 = (sel_transform[...,:2,:] * fact2).sum(axis=-1) #(pt, bary[:2])

    # add third barycentric coordinate
    bary = np.concatenate((bary0, (1.-bary0.sum(axis=-1))[...,ax_]), axis=-1)
    return bary

class TriRegridder:
    """Universal regridder which uses Delaunay triangulation and barycentric
    coordinate interpolation.

    Works on any grid - triangular, rectangular or unstructured. The only
    requirements are the arrays of latitude and longitude coordinates. The grid
    arrays may be organized as 1- or 2-dimensional.
    """
    def __init__(self, clat, clon, ylat, xlon, buffer):
        #ylat = ylat[0,:5]
        #xlon = xlon[0,:5]
        # Simple Mercator-like stretching for equidistant lat/lon coords
        self.lon_coef = np.cos(ylat.mean()*rad)

        deg_range = buffer / (PalmPhysics.radius*rad)
        lat0 = ylat.min() - deg_range
        lat1 = ylat.max() + deg_range
        deg_range /= self.lon_coef
        lon0 = xlon.min() - deg_range
        lon1 = xlon.max() + deg_range
        verbose(f'Using range lat = {lat0} .. {lat1}, lon = {lon0} .. {lon1}.')

        verbose('Loading coords')
        self.ptmask = (lat0 <= clat) & (clat <= lat1) & (lon0 <= clon) & (clon <= lon1)
        self.npt = self.ptmask.sum()

        verbose(f'Selected {self.npt} out of {len(clat)} points for triangulation.')
        sclat = clat[self.ptmask]
        sclon = clon[self.ptmask]
        sclonx = sclon * self.lon_coef
        tri = Delaunay(np.transpose([sclat, sclonx]))

        # identify simplices
        xlonx = xlon * self.lon_coef
        pt = np.concatenate((ylat[:,:,ax_], xlonx[:,:,ax_]), axis=2)
        isimp = tri.find_simplex(pt)
        assert (isimp >= 0).all()

        self.bary = barycentric(tri, pt, isimp)

        self.simp = tri.simplices[isimp] #(pt,bary)

        # find global coords
        nz = np.nonzero(self.ptmask)[0]
        self.iglob = nz[self.simp]

    def loader(self, obj):
        """Prepares a slicing object which automatically adds selector indices
        for this regridder.
        """
        return SliceExtender(obj, self.ptmask)

    def regrid(self, data):
        """Regrid from point set selected using ptmask"""

        sel_data = data[...,self.simp] #(pt,bary)
        return (sel_data * self.bary).sum(axis=-1)

    def regrid_full(self, data):
        """Regrid from full source point set"""

        sel_data = data[...,self.iglob] #(pt,bary)
        return (sel_data * self.bary).sum(axis=-1)

def verify_palm_hinterp(regridder, lats, lons):
    """Regrids source lat+lon coordinates to PALM coordinates using the regridder and verifies the result."""

    diff = regridder.regrid(lats) - rt.palm_grid_lat
    log('Regridder verification for latitudes:  Error: {:10.6f} .. {:10.6f} '
        '(bias={:10.6f}, MAE={:10.6f}, RMSE={:10.6f}).',
        diff.min(), diff.max(), diff.mean(), np.abs(diff).mean(),
        np.sqrt(np.square(diff).mean()))

    diff = regridder.regrid(lons) - rt.palm_grid_lon
    log('Regridder verification for longitudes: Error: {:10.6f} .. {:10.6f} '
        '(bias={:10.6f}, MAE={:10.6f}, RMSE={:10.6f}).',
        diff.min(), diff.max(), diff.mean(), np.abs(diff).mean(),
        np.sqrt(np.square(diff).mean()))
