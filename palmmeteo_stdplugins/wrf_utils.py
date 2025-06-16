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

'''WRF (+CAMx) utility module for PALM dynamic driver generator'''

import math
import numpy as np
import pyproj
import scipy.ndimage as ndimage
import netCDF4
from palmmeteo.utils import SliceExtender
from palmmeteo.library import PalmPhysics

import metpy
metpy_version_master = int(metpy.__version__.split('.', 1)[0])
import metpy.calc as mpcalc
from metpy.interpolate import log_interpolate_1d
from metpy.units import units

# User-selectable values FIXME: move to config

# Settings for geostrophic wind
gw_gfs_margin_deg = 5. #smoothing area in degrees lat/lon
gw_wrf_margin_km = 10. #smoothing area in km
#gw_alpha = .143 #GW vertical interpolation by power law
gw_alpha = 1. #ignore wind power law, interpolate linearly

ax_ = np.newaxis

class WrfPhysics(PalmPhysics):
    # Constants directly equivalent to WRF code
    radius = 6370000.0
    base_temp = 300.        # NOT the variable T00 from WRFOUT!
    rd_d_cp = 2./7.         # from WRF v4 technote (R_d / c_p)


class WRFCoordTransform(object):
    'Coordinate transformer for WRFOUT files'

    def __init__(self, ncf):
        attr = lambda a: getattr(ncf, a)

        # Define grids

        # see http://www.pkrc.net/wrf-lambert.html
        #latlon_wgs84 = pyproj.Proj(proj='latlong',
        #    ellps='WGS84', datum='WGS84',
        #    no_defs=True) #don't use - WRF datum misuse

        self.latlon_sphere = pyproj.Proj(proj='latlong',
            a=WrfPhysics.radius, b=WrfPhysics.radius,
            towgs84='0,0,0', no_defs=True)

        self.lambert_grid = pyproj.Proj(proj='lcc',
            lat_1=attr('TRUELAT1'),
            lat_2=attr('TRUELAT2'),
            lat_0=attr('MOAD_CEN_LAT'),
            lon_0=attr('STAND_LON'),
            a=WrfPhysics.radius, b=WrfPhysics.radius,
            towgs84='0,0,0', no_defs=True)

        # resoltion in m
        self.dx = dx = attr('DX')
        self.dy = dy = attr('DY')

        # number of mass grid points
        self.nx = nx = attr('WEST-EAST_GRID_DIMENSION') - 1
        self.ny = ny = attr('SOUTH-NORTH_GRID_DIMENSION') - 1

        # distance between centers of mass grid points at edges
        extent_x = (nx - 1) * dx
        extent_y = (ny - 1) * dy

        # grid center in lambert
        center_x, center_y = pyproj.transform(self.latlon_sphere, self.lambert_grid,
            attr('CEN_LON'), attr('CEN_LAT'))

        # grid origin coordinates in lambert
        self.i0_x = center_x - extent_x*.5
        self.j0_y = center_y - extent_y*.5

    def latlon_to_ji(self, lat, lon):
        x, y = pyproj.transform(self.latlon_sphere, self.lambert_grid,
                lon, lat)
        return (y-self.j0_y)/self.dy, (x-self.i0_x)/self.dx

    def ji_to_latlon(self, j, i):
        lon, lat = pyproj.transform(self.lambert_grid, self.latlon_sphere,
            i*self.dx+self.i0_x, j*self.dy+self.j0_y)
        return lat, lon

    def verify(self, ncf):
        lat = ncf.variables['XLAT'][0]
        lon = ncf.variables['XLONG'][0]
        j, i = np.mgrid[0:self.ny, 0:self.nx]

        jj, ii = self.latlon_to_ji(lat, lon)
        d = np.hypot(jj-j, ii-i)
        print('error for ll->ji: max {0} m, avg {1} m.'.format(d.max(), d.mean()))

        llat, llon = self.ji_to_latlon(j, i)
        d = np.hypot(llat - lat, llon - lon)
        print('error for ji->ll: max {0} deg, avg {1} deg.'.format(d.max(), d.mean()))

        lat = ncf.variables['XLAT_U'][0]
        lon = ncf.variables['XLONG_U'][0]
        j, i = np.mgrid[0:self.ny, 0:self.nx+1]
        jj, ii = self.latlon_to_ji(lat, lon)
        ii = ii + .5
        d = np.hypot(jj-j, ii-i)
        print('error for U-staggered ll->ji: max {0} m, avg {1} m.'.format(d.max(), d.mean()))

class CAMxCoordTransform(WRFCoordTransform):
    'Coordinate transformer for CAMx files running from WRF'

    def __init__(self, ncf):
        attr = lambda a: getattr(ncf, a)

        # Define grids

        self.latlon_sphere = pyproj.Proj(proj='latlong',
            a=WrfPhysics.radius, b=WrfPhysics.radius,
            towgs84='0,0,0', no_defs=True)

        self.lambert_grid = pyproj.Proj(proj='lcc',
            lat_1=attr('P_ALP'),
            lat_2=attr('P_BET'),
            lat_0=attr('YCENT'),
            lon_0=attr('P_GAM'),
            a=WrfPhysics.radius, b=WrfPhysics.radius,
            towgs84='0,0,0', no_defs=True)

        # resoltion in m
        self.dx = attr('XCELL')
        self.dy = attr('YCELL')

        # number of mass grid points
        self.nx = attr('NCOLS')
        self.ny = attr('NROWS')

        # grid origin coordinates in lambert
        self.i0_x = attr('XORIG')
        self.j0_y = attr('YORIG')

    def verify(self, ncf):
        lat = ncf.variables['latitude'][:]
        lon = ncf.variables['longitude'][:]
        j, i = np.mgrid[0:self.ny, 0:self.nx]

        jj, ii = self.latlon_to_ji(lat, lon)
        d = np.hypot(jj-j, ii-i)
        print('error for ll->ji: max {0} m, avg {1} m.'.format(d.max(), d.mean()))

        llat, llon = self.ji_to_latlon(j, i)
        d = np.hypot(llat - lat, llon - lon)
        print('error for ji->ll: max {0} deg, avg {1} deg.'.format(d.max(), d.mean()))

class BilinearRegridder(object):
    '''Bilinear regridder for multidimensional data.

    By standard, the last two dimensions are always Y,X in that order.
    '''
    def __init__(self, projected_x, projected_y, preloaded=False):
        projected_x = np.asanyarray(projected_x)
        projected_y = np.asanyarray(projected_y)
        self.shape = projected_x.shape
        self.rank = len(self.shape)
        assert self.shape == projected_y.shape

        y0 = np.floor(projected_y)
        yd = projected_y - y0
        ydd = 1. - yd
        self.y0 = y0.astype('i8')

        x0 = np.floor(projected_x)
        xd = projected_x - x0
        xdd = 1. - xd
        self.x0 = x0.astype('i8')

        if preloaded:
            # Prepare slices for preloading from NetCDF files (for cases where
            # the range of loaded Y, X coordinates is much less than total
            # size. The regrid method then expects preloaded data.

            ybase = self.y0.min()
            ytop = self.y0.max()+2
            assert 0 <= ybase <= ytop, (0, ybase, ytop)
            self.ys = slice(ybase, ytop)
            self.y0 -= ybase
            self.ylen = ytop - ybase

            xbase = self.x0.min()
            xtop = self.x0.max()+2
            assert 0 <= xbase <= xtop, (0, xbase, xtop)
            self.xs = slice(xbase, xtop)
            self.x0 -= xbase
            self.xlen = xtop - xbase

        self.y1 = self.y0 + 1
        self.x1 = self.x0 + 1

        self.weights = np.array([
            ydd * xdd, #wy0x0
            ydd * xd , #wy0x1
            yd  * xdd, #wy1x0
            yd  * xd , #wy1x1
            ])

    def loader(self, obj):
        """Prepares a slicing object which automatically adds selector indices
        for this regridder.
        """
        return SliceExtender(obj, self.ys, self.xs)

    def regrid(self, data):
        # data may contain additional dimensions (before Y,X)
        dshape = data.shape[:-2]
        drank = len(dshape)

        # Prepare array for selected data
        sel_shape = (4,) + dshape + self.shape
        selection = np.empty(sel_shape, dtype=data.dtype)

        selection[0, ...] = data[..., self.y0, self.x0]
        selection[1, ...] = data[..., self.y0, self.x1]
        selection[2, ...] = data[..., self.y1, self.x0]
        selection[3, ...] = data[..., self.y1, self.x1]

        # Slice weights to match the extra dimensions
        wslice = ((slice(None),) +      #weights
            (ax_,) * drank +            #data minus Y,X
            (slice(None),) * self.rank) #regridded shape

        w = selection * self.weights[wslice]
        return w.sum(axis=0)

def calc_ph_hybrid(f, it, mu):
    pht = f.variables['P_TOP'][it]
    c3f = f.variables['C3F'][it]
    c4f = f.variables['C4F'][it]
    c3h = f.variables['C3H'][it]
    c4h = f.variables['C4H'][it]
    return (c3f[:,ax_,ax_]*mu[ax_,:,:] + (c4f[:,ax_,ax_] + pht),
            c3h[:,ax_,ax_]*mu[ax_,:,:] + (c4h[:,ax_,ax_] + pht))

def calc_ph_sigma(f, it, mu):
    pht = f.variables['P_TOP'][it]
    eta_f = f.variables['ZNW'][it]
    eta_h = f.variables['ZNU'][it]
    return (eta_f[:,ax_,ax_]*mu[ax_,:,:] + pht,
            eta_h[:,ax_,ax_]*mu[ax_,:,:] + pht)

def wrf_t(f, it):
    p = f.variables['P'][it,:,:,:] + f.variables['PB'][it,:,:,:]
    return (f.variables['T'][it,:,:,:] + WrfPhysics.base_temp) * WrfPhysics.exner(p)

def calc_gp(f, it, ph):
    terr = f.variables['HGT'][it,:,:]
    gp0 = terr * WrfPhysics.g
    gp = [gp0]
    t = wrf_t(f, it)
    for lev in range(1, ph.shape[0]):
        gp.append(WrfPhysics.barom_lapse0_gp(gp[-1], ph[lev,:,:], ph[lev-1,:,:], t[lev-1,:,:]))
    return np.array(gp)

def palm_wrf_gw(f, lon, lat, levels, tidx=0):
    '''Calculate geostrophic wind from WRF using metpy'''

    hgts, ug, vg = calcgw_wrf(f, lat, lon, levels, tidx)

    # extrapolate at the bottom
    hgts = np.r_[np.array([0.]), hgts]
    ug = np.r_[ug[0], ug]
    vg = np.r_[vg[0], vg]

    return minterp(levels, hgts, ug, vg)


def minterp(interp_heights, data_heights, u, v):
    '''Interpolate wind using power law for agl levels'''

    pdata = data_heights ** gw_alpha
    pinterp = interp_heights ** gw_alpha
    hindex = np.searchsorted(data_heights, interp_heights, side='right')
    lindex = hindex - 1
    assert lindex[0] >= 0
    assert hindex[-1] < len(data_heights)
    lbound = pdata[lindex]
    hcoef = (pinterp - lbound) / (pdata[hindex] - lbound)
    #print(data_heights)
    #print(lindex)
    #print(hcoef)
    lcoef = 1. - hcoef
    iu = u[lindex] * lcoef + u[hindex] * hcoef
    iv = v[lindex] * lcoef + v[hindex] * hcoef
    return iu, iv

def get_wrf_dims(f, lat, lon, xlat, xlong):
    '''A crude method, yet satisfactory for approximate WRF surroundings'''

    sqdist = (xlat - lat)**2 + (xlong - lon)**2
    coords = np.unravel_index(sqdist.argmin(), sqdist.shape)

    xmargin = int(math.ceil(gw_wrf_margin_km * 1000 / f.DX)) #py2 ceil produces float
    ymargin = int(math.ceil(gw_wrf_margin_km * 1000 / f.DY))
    y0, y1 = coords[0] - ymargin, coords[0] + ymargin
    x0, x1 = coords[1] - xmargin, coords[1] + xmargin
    assert 0 <= y0 < y1 < sqdist.shape[0], "Point {0} + surroundings not inside domain".format(coords[0])
    assert 0 <= x0 < x1 < sqdist.shape[1], "Point {0} + surroundings not inside domain".format(coords[1])

    return coords, (slice(y0, y1+1), slice(x0, x1+1)), (ymargin, xmargin)

def calcgw_wrf(f, lat, lon, levels, tidx=0):
    # MFDataset removes the time dimension from XLAT, XLONG
    xlat = f.variables['XLAT']
    xlslice = (0,) * (len(xlat.shape)-2) + (slice(None), slice(None))
    xlat = xlat[xlslice]
    xlong = f.variables['XLONG'][xlslice]

    (iy, ix), area, (iby, ibx) = get_wrf_dims(f, lat, lon, xlat, xlong)
    areat = (tidx,) + area
    areatz = (tidx, slice(None)) + area
    #print('wrf coords', lat, lon, xlat[iy,ix], xlong[iy,ix])
    #print(xlat[area][iby,ibx], xlong[area][iby,ibx], areat)

    # load area
    hgt = (f.variables['PH'][areatz] + f.variables['PHB'][areatz]) / 9.81
    hgtu = (hgt[:-1] + hgt[1:]) * .5
    pres = f.variables['P'][areatz] + f.variables['PB'][areatz]
    terrain = f.variables['HGT'][areat]

    # find suitable pressure levels
    yminpres, xminpres = np.unravel_index(pres[0].argmin(), pres[0].shape)
    pres1 = pres[0, yminpres, xminpres] - 1.

    aglpt = hgtu[:,iby,ibx] - terrain[iby,ibx]
    pres0 = pres[np.searchsorted(aglpt, levels[-1]), iby, ibx]
    plevels = np.arange(pres1, min(pres0, pres1)-1, -1000.)

    # interpolate wrf into pressure levels
    phgt = log_interpolate_1d(plevels, pres, hgtu, axis=0)

    # lat_lon_grid_deltas doesn't work under py2, but for WRF grid it is still
    # not very accurate, better use direct values.
    #dx, dy = mpcalc.lat_lon_grid_deltas(xlong[area], xlat[area])
    dx = f.DX * units.m
    dy = f.DY * units.m

    if metpy_version_master >= 1:
        mylat = np.deg2rad(xlat[area])
        my_geostrophic_wind = lambda sh: mpcalc.geostrophic_wind(sh, dx=dx,
                dy=dy, latitude=mylat)
    else:
        coriol = mpcalc.coriolis_parameter(np.deg2rad(xlat[area])).to('1/s')
        my_geostrophic_wind = lambda sh: mpcalc.geostrophic_wind(sh, coriol,
                dx, dy)

    # Smooth height data. Sigma=1.5 for gfs 0.5deg
    res_km = f.DX / 1000.

    ug = np.zeros(plevels.shape, 'f8')
    vg = np.zeros(plevels.shape, 'f8')
    for i in range(len(plevels)):
        sh = ndimage.gaussian_filter(phgt[i,:,:], sigma=1.5*50/res_km, order=0)
        ugl, vgl = my_geostrophic_wind(sh * units.m)
        ug[i] = ugl[iby, ibx].magnitude
        vg[i] = vgl[iby, ibx].magnitude

    return phgt[:,iby,ibx], ug, vg

# The following two functions calculate GW from GFS files, although this
# function is currently not implemented in PALM dynamic driver generation
# script

def calcgw_gfs(v, lat, lon):
    height, lats, lons = v.data(lat1=lat-gw_gfs_margin_deg ,lat2=lat+gw_gfs_margin_deg,
            lon1=lon-gw_gfs_margin_deg, lon2=lon+gw_gfs_margin_deg)
    i = np.searchsorted(lats[:,0], lat)
    if abs(lats[i+1,0] - lat) < abs(lats[i,0] - lat):
        i = i+1
    j = np.searchsorted(lons[0,:], lon)
    if abs(lons[0,i+1] - lon) < abs(lons[0,i] - lon):
        j = j+1
    #print('level', v.level, 'height', height[i,j], lats[i,j], lons[i,j])

    # Set up some constants based on our projection, including the Coriolis
    # parameter and grid spacing, converting lon/lat spacing to Cartesian
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    res_km = (dx[i,j]+dy[i,j]).magnitude / 2000.

    # Smooth height data. Sigma=1.5 for gfs 0.5deg
    height = ndimage.gaussian_filter(height, sigma=1.5*50/res_km, order=0)

    if metpy_version_master >= 1:
        geo_wind_u, geo_wind_v = mpcalc.geostrophic_wind(height * units.m,
                dx=dx, dy=dy, latitude=np.deg2rad(lats))
    else:
        f = mpcalc.coriolis_parameter(np.deg2rad(lats)).to('1/s')

        # In MetPy 0.5, geostrophic_wind() assumes the order of the dimensions
        # is (X, Y), so we need to transpose from the input data, which are
        # ordered lat (y), lon (x).  Once we get the components,transpose again
        # so they match our original data.
        geo_wind_u, geo_wind_v = mpcalc.geostrophic_wind(height * units.m, f, dx, dy)

    return height[i,j], geo_wind_u[i,j], geo_wind_v[i,j]

def combinegw_gfs(grbs, levels, lat, lon):
    heights = []
    us = []
    vs = []
    for grb in grbs:
        h, u, v = calcgw_gfs(grb, lat, lon)
        heights.append(h)
        us.append(u.magnitude)
        vs.append(v.magnitude)
    heights = np.array(heights)
    us = np.array(us)
    vs = np.array(vs)

    ug, vg = minterp(np.asanyarray(levels), heights[::-1], us[::-1], vs[::-1])
    return ug, vg

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-w', '--wrfout', help='verify wrfout file')
    parser.add_argument('-c', '--camx', help='verify camx file')
    args = parser.parse_args()

    if args.wrfout:
        f = netCDF4.Dataset(args.wrfout)

        print('Verifying coord transform:')
        t = WRFCoordTransform(f)
        t.verify(f)

        print('\nVerifying vertical levels:')
        mu = f.variables['MUB'][0,:,:] + f.variables['MU'][0,:,:]
        gp = f.variables['PH'][0,:,:,:] + f.variables['PHB'][0,:,:,:]

        print('\nUsing sigma:')
        phf, phh = calc_ph_sigma(f, 0, mu)
        gp_calc = calc_gp(f, 0, phf)
        delta = gp_calc - gp
        for lev in range(delta.shape[0]):
            print_dstat(lev, delta[lev])

        print('\nUsing hybrid:')
        phf, phh = calc_ph_hybrid(f, 0, mu)
        gp_calc = calc_gp(f, 0, phf)
        delta = gp_calc - gp
        for lev in range(delta.shape[0]):
            print_dstat(lev, delta[lev])

        f.close()

    if args.camx:
        f = netCDF4.Dataset(args.camx)
        t = CAMxCoordTransform(f)
        t.verify(f)
        f.close()
