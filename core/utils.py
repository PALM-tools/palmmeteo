#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A collection of simple, technical utilities.

These utilities need to be stateless, i.e. they must not depend on config or
runtime.
"""

import os
import re
from datetime import timedelta
import numpy as np
from .logging import die, warn, log, verbose

ax_ = np.newaxis
rad = np.pi / 180.
td0 = timedelta(hours=0)

fext_re = re.compile(r'\.(\d{3})$')

def find_free_fname(fpath, overwrite=False):
    if not os.path.exists(fpath):
        return fpath

    if overwrite:
        log('Existing file {} will be overwritten.', fpath)
        return fpath

    # Try to find free fpath.###
    path, base = os.path.split(fpath)
    nbase = len(base)
    maxnum = -1
    for fn in os.listdir(path):
        if not fn.startswith(base):
            continue
        m = fext_re.match(fn[nbase:])
        if not m:
            continue
        maxnum = max(maxnum, int(m.group(1)))
    if maxnum >= 999:
        raise RuntimeError('Cannot find free filename starting with ' + fpath)

    newpath = '{}.{:03d}'.format(fpath, maxnum+1)
    log('Filename {} exists, using {}.', fpath, newpath)
    return newpath

def tstep(td, step):
    d, m = divmod(td, step)
    if m != td0:
        raise ValueError('Not a whole timestep!')
    return d

def ensure_dimension(f, dimname, dimsize):
    """Creates a dimension in a netCDF file or verifies its size if it already
    exists.
    """
    try:
        d = f.dimensions[dimname]
    except KeyError:
        # Dimension is missing - create it and return
        return f.createDimension(dimname, dimsize)

    # Dimension is present
    if dimsize is None:
        # Wanted unlimited dim, check that it is
        if not d.isunlimited():
            raise RuntimeError('Dimension {} is already present and it is '
                    'not unlimited as requested.'.format(dimname))
    else:
        # Fixed size dim - compare sizes
        if len(d) != dimsize:
            raise RuntimeError('Dimension {} is already present and its '
                    'size {} differs from requested {}.'.format(dimname,
                        len(d), dimsize))
    return d

def getvar(f, varname, *args, **kwargs):
    """Creates a variable in a netCDF file or returns it if it already exists.
    Does NOT verify its parameters.
    """
    try:
        v = f.variables[varname]
    except KeyError:
        return f.createVariable(varname, *args, **kwargs)
    return v

class SliceExtender:
    __slots__ = ['slice_obj', 'slices']

    def __init__(self, slice_obj, *slices):
        self.slice_obj = slice_obj
        self.slices = slices

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self.slice_obj[key+self.slices]
        else:
            return self.slice_obj[(key,)+self.slices]
