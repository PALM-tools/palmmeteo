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

"""
A collection of simple, technical utilities.

These utilities need to be stateless, i.e. they must not depend on config or
runtime.
"""

import os
import re
from datetime import datetime, timedelta
from dataclasses import dataclass
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
    """Fully divide datetime td by timedelta step."""
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

def assert_dir(filepath):
    """Creates a directory for an output file if it doesn't exist already."""

    dn = os.path.dirname(filepath)
    if not os.path.isdir(dn):
        os.makedirs(dn)

@dataclass
class DTIndexer:
    """
    Calculates integral time index from start and origin. Avoids
    using the unpicklable lambdas.
    """
    origin: datetime
    timestep: timedelta

    def __call__(self, dt):
        return tstep(dt-self.origin, self.timestep)

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

class SliceBoolExtender:
    __slots__ = ['slice_obj', 'slices', 'boolindex']

    def __init__(self, slice_obj, slices, boolindex):
        self.slice_obj = slice_obj
        self.slices = slices
        self.boolindex = boolindex

    def __getitem__(self, key):
        if isinstance(key, tuple):
            v = self.slice_obj[key+self.slices]
        else:
            v = self.slice_obj[(key,)+self.slices]
        return v[...,self.boolindex]

class Workflow:
    """Indexes and maintains the workflow as a sequence of named stages"""

    def __init__(self, default_stages):
        self.default_stages = default_stages
        self.idx = {s:i for i, s in enumerate(default_stages)}
        self.snapshot_from = None

    def stage_idx(self, s):
        try:
            return self.idx[s]
        except KeyError:
            raise ValueError(f'Unknown workflow stage: "{s}". '
                             f'Valid workflow stages are: {self.default_stages}.')

    def assign_all(self):
        self.workflow = self.default_stages

    def assign_fromto(self, stage_from, stage_to):
        wf1 = self.stage_idx(stage_from) if stage_from else 0
        wf2 = self.stage_idx(stage_to)   if stage_to   else -1
        self.workflow = self.default_stages[wf1:wf2+1]
        if wf1 > 0:
            self.snapshot_from = self.default_stages[wf1-1]

    def assign_list(self, stages):
        workflow = [self.stage_idx(s) for s in stages]
        wfnext = [si+1 for si in workflow[:-1]]
        if wfnext != workflow[1:]:
            raise ValueError(f'Workflow is not contiguous! Selected: {stages}. '
                             f'Complete workflow: {self.default_stages}.')
        self.workflow = [self.default_stages[si] for si in workflow]
        if workflow[0] > 0:
            self.snapshot_from = self.default_stages[workflow[0]-1]

    def __iter__(self):
        return iter(self.workflow)
