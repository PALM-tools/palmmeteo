#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PALM-meteo: processor of meteorological input data for the PALM model system.

Creates PALM dynamic driver from various sources.
"""

try:
    from ._version import __version__, __version_tuple__
except:
    __version__ = 'undetermined'

from .config import cfg
from .runtime import rt
from .dispatch import run, main
