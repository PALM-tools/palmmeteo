#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PALM-meteo: processor of meteorological input data for the PALM model system.

Creates PALM dynamic driver from various sources.
"""

try:
    # Use generated version file first (wheel or editable install)
    from ._version import __version__, __version_tuple__
except ModuleNotFoundError:
    # If generated file is not present, this could be a git repo. Try to
    # determine version dynamically using setuptools_scm explicitly.
    try:
        from setuptools_scm import get_version
        __version__ = get_version(root='..', relative_to=__file__)
    except:
        # Cannot determine, but allow running pmeteo anyway
        __version__ = 'undetermined'

signature = f'PALM-meteo version {__version__}'

from .config import cfg
from .runtime import rt
from .dispatch import run, main
