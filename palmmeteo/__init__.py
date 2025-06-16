#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""PALM-meteo: processor of meteorological input data for the PALM model system

Creates PALM dynamic driver from various sources.
"""
from ._version import __version__, __version_tuple__
from .config import cfg
from .runtime import rt
from .dispatch import run, main
