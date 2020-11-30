#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""PALM meteo input processor.

Creates PALM dynamic driver from multiple sources.
"""

from argparse import ArgumentParser
from core import run


argp = ArgumentParser(description=__doc__)
argp.add_argument('-c', '--config', help='configuration file')
verbosity = argp.add_mutually_exclusive_group()
verbosity.add_argument('-v', '--verbose', action='store_true', help='increase verbosity')
verbosity.add_argument('-s', '--silent', action='store_true', help='print only errors')

argv = argp.parse_args()

run(argv)
