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

import sys

__all__ = ['die', 'warn', 'log', 'verbose', 'configure_log']


def die(s, *args, **kwargs):
    """Write message to error output and exit with status 1."""

    if args or kwargs:
        error_output(s.format(*args, **kwargs) + '\n')
    else:
        error_output(s + '\n')
    sys.exit(1)


def warn(s, *args, **kwargs):
    """Write message to error output."""

    if args or kwargs:
        error_output(s.format(*args, **kwargs) + '\n')
    else:
        error_output(s + '\n')


def log_on(s, *args, **kwargs):
    """Write logging or debugging message to standard output (logging is enabled)."""

    if args or kwargs:
        log_output(s.format(*args, **kwargs) + '\n')
    else:
        log_output(s + '\n')


def log_off(s, *args, **kwargs):
    """Do nothing (logging is disabled)."""

    pass


log_output = sys.stdout.write
error_output = sys.stderr.write
log = log_on
verbose = log_off

def configure_log(cfg):
    global log, verbose

    log = log_on if cfg.verbosity >= 1 else log_off
    verbose = log_on if cfg.verbosity >= 2 else log_off
