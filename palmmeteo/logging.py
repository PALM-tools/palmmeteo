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
from time import strftime

__all__ = ['die', 'warn', 'log', 'verbose', 'configure_log']

dtf = '%Y-%m-%d %H:%M:%S '
log_output = sys.stdout.write
error_output = sys.stderr.write

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

class LoggingLevel:
    def __init__(self, is_on, use_dt=False):
        self.is_on = is_on
        self.use_dt = use_dt

    def __call__(self, s, *args, **kwargs):
        """Write logging or debugging message with optional datetime if configured to do so."""

        if not self.is_on:
            return

        if args or kwargs:
            ss = s.format(*args, **kwargs)
        else:
            ss = s

        if self.use_dt:
            log_output(strftime(dtf) + ss + '\n')
        else:
            log_output(ss + '\n')

    def __bool__(self):
        return self.is_on

log = LoggingLevel(True)
verbose = LoggingLevel(False)

def configure_log(verbosity, log_datetime=False):
    log.__init__(verbosity >= 1, log_datetime)
    verbose.__init__(verbosity >= 2, log_datetime)
