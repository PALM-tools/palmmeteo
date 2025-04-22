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

from palmmeteo.plugins import Plugin, ImportPluginMixin
from palmmeteo.logging import die, warn, log, verbose

available_meteo_vars = {
    'tas':  {'desc': 'temperature at surface', 'units': 'K'},
    'ta':   {'desc': '3D temperature', 'units': 'K'},
    'qas':  {'desc': 'specific humidity at the surface', 'units': 'kg/kg'},
    'qa':   {'desc': '3D specific humidity', 'units': '1'},
    'rsds': {'desc': 'surface incident SW radiation for BVOC', 'units': 'W/m2'},
    'pa':   {'desc': '3D pressure', 'units': 'Pa'},
    'zf':   {'desc': 'layer interface heights', 'units': 'm'},
    'ua':   {'desc': 'U-wind', 'units': 'm/s'},
    'va':   {'desc': 'V-wind', 'units': 'm/s'}
}

required_variables = set()


class RequiresMeteoPluginMixin(Plugin):
    """
    Set a list of required meteorological variables in plugin metainformation:

    class Requires:
        meteo_vars = [ ... ]

    Global available_meteo_vars holds names of all variables known to the
    processor.
    """

    def __init__(self, *args, **kwargs):
        if hasattr(self, 'Requires') and hasattr(self.Requires, 'meteo_vars'):
            my_required_vars = set()
            for v in self.Requires.meteo_vars:
                if v not in available_meteo_vars:
                    raise ValueError(
                        'Unknown meteorological variable required by plugin {}.'
                        .format(self))
                else:
                    my_required_vars.add(v)

            required_variables.update(my_required_vars)
        else:
            raise AttributeError(
                'Missing Requires.meteo_vars in plugin {} derived from '
                'RequiresMeteoPluginMixin.'.format(self))

#TODO WRF and CAMx plugins will be ported to provides/requires architecture later

class SomeMeteoPlugin(ImportPluginMixin):
    def import_data(self, *args, **kwargs):
        log('Importing meteo data...')

    class Provides:
        meteo_vars = ['tas', 'pa']


class EmisPlugin(RequiresMeteoPluginMixin, ImportPluginMixin):
    def import_data(self, *args, **kwargs):
        log('Import emission data')

    class Requires:
        meteo_vars = ['tas', 'pa']
