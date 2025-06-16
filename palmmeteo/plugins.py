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

import importlib
from abc import ABCMeta, abstractmethod

event_hooks = {
}

plugins = []


def eventhandler(event):
    """
    Decorator function to register method as an event handler.

    """
    def wrap(f):
        def wrapped_f(*args, **kwargs):
            f(*args, **kwargs)

        wrapped_f._event = event
        return wrapped_f
    return wrap


class PluginMeta(ABCMeta):
    """
    Meta class for plugin classes

    Inherits from ABCMeta so we can use @abstractmethod decorator in plugin
    mixins.

    Checks for event handler methods in plugins and fills in event_hooks list.
    Allows only one method name to be the handler of a specific event.
    E.g. ImportPluginMixin registers event 'import' with 'import_data' method as
    its handler. Import plugins inherited from ImportPluginMixin then must
    implement 'import_data' method to handle the 'import' event.
    """
    def __new__(cls, name, bases, dct):
        inst = super().__new__(cls, name, bases, dct)
        for n, o in dct.items():
            if callable(o) and hasattr(o, '_event'):
                if o._event in event_hooks:
                    raise ValueError(
                        'Hook already defined for event {}'.format(o._event))

                event_hooks[o._event] = {'class': name, 'method': n}

        return inst


class Plugin(metaclass=PluginMeta):
    """
    Base class for plugins

    Expects to receive configuration object in cfg argument.
    """
    def __init__(self, *args, **kwargs):
        if 'cfg' in kwargs:
            self.cfg = kwargs['cfg']

        if 'rt' in kwargs:
            self.rt = kwargs['rt']

    @eventhandler('check_config')
    def check_config(self, *args, **kwargs):
        """Validate plugin-related configuration.

        Any plugin can optinally implement the check_config method for
        validating configuration. It is not required, so the method is not
        abstract. It has to be stateless, runtime data are not stored to disk.
        """
        pass


class ImportPluginMixin(Plugin):
    """
    Base class mixin for plugins importing data.
    Registers 'import_data' method as a handler for event 'import_data'.

    Abstract methods required to be implemented by derived classes:
        import_data
    """
    @abstractmethod
    @eventhandler('import_data')
    def import_data(self, *args, **kwargs):
        pass


class HInterpPluginMixin(Plugin):
    """
    Base class mixin for plugins importing data.
    Registers 'interpolate_horiz' method as a handler for event 'hinterp'.

    Abstract methods required to be implemented by derived classes:
        interpolate_horiz
    """
    @abstractmethod
    @eventhandler('hinterp')
    def interpolate_horiz(self, *args, **kwargs):
        pass


class VInterpPluginMixin(Plugin):
    """
    Base class mixin for plugins importing data.
    Registers 'interpolate_vert' method as a handler for event 'vinterp'.

    Abstract methods required to be implemented by derived classes:
        interpolate_vert
    """
    @abstractmethod
    @eventhandler('vinterp')
    def interpolate_vert(self, *args, **kwargs):
        pass


class SetupPluginMixin(Plugin):
    """
    Base class mixin for setup plugins.
    Registers 'setup_model' method as a handler for event 'setup_model'.

    Abstract methods required to be implemented by derived classes:
        setup_model
    """
    @abstractmethod
    @eventhandler('setup_model')
    def setup_model(self, *args, **kwargs):
        pass


class WritePluginMixin(Plugin):
    """
    Base class mixin for writer plugins.
    Registers 'import_data' method as a handler for event 'write'.

    Abstract methods required to be implemented by derived classes:
        write_data
    """
    @abstractmethod
    @eventhandler('write')
    def write_data(self, *args, **kwargs):
        pass


def plugin_factory(plugin, *args, **kwargs):
    try:
        mod_name, cls_name = plugin.rsplit('.', 1)
        mod_obj = importlib.import_module(mod_name)
        cls_obj = getattr(mod_obj, cls_name)
    except ValueError:
        cls_obj = globals()[plugin]

    plugin_instance = cls_obj(*args, **kwargs)

    return plugin_instance
