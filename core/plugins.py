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

        if 'rt_cfg' in kwargs:
            self.rt_cfg = kwargs['rt_cfg']


class ImportPluginMixin(Plugin):
    """
    Base class mixin for plugins importing data.
    Registers 'import_data' method as a handler for event 'import'.

    Abstract methods required to be implemented by derived classes:
        import_data
    """
    @abstractmethod
    @eventhandler('import')
    def import_data(self, *args, **kwargs):
        pass


class SetupPluginMixin(Plugin):
    """
    Base class mixin for setup plugins.
    Registers 'import_data' method as a handler for event 'setup'.

    Abstract methods required to be implemented by derived classes:
        setup_model
    """
    @abstractmethod
    @eventhandler('setup')
    def setup_model(self, *args, **kwargs):
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
