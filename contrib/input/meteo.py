from core.plugins import Plugin, ImportPluginMixin
from core.logging import die, warn, log, verbose

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
