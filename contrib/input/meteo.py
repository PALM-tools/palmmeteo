import os
import glob
from datetime import datetime, timedelta
import netCDF4

from core.plugins import Plugin, ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin
from core.logging import die, warn, log, verbose
from core.config import cfg
from core.runtime import rt

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


class WRFPlugin(ImportPluginMixin, HInterpPluginMixin, VInterpPluginMixin):
    def import_data(self, *args, **kwargs):
        log('Importing WRF data...')
        ######################################
        # get time extent of the PALM simulation
        #####################################
        # get complete list of wrf files
        wrf_file_list = glob.glob(os.path.join(cfg.paths.wrf_output,
            cfg.paths.wrf_file_mask))
        # get simulation origin and final time as datetime
        rt.end_time = rt.start_time + timedelta(hours=cfg.simulation.length_hours)
        rt.end_time_rad = rt.end_time
        verbose('PALM simulation extent {} - {} ({} hours).', rt.start_time,
                rt.end_time, cfg.simulation.length_hours)
        if rt.nested_domain:
            log('Nested domain - process only initialization. '
                    'Set end_time = start_time')
            rt.end_time = rt.start_time

        # get wrf times and sort wrf files by time
        print('Analyse WRF files dates:')
        file_times = []
        for wrf_file in wrf_file_list:
            # get real time from wrf file
            with netCDF4.Dataset(wrf_file, 'r') as nc_wrf:
                ta = nc_wrf.variables['Times'][:]
                t = ta.tobytes().decode("utf-8")
                td = datetime.strptime(t, '%Y-%m-%d_%H:%M:%S')
                verbose('{}: {}', os.path.basename(wrf_file), td)
                if rt.start_time <= td <= rt.end_time:
                    file_times.append((td, wrf_file))
        if not file_times:
            die('No suitable WRF files found!')

        file_times.sort()
        times, wrf_files = zip(*file_times) #unzip
        rt.times = list(times)
        rt.wrf_files = list(wrf_files)
        verbose('PALM output times: {}', ', '.join(map(str, times)))

        if times[0] != rt.start_time:
            die('WRF files do not contain PALM origin_time timestep - cannot process!')
        if times[-1] != rt.end_time:
            die('WRF files do not contain PALM end_time timestep - cannot process!')

        if not rt.nested_domain and len(times) != simulation_hours+1:
            die('Number of WRF files does not aggre with number of simulation hours')

    def interpolate_horiz(self, *args, **kwargs):
        pass

    def interpolate_vert(self, *args, **kwargs):
        pass

    class Provides:
        meteo_vars = ['tas', 'pa']


class EmisPlugin(RequiresMeteoPluginMixin, ImportPluginMixin):
    def import_data(self, *args, **kwargs):
        print('Import emission data')

    class Requires:
        meteo_vars = ['tas', 'pa']
