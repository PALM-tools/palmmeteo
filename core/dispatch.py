import netCDF4

import core.plugins as plg
from .logging import die, warn, log, verbose, configure_log
from .config import load_config, cfg
from .runtime import rt, basic_init


def build_exec_queue(event, from_plugins):
    # logika vytvareni fronty muze byt slozitejsi nez jen prosty seznam (strom, mozna paralelizace....)
    queue = []
    for plugin in from_plugins:
        if isinstance(plugin, getattr(plg, plg.event_hooks[event]['class'])):
            queue.append(plugin)

    return queue


def execute_event(event):
    queue = build_exec_queue(event, rt.plugins)

    kwargs = {}
    common_files = []
    try:
        # Prepare common files or other common processing for specific events
        if event == 'import':
            f = netCDF4.Dataset(rt.paths.imported, 'w', format='NETCDF4')
            common_files.append(f)
            kwargs['fout'] = f
        elif event == 'hinterp':
            f = netCDF4.Dataset(rt.paths.hinterp, 'w', format='NETCDF4')
            common_files.append(f)
            kwargs['fout'] = f
        elif event == 'vinterp':
            f = netCDF4.Dataset(rt.paths.vinterp, 'w', format='NETCDF4')
            common_files.append(f)
            kwargs['fout'] = f

        # Execute each plugin in a queue
        for plugin in queue:
            getattr(plugin, plg.event_hooks[event]['method'])(**kwargs)
    finally:
        for f in common_files:
            try:
                f.close()
            except:
                warn('Error closing file {}!', f)


def run(argv):
    load_config(argv)
    basic_init(rt)
    configure_log(cfg)
    rt.plugins = [plg.plugin_factory(p, cfg=cfg, rt=rt)
                      for p in cfg.plugins]

    for event in cfg.workflow:
        execute_event(event)
