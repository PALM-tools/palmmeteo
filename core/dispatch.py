import core.plugins as plg
from .config import load_config, cfg, rt_cfg


def build_exec_queue(event, from_plugins):
    # logika vytvareni fronty muze byt slozitejsi nez jen prosty seznam (strom, mozna paralelizace....)
    queue = []
    for plugin in from_plugins:
        if isinstance(plugin, getattr(plg, plg.event_hooks[event]['class'])):
            queue.append(plugin)

    return queue


def execute_event(event):
    queue = build_exec_queue(event, rt_cfg.plugins)
    for plugin in queue:
        getattr(plugin, plg.event_hooks[event]['method'])()


def run(argv):
    load_config(argv.config)
    rt_cfg.plugins = [plg.plugin_factory(p, cfg=cfg, rt_cfg=rt_cfg)
                      for p in cfg.plugins]

    for event in cfg.workflow:
        execute_event(event)
