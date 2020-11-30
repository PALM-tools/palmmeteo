import core.plugins as plg
import core.logging
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
    for plugin in queue:
        getattr(plugin, plg.event_hooks[event]['method'])()


def run(argv):
    load_config(argv)
    basic_init(cfg)
    core.logging.configure(cfg)
    rt.plugins = [plg.plugin_factory(p, cfg=cfg, rt=rt)
                      for p in cfg.plugins]

    for event in cfg.workflow:
        execute_event(event)
