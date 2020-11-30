import os

class RuntimeObj(object):
    """An object for holding runtime-related values.

    May be nested.
    """
    pass

def basic_init(cfg):
    """Performs initializaiton of basic values from config."""

    # Paths
    rt.paths = RuntimeObj()
    path_exp = dict(scenario=cfg.scenario)
    rt.paths.base = cfg.paths.base.format(**path_exp)
    rt.paths.palm_input = os.path.join(rt.paths.base, cfg.paths.palm_input.format(**path_exp))

    # Domain
    rt.stretching = (cfg.domain.dz_stretch_factor != 1.0)

rt = RuntimeObj()
