import os
import yaml
from .logging import die, warn, log, verbose


class ConfigObj(object):
    def ingest_dict(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                # For a dictionary, we recurse
                try:
                    vl = getattr(self, k)
                except AttributeError:
                    # not yet present: create a new empty child node
                    vl = ConfigObj()
                    setattr(self, k, vl)
                try:
                    vl.ingest_dict(v)
                except AttributeError:
                    die('Trying to replace a non-dictionary setting "{}" with '
                            'a dictionary', k)
            else:
                # Lists and other objects are copied as-is
                setattr(self, k, v)


def load_config(argv):
    """Loads all configuration.

    Configuration is loaded in this order:
    1) default configuration
    2) configfile
    3) command-line options
    Each step may overwrite values from previous steps.
    """
    global cfg

    # load default config
    default_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
            'default_config.yaml')
    with open(default_cfg_path, 'r') as config_file:
        cfg.ingest_dict(yaml.load(config_file))

    # load settings from selected configfile (if available)
    if argv.config:
        with open(argv.config, 'r') as config_file:
            cfg.ingest_dict(yaml.load(config_file))

    # load extra settings from commandline
    cfg.ingest_dict(vars(argv))

    # Basic verification
    if not cfg.scenario:
        die('"scenario" (scenario name) must be configured')


cfg = ConfigObj()
