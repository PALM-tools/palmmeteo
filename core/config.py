import os
import yaml
from .logging import die, warn, log, verbose


class ConfigObj(object):
    def ingest_dict(self, d, overwrite=True):
        for k, v in d.items():
            if isinstance(v, dict):
                # For a dictionary, we recurse
                try:
                    vl = vars(self)[k]
                except KeyError:
                    # not yet present: create a new empty child node
                    vl = ConfigObj()
                    vars(self)[k] = vl
                try:
                    vl.ingest_dict(v, overwrite)
                except KeyError:
                    die('Trying to replace a non-dictionary setting "{}" with '
                            'a dictionary', k)
            else:
                # Lists and other objects are copied as-is. Non-null values
                # are overwritten only if overwrite = True
                if overwrite or vars(self).get(k, None) is None:
                    vars(self)[k] = v


def load_config(argv):
    """Loads all configuration.

    Configuration is loaded in this order:
    1) initial configuration values
    2) configfile
    3) command-line options
    Each step may overwrite values from previous steps.
    """
    global cfg

    # load initial configuration segments
    for segment in ['core', 'workflow']:
        cfg_segment_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                'config_init_{}.yaml'.format(segment))
        with open(cfg_segment_path, 'r') as f:
            cfg.ingest_dict(yaml.load(f))

    # load settings from selected configfile (if available)
    if argv.config:
        with open(argv.config, 'r') as config_file:
            cfg.ingest_dict(yaml.load(config_file))

    # apply settings for selected operation mode (preserving user-set values)
    try:
        mode_settings = vars(cfg.operation_modes)[cfg.mode]
    except KeyError:
        die('Unknown operation mode: "{}". Available modes are: {}', cfg.mode,
                ', '.join(vars(cfg.operation_modes).keys()))
    cfg.ingest_dict(vars(mode_settings), overwrite=False)

    # load extra settings from commandline
    cfg.ingest_dict(vars(argv))

    # Basic verification
    if not cfg.case:
        die('"case" (case name) must be configured')


cfg = ConfigObj()
