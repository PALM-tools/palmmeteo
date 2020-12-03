import os
import yaml
import datetime
from .logging import die, warn, log, verbose


class ConfigObj(object):
    def ingest_dict(self, d, overwrite=True, extend=False):
        for k, v in d.items():
            if isinstance(v, dict):
                # For a dictionary (top-level or with only dictionaries above,
                # i.e. a subsection), we recurse
                try:
                    vl = vars(self)[k]
                except KeyError:
                    # not yet present: create a new empty child node
                    vl = ConfigObj()
                    vars(self)[k] = vl
                try:
                    vl.ingest_dict(v, overwrite, extend)
                except KeyError:
                    die('Trying to replace a non-dictionary setting "{}" with '
                            'a dictionary', k)
            elif extend and isinstance(v, list):
                # We extend lists if requested
                vl = vars(self).setdefault(k, [])
                try:
                    vl.extend(v)
                except AttributeError:
                    die('Trying to extend a non-list setting "{}" with '
                            'a list', k)
            else:
                # Non-extended lists and all other objects are considered as
                # values and they are copied as-is (including their subtrees if
                # present). Non-null values are overwritten only if
                # overwrite=True
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
    for segment in ['core', 'workflow', 'plugins']:
        cfg_segment_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                'config_init_{}.yaml'.format(segment))
        with open(cfg_segment_path, 'r') as f:
            cfg.ingest_dict(yaml.load(f))

    # load settings from selected configfile (if available)
    if argv.config:
        with open(argv.config, 'r') as config_file:
            cfg.ingest_dict(yaml.load(config_file))

    # apply settings for selected tasks
    for task in cfg.tasks:
        try:
            task_set = vars(cfg.task_config)[task]
        except KeyError:
            die('Unknown task: "{}". Available tasks are: {}', task,
                    ', '.join(vars(cfg.task_config).keys()))

        if task_set.set:
            cfg.ingest_dict(vars(task_set.set), overwrite=False, extend=False)
        if task_set.extend:
            cfg.ingest_dict(vars(task_set.extend), overwrite=False, extend=True)

    # load extra settings from commandline
    cfg.ingest_dict(vars(argv))
    if argv.verbosity_arg is not None:
        cfg.verbosity = argv.verbosity_arg

    # Basic verification
    if not cfg.case:
        die('"case" (case name) must be configured')
    if cfg.simulation.origin_time.tzinfo is None:
        cfg.simulation.origin_time = cfg.simulation.origin_time.replace(
                tzinfo=datetime.timezone.utc)


cfg = ConfigObj()
