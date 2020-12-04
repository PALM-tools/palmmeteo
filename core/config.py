import os
import yaml
import datetime
from collections import defaultdict
from .logging import die, warn, log, verbose

class ConfigError(Exception):
    def __init__(self, desc, section=None, key=None):
        self.desc = desc
        self.section = section
        self.key = key

        # Build message
        s = ['Configuration error: ', desc]
        if section:
            s.extend([', item: ', ':'.join(section.get_path()+[key])])
            try:
                v = section.d[key]
            except KeyError:
                s.append(', missing value')
            else:
                s.extend([', value=', str(v)])
        s.append('.')
        self.msg = ''.join(s)

    def __str__(self):
        return self.msg

class ConfigObj(object):
    def __init__(self, parent=None, name=None):
        self.parent = parent
        self.name = name
        self.d = {}

    def __getattr__(self, name):
        try:
            return self.d[name]
        except KeyError:
            raise AttributeError('Attribute {} not found. Possibly a missing '
                    'configuration setting in section {}.'.format(name,
                        ':'.join(self.get_path())))

    def ingest_dict(self, d, overwrite=True, extend=False, check_exist=False):
        for k, v in d.items():
            if isinstance(v, ConfigObj):
                # we are actually ingesting a subtree - replace by its dict
                v = v.d

            if isinstance(v, dict):
                # For a dictionary (top-level or with only dictionaries above,
                # i.e. a subsection), we recurse
                try:
                    vl = self.d[k]
                except KeyError:
                    # not yet present: create a new empty child node
                    vl = ConfigObj(self, k)
                    self.d[k] = vl
                try:
                    vl.ingest_dict(v, overwrite, extend, check_exist)
                except AttributeError:
                    raise ConfigError('Trying to replace a non-dictionary '
                            'setting with a dictionary', self, k)
            elif extend and isinstance(v, list):
                # We extend lists if requested
                vl = self.d.setdefault(k, [])
                try:
                    vl.extend(v)
                except AttributeError:
                    raise ConfigError('Trying to extend a non-list setting with '
                            'a list', self, k)
            elif v is None and isinstance(self.d.get(k), ConfigObj):
                # This is a special case: we are replacing an existing section
                # with None. That most probably means that the user has
                # presented an empty section (possibly with all values
                # commented out). In that case, we do not want to erase that
                # section. To actually erase a whole section, the user can
                # still present empty dictionary using the following syntax:
                # section_name: {}
                pass
            else:
                # Non-extended lists and all other objects are considered as
                # values and they are copied as-is (including their subtrees if
                # present). Non-null values are overwritten only if
                # overwrite=True.
                if overwrite:
                    if check_exist and k not in self.d:
                        warn('WARNING: ignoring an unknown setting {}={}.',
                                ':'.join(self.get_path()+[k]), v)
                    self.d[k] = v
                else:
                    if self.d.get(k, None) is None:
                        self.d[k] = v

    def get_path(self):
        if self.parent is None:
            return []
        path = self.parent.get_path()
        path.append(self.name)
        return path


duration_units = {
        'd': 'days',
        'h': 'hours',
        'm': 'minutes',
        's': 'seconds',
        }

def parse_duration(section, item):
    def err():
        raise ConfigError('Bad specification of duration. The correct format is '
                '{num} {unit} [{num} {unit} ...], where {unit} is one of d, h, '
                'm, s. Example: "1 m 3.2 s".', section, item)

    try:
        s = section.d[item]
    except KeyError:
        err()

    words = s.split()
    n = len(words)
    if n % 2:
        err()

    d = defaultdict(int)
    for i in range(0, n, 2):
        ns, unit = words[i:i+2]
        try:
            num = int(ns)
        except ValueError:
            try:
                num = float(ns)
            except ValueError:
                err()
        try:
            u = duration_units[unit]
        except KeyError:
            err()
        d[u] += num

    return datetime.timedelta(**d)


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
            cfg.ingest_dict(yaml.load(config_file), check_exist=True)

    # apply settings for selected tasks
    for task in cfg.tasks:
        try:
            task_set = cfg.task_config.d[task]
        except KeyError:
            die('Unknown task: "{}". Available tasks are: {}', task,
                    ', '.join(cfg.task_config.d.keys()))

        if task_set.set:
            cfg.ingest_dict(task_set.set.d, overwrite=False, extend=False)
        if task_set.extend:
            cfg.ingest_dict(task_set.extend.d, overwrite=False, extend=True)

    # load extra settings from commandline
    cfg.ingest_dict(vars(argv))
    if argv.verbosity_arg is not None:
        cfg.verbosity = argv.verbosity_arg

    # Basic verification
    if not cfg.case:
        raise ConfigError('Case name must be specified', cfg, 'case')
    if cfg.simulation.origin_time.tzinfo is None:
        cfg.simulation.origin_time = cfg.simulation.origin_time.replace(
                tzinfo=datetime.timezone.utc)


cfg = ConfigObj()
