import yaml


class ConfigObj:
    def ingest_yaml(self, y):
        for k, v in y.items():
            setattr(cfg, k, v)


def load_config(file_name):
    global cfg
    with open(file_name, 'r') as config_file:
        cfg.ingest_yaml(yaml.load(config_file))


cfg = ConfigObj()
rt_cfg = ConfigObj()
