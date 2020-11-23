from .plugins import SetupPluginMixin, ImportPluginMixin


class SetupPlugin(SetupPluginMixin):
    def setup_model(self, *args, **kwargs):
        print('Model setup...')


class InitStaticPlugin(ImportPluginMixin):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)

    def import_data(self, *args, **kwargs):
        print('Importing static data')

    class Provides:
        pass
