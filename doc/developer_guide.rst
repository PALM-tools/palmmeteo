=================================================
PALM Meteorogical input processor developer guide
=================================================

Directory structure
===================

* core: general purpose builtin classes and utility functions
* contrib: model specific plugins
* doc: documentation, example configuration, etc.


Configuration
=============

Configuration is stored in YAML files. User configuration file is
specified by the ``-c/--config`` option to the main program. Some configuration
options may be specified on the command line as main program options.

The configuration system is designed to provide as many sane defaults
as possible in the ``default_config.yaml`` file. Users can override
the defaults in the user configuration file or by using command line options.


Runtime
=======

Shared data, including runtime configuration, are stored in the 
configuration-like ``runtime.rt`` object.


Workflow
========

Program workflow is controlled by the ``workflow`` list of events in the main
configuration section. For the time being, the events are triggered in the
``dispatch.run`` function one by one according to the order of the list.
The ``dispatch.build_exec_queue`` function creates a list of *plugins*
(see below) to be executed for a specific event. At the moment it only takes
the plain list named ``plugins`` from the main configuration section.

The plan for later is to implement a more complex execution queue logic, eg.
tree, parallel execution...


Plugins
=======

Tasks are implemented as plugins. A plugin is a class which implements
event-specific methods. The ``event_hooks`` dictionary is a **1:1 mapping**
of method names to events constructed from registered plugins by searching for
methods annotated with an ``@eventhandler('event_name')`` decorator, e.g.::

    class CleanUpPlugin(Plugin):
        @eventhandler('cleanup')
        def do_clean_up(self, *args, **kwargs):
            ...

maps the ``do_clean_up`` method name to the ``cleanup`` event. At the moment,
only one method name is allowed for one event. To allow more than one class
to implement a specific event handler, mixin classes are created first
for each event, e.g.::

    class ImportPluginMixin(Plugin):
        @abstractclass
        @eventhandler('import')
        def import_data(self, *args, **kwargs):
            ...

declares an ``import_data`` method name as an event handler for the ``import``
event. Mixin (interface) methods should be declared as abstract to prevent
instantiation of the mixin class. Classes meant to actually implement some
event handling should be inherited from one or more mixins, e.g.::

    class WRFPlugin(ImportPluginMixin, CleanUpPlugin):
        def import_data(self, *args, **kwargs):
            ...

        def do_clean_up(self, *args, **kwargs):
            ...


Logging, debugging
==================

The ``logging`` module provides four functions:

* die: print message to stderr and stop execution
* warn: print message to stderr
* log: print debugging message to stderr unless cfg.silent is set to True
* verbose: print debugging to stderr if cfg.verbose is set to True
