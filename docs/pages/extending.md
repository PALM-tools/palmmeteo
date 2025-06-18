# PALM-meteo developer guide

The core of the PALM-meteo program is implemented within the Python package
#palmmeteo. All plugins are placed in separate Python packages with names
starting with with `palmmeteo_`. When PALM-meteo starts, it locates all
available packages of such names and treats them as plugins. The included
Python package #palmmeteo_stdplugins contains the built-in plugins for the
standard PALM-meteo functionality.

## Creating new plugins

In order to create a new PALM-meteo plugin, you need to create a Python package
whose name starts with `palmmeteo_`. The package needs to be available when
PALM-meteo is started. For that, you can either put in a subdirectory of
the current directory from which PALM-meteo is started, or install it in one
of your system's Python library paths, or add its parent directory into
`sys.path`, e.g. using the environment variable `PYTHONPATH`.

### Plugin code

Suppose you want to enable PALM-meteo to process inputs from a new
meteorological model called Foo. You can place the code within Python package
`palmmeteo_foo`. Inside that, you will create a single Python module `foo.py`
which contains the class FooPlugin. This class will implement data import from
the Foo model and use PALM-meteo's library functions for the next stages
(horizontal and vertical interpolation).

If you are unsure about the implementation, you may also check the class
#palmmeteo_stdplugins.wrf.WRFPlugin
as your reference.

### Basic configuration

Each plugin package must contain the file `config_init.yaml` placed in its
directory. This configuration file contains all plugin-relevant configuration
options with their default values. The user configuration files can then
replace the default values with their own.

From the user perspective, a plugin's functionality is typically enabled by
selecting a `task`. Each task represents a configuration subtree with options
which are set or list items which are appended when the user enables the task.

Without using a dedicated task, the user would need to enable the plugin by
manually adding its fully qualified class name to the list of enabled plugins,
which resides in the configuration section `plugins:`, but that would not be
user friendly. Moreover, you may also implement multiple different tasks which
enable the same plugin but specify different configuration options that are set
with each task.

Continuing with our Foo example, you want to create a task named `foo` which
enables `FooPlugin` by adding it to the list of enabled plugins. To do that,
include this in your `config_init.yaml` file:

~~~~~~~~~{.yaml}
# Configuration items for selected tasks
task_config:
    foo:
        # Configuration items to be set (if unset by user)
        set:

        # Configuration items to be extended (added) to existing lists (or new
        # lists are created).
        extend:
            plugins:
                - palmmeteo_foo.foo.FooPlugin
~~~~~~~~~

## Implementing workflow items

In order to implement actual functionality, it is best to add one or more of the
supplied abstract base classes as a base class for your plugin:

 - #palmmeteo.plugins.SetupPluginMixin
 - #palmmeteo.plugins.ImportPluginMixin
 - #palmmeteo.plugins.HInterpPluginMixin
 - #palmmeteo.plugins.VInterpPluginMixin
 - #palmmeteo.plugins.WritePluginMixin

### Accessing configuration

The final configuration based on default and user-configured values is
accessible in the `cfg` object, which is best imported as

    from palmmeteo.config import cfg

This is an instance of the class #palmmeteo.config.ConfigObj and its
configuration items may be accessed both with the attribute access (dot
notation, e.g. `cfg.mysection.myoption`) and with the item access (bracket
notation, e.g. `cfg.mysection["myoption"]`). It also supports iteration of
key-value pairs and testing of item presence (e.g. `if "myoption" in
cfg.mysection`).

After the configuration is loaded, the `cfg` object is set as read-only.

### Logging

The #palmmeteo.logging module provides these logging functions:

 - [verbose](#palmmeteo.logging.verbose):
   Prints a verbose message into `sys.stdout` only when `--verbose` is enabled.
 - [log    ](#palmmeteo.logging.log    ):
   Prints a standard logging message into `sys.stdout` unless `--silent` is used.
 - [warn   ](#palmmeteo.logging.warn   ):
   Prints a non-fatal warning message into `sys.stderr`.
 - [die    ](#palmmeteo.logging.die    ):
   Prints an error message into `sys.stderr` and terminates the progam immediately.

All of these functions accept the message as their first argument and optional
positional and/or keyword arguments which are used as formatting arguments.
A newline is always added at the end.

### Data persistence

PALM-meteo supports restarting individual workflow stages by storing the
intermediate meteorological data in NetCDF files and other auxiliary data in
pickle files. Any auxiliary data which need to persist from one stage to the
following stages must be stored in the runtime object `rt` (instance of
#palmmeteo.runtime.RuntimeObj), which is imported as

    from palmmeteo.runtime import rt

The object supports attribute access (e.g. `rt.myvariable`).

Although the plugin objects themselves do persist between the stages of
a single PALM-meteo execution, they are not stored for restarting of
stages.

The `rt` object is shared among all plugins and it is often used to exchange
data among plugins. For other cases, plugin developers are encouraged to use
names which are not likely to conflict with other plugins.

Because the configured paths undergo processing after the configuration is
loaded, the processed paths are accessible under `rt.paths`.

### Using the PALM-meteo library

Most common tasks are already implemented and they are provided within the
module [palmmeteo.library](#palmmeteo.library). See the code of standard
plugins as an example.
