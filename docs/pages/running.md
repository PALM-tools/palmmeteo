# Running PALM-meteo

When PALM-meteo is installed using Method 2 or Method 3, the main running
script `pmeteo` is created in the installation directory. This script is
executed as

    ./pmeteo {OPTIONS}

When using Method 1 or Method 2, this command is also added to the system's
`PATH` variable, which means that PALM-meteo can be also run as

    pmeteo {OPTIONS}

as long as the virtual environment (if used) is enabled.

Finally, as a backup option if neither of the scripts is installed, the
`palmmeteo` Python module can be executed using

    python3 -m palmmeteo {OPTIONS}

Here is the complete list of command line options:

\verbinclude commandline.txt

## Running a subset of stages

To generate the PALM dynamic driver, all stages of the workflow need to be
executed sequentially, which is the default. PALM-meteo also fully supports
stopping the execution after any stage and restarting it from the next stage
later (even on a different machine, if you transfer the intermediate and
snapshot files, which are placed in the `METEO` directory by default).

If the execution of any stage crashes for some reason (e.g. power loss,
insufficient disk space, bug in the code for which you have a fix), you may
also restart it safely. It is also fully supported to restart from any stage
which was already performed, unless the snapshot from the previous stage was
deleted (configuration option `intermediate_files:delete_after_success`, off by
default).

Advanced users may also attempt other usage patterns, such as restarting from
a later stage after changing the configuration. However, the success of this
approach depends on the support by the involved plugins. If the configuration
has been changed, it is advisable to include the first stage `check_config` and
potentially also the second stage `setup_model` at the beginning. The remaining
stages should always form a continuous string because the main model data can
only be transferred between adjacent stages.
