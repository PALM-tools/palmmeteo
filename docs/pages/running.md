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
