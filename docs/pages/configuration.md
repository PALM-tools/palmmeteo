# PALM-meteo configuration

Each case for PALM-meteo is configured using one or more *YAML* configuration
files. If more than one files is used, the latter files extend and overwrite
the prior files.

## Configuration options

The list of all options which are considered
user-configurable is supplied in the file `template.yaml`. This file
contains the options with their default values and documentation as comments:

\include template.yaml

This file may be used as a template for a new PALM-meteo configuration by
uncommenting the values that the user wants to change.

Any valid options not listed in the file `template.yaml` are intended for
developers only.
