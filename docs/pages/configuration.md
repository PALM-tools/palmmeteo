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

### Specifying paths

The option `paths.base` specifies the base path for all other paths. It may be
either an absolute path or a path relative to the current directory (the
directory from which PALM-meteo was started). All other paths are either
absolute or relative to `paths.base`.

Any path option may contain **replacement strings** in curly braces. These are
used to simplify the specification of typical paths. For example, the string
`{domain}` is (by default) replaced either with `_N02`, `_N03` etc. or with an
empty string for domain number 1, so that it represents a typical file suffix
for PALM input and output files (such as the dynamic or static driver).

Even the replacement strings may be customized using the configuration section
`path_strings:`; each item specifies a Python expression which gets evaluated
as a replacement string and it may reference other configuration options.

The default configuration of paths expects a typical PALM setup: the `JOBS`
directory placed next to the PALM-meteo instalation directory, within which
there are input and output files for individual jobs. Apart from PALM's
standard per-job subdirectories such as `INPUT`, PALM-meteo adds the
subdirectory `METEO` with its intermediate files. In addition to that, some
standard plugins expect other directories with their inputs (such as the `WRF`
directory with WRFOUT files, unless configured otherwise).

### Specifying time durations

The configuration options which expect duration (not absolute time), such as
`simulation.length` and `simulation.timestep`, are specified as a number and
unit, which is one of: `d` (days), `h` (hours), `m` (minutes) or `s` (seconds).

The number and unit are separated by space(s). You may also use decimal numbers
or combine multiple units, so the string `1 d 3.5 h` translates to 27.5
hours.
