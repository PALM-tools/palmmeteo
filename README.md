# PALM-meteo: processor of meteorological input data for the PALM model system

PALM-meteo is an advanced and modular tool to create PALM's *dynamic driver*
with initial and boundary conditions (IBC) and other time-varying data,
typically using (but not limited to) outputs from mesoscale models.

## Functionality

The PALM-meteo workflow consists of these _stages_:

1. **Model setup**: setting up basic items such as the PALM model domain.
   Currently this requireds providing the already prepared PALM *static
   driver*.

2. **Input loading**: selection of requested variables, area selection,
   transformation and/or unit conversion where required. Certain input
   variables are also temporally disagregated.

3. **Horizontal interpolation** from the input model grid to PALM grid.
   Includes geographic projection conversion where required. Models with
   traditional rectangular grid are regridded using bilinear interpolation, the
   ICON model with the icosahedral grid uses Delaunay triangulation and
   barycentric interpolation.

4. **Vertical interpolation** from input model levels (which may be
   terrain-following, isobaric, eta, hybrid etc.) to PALM model levels
   (altitude-based). Part of this process is terrain matching, as the
   high-resolution PALM terrain may differ, even significantly, from the input
   model terrain. This process includes configurable vertical adaptation, where
   the lowest input layer is shifted to match the PALM terrain while in the
   higher layers, the vertical shifts are progressively smaller until they
   reach the _transition level_, above which they are not shifted at all.

5. **Output generation** creates the final PALM dynamic driver. Final
   adjustments are performed here, notably the mass balancing which is
   performed on all boundaries (respecting terrain on lateral boundaries), so
   that the PALM's internal mass balancing (which is performed only on the top
   boundary as a last resort) is not overused.

Currently PALM-meteo supports these meteorological inputs:

### Meteorological IBC
- WRF
- ICON
- Aladin
- Synthetic inputs (detailed profile specification etc.)

### Radiation inputs (optional)
- WRF
- ICON
- Aladin

### Chemical IBC (optional)
- CAMx
- CAMS

PALM-meteo is higly modular and more input sources will be likely added in the
future. A detailed technical description will be made available in the upcoming
scientific paper.

## Installation

There are three basic ways to install PALM-meteo.

### 1. Simple minimal installation using PIP

PALM-meteo is available in PyPI and it can be installed with the simple
command:

```
pip3 install palmmeteo
```

This will also install the `pmeteo` command for running PALM-meteo.
Depending on your operating system and environment, you may want to create
a virtual environment first, or use the `--user` option, or use
[pipx](https://pipx.pypa.io/) instead of pip.

However, this method will install only the bare minimum to run PALM-meteo,
without documentation and tests, so it is only recommended for experienced
users who want a quick installation.

### 2. Full in-place instalation with a virtual environment

This is the recommended method for most users.

1. Make sure that you have python3, pip, and python-venv installed and available.
   This is operating system dependent.
2. Download and extract PALM-meteo into a directory of your choice. If you have
   already installed the PALM model system, then PALM-meteo is already
   extracted under `packages/dynamic_driver/palm_meteo`, albeit probably not in
   the newest version. Or you can download PALM-meteo from the
   [PALM gitlab server](https://gitlab.palm-model.org/dynamic_driver/palm_meteo)
   or from
   [GitHub](https://github.com/rolicot/palmmeteo/releases).
3. From within the PALM-meteo install directory, run `./setup`. This will
   create a new virtual environment, install PALM-meteo with all dependencies
   and create a symlink to the `pmeteo` command (which enables the virtual
   environment automatically).

### 3. Advanced installation for developers

1. Clone the PALM-meteo git repository.
2. Create and enable your own virtual environment if you prefer.
3. From within the PALM-meteo install directory, run `./setup_novenv`. This will
   install PALM-meteo dependencies and create the `pmeteo` running script.

### Testing

When installed using method 2, the install script performs the test at the end
of a successful instalation. You may also exectute them at any time from the
PALM-meteo install directory using the command

```
tests/integration_tests/all_tests.sh
```

or for indivudal tests:

```
tests/integration_tests/test_XX_NAMEOFTEST.sh
```

You may also use the directory `tests/integration_tests/simple_wrf` as
a tutorial with an example WRF case for PALM-meteo.

## Usage

For each dynamic driver, a *YAML* configuration file needs to be prepared. This
file uses sensible defaults for most options, so it does not need to be very
long, as demonstrated by the `example.yaml` file. However for the beginners it
is best to start by making a copy of the `template.yaml` file, which contains
all possible options with their defaults and documentation, and modifying it
according to your needs.

### Basic model configuration

The main part of configuration is selecting a single or multiple *tasks* by
adding a list item in the `tasks:` configuration section.  Selecting a task
means just telling what PALM-meteo what it has to do, which typically involves
creating IBC and/or other PALM inputs using the selected method, such as using
a specific input model.

These are the currently supported tasks (obviously many of them are mutually
 exclusive):
 
- `wrf`: Create IBC from WRF model outputs.
- `wrf_rad`: Create PALM radiation inputs from WRF model outputs (typically
  *AUXHIST* outputs with potentially different time step from standard
  *WRFOUT*).
- `icon2`:   Create IBC from ICON outputs in the *NetCDF* format.
- `aladin`:  Create IBC from Aladin outputs in the *grib* format.
- `camx`:    Create chemistry IBC from CAMx model outputs.
- `cams`:    Create chemistry IBC from CAMS model outputs.

When the specified task(s) are selected, the task configuration mechanism
enables the required plugins and pulls in the respective task-specific
configuration defauls, which may be overwritten within the configuration file.

### Running the model

With a prepared configuration file such as `myconfig.yaml`, simply run

```
./pmeteo -c myconfig.yaml
```

in the project directory. See also the output of `./pmeteo -h`.

## License and authors

PALM-meteo is distributed under the GNU GPL v3+ license (see the `LICENSE`
file).  It was created by the Institute of Computer Science of the Czech
Academy of Sciences (ICS CAS) with contributions by the Charles University
in Prague (MFF UK), the Czech Hydrometeorological Institute (CHMI) and the Deutsche
Wetterdienst (DWD).
