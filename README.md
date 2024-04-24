# PALM-METEO

PALM-METEO is an advanced and modular tool to create PALM's *dynamic driver*
with initial and boundary conditions (IBC) and other time-varying data,
typically using (but not limited to) outputs from mesoscale models.

## Functionality

The PALM-METEO workflow consists of these items:

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
   model terrain. This process includes configurable stretching, where the
   lowest layer matches the terrain while the vertical shifts are progressively
   smaller in the upper layers.

5. **Output generation** creates the final PALM dynamic driver. Final
   adjustments are performed here, notably the mass balancing which is
   performed on all boundaries (respecting terrain on lateral boundaries), so
   that the PALM's internal mass balancing (which is performed only on the top
   boundary as a last resort) is not overused.

Currently PALM-METEO supports these meteorological inputs:

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

PALM-METEO is higly modular and more input sources will be likely added in the
future. A detailed technical description will be made available in the upcoming
scientific paper.

## Installation

PALM-METEO may be used out-of-the box with the project directory as long as the
all required libraries are available. The easiest way to install them is using

```
pip3 install -r requirements.txt
```

however if you prefer slightly different versions of the libraries specified in
`requirements.txt`, potentially from your operating system's distribution, you
may try them as well.

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
means just telling what PALM-METEO what it has to do, which typically involves
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
./main.py -c myconfig.yaml
```

in the project directory. See also the output of `./main.py -h`.

## License and authors

PALM-METEO is distributed under the GNU GPL v3+ license (see the `LICENSE`
file).  It was created by the Institute of Computer Science of the Czech
Academy of Sciences (ICS CAS) with contributions by the Deutsche Wetterdienst
(DWD).
