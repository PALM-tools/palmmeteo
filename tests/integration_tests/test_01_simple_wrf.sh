#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

if do_test simple_wrf -c "$basedir/simple_wrf/simple_wrf.yaml"; then
    compare_file "$basedir/simple_wrf/INPUT/simple_wrf_dynamic"
    compare_file "$basedir/simple_wrf/METEO/import.nc"
    compare_file "$basedir/simple_wrf/METEO/hinterp.nc"
    compare_file "$basedir/simple_wrf/METEO/vinterp.nc"
fi

do_summary
