#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

if do_test simple_wrf_vinterp -c "$basedir/simple_wrf/simple_wrf.yaml" -w vinterp -W vinterp; then
    compare_file "$basedir/simple_wrf/METEO/vinterp.nc"
fi

do_summary
