#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

if do_test simple_wrf_hinterp -c "$basedir/simple_wrf/simple_wrf.yaml" -w hinterp -W hinterp; then
    compare_file "$basedir/simple_wrf/METEO/hinterp.nc"
fi

do_summary
