#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

if do_test simple_wrf_import -c "$basedir/simple_wrf/simple_wrf.yaml" -w import_data -W import_data; then
    compare_file "$basedir/simple_wrf/METEO/import.nc"
fi

do_summary
