#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

if do_test simple_wrf_write -c "$basedir/simple_wrf/simple_wrf.yaml" -w write -W write; then
    compare_file "$basedir/simple_wrf/INPUT/simple_wrf_dynamic"
fi

do_summary
