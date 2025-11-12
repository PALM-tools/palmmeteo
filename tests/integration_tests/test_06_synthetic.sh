#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

if do_test synthetic -c "$basedir/synthetic/synthetic.yaml"; then
    compare_file "$basedir/synthetic/INPUT/synthetic_dynamic"
fi

do_summary
