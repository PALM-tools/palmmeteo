#!/bin/sh

. "$(dirname "$0")/.test_lib.sh"

# Do all test with a single summary
isgroup=y
for testname in "$basedir"/test_*.sh; do
    . "$testname"
done
unset isgroup

do_summary
