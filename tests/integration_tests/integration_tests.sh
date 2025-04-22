#!/bin/sh

# Copyright 2018-202% Institute of Computer Science of the Czech Academy of
# Sciences, Prague, Czech Republic. Authors: Pavel Krc
#
# This file is part of PALM-METEO.
#
# PALM-METEO is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM-METEO is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM-METEO. If not, see <https://www.gnu.org/licenses/>.

set -e

basedir="tests/integration_tests"
pmeteo="./pmeteo"
ncdiffp="$basedir/ncdiffp"

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

compare_file() {
    generated="$1"
    reference="$1.ref"
    echo
    echo ==============================
    echo "= Comparing generated $generated to $reference"
    echo ==============================
    echo
    if "$ncdiffp" "$reference" "$generated" -c -p -s -a ; then
        printf "${GREEN}File $generated matched successully.${NC}\n"
        nfok=$(( $nfok + 1 ))
    else
        printf "${RED}File $generated does not match!${NC}\n"
        nfbad=$(( $nfbad + 1 ))
        retval=1
    fi
}

ntok=0
ntbad=0
nfok=0
nfbad=0
retval=0
for testcase in "$@"; do
    casedir="$basedir/$testcase"
    cfgfile="$casedir/$testcase.yml"
    echo
    echo ==============================
    echo "= Running integration test $testcase ($cfgfile)"
    echo ==============================
    echo
    if "$pmeteo" -c "$cfgfile"; then
        printf "${GREEN}Integration test $testcase executed successfully.${NC}\n"
        ntok=$(( $ntok + 1 ))

        compare_file "$casedir/INPUT/${testcase}_dynamic"
        compare_file "$casedir/METEO/import.nc"
        compare_file "$casedir/METEO/hinterp.nc"
        compare_file "$casedir/METEO/vinterp.nc"
    else
        printf "${RED}Integration test $testcase failed!${NC}\n"
        ntbad=$(( $ntbad + 1 ))
        retval=1
    fi
done

echo
echo ==============================
echo = Summary
echo ==============================
echo
[ $ntok  -le 0 ] || printf "Tests finished:   $GREEN$ntok$NC\n"
[ $ntbad -le 0 ] || printf "Tests failed:     $RED$ntbad$NC\n"
[ $nfok  -le 0 ] || printf "Files matched:    $GREEN$nfok$NC\n"
[ $nfbad -le 0 ] || printf "Files mismatched: $RED$nfbad$NC\n"

exit $retval
