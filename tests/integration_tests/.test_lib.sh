#!/bin/sh

# Copyright 2018-2025 Institute of Computer Science of the Czech Academy of
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

init_tests() {
    set -e

    # Paths
    basedir="tests/integration_tests"
    pmeteo="./pmeteo"
    ncdiffp="$basedir/ncdiffp"

    # Test results
    ntok=0
    ntbad=0
    nfok=0
    nfbad=0
    retval=0
    unset isgroup

    RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'

    pmeteo_integration_tests_initialized=y
}

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

do_test() {
    testname="$1"
    shift 1

    echo
    echo ==============================
    echo "= Running integration test $testname"
    echo ==============================
    echo
    if "$pmeteo" "$@"; then
        printf "${GREEN}Integration test $testname executed successfully.${NC}\n"
        ntok=$(( $ntok + 1 ))
        return 0
    else
        printf "${RED}Integration test $testname failed!${NC}\n"
        ntbad=$(( $ntbad + 1 ))
        retval=1
        return 1
    fi
}

do_summary() {
    # Skip if we are doing a group
    ! [ $isgroup ] || return 0

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
}

# Init unless already initialized
[ $pmeteo_integration_tests_initialized ] || init_tests
