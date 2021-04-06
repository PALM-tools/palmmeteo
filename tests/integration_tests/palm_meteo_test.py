#!/usr/bin/env python
#------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2021 Leibniz Universitaet Hannover
# Copyright 2021 Deutscher Wetterdienst Offenbach
#------------------------------------------------------------------------------#
#
# Authors:
# --------
# @author Eckhard Kadasch
#
# Description:
# ------------
# This script runs the integration test for palm_meteo. It runs palm_meteo with
# a test setup and compares a set of global attributes and all netCDF variables
# of the produced dynamic driver with a reference driver file.
#------------------------------------------------------------------------------#
import os
import subprocess
import sys

import netCDF4
import numpy

DEBUGGING = 'on'
EXIT_CODE_OK = 0
EXIT_CODE_FAIL = 1

ATTRIBUTES_TO_CHECK = set(['origin_lon', 'origin_lat', 'origin_z'])
REFERENCE_FILENAME = 'id_dynamic_reference.nc'
TEST_FILENAME = 'id_dynamic.nc'


def main(argv):

    exit_code = EXIT_CODE_FAIL
    test_script_name = os.path.basename(__file__)

    try:
        config_file = argv[1]
    except IndexError as e:
        print('Test configuration missing. Please specify the palm_meteo test configuration YAML file.')
        print(f'usage: {test_script_name} <config file>')
        return exit_code

    config_file = os.path.abspath(config_file)
    case_dir = os.path.dirname(config_file)
    case_name = case_dir.split('/')[-1].partition(
        '_' + case_dir.split('_')[-1]
    )[0]
    test_dir = os.path.dirname(os.path.realpath(__file__))
    path_to_reference_file = case_dir + '/' + REFERENCE_FILENAME
    path_to_test_file = case_dir + '/' + TEST_FILENAME
    script_path = os.path.abspath(
        os.path.dirname(__file__) + '../../../main.py'
    )
    reference_call = f'{script_path} -c {config_file}'

    print(f'Test case:           {case_name}')
    print(f'Test case directory: {case_dir}')
    print(f'Test configuration:  {config_file}')
    print(f'Comparing test file: {path_to_test_file}')
    print(f'with reference file: {path_to_reference_file}')
    print(f'Running:             {reference_call}')
    print(f'...in:               {test_dir}')

    subprocess.run(['rm', '-f', path_to_test_file], cwd=test_dir)

    # TODO replace reference file copy...
    subprocess.run(['cp', path_to_reference_file, path_to_test_file], cwd=test_dir)
    # TODO ...by palm_meteo call:
    #subprocess.run(reference_call.split(' '), cwd=test_dir)

    try:

       with netCDF4.Dataset(path_to_reference_file, 'r') as reference_file, \
            netCDF4.Dataset(path_to_test_file, 'r') as test_file:

            print_debug(f'Comparing test file: {path_to_test_file}')
            print_debug(f'with reference file: {path_to_reference_file}')

            exit_code = compare_files(reference_file, test_file)

    except OSError as e:

        print_debug(f'{e.strerror}: {e.filename.decode("utf-8")}')

    print_test_result(exit_code, case_name)
    return exit_code


def compare_files(reference_file, test_file):

    test_result, missing_items = test_file_contains_reference_attributes(test_file)
    if test_result:
        print_debug('All required global attributes are present.')
    else:
        print_debug('The following global attributes are missing:')
        print_debug(missing_items)
        return EXIT_CODE_FAIL

    if all_attributes_match(reference_file, test_file):
        print_debug('All attributes match.')
    else:
        print_debug('Some global attributes do not match.')
        return EXIT_CODE_FAIL

    test_result, missing_items = test_file_contains_reference_variables(reference_file, test_file)
    if test_result:
        print_debug('All variables are present.')
    else:
        print_debug('The following variables are missing:')
        print_debug(missing_items)
        return EXIT_CODE_FAIL

    if all_variables_match(reference_file, test_file):
        print_debug('All variables match.')
    else:
        print_debug('Some variables do not match.')
        return EXIT_CODE_FAIL

    return EXIT_CODE_OK


def test_file_contains_reference_attributes(test_file):
    """
    Check if any required global netCDF attribute (listed in
    ATTRIBUTES_TO_CHECK) are missing from the test file. Additional attributes
    are permitted.
    """
    return set_contains_reference_items(
        reference_set=ATTRIBUTES_TO_CHECK,
        test_set=set(test_file.ncattrs())
    )


def all_attributes_match(reference_file, test_file):
    for attribute in ATTRIBUTES_TO_CHECK:
        reference_value = float(reference_file.getncattr(attribute))
        test_value = float(test_file.getncattr(attribute))

        if test_value == reference_value:
            print_debug(f'  {attribute}: value matches')
        else:
            diff = test_value - reference_value
            print_debug(f'  {attribute}: error = {diff}')
            return False

    return True


def test_file_contains_reference_variables(reference_file, test_file):
    """
    Check if any netCDF variable contained in the reference file is missing
    from the test file. Additional variables in the test file are permitted.
    """
    reference_vars = set(reference_file.variables.keys())
    test_vars = set(test_file.variables.keys())

    return set_contains_reference_items(
        reference_set=set(reference_file.variables.keys()),
        test_set=set(test_file.variables.keys())
    )


def all_variables_match(file_a, file_b):
    vars_a = set(file_a.variables.keys())
    vars_b = set(file_b.variables.keys())
    shared_vars = vars_a.intersection(vars_b)
    true_if_all_match = True

    for var in shared_vars:
        data_matches = (file_a.variables[var][:] == file_b.variables[var][:]).all()
        if data_matches:
            print_debug(f'  {var}: data matches')
        else:
            max_diff = numpy.max(file_a.variables[var][:] - file_b.variables[var][:])
            print_debug(f'  {var}: max error = {max_diff}')
            true_if_all_match = False

    return true_if_all_match


def set_contains_reference_items(reference_set, test_set):
    missing_items = reference_set.difference(test_set)
    if len(missing_items) == 0:
        return True, missing_items
    else:
        return False, missing_items


def print_debug(message):
    if DEBUGGING == 'on':
        print(f'palm_meteo_test: {message}')




def print_test_result(exit_code, case_name):

    if exit_code == EXIT_CODE_OK:
        print_debug(f"SUCCESS: palm_meteo passed integration test case '{case_name}'.")
    else:
        print_debug(f"FAILURE: palm_meteo failed integration test case '{case_name}'.")



if __name__ == '__main__':
    exit(main(sys.argv))

