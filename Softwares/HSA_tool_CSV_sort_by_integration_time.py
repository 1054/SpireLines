#!/usr/bin/env python
# 
# Query the Herschel Science Archive (HSA) by RA Dec
# 



# See HSA documentation:
# http://archives.esac.esa.int/hsa/archive-help/cl/
# see "TAP resources"






# Import python packages

import os, sys, time, numpy, astropy
import astropy.io.ascii as asciitable
from astropy.time import Time



# Read user input

input_csv_file = ''
output_csv_file = ''

i = 1
while i < len(sys.argv):
    str_arg = sys.argv[i].lower().replace('--','-').replace('_','-')
    if str_arg == "-out" or str_arg == "-output":
        if i+1 < len(sys.argv):
            i = i+1
            output_csv_file = sys.argv[i]
    else:
        if input_csv_file == '':
            input_csv_file = sys.argv[i]
    i = i+1

if input_csv_file == '':
    print('Usage: ')
    print('    HSA_tool_CSV_sort_by_integration_time.py aaa.csv')
    print('')
    sys.exit()

if output_csv_file == '':
    output_csv_file, output_csv_ext = os.path.splitext(input_csv_file)
    output_csv_ext = 'csv'
    output_csv_file = output_csv_file + '_sorted'
else:
    output_csv_file, output_csv_ext = os.path.splitext(output_csv_file)

# 
# Sort by integration time
# 
if os.path.isfile(input_csv_file):
    input_data = asciitable.read(input_csv_file)
    if 'start_time' in input_data.colnames and \
        'end_time' in input_data.colnames:
        start_time_var = Time(input_data['start_time'], format='isot', scale='utc')
        end_time_var = Time(input_data['end_time'], format='isot', scale='utc')
        integration_time = end_time_var - start_time_var
        #print(start_time_var, end_time_var, integration_time.value * 60.0)
        input_data['obs_duration_minutes'] = integration_time.value * 60.0
        input_data.sort('obs_duration_minutes')
        input_data.reverse()
        asciitable.write(input_data, output_csv_file + '.' + output_csv_ext, Writer=asciitable.Csv, overwrite=True)
        print('Output to "%s"!'%(output_csv_file + '.' + output_csv_ext))
        asciitable.write(input_data[['observation_id','obs_duration_minutes','observer','aor']], output_csv_file + '_observation_id.' + output_csv_ext, Writer=asciitable.Csv, overwrite=True)
        print('Output to "%s"!'%(output_csv_file + '_observation_id.' + output_csv_ext))
        






