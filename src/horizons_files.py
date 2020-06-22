"""
JPL Horizons Files
Load data in files downloaded from Horizons
Download settings:
Ephemeris Types: Vector
Candidate Origin: Solar System Barycenter
Time Span: Start JD2440400.5, Stop JD2477600.5, Step 5 days
Table Settings: quantities code=2, delta-T Yes, CSV yes
"""

import numpy as np
import pandas as pd
import re

# ********************************************************************************************************************
def hrzn_target_body(lines):
    """Search horizons download file for name of target body"""
    for line in lines:
        if line[0:18] == 'Target body name: ':
            target_body = line[18:]
            break
    # regex search pattern
    p = re.compile('(^[\w]+) \(([\d]+)\)[\s]+\{source: ([\w]+)\}\n$')
    # extract name, number, and source from the target body string
    srch = p.search(target_body)
    body_name = srch.group(1)
    body_number = srch.group(2)
    integration_source = srch.group(3)
    return (body_name, body_number, integration_source)

# ********************************************************************************************************************
def hrzn_find_start(lines):
    """Search horizons download file for first line of data"""
    for i, line in enumerate(lines):
        if line == '$$SOE\n':
            return i+1

# ********************************************************************************************************************
def hrzn_find_end(lines):
    """Search horizons download file for last line of data"""
    n = len(lines)
    for i, line in enumerate(lines[::-1]):
        if line == '$$EOE\n':
            return n-i-1

# ********************************************************************************************************************
def write_data_line(fh, line, prefix):
    """Write out one line of data from the text file to the CSV"""
    entries_raw = line.split(',')
    entries = [x.lstrip() for x in entries_raw[0:-1]]
    entries[1] = entries[1].replace('A.D. ', '')
    line_out = f'{prefix},' + ','.join(entries) + '\n'
    fh.write(line_out)

# ********************************************************************************************************************
def hrzn_txt2csv(fname_txt):
    """Convert a text file as downloaded into a CSV ready to be imported to Pandas or a database"""

    # Read the file
    with open(fname_txt) as fh:
        lines = fh.readlines()

    # Get body name, body number, and source of the integration
    body_name, body_number, integration_source = hrzn_target_body(lines)
    # Prefix for each row of output
    prefix = f'{body_number},{body_name},{integration_source}'

    # Get the start and end line for data in the file
    i0 = hrzn_find_start(lines)
    i1 = hrzn_find_end(lines)

    # The header row is always 3 before the first data line
    header_row = i0 - 3

    # The column names in this file
    col_names_txt = lines[header_row]
    col_names_file = [x.strip() for x in col_names_txt.split(',')][0:-1]

    # Convert these to MSE column names
    name_map = {
        'JDTDB': 'JD',
        'Calendar Date (TDB)': 'CalendarDateTime',
        'delta-T': 'delta_T',
        'X': 'qx',
        'Y': 'qy',
        'Z': 'qz',
        'VX': 'vx',
        'VY': 'vy',
        'VZ': 'vz',    
    }

    # The header for the output CSV file; includes BodyNumber and BodyName
    col_names = ['BodyNumber', 'BodyName','IntegrationSource'] + [name_map[x] for x in col_names_file]
    header = ','.join(col_names) + '\n'

    # Name of the CSV output file
    fname_csv = fname_txt.replace('.txt', '.csv')

    # Write out the CSV file
    with open(fname_csv, 'w') as fh:
        fh.write(header)
        for line in lines[i0:i1]:
            write_data_line(fh, line, prefix)

    # Return the name of the CSV file
    return fname_csv

# ********************************************************************************************************************
def hrzn_csv2df(fname_csv):
    """Convert a horizons CSV file to a Pandas DataFrame"""

    # Column names and their data types
    dtype = {
        'BodyNumber': np.int32,
        'BodyName': str,
        'IntegrationSource': str,
        'JD': np.float64,
        # 'CalendarDateTime': pd.to_datetime,
        'delta_T': np.float64,
        'qx': np.float64,
        'qy': np.float64,
        'qz': np.float64,
        'vx': np.float64,
        'vy': np.float64,
        'vz': np.float64,
    }

    # Read the CSV
    df = pd.read_csv(fname_csv, sep=',', header=0, dtype=dtype, parse_dates=['CalendarDateTime'])

    return df
