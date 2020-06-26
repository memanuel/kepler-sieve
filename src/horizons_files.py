"""
JPL Horizons Files
Load data in files downloaded from Horizons
Download settings:
Ephemeris Types: Vector
Candidate Origin: Solar System Barycenter
Time Span: Start JD2440400.5, Stop JD2477600.5, Step 5 days
Table Settings: quantities code=2, delta-T Yes, CSV yes
"""

# Core
import numpy as np
import pandas as pd

# Utility
import os
import glob
import re
from tqdm.auto import tqdm

# Database
from sqlalchemy import create_engine

# MSE
import db_config

# Create database engine - once for the whole module
db_url = f'mysql+pymysql://{db_config.username}:{db_config.password}@{db_config.hostname}/JPL'
engine = create_engine(db_url)

# ********************************************************************************************************************
def hrzn_target_body(lines, is_ast: bool):
    """Search Horizons download file for name of target body"""

    # find line in the file header describing the target body
    for line in lines:
        if line[0:18] == 'Target body name: ':
            target_body = line[18:]
            break

    # Format of the name varies for major bodies and small bodies
    if not is_ast:
        # regex search pattern
        p = re.compile('(^.+) \(([\d]+)\)[\s]+\{source: (.+)\}\n$')
        # extract name, number, and source from the target body string
        srch = p.search(target_body)
        body_name = srch.group(1)
        body_number = srch.group(2)
        integration_source = srch.group(3)
    else:
        # regex search pattern
        p = re.compile('(^[0-9]+) (.+) \(.*\)[\s]+\{source: (.+)\}\n$')
        # extract name, number, and source from the target body string
        srch = p.search(target_body)
        body_number = srch.group(1)
        body_name = srch.group(2)
        integration_source = srch.group(3)

    # Return tuple with the name, number, and integration source
    return (body_name, body_number, integration_source)

# ********************************************************************************************************************
def hrzn_find_start(lines):
    """Search Horizons download file for first line of data"""
    for i, line in enumerate(lines):
        if line == '$$SOE\n':
            return i+1

# ********************************************************************************************************************
def hrzn_find_end(lines):
    """Search Horizons download file for last line of data"""
    n = len(lines)
    for i, line in enumerate(lines[::-1]):
        if line == '$$EOE\n':
            return n-i-1

# ********************************************************************************************************************
month_tbl = {
    '-Jan-': '-01-',
    '-Feb-': '-02-',
    '-Mar-': '-03-',
    '-Apr-': '-04-',
    '-May-': '-05-',
    '-Jun-': '-06-',
    '-Jul-': '-07-',
    '-Aug-': '-08-',
    '-Sep-': '-09-',
    '-Oct-': '-10-',
    '-Nov-': '-11-',
    '-Dec-': '-12-',
}

def write_data_line(fh, line, prefix):
    """Write out one line of data from the text file to the CSV"""
    entries_raw = line.split(',')
    entries = [x.lstrip() for x in entries_raw[0:-1]]
    date_time_str = entries[1].replace('A.D. ', '')
    month_mmm = date_time_str[4:9]
    month_nn = month_tbl[month_mmm]
    entries[1] = date_time_str.replace(month_mmm, month_nn)
    line_out = f'{prefix},' + ','.join(entries) + '\n'
    fh.write(line_out)

# ********************************************************************************************************************
def hrzn_txt2csv(fname_txt):
    """Convert a Horizons text file as downloaded into a CSV ready to be imported to Pandas or a database"""

    # Determine if this is an asteroid file from the file name
    fname_base = os.path.basename(fname_txt)
    is_ast = fname_base[0:4] == 'ast_'

    # Read the file
    with open(fname_txt) as fh:
        lines = fh.readlines()

    # Get body name, body number, and source of the integration
    body_name, body_number, integration_source = hrzn_target_body(lines, is_ast)
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
    """Convert a Horizons CSV file to a Pandas DataFrame"""

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

# ********************************************************************************************************************
def hrzn_df2db(df):
    """Load a Pandas DataFrame into DB table JPL.HorizonsImport"""
    with engine.connect() as con:
        df.to_sql(name='HorizonsImport', con=con, if_exists='append', index=False)

# ********************************************************************************************************************
def hrzn_txt2db(fname_txt):
    """Load a Horizons TXT file into DB table JPL.HorizonsImport via a DataFrame"""
    
    # Convert text file to a CSV file
    fname_csv = hrzn_txt2csv(fname_txt)

    # Load CSV file to a DataFrame
    df = hrzn_csv2df(fname_csv)

    # Load the DataFrame into the DB Table
    hrzn_df2db(df)

# ********************************************************************************************************************
def hrzn_load():
    """Load all the Horizons TXT files into DB table JPL.HorizonsImport"""
    # Truncate the JPL.HorizonsImport table
    sql = "truncate table JPL.HorizonsImport"
    with engine.connect() as conn:
        conn.execute(sql)

    # Get a list of all the TXT files with Horizons data to be loaded
    runs = ['planets', 'moons', 'asteroids']
    files = []
    for run in runs:
        files += sorted(glob.glob(f'../data/jpl/horizons/{run}/*.txt'))

    # Iterate through all the text files in order; load each one
    for fname_txt in tqdm(files):
        hrzn_txt2db(fname_txt)
