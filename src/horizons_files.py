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
# import os
from pathlib import Path
import glob
import re
from tqdm.auto import tqdm

# Database
import sqlalchemy

# MSE Imports
import db_config
from asteroid_element import load_ast_elt

# Typing
from typing import List, Tuple, Dict, TextIO

# Create database engine - once for the whole module
# db_url: str = f'mysql+pymysql://{db_config.username}:{db_config.password}@{db_config.hostname}/JPL'
db_url: str = f'mysql+pymysql://{db_config.username}:{db_config.password}@{db_config.hostname}'
engine: sqlalchemy.engine = sqlalchemy.create_engine(db_url)

# ********************************************************************************************************************
def hrzn_target_body(lines: List[str], is_ast: bool) -> Tuple[str, int, str,]:
    """Search Horizons download file for name of target body"""

    # Find line in the file header describing the target body
    line: str
    for line in lines:
        if line[0:18] == 'Target body name: ':
            target_body: str = line[18:]
            break

    # Format of the name varies for major bodies and small bodies
    body_name: str
    body_number_str: str
    integration_source: str
    if not is_ast:
        # regex search pattern
        p = re.compile('(^.+) \(([\d]+)\)[\s]+\{source: (.+)\}\n$')
        # extract name, number, and source from the target body string
        srch = p.search(target_body)
        body_name = srch.group(1)
        body_number_str = srch.group(2)
        integration_source = srch.group(3)
    else:
        # regex search pattern
        p = re.compile('(^[0-9]+) (.+) \(.*\)[\s]+\{source: (.+)\}\n$')
        # Extract name, number, and source from the target body string
        srch = p.search(target_body)
        body_number_str = srch.group(1)
        body_name = srch.group(2)
        integration_source = srch.group(3)

    # Convert body number to an integer
    body_number: int = int(body_number_str)

    # Return tuple with the name, number, and integration source
    return (body_name, body_number, integration_source)

# ********************************************************************************************************************
def hrzn_find_start(lines: List[str]) -> int:
    """Search Horizons download file for first line of data"""
    i: int
    line: str
    for i, line in enumerate(lines):
        if line == '$$SOE\n':
            return i+1

# ********************************************************************************************************************
def hrzn_find_end(lines: List[str]) -> int:
    """Search Horizons download file for last line of data"""
    n: int = len(lines)
    i: int
    line: str
    for i, line in enumerate(lines[::-1]):
        if line == '$$EOE\n':
            return n-i-1

# ********************************************************************************************************************
month_tbl: Dict[str, str] = {
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

def write_data_line(fh: TextIO, line: str, prefix: str):
    """Write out one line of data from the text file to the CSV"""
    entries_raw: List[str] = line.split(',')
    entries: List[str] = [x.lstrip() for x in entries_raw[0:-1]]
    date_time_str: str = entries[1].replace('A.D. ', '')
    month_mmm: str = date_time_str[4:9]
    month_nn: str = month_tbl[month_mmm]
    entries[1]: str = date_time_str.replace(month_mmm, month_nn)
    line_out: str = f'{prefix},' + ','.join(entries) + '\n'
    fh.write(line_out)

# ********************************************************************************************************************
# Table with the type of each body given its name
# Typer are 'S' (star), 'PS' (planetary system barycenter), 'PB' (planet body), M' (moon), 'A' (asteroid)
planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
body_type_tbl: Dict[str, str] = {
    'Sun': 'S',
}
for planet in planets:
    body_type_tbl[f'{planet} Barycenter'] = 'PS'
# The name of Earth-Moon Barycenter doesn't follow the usual pattern
body_type_tbl['Earth-Moon Barycenter'] = 'PS'
for planet in planets:
    body_type_tbl[f'{planet}'] = 'PB'

# ********************************************************************************************************************
def hrzn_txt2csv(fname_txt: str):
    """Convert a Horizons text file as downloaded into a CSV ready to be imported to Pandas or a database"""

    # Determine if this is an asteroid file from the file name
    # fname_base: str = os.path.basename(fname_txt).replace('.txt', '')
    path_txt = Path(fname_txt)
    fname_base: str = path_txt.stem
    is_ast: bool = fname_base[0:4] == 'ast_'

    # Read the file
    with open(fname_txt) as fh:
        lines: List[str] = fh.readlines()

    # Get body name, body number, and source of the integration
    body_name: str
    body_number: int
    integration_source: str
    body_name, body_number, integration_source = hrzn_target_body(lines, is_ast)

    # The type of this body: one of 'S' (star), 'P' (planet), 'M' (moon), 'A' (asteroid)    
    if not is_ast:
        # Look up this body name on the table
        # The planet barycenters and bodies are enumerated; anything left that isn't an asteroid will be a moon.
        body_type: str = body_type_tbl.get(body_name, 'M')           
        # Mercury and Venus are special cases; there is only one planet in these systems
        # and JPL describes e.g. Mercury Barcenter with ID 1 instead of 199
        if body_name == 'Mercury Barycenter':
            body_number = 1
        if body_name == 'Venus Barycenter':
            body_number = 2
    # We already determined this body is an asteroid from its file name
    else:
        body_type = 'A'

    # Prefix for each row of output
    prefix = f'{body_type},{body_number},{body_name},{integration_source}'

    # Get the start and end line for data in the file
    i0: int = hrzn_find_start(lines)
    i1: int = hrzn_find_end(lines)

    # The header row is always 3 before the first data line
    header_row: int = i0 - 3

    # The column names in this file
    col_names_txt: str = lines[header_row]
    col_names_file: List[str] = [x.strip() for x in col_names_txt.split(',')][0:-1]

    # Convert these to MSE column names
    name_map: Dict[str, str] = {
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
    col_names: List[str] = ['BodyTypeCD','BodyNumber', 'BodyName','IntegrationSource'] + \
                            [name_map[x] for x in col_names_file]
    header: str = ','.join(col_names) + '\n'

    # Name of the CSV output file
    # fname_csv: str = fname_txt.replace('.txt', '.csv')
    path_csv = Path(db_config.directory_csv, 'jpl', fname_base).with_suffix('.csv')
    fname_csv: str = path_csv.as_posix()

    # Write out the CSV file
    fh: TextIO
    with open(fname_csv, 'w') as fh:
        fh.write(header)
        for line in lines[i0:i1]:
            write_data_line(fh, line, prefix)

    # Return the name of the CSV file
    return fname_csv

# ********************************************************************************************************************
def hrzn_csv2df(fname_csv: str):
    """Convert a Horizons CSV file to a Pandas DataFrame"""

    # Column names and their data types
    dtype: Dict[str, dtype] = {
        'BodyTypeCD': str,
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
    df: pd.DataFrame = pd.read_csv(fname_csv, sep=',', header=0, dtype=dtype, parse_dates=['CalendarDateTime'])

    return df

# ********************************************************************************************************************
def hrzn_csv2db(fname_csv: str, conn: sqlalchemy.engine.Connection):
    """Load a CSV file into DB table JPL.HorizonsImport"""
    # Build the sql string to load the CSV data file
    sql = \
    f"""
    load data infile '{fname_csv}'
    into table JPL.HorizonsImport 
    fields terminated by ','
    lines terminated by '\n'
    ignore 1 lines
    (BodyTypeCD, BodyNumber, BodyName, IntegrationSource, JD, CalendarDateTime, delta_T, qx, qy, qz, vx, vy, vz)
    set HorizonsImportID = NULL;
    """
    
    # Execute the load data infile command
    # This requires that the kepler user has privileges on both the JPL database AND the global FILE privilege
    # See sql folder in project for user configuration.
    conn.execute(sql)

# ********************************************************************************************************************
def hrzn_df2db(df: pd.DataFrame, conn: sqlalchemy.engine.Connection):
    """Load a Pandas DataFrame into DB table JPL.HorizonsImport"""
    df.to_sql(name='HorizonsImport', con=conn, schema='JPL', if_exists='append', index=False)

# ********************************************************************************************************************
def hrzn_txt2db(fname_txt: str, conn: sqlalchemy.engine.Connection):
    """Load a Horizons TXT file into DB table JPL.HorizonsImport via a CSV file"""   
    # Convert text file to a CSV file
    fname_csv: str = hrzn_txt2csv(fname_txt)

    # Load the DataFrame into the DB Table
    hrzn_csv2db(fname_csv, conn)

# ********************************************************************************************************************
def hrzn_load_db():
    """Load all the Horizons TXT files into DB table JPL.HorizonsImport"""
    # Get a list of all the TXT files with Horizons data to be loaded
    runs = ['planets/daily', 'moons/daily', 'moons/weekly', 'asteroids/weekly']
    files = []
    for run in runs:
        files += sorted(glob.glob(f'../data/jpl/horizons/{run}/*.txt'))

    # Do all the database operations using just one DB connection
    with engine.connect() as conn:
        # Truncate the JPL.HorizonsImport table
        sql = "truncate table JPL.HorizonsImport;"
        conn.execute(sql)

        # Iterate through all the text files in order; load each one
        for fname_txt in tqdm(files):
            hrzn_txt2db(fname_txt, conn)

# ********************************************************************************************************************
def ast_elt_load_db():
    """Load JPL files with orbital elements of numbered and unnumbered asteroids into DB table JPL.AsteroidElement"""

    # Load the orbital elements as a DataFrame
    ast_elt = load_ast_elt()

    # Change column names to match database conventions.
    # These are different because MariaDB has case insensitive column names, 
    # causing a clash between Omega and omega. ugh!
    mapper = {
        'Omega' : 'Omega_node',
        'omega': 'omega_peri',
    }
    ast_elt.rename(columns=mapper, inplace=True)

    with engine.connect() as conn:
        # Truncate the JPL.AsteroidElement database table
        sql = "truncate table JPL.AsteroidElement;"
        # Insert the data from the DataFrame to the database
        conn.execute(sql)
        ast_elt.to_sql(name='AsteroidElement', con=conn, schema='JPL', if_exists='append', index=False)
