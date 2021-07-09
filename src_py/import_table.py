"""
Perform a batch import of CSVs to a DB table.
Example call:
$ python import_table AsteroidVectors
$ python import_table AsteroidVectors --schema KS

Michael S. Emanuel
2021-02-18
"""

# File system
import os
import sys
import glob
from pathlib import Path

# Running programs
import argparse
import subprocess
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from config import ks_root, mdbi_opt
from db_utils import find_fnames_csv, clean_empty_dirs, get_columns

# Types
from typing import List

# Root directory where CSV files are found 
dir_csv: str = os.path.join(ks_root, 'data', 'df2db')

# ********************************************************************************************************************* 
def report_csv_files(fnames_csv, verbose: bool):
    """Report CSV file names"""
    if verbose:
        nf = len(fnames_csv)
        print(f'CSV files: {nf}')
        if nf > 0:
            print(fnames_csv[0])

# ********************************************************************************************************************* 
def main():
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Load CSVs into the named database table.')
    parser.add_argument('table', nargs='?', metavar='TBL', type=str, 
                         help='The name of the table to load, e.g AsteroidVectors')
    parser.add_argument('--schema', nargs='?', metavar='SCH', type=str, default='KS',
                        help='The name of the schema, e.g KS')
    parser.add_argument('--columns', const=True, default=False, action='store_const',
                        help='specify only a subset of columns; for tables with generated columns')
    parser.add_argument('--clean', const=True, default=False, action='store_const',
                        help='just clean up empty directories and quit')

    # Unpack command line arguments
    args = parser.parse_args()
    schema = args.schema
    table = args.table
    use_cols = args.columns
    clean = args.clean

    # File name for loading
    Path(os.path.join(dir_csv, 'mariadb-import')).mkdir(parents=True, exist_ok=True)
    fname_load = os.path.join(dir_csv, 'mariadb-import', f'{table}.csv')

    # Get CSV file names
    fnames_csv = find_fnames_csv(table=table, verbose=True)

    # Do an initial pass of cleaning empty directories
    clean_empty_dirs(table=table)
    # Quit early if clean was requested
    if clean:
        sys.exit()

    # Add a progress bar
    fnames_it = tqdm_auto(fnames_csv, smoothing=0.0)

    # Get the real column names for this DB table, skipping derived columns
    # if use_cols:
    columns = get_columns(schema=schema, table=table)
    col_list = ','.join(columns)

    # Iterate through the CSVs and load them with mariadb-import
    i: int = 0
    for fname_csv in fnames_it:
        # Rename this chunk file to AsteroidVectors.csv
        # This is a limitation of mariadb-import; the CSV file name must match the table name exactly
        os.rename(fname_csv, fname_load)

        # Arguments to run mariadb-import from subprocess
        args = [
            'mariadb-import',
            f'--defaults-file={mdbi_opt}',
            '--replace',
            # '--use-threads=8',
            '--silent',
            schema, 
            fname_load
        ]
        if use_cols:
            args.insert(3, f'--columns={col_list}')

        # Run mariadb-import
        # print('\n', fname_csv)
        subprocess.run(args)

        # Delete the file; once it's loaded we don't need it anymore
        os.remove(fname_load)

        # Loop counter
        i += 1
        if (i % 100) == 0:
            clean_empty_dirs(table=table)

    # If we get here, files were loaded successfully.  Clean out any empty folders
    clean_empty_dirs(table=table)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()