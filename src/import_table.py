"""
Perform a batch import of CSVs to a DB table.
Example call:
$ python import_table AsteroidVectors

Michael S. Emanuel
2021-02-18
"""

# File system
import os
import sys
import subprocess
import glob

# MSE imports
from config import ks_root
from db_utils import find_fnames_csv

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
    parser.add_argument('--schema', nargs='?', metavar='SCH', type=str, default='KS'
                        help='The name of the schema, e.g KS')

    # Unpack command line arguments
    args = parser.parse_args()

    # Table name
    schema = args.schema
    table = args.table
    # File name for loading
    fname_load = os.path.join(dir_csv, 'mariadb-import', f'{table}.csv')
    Path(fname_load).mkdir(parents=True, exist_ok=True)

    # Get CSV file names
    fnames_csv = find_fnames_csv(table=table, verbose=True)
    fnames_csv = fnames_csv[100:110]

    # Iterate through the CSVs and load them with mariadb-import
    for fname_csv in fnames_csv:
        # Rename this chunk file to AsteroidVectors.csv
        # This is a limitation of mariadb-import; the CSV file name must match the table name exactly
        os.rename(fname_csv, fname_load)

        # Arguments to run mariadb-import from subprocess
        args = [
            'mariadb-import',
            '--defaults-file=mariadb-import-options.cnf',
            '--replace',
            schema, 
            fname_load,
        ]
        # Run mariadb-import
        print('\n', fname_csv)
        subprocess.run(args)

        # # Restore the file name
        os.rename(fname_load, fname_csv)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
