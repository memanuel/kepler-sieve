import os
import sys
import subprocess
import glob

import kepler_sieve
from db_utils import find_fnames_csv

from typing import List

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
    # Directory name
    dir_csv = '/ssd1/Harvard/kepler-sieve/data/df2db'
    # Table name
    schema = 'KS'
    table = 'AsteroidVectors'
    # File name for loading
    fname_load = os.path.join(dir_csv, f'{table}.csv')

    # Get CSV file names
    # fnames_csv = ['/ssd1/Harvard/kepler-sieve/data/df2db/AsteroidVectors-chunk-0.csv']
    # fnames_csv = fnames_csv[100:110]
    fnames_csv = find_fnames_csv(table=table, verbose=True)

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

        # Restore the file name
        os.rename(fname_load, fname_csv)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
