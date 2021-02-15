# Libraries for searching file system
import os
import glob

# The root directory where CSVs are saved
dir_csv = '/ssd1/Harvard/kepler-sieve/data/df2db'

# *****************************************************************************
def find_csvs(table):
	"""Find CSVs to import for the named table"""
	search_path = os.path.join(dir_csv, table, f'pid_*', f'{table}*.csv')
	fnames = glob.glob(search_path)
	fnames.sort()
	return fnames

# *****************************************************************************
def import_table(mysql_import, table, threads, columns=None):
	"""
	Import CSVs for the named table into the database
	INPUTS;
		mysql_import: 	the table importing function, utils.import_table
		table:			name of the DB table to be imported
		threads:		number of threads to use
		columns:		columns of the destination DB table
	"""
	# Find the CSV files for this table
	fnames = find_csvs(table=table)

	# Options for mysqlsh import table utility
	options = {
		"schema": "KS", 
		"table": table, 
		"replaceDuplicates": True,
		"dialect": "csv-unix", 
		"skipRows": 1, 
		"threads": threads,
		"showProgress": True, 
	}
	# Add the columns option if it was specified
	if columns is not None:
		options['columns'] = columns		

	# Status
	print(f'Found {len(fnames)} CSV files to import to tabld {table}.')

	# Import the CSVs into the DB
	mysql_import(fnames, options)

