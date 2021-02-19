# Libraries for searching file system
import os
import glob
import shutil

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
def merge_csvs_chunk(table, fnames_in, chunk_size, i):
	"""Merge the CSVs matching a file into one big CSV"""
	n0 = chunk_size*i
	n1 = n0 + chunk_size
	fname_out = f'{table}-megachunk-{i:03d}.csv'
	path_out = os.path.join(dir_csv, table, fname_out)
	with open(path_out, 'wb') as fh_out:
		for fname in fnames_in[n0:n1]:
			with open(fname, 'rb') as fh_in:
				shutil.copyfileobj(fh_in, fh_out)
				fh_out.write(b"\n")

# *****************************************************************************
def merge_csvs(table, chunk_size):
	"""Merge the CSVs matching a file into a batch of big CSV"""
	fnames_in = find_csvs(table)
	# i_max = len(fnames_in) // chunk_size
	i_max = 1
	for i in range(i_max):
		merge_csvs_chunk(table=table, fnames_in=fnames_in, chunk_size=chunk_size, i=i)

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

	# Delete these files after they are imported

