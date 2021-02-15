# MSE utility for importing CSV tables to database
from mysqlsh_util import find_csvs

# Libraries for searching file system
import os
import glob

# The root directory where CSVs are saved
dir_csv = '/ssd1/Harvard/kepler-sieve/data/df2db'

# Columns for orbital elements
columns_elt = ['TimeID', 'AsteroidID', 'MJD', 'a', 'e', 'inc', 'Omega_node', 'omega_peri', 'f', 'M']

def import_table(table, threads, columns=None):
	"""Import CSVs for the named table into the database"""
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

	# Import the CSVs into the DB
	util.import_table(fnames, options)

# Import state vectors and orbital elements
import_table(table='AsteroidVectors', threads=40)
import_table(table='AsteroidElements', threads=40)

