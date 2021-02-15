# Libraries for searching file system
import os
import glob

# The root directory where CSVs are saved
dir_csv = '/ssd1/Harvard/kepler-sieve/data/df2db'

def find_csvs(table):
	"""Import CSVs for the named table into the database"""
	search_path = os.path.join(dir_csv, table, f'pid_*', f'{table}*.csv')
	fnames = glob.glob(search_path)
	fnames.sort()
	return fnames

