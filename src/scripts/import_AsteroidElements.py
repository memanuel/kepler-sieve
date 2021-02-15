# MSE utility for importing CSV tables to database
from mysqlsh_util import import_table

# Libraries for searching file system
import os
import glob

# The root directory where CSVs are saved
dir_csv = '/ssd1/Harvard/kepler-sieve/data/df2db'

# Columns for orbital elements
columns = ['TimeID', 'AsteroidID', 'MJD', 'a', 'e', 'inc', 'Omega_node', 'omega_peri', 'f', 'M']

# Import orbital elements
import_table(mysql_import=util.import_table, table='AsteroidElements', threads=40, columns=columns)

