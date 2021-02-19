# MSE utility for importing CSV tables to database
from mysqlsh_util import import_table

# Libraries for searching file system
import os
import glob

# The root directory where CSVs are saved
dir_csv = '/ssd1/Harvard/kepler-sieve/data/df2db'

# Import state vectors
table='AsteroidVectors'
import_table(mysql_import=util.import_table, table=table, threads=20)



