"""
Global configuration options for KS project

Michael S. Emanuel
2021-02-19
"""

# Libraries
import os

# ********************************************************************************************************************* 
# The base installation directory of the project
ks_root: str = '/home/michael/Harvard/kepler-sieve'

# The mariadb-import configuration file
mdbi_opt: str = os.path.join(ks_root, 'cfg', 'mariadb-import-options.cnf')
