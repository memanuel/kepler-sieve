"""
Database configuration settings
Includes one shared sqlalchemy database engine variable that can be shared by multiple consumers.
"""

import sqlalchemy

# Credential to connect to the kepler database instance
hostname = 'Thor.elliptic.loc'
username = 'kepler'
password = 'kepler'

# Folder on local filesystem used to load data files
directory_csv = '/ssd1/tmp/mysql'

# Create database engine - once for the whole module
db_url: str = f'mysql+pymysql://{username}:{password}@{hostname}'
db_engine: sqlalchemy.engine = sqlalchemy.create_engine(db_url)
