"""
Database configuration settings
Create one database engine pool with 64 engines; each has a connection pool of size 1.
"""

import sqlalchemy
from db_config import db_url

# Create a single shared collection of database engines
engine_count: int = 64
db_engines = \
    tuple([sqlalchemy.create_engine(db_url, pool_size=1) for i in range(engine_count)])
