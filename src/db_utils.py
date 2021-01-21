# Core
import pandas as pd
import sqlalchemy

# Utilities
from pathlib import Path
import os
import time

# MSE imports
from db_config import db_engine

# ********************************************************************************************************************* 
# Directory for inserting CSV data into database; make if missing
dir_csv: str = '../data/df2db'
Path(dir_csv).mkdir(parents=True, exist_ok=True)

# ********************************************************************************************************************* 
def df2db(df: pd.DataFrame, schema: str, table: str, truncate: bool=False, report_time: bool=False):
    """
    Insert the contents of a Pandas DataFrame into a SQL table.
    INPUTS:
        df:     The DataFrame to insert
        schema: The schema of the destination DB table
        table:  The name of the destination DB table 
        truncate: Flag indicating whether to first truncate the destination table
    OUTPUTS:
        None.  Modifies the DB table on the server.
    """

    # Delimited table name
    schema_table: str = f'{schema}.{table}'

    # Get DB metadata
    with db_engine.connect() as conn:
        md = sqlalchemy.MetaData(bind=conn, schema=schema)
        md.reflect()
        tbl_md = md.tables[schema_table]
        columns = [c.name for c in tbl_md.columns]

    # Filter DataFrame to the desired columns
    df = df[columns]

    # SQL to truncate the table
    sql_truncate: str = f'TRUNCATE {schema_table};'

    
    # File name of CSV
    fname_csv = os.path.join(dir_csv, f'{table}.csv')
    if report_time:
        print(f'CSV file name: {fname_csv}')
    
    # Convert file to CSV
    t0 = time.time()
    df.to_csv(fname_csv, columns=columns, index=False)   
    t1 = time.time()
    # Report elapsed time if requested
    print(f'Elapsed Time for CSV conversion: {(t1-t0):5.2f} seconds.')

    # List of column names
    col_list = '(' + ','.join(columns) + ')'
    # print(col_list)

    # SQL to Load CSV into database
    sql_load_csv = \
f"""
LOAD DATA LOCAL INFILE 
'{fname_csv}'
IGNORE
INTO TABLE {schema_table}
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n'
IGNORE 1 LINES
{col_list}
"""
    
    # Review the SQL statement
    # print(sql_load_csv)
        
    # Truncate the table if requested, then insert the contents of the Pandas DF
    with db_engine.connect() as conn:
        if truncate:
            conn.execute(sql_truncate)
        conn.execute(sql_load_csv)

    # Report elapsed time if requested
    t2 = time.time()
    if report_time:
        print(f'Elapsed Time for DB insertion: {(t2-t1):5.2f} seconds.')
