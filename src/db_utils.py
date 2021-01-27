# Core
import pandas as pd
import sqlalchemy

# Utilities
from pathlib import Path
import os
import time
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from utils import hash_id_crc32
from db_config import db_engine

# Types
from typing import Optional, Dict

# ********************************************************************************************************************* 
# Directory for inserting CSV data into database; make if missing
dir_csv: str = '../data/df2db'
Path(dir_csv).mkdir(parents=True, exist_ok=True)

# ********************************************************************************************************************* 
def sp_bind_args(sp_name: str, params: Optional[Dict]):
    """Bind arguments to a SQL stored procedure.  Return a sqlalchemy text object."""
    # Combine the arguments into a string formatted as expected by SQL alchemy
    arg_str = ', '.join(f':{k}' for k in params.keys())
    # Assemble the SQL string
    sql_str = f'CALL {sp_name}({arg_str});'
    # Bind the parameters into a SQL ALchemy text object
    sql_stmt = sqlalchemy.text(sql_str).bindparams(**params)
    return sql_stmt

# ********************************************************************************************************************* 
def sp2df(sp_name: str, params: Optional[Dict]=None):
    """
    Execute a SQL stored procedure and return a DataFrame.
    INPUTS:
        sp_name:  Name of the stored procedure
        params:   Dictionary of parameters for stored procedure.  key = parameter name, value = parameter value
    OUTPUT:
        df:  Pandas DataFrame with resultset        
    """

    # Bind arguments
    sql_stmt = sp_bind_args(sp_name=sp_name, params=params)

    # Execute the bound SQL and return as a DataFrame
    with db_engine.connect() as conn:
        try:
            df = pd.read_sql(sql_stmt, conn)
        # OK to run an SP with no results; just return None instead
        except sqlalchemy.exc.ResourceClosedError:
            df = None
    return df

# ********************************************************************************************************************* 
def sp_run(sp_name: str, params: Optional[Dict]=None):
    """
    Execute a SQL stored procedure.
    INPUTS:
        sp_name:  Name of the stored procedure
        params:   Dictionary of parameters for stored procedure.  key = parameter name, value = parameter value
    OUTPUT:
        None
    """

    # Bind arguments
    sql_stmt = sp_bind_args(sp_name=sp_name, params=params)

    # Execute the bound SQL and return as a DataFrame
    with db_engine.connect() as conn:
        conn.execute(sql_stmt)
    return df

# ********************************************************************************************************************* 
def df2db(df: pd.DataFrame, schema: str, table: str, truncate: bool=False, verbose: bool=False):
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
    row_count: int = df.shape[0]

    # SQL to truncate the table
    sql_truncate: str = f'TRUNCATE {schema_table};'
    
    # File name of CSV
    hash_id: int = hash_id_crc32(df)
    fname_csv = os.path.join(dir_csv, f'{table}_{hash_id}.csv')
    if verbose:
        print(f'Inserting {row_count} records from dataframe into {schema_table}...')
        print(f'CSV file name: {fname_csv}')
    
    # Convert file to CSV in chunks    
    t0 = time.time()   
    df.to_csv(fname_csv, columns=columns, index=False)   
    t1 = time.time()
    # Report elapsed time if requested
    if verbose:
        print(f'Elapsed Time for CSV conversion: {(t1-t0):5.2f} seconds.')

    # List of column names
    col_list = '(' + ','.join(columns) + ')'
    # print(col_list)

    # SQL to Load CSV into database
    sql_load_csv = \
f"""
LOAD DATA LOCAL INFILE 
'{fname_csv}'
REPLACE
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
    if verbose:
        print(f'Elapsed Time for DB insertion: {(t2-t1):5.2f} seconds.')

# ********************************************************************************************************************* 
def df2db_chunked(df: pd.DataFrame, schema: str, table: str, chunk_size: int, 
                  truncate: bool=False, progbar: bool=False):
    """
    Insert the contents of a Pandas DataFrame into a SQL table.
    INPUTS:
        df:         The DataFrame to insert
        schema:     The schema of the destination DB table
        chunk_size: The number of rows in each chunk
        table:  The name of the destination DB table 
        truncate: Flag indicating whether to first truncate the destination table
    OUTPUTS:
        None.  Modifies the DB table on the server.
    """
    # Set up chunks
    row_count: int = df.shape[0]
    chunk_count: int = row_count // chunk_size + 1 
    ii = list(range(chunk_count))
    if progbar:
        ii = tqdm_auto(ii)

    # Iterate over chunks of the DataFrame
    for i in ii:
        # Current chunk of the df
        i0: int = i * chunk_size
        i1: int = i0 + chunk_size
        df_i: pd.DataFrame = df.iloc[i0:i1]
        # Delegate to df2db
        truncate_i: bool = truncate if (i==0) else False
        df2db(df=df_i, schema=schema, table=table, truncate=truncate_i, verbose=False)
