# Core
import pandas as pd
import sqlalchemy

# Utilities
import multiprocessing
from pathlib import Path
import os
import time
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from utils import hash_id_crc32
from db_config import db_engine, db_url

# Types
from typing import Optional, List, Dict

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
def sql_run(sql_str: str, params: Optional[Dict]=None) -> None:
    """
    Execute a SQL string.
    INPUTS:
        sql_str:  String with generic SQL 
        params:   Dictionary of parameters for sql statement.  key = parameter name, value = parameter value
    OUTPUT:
        None
    """

    # Bind the parameters into a SQL ALchemy text object
    sql_stmt = sqlalchemy.text(sql_str).bindparams(**params)
    # Execute the bound SQL and return as a DataFrame
    with db_engine.connect() as conn:
        conn.execute(sql_stmt)    

# ********************************************************************************************************************* 
def truncate_table(schema: str, table: str) -> None:
    """Truncate the named table"""

    # Delimited table name
    schema_table: str = f'{schema}.{table}'    
    # SQL to truncate the table
    sql_truncate: str = f'TRUNCATE TABLE {schema_table};'
    # Execute the truncation
    with db_engine.connect() as conn:
        conn.execute(sql_truncate)    


# ********************************************************************************************************************* 
def get_columns(schema: str, table: str) -> List[str]:
    """Get list of column names for DB table; skip derived columns"""
    # Delimited table name
    schema_table: str = f'{schema}.{table}'

    # Get DB metadata
    with db_engine.connect() as conn:
        md = sqlalchemy.MetaData(bind=conn, schema=schema)
        md.reflect()
        tbl_md = md.tables[schema_table]
        # Get list of all available DB columns; exclude computed columns!
        columns = [c.name for c in tbl_md.columns if not c.computed]
    return columns

# ********************************************************************************************************************* 
def df2db(df: pd.DataFrame, schema: str, table: str, columns: List[str]=None,
          truncate: bool=False, verbose: bool=False, conn=None):
    """
    Insert the contents of a Pandas DataFrame into a SQL table.
    INPUTS:
        df:       The DataFrame to insert
        schema:   The schema of the destination DB table
        table:    The name of the destination DB table
        columns:  List of columns to insert
        truncate: Flag indicating whether to first truncate the destination table
    OUTPUTS:
        None.  Modifies the DB table on the server.
    """

    # Delimited table name
    schema_table: str = f'{schema}.{table}'

    # Get columns from DB metadata if they were not provided by caller
    if columns is None:
        columns = get_columns(schema=schema, table=table)

    # Filter DataFrame to the desired columns
    df = df[columns]
    row_count: int = df.shape[0]

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
    if truncate:
        truncate_table(schema=schema, table=table)
    if conn is None:
        with db_engine.connect() as conn:
            conn.execute(sql_load_csv)
    else:
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
    # Get DB columns
    columns = get_columns(schema=schema, table=table)

    # Truncate the table if requested
    if truncate:
        truncate_table(schema=schema, table=table)

    # Subfunction to process one chunk
    def process_chunk(i: int, conn):
        """Process chunk i"""
        # print(f'i={i}. conn={conn}.')
        # Current chunk of the df
        i0: int = i * chunk_size
        i1: int = i0 + chunk_size
        df_i: pd.DataFrame = df.iloc[i0:i1]
        # Delegate to df2db
        df2db(df=df_i, schema=schema, table=table, columns=columns, 
            truncate=False, verbose=False, conn=conn)

    # Subfunction for a core to process a batch of chunks
    def process_batch(ii: List[int], conn, progbar: bool):
        """Process batch of chunks on one CPU core"""
        # if len(ii) > 0:
        #    print(f'Started core {ii[0]}. conn={conn}.')
        if progbar:
            ii = tqdm_auto(ii)
        for i in ii:
            process_chunk(i=i, conn=conn)

    # Set up chunks
    row_count: int = df.shape[0]
    chunk_count: int = row_count // chunk_size + 1
    cpu_max: int = 2
    cpu_count: int = min(multiprocessing.cpu_count() // 2, cpu_max)
    
    # Iterate over chunks of the DataFrame on one CPU core
    ii = list(range(chunk_count))
    if progbar:
       ii = tqdm_auto(ii)
    with db_engine.connect() as conn:
        for i in ii:
            process_chunk(i=i, conn=conn)
    
    # # Create list of engines and DB connections
    # engines= [sqlalchemy.create_engine(db_url) for c in range(cpu_count)]
    # conns = [engine.connect() for engine in engines]

    # # Iterate over cpu cores
    # for c in range(cpu_count):
    #     # The connection used by this core
    #     conn = conns[c]
    #     # The indices i handled by cpu c
    #     ii = list(range(c, chunk_count, cpu_count))
    #     progbar_c = progbar and (c==0)
    #     # Dispatch to cpu c
    #     kwargs={
    #         'ii' : ii, 
    #         'conn' : conn, 
    #         'progbar': progbar_c}
    #     args = (ii, conn, progbar_c)
    #     p = multiprocessing.Process(target=process_batch, args=args)
    #     p.start()

    # # Close all the connections
    # for conn in conns:
    #     conn.close()
