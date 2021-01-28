# Core
import pandas as pd
import dask
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
conn_type = sqlalchemy.engine.base.Connection

# ********************************************************************************************************************* 
# Directory for inserting CSV data into database; make if missing
dir_csv: str = '../data/df2db'
Path(dir_csv).mkdir(parents=True, exist_ok=True)

# Create a single shared collection of database engines
engine_count: int = 16
db_engines = \
        tuple([sqlalchemy.create_engine(db_url, pool_size=2) for i in range(engine_count)])
        # tuple([sqlalchemy.create_engine(db_url, pool_size=1, poolclass=NullPool) for i in range(engine_count)])

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
def df2csv(df: pd.DataFrame, fname_csv: str, columns: List[str], chunksize: int = 0) -> List[str]:
    """
    Convert a DataFrame to one or a collection of CSVs with a specified chunksize of rows.csv
    INPUTS:
        df:        A Pandas DataFrame
        fname_csv: Filename of output, e.g. TableName.csv for a DB export.
                   Chunk number will be automatically appended when used.
        chunksize: The number of rows in each chunk.  Use 0 for one big chunk (pandas to csv)
    OUTPUTS:
        fnames: List of output filenames
    """
    # If chunksize was passed, use Dask
    if chunksize > 0:
        # Convert the Pandas into a Dask DataFrame
        ddf = dask.dataframe.from_pandas(df, chunksize=chunksize, name='Integration')
        # Export it to CSV in chunks
        fname_chunk = fname_csv.replace('.csv', '-chunk-*.csv')
        fnames = ddf.to_csv(filename=fname_chunk, columns=columns, index=False)
    # If no chunksize was specified, dump the whole frame into one CSV file
    else:
        df.to_csv(fname_csv, columns=columns, index=False)
        fnames = [fname_csv]
    return fnames

# ********************************************************************************************************************* 
def csv2db_stage(schema: str, table: str, columns: List[str], fname_csv: str, i: int, conn: conn_type):
    """
    Load one CSV file into a staging table for the named DB table.
    INPUTS:
        schema:    Schema of the DB table
        table:     Name of the DB table
        columns:   List of columns of the DB table
        fname_csv: Name of the CSV file
        i:         Chunk number
        conn:      DB connection object
    OUTPUTS:
        None.  Modifies the database by creating a staging table
    """

    # Destination table with schema
    schema_table = f'{schema}.{table}'
    # Staging table for this chunk
    staging_table = f'{schema}.{table}_chunk_{i:03d}'
    # List of column names
    col_list = '(' + ','.join(columns) + ')'

    # SQL to create staging table
    sql_clone_table = f"CREATE OR REPLACE TABLE {staging_table} LIKE {schema_table};"

    # SQL to Load CSV into database into staging table
    sql_load_csv = \
        f"""
        LOAD DATA LOCAL INFILE 
        '{fname_csv}'
        REPLACE
        INTO TABLE {staging_table}
        FIELDS TERMINATED BY ','
        LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        {col_list}
        """

    # Clone the table and load the CSV file
    conn.execute(sql_clone_table)
    conn.execute(sql_load_csv)

# ********************************************************************************************************************* 
def csvs2db_stage(schema: str, table: str, columns: List[str], fnames_csv: List[str], 
                  ii: List[int], c: int, progbar: bool):
    """
    Load a batch of CSV files into a staging table for the named DB table.
    INPUTS:
        schema:     Schema of the DB table
        table:      Name of the DB table
        columns:    List of columns of the DB table
        fnames_csv: List of the CSV file names
        ii:         List of chunk numbers for this batch
        c:          CPU (process) number; used to choose the DB engine from the shared pool
        progbar:    Whether to show a tqdm progress bar
    OUTPUTS:
        None.  Modifies the database by creating a staging table
    """
    if progbar:
        ii = tqdm_auto(ii)
    # Use the DB engine corresponding to this CPU number
    db_engine_c = db_engines[c]
    # print(f'CPU {c}, DB engine {db_engine_c}')
    # All SQL statements on this engine use the same connection
    with db_engine_c.connect() as conn:
        for i in ii:
            csv2db_stage(schema=schema, table=table, columns=columns, fname_csv=fnames_csv[i], i=i, conn=conn)
    # print(f'Completed csv2db_stage batch on {ii}.')

# ********************************************************************************************************************* 
def csv2db(schema: str, table: str, columns: List[str], i: int, conn):
    """
    Load one CSV file into a staging table for the named DB table.
    INPUTS:
        schema:    Schema of the DB table
        table:     Name of the DB table
        i:         Chunk number
        conn:      DB connection object
    OUTPUTS:
        None.  Modifies the database by inserting from the staging table into the main table
    """

    # Destination table with schema
    schema_table = f'{schema}.{table}'
    # Staging table for this chunk
    staging_table = f'{schema}.{table}_chunk_{i:03d}'
    # List of column names
    col_list = '(' + ','.join(columns) + ')'
    # Columns to select
    select_fields = ',\n'.join(columns)

    # SQL to insert into the main table from the staging table
    sql_insert_stage = \
        f"""
        REPLACE INTO {schema_table}
        {col_list}
        SELECT
        {select_fields}
        FROM {staging_table};
        """

    # SQL to drop the staging table, which is no longer needed
    sql_drop_stage = f"DROP TABLE {staging_table};"

    # Insert from the staging table, then drop it
    conn.execute(sql_insert_stage)
    conn.execute(sql_drop_stage)

# ********************************************************************************************************************* 
def csvs2db(schema: str, table: str, fnames_csv: List[str], columns: List[str], progbar: bool):
    """
    Load a batch of CSVs into the named DB table.
    First delegates to csv2db to populate the staging tables in parallel.
    Then rolls up from all of the staging tables into the big table at the end.
    """
    # List of indices
    chunk_count: int = len(fnames_csv)
    ii = list(range(chunk_count))
    if progbar:
        ii = tqdm_auto(ii)

    # Multiprocessing
    cpu_max: int = 16
    cpu_count: int = min(multiprocessing.cpu_count() // 2, chunk_count, cpu_max, engine_count)
    
    # Prepare inputs for CSV staging function
    stage_inputs = []
    for c in range(cpu_count):
        # The indices i handled by cpu c
        ii_c = tuple(range(c, chunk_count, cpu_count))
        progbar_c = progbar and (c==0)
        # Arguments to csv2db_stage
        args = (schema, table, columns, fnames_csv, ii_c, c, progbar_c)
        # List of argument tuples for the staging function
        stage_inputs.append(args)
    # Convert from list to tuple
    stage_inputs = tuple(stage_inputs)
    
    # Run in CSV staging in parallel
    pool = multiprocessing.Pool(processes=cpu_count)
    pool.starmap(csvs2db_stage, stage_inputs)

    # Roll up from staging table to schema on one CPU core
    with db_engine.connect() as conn:
        for i in ii:        
            csv2db(schema=schema, table=table, columns=columns, i=i, conn=conn)

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
    fname_csv = os.path.join(dir_csv, f'{table}.csv')
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
