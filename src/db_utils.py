# Core
import pandas as pd
import dask.dataframe
import sqlalchemy

# Utilities
import os
import multiprocessing
from pathlib import Path
import glob
import time
import traceback
from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from utils import print_time
from db_config import db_engine, db_url

# Types
from typing import Optional, List, Dict
conn_type = sqlalchemy.engine.base.Connection

# ********************************************************************************************************************* 
# Directory for inserting CSV data into database; make if missing
dir_csv: str = '../data/df2db'
Path(dir_csv).mkdir(parents=True, exist_ok=True)

# Create a single shared collection of database engines
engine_count: int = 32
db_engines = \
    tuple([sqlalchemy.create_engine(db_url, pool_size=2) for i in range(engine_count)])

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
        except:
            print(f'sp2df(sp_name) failed!')
            traceback.print_exc()
            df = None
    return df

# ********************************************************************************************************************* 
def sp_run(sp_name: str, params: Optional[Dict]=None) -> None:
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
        try:
            conn.execute(sql_stmt)
        except:
            print(f'sp_run({sp_name}) failed!')
            print(sql_stmt)
            traceback.print_exc()

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
        try:
            conn.execute(sql_stmt)
        except:
            print(f'sql_run failed!')
            print(sql_stmt)
            traceback.print_exc()

# ********************************************************************************************************************* 
def truncate_table(schema: str, table: str) -> None:
    """Truncate the named table"""

    # Delimited table name
    schema_table: str = f'{schema}.{table}'    
    # SQL to truncate the table
    sql_truncate: str = f'TRUNCATE TABLE {schema_table};'
    # Execute the truncation
    with db_engine.connect() as conn:
        try:
            conn.execute(sql_truncate)
        except:
            print(f'truncate_table failed!')
            print(sql_truncate)
            traceback.print_exc()

# ********************************************************************************************************************* 
def get_columns(schema: str, table: str) -> List[str]:
    """Get list of column names for DB table; skip derived columns"""
    # Delimited table name
    schema_table: str = f'{schema}.{table}'

    # Get DB metadata
    with db_engine.connect() as conn:
        try:
            md = sqlalchemy.MetaData(bind=conn, schema=schema)
            md.reflect()
            tbl_md = md.tables[schema_table]
            # Get list of all available DB columns; exclude computed columns!
            columns = [c.name for c in tbl_md.columns if not c.computed]
        except:
            print(f'get_columns failed!')
            traceback.print_exc()
            columns = []
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
    # Set up dask parallel computing
    compute_kwargs = {
        'scheduler':'threads',
    }
    # Put all the interesting steps in a try / catch so calling program won't crash
    try:
        # If chunksize was passed, use Dask
        if chunksize > 0:
            # Convert the Pandas into a Dask DataFrame
            ddf = dask.dataframe.from_pandas(df, chunksize=chunksize)
            # Export it to CSV in chunks
            fname_chunk = fname_csv.replace('.csv', '-chunk-*.csv')
            fnames = ddf.to_csv(filename=fname_chunk, columns=columns, index=False, compute_kwargs=compute_kwargs)
        # If no chunksize was specified, dump the whole frame into one CSV file
        else:
            df.to_csv(fname_csv, columns=columns, index=False)
            fnames = [fname_csv]
    # Catch exception and return an empty list of CSV file names
    except:
        print(f'df2csv() to fname_csv failed!')
        traceback.print_exc()
        fnames = []
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
    # Display progress bar if requested
    if progbar:
        ii = tqdm_auto(ii)
    # Use the DB engine corresponding to this CPU number
    db_engine_c = db_engines[c]
    # All SQL statements on this engine use the same connection
    with db_engine_c.connect() as conn:
        for i in ii:
            csv2db_stage(schema=schema, table=table, columns=columns, fname_csv=fnames_csv[i], i=i, conn=conn)

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
def csvs2db(schema: str, table: str, 
            fnames_csv: Optional[List[str]]=None, columns: Optional[List[str]]=None, progbar: bool=True):
    """
    Load a batch of CSVs into the named DB table.
    First delegates to csv2db to populate the staging tables in parallel.
    Then rolls up from all of the staging tables into the big table at the end.
    """
    # Get list of CSV files if they were note provided
    if fnames_csv is None:
        search_path = os.path.join(dir_csv, table, f'{table}-chunk*.csv')
        fnames_csv = glob.glob(search_path)
        fnames_csv.sort()    

    # Get columns from DB metadata if they were not provided by caller
    if columns is None:
        columns = get_columns(schema=schema, table=table)

    # List of indices
    chunk_count: int = len(fnames_csv)
    ii = list(range(chunk_count))

    # Multiprocessing
    cpu_max: int = 32
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
            try:
                csv2db(schema=schema, table=table, columns=columns, i=i, conn=conn)
            except:
                print('csvs2db failed! Table {schema}.{table}, chunk i={i}.')
                traceback.print_ext()

# ********************************************************************************************************************* 
def df2db(df: pd.DataFrame, schema: str, table: str, columns: List[str]=None, 
          truncate: bool=False, chunksize: int=0, verbose: bool=False, progbar: bool=True):
    """
    Insert the contents of a Pandas DataFrame into a SQL table.
    INPUTS:
        df:         The DataFrame to insert
        schema:     The schema of the destination DB table
        table:      The name of the destination DB table
        columns:    List of columns to insert; read from DB metadata if omitted
        truncate:   Flag indicating whether to first truncate the destination table
        chunksize:  Number of rows in each chunks
        verbose:    Verbosity of output (true / false)
    OUTPUTS:
        None.  Modifies the DB table on the server.
    """
    
    # Get columns from DB metadata if they were not provided by caller
    if columns is None:
        columns = get_columns(schema=schema, table=table)

    # Filter DataFrame to the desired columns
    df = df[columns]
    row_count: int = df.shape[0]

    # File name of CSV
    fname_csv = os.path.join(dir_csv, f'{table}', f'{table}.csv')

    # Build CSV in parallel with dask
    if verbose:
        print(f'Extracting {row_count} records from DataFrame into CSV files in chunks of {chunksize} rows...')        
        print(f'CSV file name: {fname_csv}')
    t0 = time.time()
    fnames_csv = df2csv(df=df, fname_csv=fname_csv, columns=columns, chunksize=chunksize)
    t1 = time.time()

    # Report elapsed time if requested
    if verbose:
        print_time(time=(t1-t0), msg='Elapsed Time for CSV Conversion')

    # Truncate the table if requested
    if truncate:
        truncate_table(schema=schema, table=table)

    # Insert from the CSVs into the DB table using multiprocessing
    try:   
        csvs2db(schema=schema, table=table, fnames_csv=fnames_csv, columns=columns, progbar=progbar)
    except:
        print('df2db failed!')
        print(f'Table {schema}.{table}.')
        print('Columns:\n', columns)
        
    # Report elapsed time if requested
    t2 = time.time()
    if verbose:
        print_time(time=(t2-t1), msg='Elapsed Time for DB insertion')
