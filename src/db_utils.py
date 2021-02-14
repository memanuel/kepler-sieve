# Data
import pandas as pd
import dask.dataframe
import sqlalchemy

# Algorithms
import multiprocessing
import itertools

# File system
from pathlib import Path
import os
import glob

# UI
from tqdm.auto import tqdm as tqdm_auto
import traceback
import time

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

# Get the process_id to avoid collisions on the staging tables
pid: int = os.getpid()

# Treat pandas chained assignement as error
pd.set_option('mode.chained_assignment', 'raise')

# ********************************************************************************************************************* 
def make_db_engines(single_thread: bool) -> None:
    """Create a shared set of DB engine objects (global variable).  Used to support multithreading."""
    global db_engines
    if single_thread:
        db_engines = (db_engine,)
    else:
        from db_engine_pool import db_engines

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
def sp2df(sp_name: str, params: Dict=dict()):
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
        # Create directory including the pid
        Path(fname_csv).parent.mkdir(parents=True, exist_ok=True)
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
        # traceback.print_exc()
        # fnames = []
        raise
    return fnames

# ********************************************************************************************************************* 
def dest_table_name(schema: str, table: str):
    """Generate the fully delimited name of the destination TB table."""
    # Destination table with schema
    dest_table = f'{schema}.{table}'
    return dest_table

# ********************************************************************************************************************* 
def staging_table_name(table: str, i: int):
    """Consistently generate schema_table and staging_table names. i is the chunk number."""
    # Staging table for this chunk
    staging_table = f'temp.{table}_pid_{pid}_chunk_{i:03d}'
    return staging_table

# ********************************************************************************************************************* 
def csv2db(schema: str, table: str, columns: List[str], fname_csv: str, conn: conn_type):
    """
    Load one CSV file directly into the named DB table.
    INPUTS:
        schema:    Schema of the DB table
        table:     Name of the DB table
        columns:   List of columns of the DB table
        fname_csv: Name of the CSV file
        conn:      DB connection object
    OUTPUTS:
        None. Modifies the database table.
    """

    # Destination and staging table names
    dest_table = dest_table_name(schema=schema, table=table)
    # List of column names
    col_list = '(' + ','.join(columns) + ')'

    # SQL to Load CSV into database into staging table
    sql_load_csv = \
        f"""
        LOAD DATA LOCAL INFILE 
        '{fname_csv}'
        REPLACE INTO TABLE {dest_table}
        FIELDS TERMINATED BY ','
        LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        {col_list}
        """

    # Clone the table and load the CSV file
    conn.execute(sql_load_csv)
    # Delete the CSV file after it has been successfully loaded
    os.remove(fname_csv)

# ********************************************************************************************************************* 
def csvs2db(schema: str, table: str, columns: List[str], fnames_csv: List[str], progbar: bool):
    """
    Directly load a batch of CSV files into the named DB table.
    INPUTS:
        schema:     Schema of the DB table
        table:      Name of the DB table
        columns:    List of columns of the DB table
        fnames_csv: List of the CSV file names
        progbar:    Whether to show a tqdm progress bar
    OUTPUTS:
        staging_tables: List of staging table names  
        Modifies the database by creating a staging table
    """
    # Display progress bar if requested
    ii = list(range(len(fnames_csv)))
    if progbar:
        ii = tqdm_auto(ii)
    
    # All SQL statements on this engine use the same connection
    with db_engine.connect() as conn:
        for i in ii:
            fname_csv = fnames_csv[i]
            csv2db(schema=schema, table=table, columns=columns, fname_csv=fname_csv, conn=conn)

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
        staging_table: Name of the staging dable in the DB
        Modifies the database by creating a staging table
    """

    # Destination and staging table names
    dest_table = dest_table_name(schema=schema, table=table)
    staging_table = staging_table_name(table=table, i=i)
    # List of column names
    col_list = '(' + ','.join(columns) + ')'

    # SQL to create staging table
    sql_clone_table = f"CREATE OR REPLACE TABLE {staging_table} LIKE {dest_table};"

    # SQL to Load CSV into database into staging table
    sql_load_csv = \
        f"""
        LOAD DATA LOCAL INFILE 
        '{fname_csv}'
        REPLACE INTO TABLE {staging_table}
        FIELDS TERMINATED BY ','
        LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        {col_list}
        """

    # Clone the table and load the CSV file
    conn.execute(sql_clone_table)
    conn.execute(sql_load_csv)
    # Delete the CSV file after it has been successfully loaded
    os.remove(fname_csv)

    # Return the name of the staging table
    return staging_table

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
        staging_tables: List of staging table names  
        Modifies the database by creating a staging table
    """
    # Display progress bar if requested
    if progbar:
        ii = tqdm_auto(ii)
    # List of staging tables assembled
    staging_tables: List[str] = []
    # Use the DB engine corresponding to this CPU number
    db_engine_c = db_engines[c]
    # All SQL statements on this engine use the same connection
    with db_engine_c.connect() as conn:
        for i in ii:
            staging_table = \
                csv2db_stage(schema=schema, table=table, columns=columns, fname_csv=fnames_csv[i], i=i, conn=conn)
            staging_tables.append(staging_table)

    return staging_tables

# ********************************************************************************************************************* 
def stage2db(dest_table: str, staging_table: str, columns: List[str], conn: conn_type):
    """
    Finishing loading a CSV file by inserting from the staging table into the named DB table.
    INPUTS:
        dest_table:     Fully delimited name of the destination DB table
        staging_table:  Fully delimited name of the staging DB table (source for insert)
        conn:           DB connection object
    OUTPUTS:
        None.  Modifies the database by inserting from the staging table into the main table
    """

    # List of column names
    col_list = '(' + ','.join(columns) + ')'
    # Columns to select
    select_fields = ',\n'.join(columns)

    # SQL to insert into the main table from the staging table
    sql_insert_stage = \
        f"""
        REPLACE INTO {dest_table}
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
def csvs2db_stage(schema: str, table: str, 
                  fnames_csv: Optional[List[str]]=None, 
                  columns: Optional[List[str]]=None,
                  single_thread: bool=False, 
                  progbar: bool=True):
    """
    Load a batch of CSVs into the named DB table.
    First delegates to csv2db_stage to populate the staging tables in parallel.
    Then rolls up from all of the staging tables into the big table at the end.
    """
    # Set DB engine pool based on single_thread flag
    make_db_engines(single_thread=single_thread)

    # Get list of CSV files if they were not provided
    if fnames_csv is None:
        search_path = os.path.join(dir_csv, table, f'pid_*', f'{table}-chunk-*.csv')
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
    engine_count: int = len(db_engines)
    cpu_count_default: int = min(multiprocessing.cpu_count() // 2, chunk_count, cpu_max, engine_count)
    cpu_count: int = 1 if single_thread else cpu_count_default
    
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
    staging_tables_pool = pool.starmap(csvs2db_stage, stage_inputs)
    # Flatten the list of staging table collections into one complete 
    staging_tables = list(itertools.chain.from_iterable(staging_tables_pool))

    # The destination table name
    dest_table = dest_table_name(schema=schema, table=table)
    # Roll up from staging table to schema on one CPU core
    with db_engine.connect() as conn:
        for staging_table in staging_tables:
            stage2db(dest_table=dest_table, staging_table=staging_table, columns=columns, conn=conn)

# ********************************************************************************************************************* 
def df2db(df: pd.DataFrame, schema: str, table: str, columns: List[str]=None, 
          chunksize: int=0, single_thread: bool=False, use_stage:bool=False,
          verbose: bool=False, progbar: bool=True):
    """
    Insert the contents of a Pandas DataFrame into a SQL table.
    INPUTS:
        df:         The DataFrame to insert
        schema:     The schema of the destination DB table
        table:      The name of the destination DB table
        columns:    List of columns to insert; read from DB metadata if omitted
        chunksize:  Number of rows in each chunk
        single_thread: Whether to run single threaded
        verbose:    Verbosity of output (true / false)
    OUTPUTS:
        None.  Modifies the DB table on the server.
    """
    # When running in single_thread mode, chunksize must be zero
    if single_thread:
        chunk_size = 0

    # Get columns from DB metadata if they were not provided by caller
    if columns is None:
        columns = get_columns(schema=schema, table=table)

    # Filter DataFrame to the desired columns
    df = df[columns]
    row_count: int = df.shape[0]

    # File name of CSV - include process id to avoid collisions when running multiple program instances
    fname_csv = os.path.join(dir_csv, f'{table}', f'pid_{pid}', f'{table}.csv')

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

    # Insert from the CSVs into the DB table using multiprocessing
    try:
        if use_stage:
            csvs2db_stage(schema=schema, table=table, fnames_csv=fnames_csv, columns=columns,
                    single_thread=single_thread, progbar=progbar)
        else:
            csvs2db(schema=schema, table=table, columns=columns, fnames_csv=fnames_csv, progbar=progbar)
    except:
        print('df2db failed!')
        print(f'Table {schema}.{table}.')
        print('Columns:\n', columns)
        raise

    # If we reach this point, the whole operation was a success.  
    # Clean up the now empty directory that previously had the CSV files
    Path(fname_csv).parent.rmdir()

    # Report elapsed time if requested
    t2 = time.time()
    if verbose:
        print_time(time=(t2-t1), msg='Elapsed Time for DB insertion')
