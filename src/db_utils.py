"""
Utilities for working with MariaDB Databases.

Michael S. Emanuel
2021-02-18
"""

# Data
import pandas as pd
import dask.dataframe
import sqlalchemy

# Algorithms
import multiprocessing
import itertools
import subprocess

# File system
from pathlib import Path
import os
import glob

# UI
from tqdm.auto import tqdm as tqdm_auto
import traceback
import time

# MSE imports
from config import ks_root
from utils import print_time
from db_config import db_engine, db_url

# Types
from typing import Optional, List, Dict
conn_type = sqlalchemy.engine.base.Connection

# ********************************************************************************************************************* 
# Directory for inserting CSV data into database; make if missing
dir_csv: str = os.path.join(ks_root, 'data', 'df2db')
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
        md = sqlalchemy.MetaData(bind=conn, schema=schema)
        md.reflect()
        tbl_md = md.tables[schema_table]
        # Get list of all available DB columns; exclude computed columns!
        columns = [c.name for c in tbl_md.columns if not c.computed]
    return columns

# ********************************************************************************************************************* 
def get_generated_columns(schema: str, table: str) -> List[str]:
    """Get list of generated column names for DB table"""
    # Delimited table name
    schema_table: str = f'{schema}.{table}'

    # Get DB metadata
    with db_engine.connect() as conn:
        md = sqlalchemy.MetaData(bind=conn, schema=schema)
        md.reflect()
        tbl_md = md.tables[schema_table]
        # Get list of all available DB columns; exclude computed columns!
        columns = [c.name for c in tbl_md.columns if c.computed]
    return columns
    
# ********************************************************************************************************************* 
def drop_columns(schema: str, table: str, cols: List[str], conn) -> None:
    """Drop the named columns"""
    # Special case - if cols is empty, do nothing
    if len(cols) == 0:
        return
    
    # Form the delimited table name
    schema_table = f'{schema}.{table}' if len(schema)>0 else table
    
    # SQL to drop the columns in just one statement
    sql = f'ALTER TABLE {schema_table}\n'
    for col in cols:
        sql += f'DROP COLUMN {col},\n'        
    sql = sql[:-2] + ';'

    # Exectue the drop columns statement
    conn.execute(sql)

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
    # Special case: if the DataFrame has zero rows, just skip it completely!
    # Otherwise would end up with a one line CSV that has the field names, which we don't want.
    if df.shape[0] == 0:
        fnames = []
        return fnames

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
def find_fnames_csv(table: str, verbose: bool) -> List[str]:
    """Generate a list of CSV file names for the given file name"""
    search_path = os.path.join(dir_csv, table, f'pid_*', f'{table}-chunk-*.csv')
    fnames_csv = sorted(glob.glob(search_path))

    # Report results
    if verbose:
        nf = len(fnames_csv)
        print(f'Found {nf} CSV files for table {table}.')
        if nf > 0:
            print(fnames_csv[0])

    return fnames_csv

# ********************************************************************************************************************* 
def csv2db(schema: str, table: str, columns: List[str], fname_csv: str):
    """
    Load one CSV file directly into the named DB table by invoking mariadb-import utility
    INPUTS:
        schema:    Schema of the DB table
        table:     Name of the DB table
        columns:   List of columns of the DB table
        fname_csv: Name of the CSV file
    OUTPUTS:
        None. Modifies the database table.
    """
    # List of column names
    col_list = ','.join(columns)

    # File name for loading
    Path(os.path.join(dir_csv, 'mariadb-import')).mkdir(parents=True, exist_ok=True)
    fname_load = os.path.join(dir_csv, 'mariadb-import', f'{table}.csv')

    # Rename this chunk file to AsteroidVectors.csv
    # This is a limitation of mariadb-import; the CSV file name must match the table name exactly
    os.rename(fname_csv, fname_load)

    # Arguments to run mariadb-import from subprocess
    args = [
        'mariadb-import',
        f'--defaults-file={mdbi_opt}',
        '--replace',
        # f'--columns={col_list}',
        # '--use-threads=1',
        '--silent',
        schema, 
        fname_load,
    ]

    # Run mariadb-import
    # print('\n', fname_csv)
    subprocess.run(args)

    # Remove the file that was loaded 
    os.remove(fname_load)

# ********************************************************************************************************************* 
def csv2db_ldi(schema: str, table: str, columns: List[str], fname_csv: str, conn: conn_type):
    """
    Load one CSV file directly into the named DB table using the Load Data Infile command.
    This is MUCH SLOWER than using mariadb-import
    INPUTS:
        schema:    Schema of the DB table
        table:     Name of the DB table
        columns:   List of columns of the DB table
        fname_csv: Name of the CSV file
        conn:      DB connection object
    OUTPUTS:
        None. Modifies the database table.
    """

    # Destination table name including schema
    dest_table = dest_table_name(schema=schema, table=table)
    # List of column names
    col_list = '(' + ','.join(columns) + ')'

    # SQL to Load CSV into database into staging table
    sql_load_csv = \
        f"""
        LOAD DATA LOCAL INFILE 
        '{fname_csv}'
        REPLACE
        INTO TABLE {dest_table}
        FIELDS TERMINATED BY ','
        LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        {col_list}
        """

    # Load the CSV file
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
    
    # Load the CSVs into the database
    for i in ii:
        fname_csv = fnames_csv[i]
        csv2db(schema=schema, table=table, columns=columns, fname_csv=fname_csv)

# ********************************************************************************************************************* 
def clean_empty_dirs(fnames_csv: List[str]):
    """Loop through a list of file names; clean out any folders that are now empty"""
    # Get set of distinct folders on this list of files
    folders = [str(Path(fname_csv).parent) for fname_csv in fnames_csv]
    folders_unq = sorted(set(folders))

    # Delete these folders if empty
    for folder in folders_unq:
        # Better to ask foregiveness
        try:
            Path(folder).rmdir()
        except OSError:
            pass

# ********************************************************************************************************************* 
def df2db(df: pd.DataFrame, schema: str, table: str, columns: List[str]=None, 
          chunksize: int=0, single_thread: bool=False, 
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

    # Insert from the CSVs into the DB table; single threaded, but mariadb-import is at least fast
    try:
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
