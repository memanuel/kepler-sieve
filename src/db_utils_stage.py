"""
Attempt to load DB tables in parallel using staging tables.
This was a FAILURE.  
Keeping code here for reference in case fragments prove useful in the future.
In particular, it has a use of parallel processing with a thread pool and starmap is syntactically correct.
It's just not fast for inserting records into MariaDB tables with LOAD DATA INSERT statements.

Michael S. Emanuel
2021-02-18
"""

# Algorithms
import itertools
import multiprocessing


# ********************************************************************************************************************* 
def make_db_engines(single_thread: bool) -> None:
    """Create a shared set of DB engine objects (global variable).  Used to support multithreading."""
    global db_engines
    if single_thread:
        db_engines = (db_engine,)
    else:
        from db_engine_pool import db_engines

# ********************************************************************************************************************* 
def staging_table_name(table: str, i: int):
    """Consistently generate schema_table and staging_table names. i is the chunk number."""
    # Staging table for this chunk
    staging_table = f'temp.{table}_pid_{pid}_chunk_{i:03d}'
    return staging_table

# ********************************************************************************************************************* 
def csv2db_stage(schema: str, table: str, columns: List[str], cols_to_drop: List[str],
                 fname_csv: str, i: int, conn: conn_type):
    """
    Load one CSV file into a staging table for the named DB table.
    INPUTS:
        schema:         Schema of the DB table
        table:          Name of the DB table
        columns:        List of columns of the DB table
        cols_to_drop:   Columns to drop after cloning tabl
        fname_csv:      Name of the CSV file
        i:              Chunk number
        conn:           DB connection object
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
    sql_clone_table_1 = f"CREATE OR REPLACE TABLE {staging_table} LIKE {dest_table};"
    sql_clone_table_2 = f"ALTER table {staging_table} transactional=default;"
    sql_clone_table_3 = f"ALTER table {staging_table} engine='Memory';"

    # SQL to Load CSV into database into staging table
    sql_load_csv = \
        f"""
        LOAD DATA LOCAL INFILE 
        '{fname_csv}'
        IGNORE 
        INTO TABLE {staging_table}
        FIELDS TERMINATED BY ','
        LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        {col_list}
        """

    # Clone the table; need to drop generated columns
    conn.execute(sql_clone_table_1)
    conn.execute(sql_clone_table_2)
    drop_columns(schema='', table=staging_table, cols=cols_to_drop, conn=conn)
    conn.execute(sql_clone_table_3)

    # Load the CSV file
    conn.execute(sql_load_csv)
    # Delete the CSV file after it has been successfully loaded
    os.remove(fname_csv)

    # Return the name of the staging table
    return staging_table

# ********************************************************************************************************************* 
def csvs2stage(schema: str, table: str, columns: List[str], fnames_csv: List[str], 
                  ii: List[int], c: int, progbar: bool):
    """
    Load a batch of CSV files into a staging table for the named DB table.
    INPUTS:
        schema:         Schema of the DB table
        table:          Name of the DB table
        columns:        List of columns of the DB table
        cols_to_drop:   List of columns to be dropped (e.g. generated columns)
        fnames_csv: List of the CSV file names
        ii:         List of chunk numbers for this batch
        c:          CPU (process) number; used to choose the DB engine from the shared pool
        progbar:    Whether to show a tqdm progress bar
    OUTPUTS:
        staging_tables: List of staging table names  
        Modifies the database by creating a staging table
    """
    # Calculate list of columns to drop
    cols_to_drop = get_generated_columns(schema=schema, table=table)

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
                csv2db_stage(schema=schema, table=table, columns=columns, cols_to_drop=cols_to_drop,
                fname_csv=fnames_csv[i], i=i, conn=conn)
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
    staging_tables_pool = pool.starmap(csvs2stage, stage_inputs)
    # Flatten the list of staging table collections into one complete 
    staging_tables = list(itertools.chain.from_iterable(staging_tables_pool))

    # The destination table name
    dest_table = dest_table_name(schema=schema, table=table)
    # Roll up from staging table to schema on one CPU core
    with db_engine.connect() as conn:
        for staging_table in staging_tables:
            stage2db(dest_table=dest_table, staging_table=staging_table, columns=columns, conn=conn)

