def sp2df(sp_name: str, params: Dict=dict()):
    """Execute a SQL stored procedure and return a DataFrame."""
    sql_stmt = sp_bind_args(sp_name=sp_name, params=params)
    with db_engine.connect() as conn:
        df = pd.read_sql(sql_stmt, conn)
    return df