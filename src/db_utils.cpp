/** @file   db_utils.cpp
 *  @brief  Implementation of db_utils.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-06-25
 * 
 */

// *****************************************************************************
// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::sql_stmt_type;
    using ks::sql_prepared_stmt_type;

// *****************************************************************************
namespace ks {

// *****************************************************************************
db_conn_type get_db_conn(bool verbose)
{
    // Username and password for DB connection
    string db_host {"Thor"};
    string schema {"KS"};
    string user {"kepler"};
    string password {"kepler"};

    // Default timeout of 30000 ms (30 seconds) too short; set to 5 minutes
    int timeout_min = 5;
    int timeout_ms = timeout_min * 60 * 1000;

    // Wrap all of this into a sql Properties object
    pair<string, string> p1("user", user);
    pair<string, string> p2("password", password);
    pair<string, string> p3("connectTimeout", std::to_string(timeout_ms));
    Properties properties({p1, p2, p3});

    // Configure DB connection
    string url = format("jdbc:mariadb://{:s}/{:s}", db_host, schema);

    // Report DB connection settings if requested
    if (verbose)
    {
        print("DB connection settings:\n");
        print("{:16s}  : {:16s}\n", "DB host", db_host);
        print("{:16s}  : {:16s}\n", "schema", schema);
        print("{:16s}  : {:16s}\n", p1.first, p1.second);
        print("{:16s}  : {:16s}\n", p2.first, p2.second);
        print("{:16s}  : {:16s}\n", p3.first, p3.second);
    }

    // Instantiate DB driver
    Driver* driver = get_driver_instance();

    // Establish DB connection
    unique_ptr<Connection> conn(driver->connect(url, properties));

    // Return the DB connection
    return conn;
}

// *****************************************************************************
int result_set_size(ResultSet *rs)
{
    // Get current row
    int current_row = rs->getRow();
    // Go to last row and get its row number
    rs->last();
    int rows = rs->getRow();
    // Return to starting row
    // rs.first();
    rs->absolute(current_row);
    return rows;
}

// *****************************************************************************
//*Join a stored proecedure name and a vector of string parameters into a single SQL statement
string sql_sp_bind_params(const string sp_name, const vector<string>& params)
{
    return format("CALL {}({});", sp_name, join(params, ", "));
}

// *****************************************************************************
ResultSet* sp_run(db_conn_type &conn, const string sp_name, const vector<string> &params)
{
    // Create a new SQL statement
    sql_stmt_type stmt(conn->createStatement());
    // Bind the vector of SQL parameters to the named SP
    string sql = sql_sp_bind_params(sp_name, params);
    // print("SP call with bound parameters:\n{:s}\n", sql);

    // Execute stored procedure into a SQL resultset object
    ResultSet* rs = stmt->executeQuery(sql);
    // std::unique_ptr<ResultSet> rs(stmt->executeQuery(sql));

    // Workaround to MariaDB behavior complaining about statements being out of sync:
    // Make sure there are no more ResultSets
    while (stmt->getMoreResults())
    {
        ResultSet* throwaway = stmt->getResultSet();
        throwaway->close();
        delete throwaway;
    }

    // Return the resultset; consumer MUST run rs->close() and delete rs when done!
    return rs;
}

// *****************************************************************************
int sp_run_int(db_conn_type &conn, const string sp_name)
{
    // Create a SQL statement
    sql_stmt_type stmt(conn->createStatement());
    // The SQL string
    string sql = format("CALL {:s}();", sp_name);
    // Execute stored procedure into a SQL resultset object
    ResultSet* rs = stmt->executeQuery(sql);
    // Workaround for MariaDB sync error
    while (stmt->getMoreResults()) {
        ResultSet* throwaway = stmt->getResultSet();
        throwaway->close();
    }
    // Get the result, which is the value in the first column of the [only] row in the output
    rs->next();
    int res = rs->getInt(1);
    // Close resultset and free memory
    rs->close();
    delete rs;
    // Now return the result
    return res;
}

// *****************************************************************************
} // Namespace ks
