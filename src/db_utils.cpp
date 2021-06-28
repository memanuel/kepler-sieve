/*****************************************************************************
 * Michael S. Emanuel
 * 2021-06-25
 * ****************************************************************************/

// *****************************************************************************
// Local dependencies
#include "db_utils.hpp"

// *****************************************************************************
// Names used
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

    // Wrap username/password into sql Properties object
    pair<string, string> p1("user", user);
    pair<string, string> p2("password", password);
    Properties properties({p1, p2});

    // Configure DB connection
    string url = format("jdbc:mariadb://{:s}/{:s}", db_host, schema);

    // Report DB connection settings if requested
    if (verbose)
    {
        print("DB connection settings:\n");
        print("DB host : {:s}\n", db_host);
        print("schema  : {:s}\n", schema);
        print("{:s}    : {:s}\n", p1.first, p1.second);
        print("{:s}: {:s}\n", p2.first, p2.second);
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
string sql_sp_bind_params(const string sp_name, const vector<string> &params)
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
    ResultSet *rs = stmt->executeQuery(sql);

    //Workaround: Makes sure there are no more ResultSets
    while (stmt->getMoreResults()) {
        ResultSet *throwaway = stmt->getResultSet();
        throwaway->close();
    }

    // Return the resultset
    return rs;
}

// *****************************************************************************
} // Namespace ks
