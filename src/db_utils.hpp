/** @file   db_utils.hpp
 *  @brief  Utilities for working with MariaDB SQL databases.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-06-25
 * 
 */

/* ****************************************************************************/
#pragma once

// *****************************************************************************
// Library dependencies
#include <string>
    using std::string;
    using std::string_view;
#include <vector>
    using std::vector;
#include <utility>
    using std::pair;
#include <memory>
    using std::unique_ptr;
#include <fmt/format.h>
    using fmt::print;
    using fmt::format;
#include <boost/algorithm/string/join.hpp>
    using boost::algorithm::join;
#include <mariadb/conncpp.hpp>
    using sql::SQLString;
    using sql::Driver;
    using sql::Properties;
    using sql::Connection;
    using sql::Statement;
    using sql::ResultSet;
    using sql::mariadb::get_driver_instance;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Type names
using db_conn_type = std::unique_ptr<sql::Connection>;
using sql_stmt_type = std::unique_ptr<sql::Statement>;
using sql_prepared_stmt_type = std::unique_ptr<sql::PreparedStatement>;

// *****************************************************************************
// Functions defined in db_utils
// *****************************************************************************

// *****************************************************************************
/** Get a connection to the Thor database server running a MariaDB instance.
 *  Connect to the default schema "KS".
 *  Use credentials for user "kepler".
 * */
db_conn_type get_db_conn(bool verbose=false);

// *****************************************************************************
/// Count the number of rows in a SQL resultset
int result_set_size(ResultSet* rs);

// *****************************************************************************
/// Helper function for SQL stored procedures- bind parameters into one SQL string.
string sql_sp_bind_params(const string sp_name, const vector<string> &params);

// *****************************************************************************
/** Execute a stored procedure and return a SQL resultset object.
 *  Important detail: consumers MUST run rs->close() and delete rs when the resultset no longer needed!
 *  Otherwise program will have a memory leak. */
ResultSet* sp_run(db_conn_type& conn, const string sp_name, const vector<string>& params);

/// Run a stored procedure statement that is known at compile time
ResultSet* sp_run(db_conn_type& conn, const char* sp_call);

// *****************************************************************************
/// Execute a stored procedure that returns a single integer
int sp_run_int(db_conn_type& conn, const string sp_name);

// *****************************************************************************
/// Wrap a string in double quotes
const string wrap_string(const string s);

// *****************************************************************************
} // Namespace ks
