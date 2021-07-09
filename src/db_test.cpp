/** @file   db_test.cpp
 *  @brief  Test harness for accessing SQL database
 *
 *  @author Michael S. Emanuel
 *  @date   2021-06-22
 * 
 * Example call:
 * ./db_test.x
 */

// *****************************************************************************
// Library dependencies
#include <string>
    using std::string;
#include <utility>
    using std::pair;
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "db_utils.hpp"
    using ks::get_db_conn;

// *****************************************************************************
/// Test database utilities.
int main()
{

    // Establish Connection
    std::unique_ptr<sql::Connection> conn = get_db_conn(true);

    // Do something with the connection object
    SQLString db_host = conn->getHostname();
    SQLString schema = conn->getSchema();
    print("\nGot DB connection to host {:s}, schema {:s}\n", db_host, schema);

    // Close Connection
    conn->close();

    // Normal program exit
    return 0;
}
