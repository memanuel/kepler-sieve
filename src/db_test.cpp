/*****************************************************************************
 * Test harness for accessing SQL database
 * 
 * Michael S. Emanuel
 * 2021-06-22
 * ****************************************************************************/

// *****************************************************************************
// Included libraries
#include <string>
#include <fmt/format.h>

// Local dependencies
#include "db_utils.h"

// *****************************************************************************
// Names used
using std::string;
using std::pair;
using fmt::print;
using ks::get_db_conn;

// *****************************************************************************
/** Test database utilities.
*/
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
