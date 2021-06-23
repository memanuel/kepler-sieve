// *****************************************************************************
// Included libraries
#include <iostream>
#include <string>
#include <boost/format.hpp>

// Local dependencies
#include "db_utils.h"

// *****************************************************************************
// Names used
using std::cout;
using std::string;
using std::pair;
using boost::format;

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
    cout << format("\nGot DB connection to host %s, schema %s\n") % db_host % schema;

    // Close Connection
    conn->close();

    // Normal program exit
    return 0;
}
