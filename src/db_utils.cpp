#include "db_utils.h"
#include "sky_patch.h"

// *****************************************************************************
/** Get a connection to the Thor database server running a MariaDB instance.
 *  Connect to the default schema "KS".
 *  Use credentials for user "kepler".
 * */
unique_ptr<Connection> get_db_conn(bool verbose=false)
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
    string url = (format("jdbc:mariadb://%1%/%2%") % db_host % schema).str();
    // SQLString url("jdbc:mariadb://Thor/KS");

    // Report DB connection settings if requested
    if (verbose)
    {
        cout << format("DB connection settings:\n");
        cout << format("DB host : %s\n") % db_host;
        cout << format("schema  : %s\n") % schema;
        cout << format("%s    : %s\n") % p1.first % p1.second;
        cout << format("%s: %s\n") % p2.first % p2.second;
    }

    // Instantiate DB driver
    Driver* driver = get_driver_instance();

    // Establish DB connection
    unique_ptr<Connection> conn(driver->connect(url, properties));

    // Return the DB connection
    return conn;
}
