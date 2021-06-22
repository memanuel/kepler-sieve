#include <iostream>
#include <string>
#include <boost/format.hpp>
#include <mariadb/conncpp.hpp>

using std::cout;
using std::string;
using std::pair;
using boost::format;

int main(int argc, char *argv[])
{

    // Username and password for DB connection
    pair<string, string> p1("user", "kepler");
    pair<string, string> p2("password", "kepler");
    cout << format("DB connection settings:");
    cout << format("%s    : %s\n") % p1.first % p1.second;
    cout << format("%s: %s\n") % p2.first % p2.second;

    // Instantiate Driver
    sql::Driver* driver = sql::mariadb::get_driver_instance();

    // Configure Connection
    sql::SQLString url("jdbc:mariadb://Thor/KS");
    sql::Properties properties({p1, p2});

    // Establish Connection
    std::unique_ptr<sql::Connection> conn(driver->connect(url, properties));

    // Normal program exit
    return 0;
}