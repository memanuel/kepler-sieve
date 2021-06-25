#pragma once

// *****************************************************************************
// Included files
#include <iostream>
#include <string>
#include <fmt/format.h>
#include <mariadb/conncpp.hpp>

// *****************************************************************************
// Standard library class names used
using std::cout;
using std::endl;
using std::string;
using std::pair;
using std::unique_ptr;
using fmt::print;
using fmt::format;

// SQL class names used (from MariaDB/Connector)
using sql::Driver;
using sql::SQLString;
using sql::Properties;
using sql::Connection;

// SQL functions used
using sql::mariadb::get_driver_instance;

// *****************************************************************************
// Functions defined in db_utils
unique_ptr<Connection> get_db_conn(bool verbose);
