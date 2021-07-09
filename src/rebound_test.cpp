/** @file   rebound_test.cpp
 *  @brief  Test harness for functions used to build and run rebound simulations.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-09
 * 
 * Example call:
 * ./test_rebound.x
 */

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;


// Local dependencies
#include "utils.hpp"
    using ks::is_close_abs;
    using ks::print_stars;
    using ks::report_test;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;

// *****************************************************************************
// Functions defined in this module
int main();
bool test_all(db_conn_type& conn);
bool test_massive_body();

// *****************************************************************************
int main()
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run all tests
    bool is_ok = test_all(conn);

    // Close DB connection
    conn->close();

    // Normal program exit; return 1 to signal test failure
    return is_ok ? 0 : 1;
}

// *****************************************************************************
bool test_all(db_conn_type& conn)
{
    // Accumulate overall test results
    bool is_ok = true;

    // Test massive body
    is_ok = is_ok && test_massive_body();

    // Report overall test results
    print("\n");
    print_stars();
    report_test("Test Suite on rebound", is_ok);
    return is_ok;
}

// *****************************************************************************
bool test_massive_body()
{
    // Load from disk
    MassiveBodyTable mbt = MassiveBodyTable();

    // Print contents
    print_stars();
    print("Massive Body:\n");
    print("{:8s} : {:8s} : {:8s}\n", "BodyID", "M", "GM");
    for (int32_t body_id: mbt.get_body_id())
    {
        // Only print the planets; don't show the heavy asteroids
        if (body_id > 1000000) {continue;}
        // Access this body and print it
        MassiveBody mb = mbt[body_id];
        print("{:8d} : {:8.2e} : {:8.2e}\n", mb.body_id, mb.M, mb.GM);
    }

    // Expected mass of Sun and Earth
    double M_sun = 1.0;
    double M_earth = 0.0000030034896145;
    // Test tolerance
    double tol = 1.0E-12;

    // Test that the Sun has mass 1.0
    bool is_ok_sun = (mbt.get_M(10) == M_sun);
    // Test that the Earth has mass close to expected result
    bool is_ok_earth = is_close_abs(mbt.get_M(399), M_earth, tol);
    // Overall test result
    bool is_ok = is_ok_sun && is_ok_earth;

    // Report results
    report_test("MassiveBody: check Sun and Earth", is_ok);
    return is_ok;
}
