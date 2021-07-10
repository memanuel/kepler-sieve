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
    using ks::print_stars;
    using ks::print_newline;
    using ks::report_test;
    using ks::is_close_abs;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
#include "StateVector.hpp"
    using ks::StateVector;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "rebound_utils.hpp"
    using ks::Particle;
    using ks::Orbit;
    using ks::Simulation;
    using ks::make_sim;

// *****************************************************************************
// Functions defined in this module
int main();
bool test_all(db_conn_type& conn);
bool test_massive_body();
bool test_make_sim();
bool test_make_sim_planets();

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
    // Current test result
    bool is_ok = true;

    // Accumulate overall test results
    bool is_ok_all = true;

    // Test massive body
    // is_ok = test_massive_body();
    // is_ok_all &= is_ok;
    // report_test("Test: Build MassiveBody", is_ok);

    // Test making an empty simulation
    print_stars(true);
    is_ok = test_make_sim();
    is_ok_all &= is_ok;
    report_test("Test: Build empty rebound Simulation", is_ok);

    // Test making a simulation with the planets
    // is_ok = test_make_sim();

    // Report overall test results
    print("\n");
    print_stars();
    report_test("Test Suite on rebound", is_ok);
    return is_ok;
}

// *****************************************************************************
/// Test loading MassiveBody class from disk.
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

// *****************************************************************************
/// Test building an empty rebound Simulation object.
bool test_make_sim()
{
    /// Build the simulation
    Simulation* sim = make_sim();

    // Status
    print("Built empty rebound simulation.\n");
    print("N: {:d}.\n", sim->N);
    print("t: {:f}.\n", sim->t);
    print("G: {:e}.\n", sim->G);

    // Test conditions
    bool is_ok = (sim->N==0) && (sim->G = G);

    // Free memory in sinulation object
    reb_free_simulation(sim);

    return is_ok;
}
