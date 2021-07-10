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
    using ks::is_close_rel;
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
    using ks::make_sim_planets;

// *****************************************************************************
// Constants used in this module

/// The test date for the simulation
constexpr double epoch = 59000.0;

/// Padding for date range used to spline planet vectors
constexpr int pad = 32;

/// First mjd to load in PlanetVector
constexpr int mjd0 = int(epoch-pad);

/// Last mjd to load in PlanetVector
constexpr int mjd1 = int(epoch+pad);

/// Time step in minutes to load PlanetVector
constexpr int dt_min = 5;

// *****************************************************************************
// Functions defined in this module
int main();
bool test_all();
bool test_massive_body();
bool test_make_sim();
bool test_make_sim_planets(const PlanetVector& pv);

// *****************************************************************************
int main()
{
    // // Establish DB connection
    // db_conn_type conn = get_db_conn();

    // Run all tests
    bool is_ok = test_all();

    // // Close DB connection
    // conn->close();

    // Normal program exit; return 1 to signal test failure
    return is_ok ? 0 : 1;
}

// *****************************************************************************
bool test_all()
{
    // Current test result
    bool is_ok = true;

    // Accumulate overall test results
    bool is_ok_all = true;

    // Build PlanetVector object used in tests
    PlanetVector pv = PlanetVector(mjd0, mjd1, dt_min);
    pv.load();
    pv.build_splines();
    print("Built PlanetVector from {:d} to {:d} with dt_min={:d}.\n", mjd0, mjd1, dt_min);

    // Test massive body
    // is_ok = test_massive_body();
    // is_ok_all &= is_ok;
    // report_test("Test: Build MassiveBody", is_ok);

    // Test making an empty simulation
    print_stars(true);
    is_ok = test_make_sim();
    is_ok_all &= is_ok;
    report_test("Test: Build empty rebound simulation", is_ok);

    // Test making a simulation with the planets
    print_stars(true);
    is_ok = test_make_sim_planets(pv);
    is_ok_all &= is_ok;
    string test_name = format("Test: Build rebound simulation with planets at epoch {:8.2f}", epoch);
    report_test(test_name, is_ok);

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
    // Reference to the simulation
    // Simulation& sim = *s;

    // Status
    print("Built empty rebound simulation.\n");
    print("N: {:d}.\n", sim->N);
    print("t: {:f}.\n", sim->t);
    print("G: {:e}.\n", sim->G);

    // Test conditions
    bool is_ok = (sim->N==0) && (sim->G == G);

    // Free memory in sinulation object
    reb_free_simulation(sim);

    // Return results of test
    return is_ok;
}

// *****************************************************************************
/// Test building an rebound Simulation with the planets collection.
bool test_make_sim_planets(const PlanetVector& pv)
{
    /// Build the simulation
    Simulation* sim = make_sim_planets(pv, epoch);
    // Reference to the simulation
    // Simulation& sim = *s;

    // Status
    print("Built rebound simulation for planets.\n");
    print("N: {:d}.\n", sim->N);
    print("t: {:f}.\n", sim->t);
    print("G: {:e}.\n", sim->G);

    // Display the particles
    print_newline();
    print("{:3s} {:10s}: {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} \n", 
          " i", " M", " qx", " qy", " qz", " vx", " vy", " vz");
    for (int i=0; i<sim->N; i++)
    {
        Particle p = sim->particles[i];
        print("{:3d} {:10.3e}: {:+10.6f} {:+10.6f} {:+10.6f} {:+10.6f} {:+10.6f} {:+10.6f} \n", 
              i, p.m, p.x, p.y, p.z, p.vx, p.vy, p.vz);    
    }
    print_newline();

    // Grab particles for Sun and Earth
    Particle p_sun = sim->particles[0];
    Particle p_earth = sim->particles[3];

    // Test conditions
    bool is_ok = (sim->N==11) && (sim->t == epoch) && (sim->G == G);
    is_ok &= (p_sun.m==1.0) & is_close_rel(p_earth.m, 3.0E-6, 0.01);

    // Free memory in simulation object
    reb_free_simulation(sim);

    // Return results of test
    return is_ok;
}
