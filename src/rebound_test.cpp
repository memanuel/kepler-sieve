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
#include <algorithm>
    using std::max;
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
    using ks::dist;
    using ks::is_close;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "rebound_utils.hpp"
    using ks::reb::Particle;
    using ks::reb::Orbit;
    using ks::reb::Simulation;
    using ks::reb::make_sim_planets;

// *****************************************************************************
// Constants used in this module

/// The test date for the simulation
constexpr double epoch = 59000.0;

/// The length of time for the integration test
// DEBUG
// constexpr double integration_test_time = 1000.0;
constexpr double integration_test_time = 1.0;

/// The first date for the integration test
constexpr double mjd0_integrate = epoch;

/// The last date for the integration test
constexpr double mjd1_integrate = epoch + integration_test_time;

/// Padding for date range used to spline planet vectors
constexpr int pad = 32;

/// First mjd to load in PlanetVector and PlanetElement
constexpr int mjd0_planet = int(mjd0_integrate - pad);

/// Last mjd to load in PlanetVector
constexpr int mjd1_planet = int(mjd1_integrate + pad);

/// Time step in minutes to load PlanetVector
// constexpr int dt_min = 5;
constexpr int dt_min = 1440;

// *****************************************************************************
// Functions defined in this module
int main();
bool test_all();
bool test_massive_body();
bool test_make_sim();
bool test_make_sim_planets(const PlanetVector& pv);
bool test_integration(Simulation& sim0, Simulation& sim1, double tol_dq, double tol_dv, bool verbose);

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
    PlanetVector pv = PlanetVector(mjd0_planet, mjd1_planet, dt_min);
    pv.load();
    pv.build_splines();
    print("Built PlanetVector from {:d} to {:d} with dt_min={:d}.\n", mjd0_planet, mjd1_planet, dt_min);

    // Build PlanetElement object used in tests
    PlanetElement pe = PlanetElement(mjd0_planet, mjd1_planet, dt_min);
    pe.load();
    pe.build_splines();
    print("Built PlanetElement from {:d} to {:d} with dt_min={:d}.\n", mjd0_planet, mjd1_planet, dt_min);

    // Test massive body
    // is_ok = test_massive_body();
    // is_ok_all &= is_ok;
    // report_test("Test: Build MassiveBody", is_ok);

    // // Test making an empty simulation
    // print_stars(true);
    // is_ok = test_make_sim();
    // is_ok_all &= is_ok;
    // report_test("Test: Build empty rebound simulation", is_ok);

    // // Test making a simulation with the planets
    // print_stars(true);
    // is_ok = test_make_sim_planets(pv);
    // is_ok_all &= is_ok;
    // string test_name = format("Test: Build rebound simulation with planets at epoch {:8.2f}", epoch);
    // report_test(test_name, is_ok);

    // Test integration consistency variables
    bool verbose = true;
    string test_name = "";
    double mjd0=0.0, mjd1=0.0;

    // Tolerance for positions - on spline nodes
    double tol_dq_node = 1.0E-11;
    double tol_dv_node = 1.0E-13;

    // Tolerance for positions - off spline nodes
    double tol_dq_spline = 1.0E-10;
    double tol_dv_spline = 1.0E-12;

    // Test integration consistency on integer dates (these are spline nodes); spline using vectors
    mjd0 = mjd0_integrate;
    mjd1 = mjd1_integrate;
    Simulation sim0_node = make_sim_planets(pv, mjd0);
    Simulation sim1_node = make_sim_planets(pv, mjd1);
    is_ok = test_integration(sim0_node, sim1_node, tol_dq_node, tol_dv_node, verbose);
    is_ok_all &= is_ok;
    print_stars(true);
    test_name = format("Test: Consistency of integration between {:8.2f} and {:8.2f} with splined vectors", mjd0, mjd1);
    report_test(test_name, is_ok);

    // Test integration consistency on non-integer start date; exercise element spline
    mjd0 = mjd0_integrate;
    mjd1 = mjd1_integrate;
    Simulation sim0_spline = make_sim_planets(pe, mjd0);
    Simulation sim1_spline = make_sim_planets(pe, mjd1);
    is_ok = test_integration(sim0_spline, sim1_spline, tol_dq_spline, tol_dv_spline, verbose);
    is_ok_all &= is_ok;
    print_stars(true);
    test_name = format("Test: Consistency of integration between {:8.2f} and {:8.2f} with splined elements", 
                        mjd0, mjd1);
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
    Simulation sim;
    // Reference to the simulation
    // Simulation& sim = *s;

    // Status
    print("Built empty rebound simulation.\n");
    print("N: {:d}.\n", sim.N());
    print("t: {:f}.\n", sim.t());

    // Test conditions
    bool is_ok = (sim.N()==0) && (sim.t() == 0.0);

    // Return results of test
    return is_ok;
}

// *****************************************************************************
/// Test building an rebound Simulation with the planets collection.
bool test_make_sim_planets(const PlanetVector& pv)
{
    /// Build the simulation
    Simulation sim = make_sim_planets(pv, epoch);

    // Status
    print("Built rebound simulation for planets.\n");
    print("N        : {:d}.\n", sim.N());
    print("N_active : {:d}.\n", sim.N_active());
    print("N_test   : {:d}.\n", sim.N_test());
    print("t        : {:f}.\n", sim.t());

    // Display the particles
    sim.print();

    // Grab particles for Sun and Earth
    Particle p_sun = sim.particle(0);
    Particle p_earth = sim.particle(3);

    // Test conditions
    bool is_ok = (sim.N()==11) && (sim.t() == epoch);
    is_ok &= (p_sun.m==1.0) & is_close_rel(p_earth.m, 3.0E-6, 0.01);

    // Return results of test
    return is_ok;
}

// *****************************************************************************
/// Test that integrating a simulation forward matches the expected result
bool test_integration(Simulation& sim0, Simulation& sim1, double tol_dq, double tol_dv, bool verbose)
{
    // Get the two simulation times
    double mjd0 = sim0.t();
    double mjd1 = sim1.t();

    // Integrate sim0 to time of sim1
    sim0.integrate(mjd1);

    // Test result
    bool is_ok = true;
    // Maximum distance
    double dq_max = 0.0;
    double dv_max = 0.0;

    // Find the maximum error
    for (int i=0; i<sim0.N(); i++)
    {
        // The two state vectors
        StateVector s0 = sim0.state_vector(i);
        StateVector s1 = sim1.state_vector(i);
        // The position difference
        double dq = dist(s0, s1);
        // The velocity distance
        double dv = dist(sv2vel(s0), sv2vel(s1));
        // The maximum difference so far
        dq_max = max(dq, dq_max);
        dv_max = max(dv, dv_max);
        // Is this test OK to the given tolerances?
        is_ok &= is_close(s0, s1, tol_dq, tol_dv);
    }

    // The test that was run
    print("Integration Test from {:8.2f} to {:8.2f}:\n", mjd0, mjd1);
    // Print the state vectors if in verbose mode
    if (verbose)
    {
        // Print simulation 0 state vectors
        print("Simulation 0: mjd {:d} integrated forward to {:d}.\n", int(mjd0), int(mjd1));
        sim0.print();
        // Print simulation 1 state vectors
        print("Simulation 1: mjd {:d} loaded from disk.\n", int(mjd1));
        sim1.print();
    }

    // Print the largest difference in position and velocity
    print("Largest difference:\n");
    print("Position: {:+9.2e} AU.\n", dq_max);
    print("Velocity: {:+9.2e} AU/day.\n", dv_max);

    // Retutn the test result
    return is_ok;
}