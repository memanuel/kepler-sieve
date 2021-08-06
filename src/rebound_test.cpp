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
#include "rebound_utils.hpp"
    using ks::reb::Particle;
    using ks::reb::Orbit;
#include "StateVector.hpp"
    using ks::StateVector;
    using ks::dist;
    using ks::dist_dq;
    using ks::dist_dv;
    using ks::is_close;
    using ks::is_equal;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "Simulation.hpp"    
    using ks::reb::Simulation;
    using ks::reb::make_sim_planets;
    using ks::reb::make_sim_planets_horizons;
    using ks::reb::make_sim_de435_horizons;

// *****************************************************************************
// Constants used in this module

/// The test date for the simulation
constexpr double epoch = 59000.0;

/// The length of time for the integration test
constexpr double integration_test_time = 1000.0;

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
bool test_all(db_conn_type& conn);
bool test_massive_body();
bool test_make_sim();
bool test_make_sim_planets(const PlanetVector& pv, bool verbose);
bool test_make_sim_planets_horizons(db_conn_type& conn, bool verbose);
bool test_integration(Simulation& sim0, Simulation& sim1, double tol_dq, double tol_dv, bool verbose);
bool test_copy(const PlanetVector& pv);

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
    // Current test result and 
    bool is_ok = true;
    string test_name = "";
    // Accumulate overall test results
    bool is_ok_all = true;
    // Set overall verbosity level
    bool verbose = false;

    // *****************************************************************************
    // Simple tests - build splined vectors and elements and empty simulations
    // *****************************************************************************

    // Build PlanetVector object used in tests
    PlanetVector pv = PlanetVector(mjd0_planet, mjd1_planet, dt_min, true);
    print("Built PlanetVector from {:d} to {:d} with dt_min={:d}.\n", mjd0_planet, mjd1_planet, dt_min);

    // Build PlanetElement object used in tests
    PlanetElement pe = PlanetElement(mjd0_planet, mjd1_planet, dt_min, true);
    print("Built PlanetElement from {:d} to {:d} with dt_min={:d}.\n", mjd0_planet, mjd1_planet, dt_min);

    // Test massive body
    is_ok = test_massive_body();
    is_ok_all &= is_ok;
    report_test("Test: Build MassiveBody", is_ok);

    // Test making an empty simulation
    is_ok = test_make_sim();
    is_ok_all &= is_ok;
    report_test("Test: Build empty rebound simulation", is_ok);

    // *****************************************************************************
    // Test building planet simulations from KS and Horizons data
    // *****************************************************************************

    // Test making a simulation with the planets
    is_ok = test_make_sim_planets(pv, verbose);
    is_ok_all &= is_ok;
    test_name = format("Test: Build rebound simulation with planets at epoch {:8.2f}", epoch);
    report_test(test_name, is_ok);

    // Test making a simulation with the planets
    is_ok = test_make_sim_planets_horizons(conn, verbose);
    is_ok_all &= is_ok;
    test_name = format("Test: Build rebound simulation with planets at epoch {:8.2f} from Horizons data", epoch);
    report_test(test_name, is_ok);

    // *****************************************************************************
    // Test consistency of integrated simulation with expected results
    // *****************************************************************************

    // Test integration consistency on integer dates (these are spline nodes); spline using vectors
    {
    // Date range for test on node dats
    constexpr double mjd0 = mjd0_integrate;
    constexpr double mjd1 = mjd1_integrate;
    // Tolerance for positions - on spline nodes
    double tol_dq = 1.0E-11;
    double tol_dv = 1.0E-13;
    // Build splines and run test
    Simulation sim0_node = make_sim_planets(pv, mjd0);
    Simulation sim1_node = make_sim_planets(pv, mjd1);
    is_ok = test_integration(sim0_node, sim1_node, tol_dq, tol_dv, verbose);
    is_ok_all &= is_ok;
    test_name = \
        format("Test: Consistency of integration between {:8.2f} and {:8.2f} with splined vectors", mjd0, mjd1);
    report_test(test_name, is_ok);
    }   // block for test on spline nodes

    // Test integration consistency on non-integer start date; exercise element spline
    {
    // Date range for test off node dates
    constexpr double mjd0 = mjd0_integrate+0.5;
    constexpr double mjd1 = mjd1_integrate;
    // Tolerance for positions - off spline nodes
    double tol_dq = 1.0E-5;
    double tol_dv = 1.0E-6;
    // Build splines and run test
    Simulation sim0_spline = make_sim_planets(pe, mjd0);
    Simulation sim1_spline = make_sim_planets(pv, mjd1);
    is_ok = test_integration(sim0_spline, sim1_spline, tol_dq, tol_dv, verbose);
    is_ok_all &= is_ok;
    test_name = \
        format("Test: Consistency of integration between {:8.2f} and {:8.2f} with splined elements", mjd0, mjd1);
    report_test(test_name, is_ok);
    }   // block for test off spline nodes

    // Test copying a simulation
    is_ok = test_copy(pv);
    is_ok_all &= is_ok;
    report_test("Test: Copy rebound simulation with planets", is_ok);

    // Report overall test results
    print_stars(true);
    report_test("Test Suite on rebound", is_ok);
    return is_ok;
}

// *****************************************************************************
/// Test loading MassiveBody class from disk.
bool test_massive_body()
{
    // Load from disk
    MassiveBodyTable mbt = MassiveBodyTable();

    // Print contents of MassiveBodyTable
    print_stars(true);
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

    // Status
    print_stars(true);
    print("Built empty rebound simulation.\n");
    print("N: {:d}.\n", sim.N());
    print("t: {:f}.\n", sim.t());

    // Test conditions
    bool is_ok = (sim.N()==0) && (sim.t() == 0.0);

    // Return results of test
    return is_ok;
}

// *****************************************************************************
/// Test building a rebound Simulation with the planets collection.
bool test_make_sim_planets(const PlanetVector& pv, bool verbose)
{
    /// Build the simulation
    Simulation sim = make_sim_planets(pv, epoch);

    // Status
    print_stars(true);
    print("Built rebound simulation for planets.\n");
    print("N        : {:d}.\n", sim.N());
    print("N_active : {:d}.\n", sim.N_active());
    print("N_test   : {:d}.\n", sim.N_test());
    print("t        : {:f}.\n", sim.t());

    // Display the particles
    if (verbose) {sim.print_vectors();}

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
/// Test building a rebound Simulation with the planets collection initialized with Horizons data.
bool test_make_sim_planets_horizons(db_conn_type& conn, bool verbose)
{
    /// Build the simulation
    Simulation sim = make_sim_planets_horizons(conn, epoch);
    // Simulation sim = make_sim_de435_horizons(conn, epoch);

    // Status
    print_stars(true);
    print("Built rebound simulation for planets with Horizons data.\n");
    sim.print_summary();
    // print("N        : {:d}.\n", sim.N());
    // print("N_active : {:d}.\n", sim.N_active());
    // print("N_test   : {:d}.\n", sim.N_test());
    // print("t        : {:f}.\n", sim.t());

    // Display the particles
    if (verbose) {sim.print_vectors();}

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
    // Index of maximum distance
    int dq_argmax=0;
    int dv_argmax=0;

    // Find the maximum error
    for (int i=0; i<sim0.N(); i++)
    {
        // The two state vectors
        StateVector s0 = sim0.state_vector(i);
        StateVector s1 = sim1.state_vector(i);
        // The position difference
        double dq = dist_dq(s0, s1);
        // The velocity distance
        double dv = dist_dv(s0, s1);
        // The maximum difference so far
        if (dq>dq_max) {dq_max=dq; dq_argmax=i;}
        if (dv>dv_max) {dv_max=dv; dv_argmax=i;}
        // Is this test OK to the given tolerances?
        is_ok &= is_close(s0, s1, tol_dq, tol_dv);
    }

    // The test that was run
    print_stars(true);
    print("Integration Test from {:8.2f} to {:8.2f}:\n", mjd0, mjd1);
    // Print the state vectors if in verbose mode
    if (verbose)
    {
        // Print simulation state vectors of both simulations
        print("Simulation 0: mjd {:8.2f} integrated forward to {:8.2f}.\n", mjd0, mjd1);
        sim0.print_vectors();
        print("Simulation 1: mjd {:8.2f} interpolated from disk.\n", mjd1);
        sim1.print_vectors();

        // Print orbital elements
        print("Simulation 0 orbital elements.\n");
        sim0.print_elements();
        print("Simulation 1 orbital elements.\n");
        sim1.print_elements();
    }

    // Print the largest difference in position and velocity
    string body_name_dq_max = get_body_name(sim0.body_ids[dq_argmax]);
    string body_name_dv_max = get_body_name(sim0.body_ids[dv_argmax]);
    print("Largest difference:\n");
    print("Position: {:9.2e} AU     ({:10s}).\n", dq_max, body_name_dq_max);
    print("Velocity: {:9.2e} AU/day ({:10s}).\n", dv_max, body_name_dv_max);

    // Return the test result
    return is_ok;
}

// *****************************************************************************
/// Test that copying a Simulation object works as expected
bool test_copy(const PlanetVector& pv)
{
    // Build rebound simulation with planets at reference time
    Simulation sim0 = make_sim_planets(pv, epoch);

    // Test orbital elements for Juno @ 59000
    // Copy / paste from KS.GetAsteroidElements(3, 4, 59000, 59000);
    constexpr OrbitalElement elt
    {
        .mjd    =  59000.0,
        .a      =  2.6682852999999986,
        .e      =  0.25693643000000027,
        .inc    =  0.22673642125828372,
        .Omega  =  2.9644675653853,
        .omega  =- 1.9536135288017569,
        .f      = 46.517431678528794,
        .M      = 46.171557086433715
    };
    int32_t candidate_id = 1000003;

    // Add asteroid to simulation with candidate elements
    sim0.add_test_particle(elt, candidate_id);

    // Original simulation
    print_stars(true);
    print("Simulation summary: planets + 1 asteroid.\n");
    sim0.print_summary();   
    sim0.print_vectors();

    // Copy the test result
    Simulation sim1 = sim0.copy();
    print("Simulation sim1 = sim0.copy()\n");
    // sim1.print_summary();
    // sim1.print_vectors();

    // Test result
    bool is_ok = (sim0.N() == sim1.N()) && (sim0.N_active() == sim1.N_active()) && (sim0.N_test() == sim1.N_test());
    is_ok &= (sim0.t() == sim1.t());
    // Check that all N particles have equal state vectors and masses
    for (int i=0; i<sim0.N(); i++)
    {
        // Check state vectors are equal
        const StateVector& s0 = sim0.state_vector(i);
        const StateVector& s1 = sim1.state_vector(i);
        is_ok &= is_equal(s0, s1);
        // Check masses are equal
        is_ok &= (sim0.particle(i).m == sim1.particle(i).m);
    }

    // Return the test result
    return is_ok;
}