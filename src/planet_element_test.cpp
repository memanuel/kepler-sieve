/** @file   planet_element_test.cpp
 *  @brief  Test harness for classes PlanetElement and PlanetVector for splined elements / vectors of planets.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-06
 * 
 * Example call:
 * ./planet_element_test.x
 */
// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "constants.hpp"
    using ks::cs::body_id_earth;
    using ks::cs::body_ids_planets;
#include "astro_utils.hpp"
    using ks::get_body_name;
#include "utils.hpp"
    using ks::print_newline;
    using ks::print_stars;
    using ks::report_test;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::sv2pos;
    using ks::sv2vel;
    using ks::is_close;
    using ks::print_state_vector_headers;
    using ks::print_state_vector;
    using ks::print_state_vector_sci;    
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::print_orbital_element;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
#include "BodyVector.hpp"
    using ks::BodyVector;
#include "PlanetVector.hpp"
    using ks::PlanetVector;    
#include "PlanetElement.hpp"
    using ks::PlanetElement;    

// *****************************************************************************
// Constants in this module
namespace {
/// The date with expected state vectors of Earth
constexpr double mjd_test = 58000.0;

/// Absolute difference between interp_pos and interp_vec outputs (position part only)
constexpr double tol_dq_ip = 1.0E-13;
/// Absolute difference between interp_vec output and expected position; for splined vectors
constexpr double tol_dq_vec = 1.0E-12;
/// Absolute difference between interp_vec output and expected velocity; for splined vectors
constexpr double tol_dv_vec = 1.0E-14;

/// Absolute difference between interp_vec output and expected position; for splined elements
constexpr double tol_dq_elt = 1.0E-11;
/// Absolute difference between interp_vec output and expected velocity; for splined elements
constexpr double tol_dv_elt = 1.0E-13;
}

// *****************************************************************************
// Functions defined in this module
int main(int argc, char* argv[]);
bool test_all(db_conn_type& conn, bool verbose);
bool test_body_vector(const BodyVector& bv, bool verbose);
bool test_planet_vector(const PlanetVector& pv, bool verbose);
bool test_planet_element(const PlanetElement& pe, bool verbose);
bool test_planet_element_vs_vector(
    const PlanetVector& pv, const PlanetElement& pe, double mjd0, double mjd1, int dt_min, 
    double tol_dq, double tol_dv, bool verbose);
bool test_vectors(StateVector& s0, StateVector& s1, Position& q1_ip, string class_name,
                  double tol_dq, double tol_dv, bool verbose);
StateVector expected_earth_58000();

// *****************************************************************************
int main(int argc, char* argv[])
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Set verbosity
    bool verbose = false;

    // Run all the tests
    bool is_ok = test_all(conn, verbose);

    // Close DB connection
    conn->close();

    // Normal program exit; return 0 for success, 1 for failure
    return is_ok ? 0 : 1;
}

// *****************************************************************************
bool test_all(db_conn_type& conn, bool verbose)
{
    // Inputs used in testing
    constexpr int width = 1000;
    constexpr int pad = 32;
    constexpr int mjd0 = mjd_test - width;
    constexpr int mjd1 = mjd_test + width;    
    int dt_min = 5;

    // Timer object
    Timer t;

    // Results of current test
    bool is_ok;
    // Overall test results
    bool is_ok_all = true;

    // *****************************************************************************
    // Build and test BodyVector for Earth
    // *****************************************************************************

    // Build BodyVector for Earth
    print_stars(true);
    t.tick();
    BodyVector bv = BodyVector("Earth");
    print("Built BodyVector for Earth.\n");
    t.tock_msg();

    // Test planet vectors
    is_ok = test_body_vector(bv, verbose);
    report_test("Test BodyVector", is_ok);
    is_ok_all = is_ok_all && is_ok;

    // *****************************************************************************
    // Build and test PlanetVector vs. Earth @ 58000
    // *****************************************************************************

    // Build PlanetVector
    print_stars(true);
    t.tick();
    PlanetVector pv(mjd0-pad, mjd1+pad, dt_min);
    pv.load();
    pv.build_splines();
    // PlanetVector pv = PlanetVector();
    print("Built PlanetVector object from mjd0 {:d} to mjd1 {:d} with time step {:d} minutes.\n", 
            pv.mjd0, pv.mjd1, pv.dt_min);
    t.tock_msg();

    // Test planet vectors
    is_ok = test_planet_vector(pv, verbose);
    report_test("Test PlanetVector", is_ok);
    is_ok_all = is_ok_all && is_ok;

    // *****************************************************************************
    // Build and test PlanetElement vs. Earth @ 58000
    // *****************************************************************************

    // Build PlanetElement
    print_stars(true);
    t.tick();
    PlanetElement pe(mjd0-pad, mjd1+pad, dt_min);
    pe.load();
    pe.build_splines();
    // PlanetElement pe = PlanetElement();
    print("Built PlanetElement object from mjd0 {:d} to mjd1 {:d} with time step {:d} minutes.\n", 
            pe.mjd0, pe.mjd1, pe.dt_min);
    t.tock_msg();

    // Test planet elements
    is_ok = test_planet_element(pe, verbose);
    report_test("Test PlanetElement", is_ok);
    is_ok_all = is_ok_all && is_ok;

    // *****************************************************************************
    // Compare PlanetElement vs. PlanetVector on many node dates
    // *****************************************************************************
    {
    constexpr int dt_min_test = 1440;
    constexpr int offset_test_min = 0;
    constexpr double mjd_offset_test = offset_test_min*dpm;
    constexpr double mjd0_test = mjd0 + mjd_offset_test;
    constexpr double mjd1_test = mjd1 + mjd_offset_test - 1.0;
    constexpr double tol_dq = 1.0E-09;
    constexpr double tol_dv = 1.0E-10;
    // Run comparison test
    is_ok = test_planet_element_vs_vector(pv, pe, mjd0_test, mjd1_test, dt_min_test, tol_dq, tol_dv, verbose);
    string test_name = format("Test PlanetElement vs. PlanetVector on Nodes ({:8.2f}-{:8.2f}, dt_min {:4d}).", 
        mjd0_test, mjd1_test, dt_min_test);
    report_test(test_name, is_ok);
    is_ok_all = is_ok_all && is_ok;
    }

    // *****************************************************************************
    // Compare PlanetElement vs. PlanetVector on non-node dates to exercise spline
    // *****************************************************************************
    {
    // Build a new spline of orbital elements sampled only DAILY
    constexpr int dt_min = 60;
    PlanetElement pe(mjd0-pad, mjd1+pad, dt_min);
    pe.load();
    pe.build_splines();

    constexpr int dt_min_test = 1440;
    constexpr int offset_test_min = dt_min / 2;
    constexpr double mjd_offset_test = offset_test_min*dpm;
    constexpr double mjd0_test = mjd0 + mjd_offset_test;
    constexpr double mjd1_test = mjd1 + mjd_offset_test - 1.0;
    constexpr double tol_dq = 1.0E-6;
    constexpr double tol_dv = 1.0E-8;
    // Run comparison test
    is_ok = test_planet_element_vs_vector(pv, pe, mjd0_test, mjd1_test, dt_min_test, tol_dq, tol_dv, verbose);
    string test_name = format("Test PlanetElement vs. PlanetVector off Nodes ({:8.2f}-{:8.2f}, offset_min {:5d}).", 
        mjd0_test, mjd1_test, offset_test_min);
    report_test(test_name, is_ok);
    is_ok_all = is_ok_all && is_ok;
    }

    // *****************************************************************************
    // Overall results
    // *****************************************************************************

    // Return overall test result
    return is_ok_all;
}

// *****************************************************************************
/// Make expected Earth state vectors at MJD 58000
StateVector expected_earth_58000()
{
    // Expected state vector components of Earth @ MJD 58000
    // Copy / paste from KS.GetStateVectors_Earth(58000, 58000, 1440);
    double qx =  0.9583774949482733;
    double qy = -0.31579917956842507;
    double qz = -0.0001266099206526;
    double vx =  0.005191609325627999;
    double vy =  0.016244467239388778;
    double vz =  0.0000000527623193;

    // Wrap into a StateVector
    return StateVector{.qx=qx, .qy=qy, .qz=qz, .vx=vx, .vy=vy, .vz=vz};
}

// *****************************************************************************
/// Test BodyVector class with splined Earth vectors
bool test_body_vector(const BodyVector& bv, bool verbose)
{
    // Test vs. expected state vectors of Earth @ 58000
    StateVector s0 = expected_earth_58000();

    // Calculate interpolated state vectors of Earth on the test date
    StateVector s1 = bv.interp_vec(mjd_test);

    // Calculate interpolated position of Earth
    Position q1 = bv.interp_pos(mjd_test);

    // Report test results
    return test_vectors(s0, s1, q1, "BodyVector", tol_dq_vec, tol_dv_vec, verbose);
}

// *****************************************************************************
/// Test PlanetVector class with splined vectors for planets
bool test_planet_vector(const PlanetVector& pv, bool verbose)
{
    // Test vs. expected state vectors of Earth @ 58000
    StateVector s0 = expected_earth_58000();

    // Calculate interpolated state vectors of Earth on the test date
    StateVector s1 = pv.interp_vec(body_id_earth, mjd_test);

    // Calculate interpolated position of Earth
    Position q1 = pv.interp_pos(body_id_earth, mjd_test);

    // Report test results
    return test_vectors(s0, s1, q1, "PlanetVector", tol_dq_vec, tol_dv_vec, verbose);
}

// *****************************************************************************
/// Test PlanetElement class with splined elements for planets
bool test_planet_element(const PlanetElement& pe, bool verbose)
{
    // Test vs. expected state vectors of Earth @ 58000
    StateVector s0 = expected_earth_58000();

    // Calculate interpolated position and state vector of Earth
    StateVector s1 = pe.interp_vec(body_id_earth, mjd_test);

    // Calculate interpolated position of Earth
    Position q1 = pe.interp_pos(body_id_earth, mjd_test);

    // Additional information if requested
    if (verbose)
    {
        // Extract position from s0 and s1
        Position q0 = sv2pos(s0);
        Position q2 = sv2pos(s1);

        // Print position from expected and interp_pos
        ks::print_position_headers(     "Component  :");
        ks::print_position(q0,          "expected   :");
        ks::print_position(q1,          "interp_pos :");
        ks::print_position_sci(q1-q0,   "difference :");

        // Print position from expected and interp_pos
        ks::print_position_headers(     "Component  :");
        ks::print_position(q0,          "expected   :");
        ks::print_position(q2,          "interp_vec :");
        ks::print_position_sci(q2-q0,   "difference :");

        // Print position from interp_pos and interp_vec
        ks::print_position_headers(     "Component  :");
        ks::print_position(q1,          "interp_pos :");
        ks::print_position(q2,          "interp_vec :");
        ks::print_position_sci(q2-q1,   "difference :");

        // Print detail of splined orbital elements
        print("Interpolated orbital elements for Earth @ MJD {:8.2f}.\n", mjd_test);
        OrbitalElement elt = pe.interp_elt(body_id_earth, mjd_test);
        ks::print_orbital_element_headers();
        ks::print_orbital_element(elt);
    }

    // Report test results
    return test_vectors(s0, s1, q1, "PlanetElement", tol_dq_elt, tol_dv_elt, verbose);
}

// *****************************************************************************
bool test_planet_element_vs_vector(
    const PlanetVector& pv, const PlanetElement& pe, double mjd0, double mjd1, int dt_min, 
    double tol_dq, double tol_dv, bool verbose)
{
    // Total error
    double dq_tot=0.0, dv_tot=0.0;
    // Number of observations
    int obs_count=0;
    // Largest error in position and velocity
    double dq_max=0.0, dv_max=0.0;
    // Date and body_id of largest error in position and velocity
    int dq_argmax_body_id=0, dv_argmax_body_id=0;
    double dq_argmax_mjd=0.0, dv_argmax_mjd=0.0;

    // Iterate over test times
    double dt = dt_min*dpm;
    for(double mjd=mjd0; mjd<=mjd1; mjd+= dt)
    {
        // Iterate over bodies in PlanetElement collection
        for (int32_t body_id: body_ids_planets)
        {
            // Interpolated state vectors according to PlanetVector
            StateVector s0 = pv.interp_vec(body_id, mjd);
            // Interpolated state vectors according to PlanetElement
            StateVector s1 = pe.interp_vec(body_id, mjd);
            // Difference in position and velocity for this body and time
            double dq = dist(sv2pos(s0), sv2pos(s1));
            double dv = dist(sv2vel(s0), sv2vel(s1));
            // Update total and count
            dq_tot += dq;
            dv_tot += dv;
            obs_count++;
            // Update max and min for dq
            if (dq>dq_max) 
            {
                dq_max = dq;
                dq_argmax_body_id = body_id;
                dq_argmax_mjd = mjd;
            }
            // Update max and min for dv
            if (dv>dv_max) 
            {
                dv_max = dv;
                dv_argmax_body_id = body_id;
                dv_argmax_mjd = mjd;
            }
        }   // for / body_id
    }   // for / mjd
    // Calculate mean error
    double dq_mean = dq_tot / obs_count;
    double dv_mean = dv_tot / obs_count;
    // Body name with largest errors
    const string dq_argmax_body_name = get_body_name(dq_argmax_body_id);
    const string dv_argmax_body_name = get_body_name(dv_argmax_body_id);

    // Report results
    print_stars(true);
    print("Compare state vectors over date range {:8.2f}-{:8.2f}, dt_min={:d}.\n", mjd0, mjd1, dt_min);
    print("Source 0: PlanetVector  (dt_min={:d}).\n", pv.dt_min);
    print("Source 1: PlanetElement (dt_min={:d}).\n", pe.dt_min);
    print("{:10s}: {:8s}: {:8s}        ({:10s} @ {:8s})\n", "Error", "Mean", "Max", "Body Name", "mjd");
    print("{:10s}: {:8.2e}: {:8.2e} AU     ({:10s} @ {:8.2f})\n", "Position", 
        dq_mean, dq_max, dq_argmax_body_name, dq_argmax_mjd);
    print("{:10s}: {:8.2e}: {:8.2e} AU/day ({:10s} @ {:8.2f})\n", "Velocity", 
        dv_mean, dv_max, dv_argmax_body_name, dv_argmax_mjd);
    // Test results
    return (dq_max<tol_dq) && (dv_max<tol_dv);
}

// *****************************************************************************
bool test_vectors(StateVector& s0, StateVector& s1, Position& q1, string class_name,
                  double tol_dq, double tol_dv, bool verbose)
{
    // Extract expected position and velocity objects from expected StateVector s0
    Position q0 = sv2pos(s0);
    Velocity v0 = sv2vel(s0);

    // Current test result
    bool is_ok;
    // Overall test result
    bool is_ok_all = true;
    // Current test name
    string test_name;

    // Report calculated state vectors and difference
    StateVector ds = s1-s0;
    if (verbose)
    {
        print("\nExpected & Splined state vectors at {:8.4f}.\n", mjd_test);
        print_state_vector_headers( "Component  :");
        print_state_vector(s0,      "Expected   :");
        print_state_vector(s1,      "Splined    :");
        print_state_vector_sci(ds,  "Difference :");
    }

    // Check consistency of position between state vectors and position
    double dq_ip = dist(s1, q1);
    is_ok = is_close(s1, q1, tol_dq_ip);
    is_ok_all = is_ok_all && is_ok;
    if (verbose) {print("\nDistance between position from interp_pos and interp_vec\n{:5.2e} AU.\n", dq_ip);}
    test_name = format("Test {:s}: interp_pos consistent with interp_vec", class_name);
    report_test(test_name, is_ok);

    // Calculate norm of position and velocity difference
    double dq = dist(q0, s1);
    double dv = dist(v0, s1);
    // Relative difference
    double dq_rel = dq / norm(q0);
    double dv_rel = dv / norm(v0);
    // Test results
    is_ok = is_close(s0, s1, tol_dq, tol_dv);
    is_ok_all = is_ok_all && is_ok;

    // Report results
    print("\nDistance between interpolated state vectors and DB values for Earth @ {:9.4f}.\n", mjd_test);
    print("dq: {:8.2e} AU       dq_rel: {:5.3e}.\n", dq, dq_rel);
    print("dv: {:8.2e} AU/day   dv_rel: {:5.3e}.\n", dv, dv_rel);
    test_name = format("Test {:s}: interp_vec matches expected value", class_name);
    report_test(test_name, is_ok);

    // Return overall test result
    return is_ok_all;    
}
