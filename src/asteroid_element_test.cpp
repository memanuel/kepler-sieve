/** @file   asteroid_element_test.cpp
 *  @brief  Test harness for class AsteroidElement for splined elements of asteroids.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-16
 * 
 * Example call:
 * ./asteroid_element_test.x
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
    using ks::SolarSystemBody_bv;
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
    using ks::norm;
    using ks::dist_dq;
    using ks::dist_dv;
    using ks::is_close;
    using ks::print_state_vector_headers;
    using ks::print_state_vector;
    using ks::print_state_vector_sci;    
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::print_orbital_element;
#include "AsteroidElement.hpp"
    using ks::AsteroidElement;

// *****************************************************************************
// Constants in this module
namespace{

// Date range for testing
constexpr double mjd0 = 58000.0;
constexpr double mjd1 = 59000.0;
constexpr int dt = 4;

// Test is on Juno with asteroid_id 3
constexpr int32_t asteroid_id = 3;

/// Absolute difference between interp_pos and interp_vec outputs (position part only)
constexpr double tol_dq_ip = 1.0E-14;
/// Absolute difference between interp_vec output and expected position, on node dates
constexpr double tol_dq_on_node = 1.0E-14;
/// Absolute difference between interp_vec output and expected velocity, on nodes dates
constexpr double tol_dv_on_node = 1.0E-16;

/// Absolute difference between interp_vec output and expected position, splined off node dates
constexpr double tol_dq_off_node = 1.0E-7;
/// Absolute difference between interp_vec output and expected velocity, splined off node dates
constexpr double tol_dv_off_node = 1.0E-9;

// Expected state vector of Juno @58000.  
// Copy / paste from KS.GetAsteroidVectors(3, 4, 58000, 58000);
constexpr StateVector s0_good
{
    .qx =  1.0693547365201785,
    .qy = -2.684939245391761,
    .qz =  0.5674675777224312,
    .vx =  0.007764282851412018,
    .vy =  0.005217549084953882,
    .vz = -0.001498976011266847
};

// Expected state vector of Juno @ 59000.  
// Copy / paste from KS.GetAsteroidVectors(3, 4, 59000, 59000);
constexpr StateVector s1_good
{
    .qx = -2.901474180996998,
    .qy = -1.1922326535693424,
    .qz =  0.39014258970008053,
    .vx =  0.001943496924910133,
    .vy = -0.008331238903688531,
    .vz =  0.0018120636681717218
};

// Test orbital elements for Juno @ 58000
// Copy / paste from KS.GetAsteroidElements(3, 4, 58000, 58000);
constexpr OrbitalElement elt0
{
    .mjd    =  58000.0,
    .a      =  2.6685312251581927,
    .e      =  0.25685345626072426,
    .inc    =  0.22671752380617316,
    .Omega  =  2.9645865407453207,
    .omega  = -1.951164916093761,
    .f      = 48.04726436026028,
    .M      = 42.2236141353505
};

// Test orbital elements for Juno @ 59000
// Copy / paste from KS.GetAsteroidElements(3, 4, 59000, 59000);
constexpr OrbitalElement elt1
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

// Prefix for state vectors and orbital elements
const string pfx_header = format("{:8s}", "Date");
const string pfx_row_0 = format("{:8.2f}", mjd0);
const string pfx_row_1 = format("{:8.2f}", mjd1);

// *****************************************************************************
}   // anonymous namespace

// *****************************************************************************
// Functions defined in this module
int main(int argc, char* argv[]);
bool test_all(db_conn_type& conn, bool verbose);
bool test_asteroid_element(const AsteroidElement& ae, bool verbose);
bool test_vectors(const StateVector& s_good, const StateVector& s_pred, const Position& p_pred, 
                  const string test_name, double tol_dq, double tol_dv, bool verbose);

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
    constexpr int pad = 32;
    constexpr int mjd0_ae = mjd0 - pad;
    constexpr int mjd1_ae = mjd1 + pad;    

    // Results of current test
    bool is_ok = true;
    // Name of current test
    string test_name = "";
    // Overall test results
    bool is_ok_all = true;

    // *****************************************************************************
    // Compare AsteroidElement vs. known state vectors on node dates
    // *****************************************************************************
    {
    // Build a spline of orbital elements sampled at frequency dt specified above
    AsteroidElement ae(asteroid_id, asteroid_id+1, mjd0_ae, mjd1_ae, dt);
    // ae.load(conn, false);

    // // Predicted state vectors on mjd0 and mjd1
    // StateVector s0_pred = ae.interp_vec(asteroid_id, mjd0);
    // // StateVector s1_pred = ae.interp_vec(asteroid_id, mjd1);
    // // Predicted position on mjd0 and mjd1
    // Position p0_pred = ae.interp_pos(asteroid_id, mjd0);
    // // Position p1_pred = ae.interp_pos(asteroid_id, mjd1);

    // // Set tolerances for on-node test
    // double tol_dq = tol_dq_on_node;
    // double tol_dv = tol_dv_on_node;

    // // Compare vectors on mjd0
    // test_name = format("Test AsteroidElement on Juno @ {:8.2f} (on node).\n", mjd0);
    // is_ok = test_vectors(s0_good, s0_pred, p0_pred, test_name, tol_dq, tol_dv, verbose);
    // is_ok_all = is_ok_all && is_ok;
    }

    // *****************************************************************************
    // Overall results
    // *****************************************************************************

    // DEBUG
    is_ok_all &= is_ok;
    // Return overall test result
    return is_ok_all;
}

// *****************************************************************************
/// Test vectors predicted by AsteroidElement class vs. expected state vectors from database
bool test_vectors(const StateVector& s_good, const StateVector& s_pred, const Position& p_pred, 
                  const string test_name, double tol_dq, double tol_dv, bool verbose)
{
    // Extract position from s_good and s_pred
    Position q_good = sv2pos(s_good);
    Position q_pred = sv2pos(s_pred);
    // Extract velocity from s_good and s_pred
    Velocity v_good = sv2vel(s_good);
    Velocity v_pred = sv2vel(s_pred);
    // Difference between two interpolation methods for position
    Position dq_ip = q_pred - q_good;
    // Difference between interpolated and expected position
    Position dq = q_pred - q_good;
    // Difference between interpolated and expected velocity
    Velocity dv = v_pred - v_good;

    // Additional information if requested
    if (verbose)
    {
        // Print predicted position with two methods
        ks::print_position_headers(     "Component  :");
        ks::print_position(p_pred,      "interp_pos :");
        ks::print_position(q_pred,      "interp_vec :");
        ks::print_position_sci(dq_ip,   "difference :");

        // Print expected and predicted position
        ks::print_position_headers(     "Component  :");
        ks::print_position(q_good,      "expected   :");
        ks::print_position(q_pred,      "interp_vec :");
        ks::print_position_sci(dq,      "difference :");

        // Print expected and predicted velocity
        ks::print_position_headers(     "Component  :");
        ks::print_position(q_good,      "expected   :");
        ks::print_position(q_pred,      "interp_vec :");
        ks::print_position_sci(dq,      "difference :");
    }

    // Calculate distances for the three tests
    double dist_dq_ip = norm(dq_ip);
    double dist_dq = norm(dq);
    double dist_dv = norm(dv);
    // Relative difference
    double dist_dq_rel = dist_dq / norm(q_good);
    double dist_dv_rel = dist_dv / norm(v_good);
    // Test results
    bool is_ok = is_close(s_good, s_pred, tol_dq, tol_dv) && (dist_dq_ip < tol_dq_ip);

    // Report results
    print("\nDistance between state vectors - {:s}.\n", test_name);
    print("dq: {:8.2e} AU       dq_rel: {:5.3e}.\n", dist_dq, dist_dq_rel);
    print("dv: {:8.2e} AU/day   dv_rel: {:5.3e}.\n", dist_dv, dist_dv_rel);
    report_test(test_name, is_ok);
    return is_ok;
}
