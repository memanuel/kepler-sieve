/** @file   candidate_elt_test.cpp
 *  @brief  Test harness for CandidateElement class.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-06
 * 
 * Example call:
 * ./candidate_elt_test.x
 */

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "utils.hpp"
    using ks::norm;
    using ks::report_test;
#include "StateVector.hpp"
    using ks::StateVector;
    using ks::Position;
    using ks::Velocity;
    using ks::print_state_vector_headers;
    using ks::print_state_vector;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::print_orbital_element_headers;
    using ks::print_orbital_element;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
    using ks::print_detection;
#include "PlanetElement.hpp"
    using ks::PlanetElement;    
#include "CandidateElement.hpp"
    using ks::CandidateElement;

// *****************************************************************************
// Constants in this module
namespace{

// Date range for testing
constexpr double mjd0 = 58000.0;
constexpr double mjd1 = 59000.0;
constexpr int dt_min = 5;

// Set candidate_id to match body_id
constexpr int32_t candidate_id = 1000003;

// Test orbital elements for Juno @ 58000
// Copy / paste from KS.GetAsteroidElements(3, 4, 58000, 58000);
constexpr OrbitalElement elt0
{
    .mjd    =  58000.0,
    .a      =  2.6685312251581927,
    .e      =  0.25685345626072426,
    .inc    =  0.22673642125828372,
    .Omega  =  2.9645865407453207,
    .omega  =- 1.951164916093761,
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

// Expected state vector of Juno @58000.  
// Copy / paste from KS.GetAsteroidVectors(3, 4, 58000, 58000);
constexpr StateVector s0
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
constexpr StateVector s1
{
    .qx = -2.901474180996998,
    .qy = -1.1922326535693424,
    .qz =  0.39014258970008053,
    .vx =  0.001943496924910133,
    .vy = -0.008331238903688531,
    .vz =  0.0018120636681717218
};

// *****************************************************************************
}   // anonymous namespace

// *****************************************************************************
// Functions defined in this module
int main(int argc, char* argv[]);
bool test_all(db_conn_type& conn);
bool test_calc_traj(CandidateElement& ce, StateVector s, int i);

// *****************************************************************************
int main(int argc, char* argv[])
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run all the tests
    bool is_ok = test_all(conn);

    // Close DB connection
    conn->close();

    // Normal program exit; return 0 for success, 1 for failure
    return is_ok ? 0 : 1;
}

// *****************************************************************************
bool test_all(db_conn_type& conn)
{
    // Inputs used in testing
    int d0 = 0;
    int d1 = 1000000;
 
    // Build BodyVectors for Sun and Earth
    BodyVector bv_sun = BodyVector("Sun");
    print("Built BodyVector for Sun.\n");
    BodyVector bv_earth = BodyVector("Earth");
    print("Built BodyVector for Earth.\n");

    print("\nOrbitalElement for Juno @ {:8.2f} and {:8.2f}.\n", mjd0, mjd1);
    print_orbital_element_headers();
    print_orbital_element(elt0);
    print_orbital_element(elt1);

    // Current test result
    bool is_ok;
    // Overall test result
    bool is_ok_all = true;

    // Timer object
    Timer t;

    // Initialize DetectionTable
    t.tick();
    DetectionTable dt = DetectionTable(d0, d1);
    dt.load();
    print("Loaded DetectionTable with detection_id in [{:d}, {:d}).\n", d0, d1);
    t.tock_msg();
    // print_detection(dt[10]);

    t.tick();

    // Expected state vectors
    string pfx_header = format("{:8s}", "Date");
    string pfx_row_0 = format("{:8.2f}", mjd0);
    string pfx_row_1 = format("{:8.2f}", mjd0);
    print("Expected state vectors:\n");
    print_state_vector_headers(pfx_header);
    print_state_vector(s0, pfx_row_0);
    print_state_vector(s1, pfx_row_1);

    // Create a test array
    constexpr int N_t = 2;
    constexpr double mjd[N_t] {mjd0, mjd1};

    // Build CandidateElement for Juno elements at mjd0
    CandidateElement ce(elt0, candidate_id, mjd, N_t);
    print("\nConstructed CandidateElement from OrbitalElement for Juno @ mjd {:8.2f}.\n", mjd0);

    // Ad hoc check trajectory
    ce.calc_trajectory(false);

    // Predicted state vector
    const StateVector s0_pred = ce.state_vector(0);
    const StateVector s1_pred = ce.state_vector(1);
    // Print the predicted state vectors
    print_state_vector_headers(pfx_header);
    print_state_vector(s0_pred, pfx_row_0);
    print_state_vector(s1_pred, pfx_row_1);

    // // Test the calculated trajectory
    // is_ok = test_calc_traj(ce0, s1, 1);
    // is_ok_all &= is_ok;
    // print("Calculated asteroid trajectory.\n");
    // t.tock_msg();

    // Return overall test result
    return is_ok_all;
}

// *****************************************************************************
bool test_calc_traj(CandidateElement& ce, StateVector s0, int i)
{
    // Extracted time array and test time
    const double* mjd = ce.get_mjd();
    const double mjd_elt = ce.elt.mjd;
    const double mjd_vec = mjd[i];

    // Status
    print("Orbital elements for Juno @ {:8.2f}.\n", mjd_elt);
    print("State vectors for Juno    @ {:8.2f}.\n", mjd_vec);

    // Wrap expected position and velocity objects
    Position pos0 = sv2pos(s0);
    Velocity vel0 = sv2vel(s0);

    // Calculate the trajectory
    ce.calc_trajectory();
    print("Calculated trajectory of CandidateElement.\n");

    // Predicted state vector
    const StateVector s1 = ce.state_vector(i);
    // Extract predicted position and velocity objects
    Position pos1 = sv2pos(s1);
    Velocity vel1 = sv2vel(s1);

    // Calculate norm of position and velocity difference
    double dq = dist(pos0, pos1);
    double dv = dist(vel0, vel1);

    // Test results
    double tol_dq = 1.0E-3;
    double tol_dv = 1.0E-5;
    bool is_ok = (dq < tol_dq) && (dv < tol_dv);

    // Report results
    print("Distance between predicted state vectors in Kepler model and DB values for Juno.\n");
    print("dq: {:8.2e} AU\n", dq);
    print("dv: {:8.2e} AU/day\n", dv);
    report_test("Juno trajectory (uncalibrated)", is_ok);
    return is_ok;
}

