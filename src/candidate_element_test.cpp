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
#include "utils.hpp"
    using ks::norm;
    using ks::report_test;
    using ks::print_stars;
#include "astro_utils.hpp"
    using ks::SolarSystemBody_bv;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "StateVector.hpp"
    using ks::StateVector;
    using ks::Position;
    using ks::Velocity;
    using ks::print_state_vector_headers;
    using ks::print_state_vector;
#include "Direction.hpp"
    using ks::Direction;
    using ks::print_direction;    
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
constexpr double mjd0 {58000.0};
constexpr double mjd1 {59000.0};
constexpr int dt_min {5};

// Set candidate_id to match body_id of Juno (asteroid_id=3)
constexpr int32_t candidate_id = 1000003;

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

// Prefix for state vectors
const string pfx_header = format("{:8s}", "Date");
const string pfx_row_0 = format("{:8.2f}", mjd0);
const string pfx_row_1 = format("{:8.2f}", mjd1);

// *****************************************************************************
}   // anonymous namespace

// *****************************************************************************
// Functions defined in this module
int main(int argc, char* argv[]);
bool test_all(db_conn_type& conn);
bool test_calc_traj(bool is_calibrated, bool verbose);
bool test_calc_direction(bool verbose);

// *****************************************************************************
int main(int argc, char* argv[])
{
    // Start timer
    Timer t;
    t.tick();

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run all the tests
    bool is_ok = test_all(conn);

    // Close DB connection
    conn->close();

    // Report total elapsed time
    print_stars();
    t.tock_msg("\nTOTAL - candidate_element_test\nElapsed Time");

    // Normal program exit; return 0 for success, 1 for failure
    return is_ok ? 0 : 1;
}

// *****************************************************************************
bool test_all(db_conn_type& conn)
{
    // Inputs used in testing
    // int d0 = 0;
    // int d1 = 1000000;
 
    // Overall test result
    bool is_ok = true;

    // Get reference elements
    print_stars(true);
    print("OrbitalElement for Juno @ {:8.2f} and {:8.2f} (expected results).\n", mjd0, mjd1);
    print_orbital_element_headers();
    print_orbital_element(elt0);
    print_orbital_element(elt1);

    // Timer object
    Timer t;

    // Test the calculated trajectory - uncalibrated
    {
        bool calibrated = false;
        bool verbose = false;
        is_ok &= test_calc_traj(calibrated, verbose);
    }

    // Test the calculated trajectory - calibrated
    {
        bool calibrated = true;
        bool verbose = false;
        is_ok &= test_calc_traj(calibrated, verbose);
    }

    // Test the calculated direction
    {
        bool verbose = true;
        is_ok &= test_calc_direction(verbose);
    }

    // Return overall test result
    return is_ok;
}

// *****************************************************************************
bool test_calc_traj(bool is_calibrated, bool verbose)
{
    // Status - name of test being run
    string cal_type = is_calibrated ? "calibrated" : "uncalibrated";
    print_stars(true);
    print("test_calc_traj - {:s}; test calculated trajectory vs. expected state vectors.\n", cal_type);

    // Expected state vectors
    if (verbose)
    {
        print("\nExpected state vectors:\n");
        print_state_vector_headers(pfx_header);
        print_state_vector(s0, pfx_row_0);
        print_state_vector(s1, pfx_row_1);
    }

    // Create a test array
    constexpr int N_t = 2;
    constexpr double mjd[N_t] {mjd0, mjd1};

    // Build CandidateElement for Juno elements at mjd0
    CandidateElement ce = CandidateElement(elt0, candidate_id, mjd, N_t);
    if (verbose)
    {
        print("\nConstructed CandidateElement from OrbitalElement for Juno @ mjd {:8.2f}.\n", mjd0);
        print("candidate_id={:d}, mjd[0]={:.1f}, mjd[1]={:.1f}, N_t={:d}.\n", candidate_id, mjd[0], mjd[1], N_t);
    }

    Timer t;
    t.tick();

    // Calibrate if requested
    if (is_calibrated) {ce.calibrate();}

    // Calculate trajectory
    ce.calc_trajectory();
    if (verbose) {print("Calculated trajectory of CandidateElement.\n");}

    // Predicted state vector
    const StateVector s0_pred = ce.state_vector(0);
    const StateVector s1_pred = ce.state_vector(1);
    // Print the predicted state vectors
    if (verbose)
    {
        print("\nPredicted state vectors with Kepler model:\n");
        print_state_vector_headers(pfx_header);
        print_state_vector(s0_pred, pfx_row_0);
        print_state_vector(s1_pred, pfx_row_1);
    }

    // Calculate norm of position and velocity difference
    double dq0 = dist_dq(s0, s0_pred);
    double dv0 = dist_dv(s0, s0_pred);
    double dq1 = dist_dq(s1, s1_pred);
    double dv1 = dist_dv(s1, s1_pred);

    // Set tolerance for the recovery of elements on the nodes
    double tol_dq_node = 1.0E-13;
    double tol_dv_node = 1.0E-15;
    // Set tolerance for the Kepler projection 1000 days away
    double tol_dq = 1.0E-2;
    double tol_dv = 1.0E-4;

    // Tighten up tolerances if calibrated to simulation
    if (is_calibrated)
    {
        tol_dq = 1.0E-10;
        tol_dv = 1.0E-12;
    }

    // Test results
    bool is_ok = (dq0 < tol_dq_node) && (dv0 < tol_dv_node) && (dq1 < tol_dq) && (dv1 < tol_dv);

    // Report results
    string cal_des = is_calibrated ? "calibrated" : "uncalibrated";
    print("\nDistance between predicted vs. expected state vectors for Juno with {:s} Kepler model.\n", cal_des);
    print("{:6s}:  {:8.2f}:  {:8.2f}\n",         "Date", mjd0,   mjd1);
    print("{:6s}:  {:8.2e}:  {:8.2e} AU\n",      "dq",   dq0,    dq1);
    print("{:6s}:  {:8.2e}:  {:8.2e} AU/day\n",  "dv",   dv0,    dv1);
    string test_name = format("Juno trajectory ({:s})", cal_des);
    report_test(test_name, is_ok);

    // Report the elapsed time
    t.tock_msg();
    return is_ok;
}

// *****************************************************************************
bool test_calc_direction(bool verbose)
{
    // Status - name of test being run
    print_stars(true);
    print("test_calc_direction - test calculated direction vs. detection on a hit for Juno.\n");
    
    Timer t;
    t.tick();

    // Build CandidateElement for Juno elements at mjd0
    CandidateElement ce = CandidateElement(elt0, candidate_id);
    // if (verbose) {print("Constructed CandidateElement from OrbitalElement for Juno @ mjd {:8.2f}.\n", mjd0);}

    // Calibrate 
    ce.calibrate();
    // if (verbose) {print("Calibrated orbit of CandidateElement.\n");}

    // Calculate trajectory
    ce.calc_trajectory();
    // if (verbose) {print("Calculated trajectory of CandidateElement.\n");}

    // Calculate direction
    ce.calc_direction();
    if (verbose) {print("Calculated directions on CandidateElement.\n");}

    // First hit on Juno is on detection_id 85817184
    int32_t detection_id {85817184};
    // The detection object
    Detection det {ce.dt[detection_id]};
    // print("Detection with detection_id {:d}:\n", detection_id);
    print_detection(det);
    // Extract the direction from this detection
    Direction u_det = det2dir(det);

    // The associated detection_time_id of this detection; corresponds to mjd 58940.4147454002
    int32_t detection_time_id = det.detection_time_id;
    // DetectionTime dt {ce.dtt[detection_time_id]};
    // print("\nDetectionTime:\n");
    // print("detection_time_id = {:d}, mjd = {:8.6f}.\n", dt.detection_time_id, dt.mjd);

    // Extract the direction on the test date
    int i = detection_time_id-1;
    Direction u_ast = ce.direction(i);
    double mjd_i = ce.mjd[i];
    print("Calculated direction at i={:d}, mjd = {:8.6f} vs. detection.\n", i, mjd_i);
    print_direction(u_ast);
    print_direction(u_det);

    // Test result
    // double du_rad = dist_rad(u_det, u_ast);
    // Distance between detection and predicted direction in arc seconds
    double du_sec = dist_sec(u_det, u_ast);
    // Tolerance in arc seconds
    double tol = 10.0;
    // Report
    bool is_ok = (du_sec < tol);
    print("Difference in directions: {:8.2f} arc seconds.\n", du_sec);
    report_test("test_calc_direction", is_ok);
    // Report the elapsed time
    t.tock_msg();
    return is_ok;

}