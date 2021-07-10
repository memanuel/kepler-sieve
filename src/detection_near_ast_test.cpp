/** @file detection_near_ast_test.cpp
 *  @brief Test harness for functions used in detection_near_ast.cpp.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-05
 * 
 * Example call:
 * ./detection_near_ast_test.x
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;
#include <gsl/gsl_spline.h>

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
#include "utils.hpp"
    using ks::is_close_abs;
    using ks::print_stars;
    using ks::print_newline;
    using ks::report_test;
    using ks::time2hms;
#include "StateVector.hpp"
    using ks::norm;
    using ks::dist;
    using ks::sv2pos;
    using ks::sv2vel;
    using ks::print_state_vector_headers;
    using ks::print_state_vector;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
#include "SkyPatchNeighbor.hpp"
    using ks::SkyPatch;
    using ks::SkyPatchNeighbor;
    using ks::sky_patch::N_sp;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
#include "AsteroidSkyPatch.hpp"
    using ks::AsteroidSkyPatch;
    using ks::AsteroidSkyPatchTable;
#include "BodyVector.hpp"
    using ks::BodyVector;
    using ks::save_vectors;
#include "AsteroidElement.hpp"
    using ks::AsteroidElement;

// *****************************************************************************
// Declare functions defined in this module
int main(int argc, char* argv[]);
bool test_all(db_conn_type& conn);
void test_detection_table(DetectionTable& dt, int detection_id);
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id);
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt);
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);
void test_body_vector(BodyVector& bv);
void test_asteroid_element(AsteroidElement& ast_elt);
void test_asteroid_element_vectors(AsteroidElement& ast_elt);

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
    // Inputs to build DetectionCandidateTable, AsteroidSkypatchTable and AsteroidElement
    int d0 = 0;
    int d1 = 100;
    int n0 = 0;
    int n1 = 100;
    int mjd0 = 58000;
    int mjd1 = 61000;
    int time_step = 4;
    bool progbar = true;

    // Current test result
    // bool is_ok;
    // Overall test result
    // TODO calculate this for real
    bool is_ok_all = true;

    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors...\n");
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table.\n");

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(d0, d1);
    dt.load();

    // Test DetectionTable
    test_detection_table(dt, d0);

    // Build AsteroidSkyPatch table
    print_newline();
    AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, n0, n1, progbar);

    // Rebuild BodyVectors from DB and save to disk
    // save_vectors("Sun");
    // save_vectors("Earth");
    // print("Saved vectors for Earth.\n");

    // Initialize BodyVector for the Sun from disk
    print("Loading BodyVector object with saved vectors for Sun.\n");
    BodyVector bv("Sun");
    print("Loaded BodyVector object with saved vectors for Sun.\n");

    // Test body vectors for Sun
    test_body_vector(bv);

    // Initialize AsteroidElement
    AsteroidElement elt(n0, n1, mjd0, mjd1, time_step);
    elt.load(conn, progbar);

    // Test asteroid elements - splining orbital elements
    test_asteroid_element(elt);

    // Test asteroid elements - state vectors from splined elements
    test_asteroid_element_vectors(elt);

    // Return overall result
    return is_ok_all;
}

// *****************************************************************************
void test_detection_table(DetectionTable& dt, int detection_id)
{
    // Print first 10 detections
    int i0=detection_id;
    int i1=std::min(dt.d1, i0+10);
    print("\nSample data: detections with IDs in [{:d}, {:d}) :\n", i0, i1);
    print("{:12s} {:10s} {:8s}  {:9s}  {:9s}  {:9s}\n", "DetectionID", "SkyPatchID", "mjd", "ux", "uy", "uz");
    for (int i=i0; i<i1;i++)
    {
        // The ith detection on the table
        Detection d = dt[i];
        // Holes in the detection vector indicated by d.detectio_id=0; skip these
        if (!d.detection_id) {continue;}
        // Print this detection
        print("{:11d} {:10d} {:8.3f} {:+9.6f} {:+9.6f} {:+9.6f}.\n", 
            d.detection_id, d.sky_patch_id, d.mjd, d.ux, d.uy, d.uz);
    }
}

// *****************************************************************************
/// Test that BodyVector class for Sun splines state vectors that match copy / pasted values for Sun @ 58400.
void test_body_vector(BodyVector& bv)
{
    // This test checks state vectors for Sun at mjd 58400
    double mjd = 58400.0;
    
    // State vectors of the sun; copy / pasted from KS.GetStateVectors_Sun(58400, 58400, 1);
    double qx = -0.0001095835967748;
    double qy =  0.007235858951602957;
    double qz = -0.0000736284237584;
    double vx = -0.0000075730407095;
    double vy =  0.0000026357733813;
    double vz =  0.0000001892676823;
    // Wrap expected results into Position and Vector objects
    Position q0 = Position {.qx=qx, .qy=qy, .qz=qz};
    Velocity v0 = Velocity {.vx=vx, .vy=vy, .vz=vz};
    // StateVector s0 = StateVector {.qx=qx, .qy=qy, .qz=qz, .vx=vx, .vy=vy, .vz=vz};

    // Tolerance for tests
    double tol_dq_ip = 1.0E-12;
    double tol_dq = 1.0E-8;
    double tol_dv = 1.0E-8;
    // Result of current test
    bool is_ok;

    // Calulate splined position and state vector
    StateVector s1 = bv.interp_vec(mjd);
    Position q1_ip = bv.interp_pos(mjd);
    // Extract position and velocity from state vector
    Position q1 = sv2pos(s1);
    Velocity v1 = sv2vel(s1);

    // Report splined orbital elements
    print("\nSplined Sun state vectors at mjd {:8.2f} :\n", mjd);
    print_state_vector_headers();
    print_state_vector(s1);

    // Check consistency of position between state vectors ans position
    double dq_ip = dist(q1, q1_ip);
    print("\nDistance between position from interp_pos and interp_vec: {:5.2e} AU.\n", dq_ip);
    is_ok = is_ok && (dq_ip < tol_dq_ip);

    // Calculate norm of position and velocity difference
    double dq = dist(q0, q1);
    double dv = dist(v0, v1);
    // Relative difference
    double dq_rel = dq / norm(q0);
    double dv_rel = dv / norm(v0);
    // Test results
    is_ok = is_ok && (dq < tol_dq) && (dv < tol_dv);

    print("\nDistance between interpolated state vectors and DB values for Sun @ {:9.4f}.\n", mjd);
    print("dq: {:8.2e} AU       dq_rel: {:5.3e}.\n", dq, dq_rel);
    print("dv: {:8.2e} AU/day   dv_rel: {:5.3e}.\n", dv, dv_rel);

    // Test that splined position matches expected results
    is_ok = dist(q0, q1) < tol_dq;
    report_test("\nTest BodyVector::interp_pos() splined position matches database", is_ok);

    // Test that splined state vectors match expected results
    is_ok = dist(v0, v1) < tol_dv;
    report_test("\nTest BodyVector::interp_vec() splined state vectors match database", is_ok);
}

// *****************************************************************************
/// Test that AsteroidElement instance splines orbital elements that match copy / pasted values for Ceres @ 58400.
void test_asteroid_element(AsteroidElement& ast_elt)
{
    // This test is designed to check Ceres (asteroid_id=1) at mjd 58400 (time_idx=100)
    int asteroid_idx = 1;
    int time_idx = 100;

    // Expected results, copy / pasted from database from KS.GetAsteroidElements(1, 2, 58400, 58400);
    double a0     = 2.7673528257126296;
    double e0     = 0.07561068735641437;
    double inc0   = 0.18489327222145555;
    double Omega0 = 1.4016429441591765;
    double omega0 = 1.2779179332926616;
    double f0     = 38.4030882294433;
    double M0     = 38.30932258771573;
    
    // Tolerance for tests
    double tol_a = 1.0E-14;
    double tol_e = 1.0E-14;
    double tol_angle = 1.0E-12;

    // Read off some asteroid elements
    print("\nAsteroidElement properties:\n");
    print("N_ast    : {:d}\n", ast_elt.N_ast);
    print("N_t      : {:d}\n", ast_elt.N_t);

    // Read the two 1D arrays
    int32_t* asteroid_ids = ast_elt.get_asteroid_id();
    double* mjds = ast_elt.get_mjd();
    // The selected asteroid_id and time
    int32_t asteroid_id = asteroid_ids[asteroid_idx];
    double mjd = mjds[time_idx];

    // Calulate splined orbital elements; these should match
    OrbitalElement elt = ast_elt.interp_elt(asteroid_id, mjd);

    // Report splined orbital elements
    print("\nSplined Asteroid index {:d} at time index {:d} :\n", asteroid_idx, time_idx);
    print("AsteroidID: {:9d}\n", asteroid_id);
    print("mjd:        {:9.2f}\n", mjd);
    print("a:          {:9.6f}\n", elt.a);
    print("e:          {:9.6f}\n", elt.e);
    print("inc:        {:9.6f}\n", elt.inc);
    print("Omega:      {:9.6f}\n", elt.Omega);
    print("omega:      {:9.6f}\n", elt.omega);
    print("f:          {:9.6f}\n", elt.f);
    print("M:          {:9.6f}\n", elt.M);

    // Test that splined orbital elements match expected results
    bool is_ok = true;
    is_ok = is_ok && is_close_abs(a0,     elt.a, tol_a);
    is_ok = is_ok && is_close_abs(e0,     elt.e, tol_e);
    is_ok = is_ok && is_close_abs(inc0,   elt.inc, tol_angle);
    is_ok = is_ok && is_close_abs(Omega0, elt.Omega, tol_angle);
    is_ok = is_ok && is_close_abs(omega0, elt.omega, tol_angle);
    is_ok = is_ok && is_close_abs(f0,     elt.f, tol_angle);
    is_ok = is_ok && is_close_abs(M0,     elt.M, tol_angle);
    report_test("\nTest AsteroidElement::interp_elt() splined elements match database", is_ok);
}

// *****************************************************************************
/// Test that AsteroidElement instance splines state vectors that match copy / pasted values for Ceres @ 58400.
void test_asteroid_element_vectors(AsteroidElement& ast_elt)
{
    // This test is designed to check Ceres (asteroid_id=1) at mjd 58400 (time_idx=100)
    int asteroid_idx = 1;
    int time_idx = 100;

    // Read the two 1D arrays
    int32_t* asteroid_ids = ast_elt.get_asteroid_id();
    double* mjds = ast_elt.get_mjd();
    // The selected asteroid_id and time
    int32_t asteroid_id = asteroid_ids[asteroid_idx];
    double mjd = mjds[time_idx];
    
    // Asteroid position; copy / pasted from database from KS.GetAsteroidElements(1, 2, 58400, 58400);
    double qx = -2.485854955060206;
    double qy = -0.6229239033462274;
    double qz =  0.4383572489736212;
    double vx =  0.0020617208493692073;
    double vy = -0.010756337836425017;
    double vz = -0.0007200628373353;

    // Wrap expected results into Position and Vector objects
    Position pos0 = Position {.qx=qx, .qy=qy, .qz=qz};
    StateVector vec0 = StateVector {.qx=qx, .qy=qy, .qz=qz, .vx=vx, .vy=vy, .vz=vz};

    // Tolerance for tests
    double tol_q = 1.0E-8;
    double tol_vec = 1.0E-8;
    // double tol_v = 1.0E-10;

    // Calulate splined position and state vector
    Position pos = ast_elt.interp_pos(asteroid_id, mjd);
    StateVector vec = ast_elt.interp_vec(asteroid_id, mjd);

    // Report splined orbital elements
    print("\nSplined Asteroid index {:d} at time index {:d} :\n", asteroid_idx, time_idx);
    print("AsteroidID: {:9d}\n",    asteroid_id);
    print("mjd:        {:9.2f}\n",  mjd);
    print("qx:         {:+9.6f}\n", vec.qx);
    print("qy:         {:+9.6f}\n", vec.qy);
    print("qz:         {:+9.6f}\n", vec.qz);
    print("vx:         {:+9.6f}\n", vec.vx);
    print("vy:         {:+9.6f}\n", vec.vy);
    print("vz:         {:+9.6f}\n", vec.vz);

    // Test that splined position from orbital elements match expected results
    {
    bool is_ok = dist(pos0, pos) < tol_q;
    report_test("\nTest AsteroidElement::interp_pos() splined state vectors match database", is_ok);
    }

    // Test that splined state vectors from orbital elements match expected results
    {
    bool is_ok = dist(vec0, vec) < tol_vec;
    report_test("\nTest AsteroidElement::interp_vec() splined state vectors match database", is_ok);
    }
}


