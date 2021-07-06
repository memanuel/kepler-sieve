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
// Tests: DetectionTable and AsteroidElement
// *****************************************************************************

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;
#include <gsl/gsl_spline.h>

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::sql_stmt_type;
    using ks::sql_prepared_stmt_type;
    using ks::get_db_conn;
    using ks::sp_run;
    using ks::sp_run_int;

#include "utils.hpp"
    using ks::is_close_abs;
    using ks::print_stars;
    using ks::print_newline;
    using ks::report_test;
    using ks::time2hms;

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

#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::norm;

#include "AsteroidElement.hpp"
    using ks::AsteroidElement;

// *****************************************************************************
// Declare functions
// void search_candidates(
//     DetectionTable& dt, AsteroidSkyPatchTable& aspt, 
//     SkyPatchNeighbor& spn, vector<DetectionNearAsteroid>& cv);
// void calculate_directions(
//     DetectionTable& dt, vector<DetectionNearAsteroid>& dv, int k0, int k1);
// void write_detections_db(db_conn_type& conn, const vector<DetectionNearAsteroid>& cv, int k0, int k1);

// Test functions
void test_detection_table(DetectionTable& dt, int detection_id);
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id);
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt);
void test_all();
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);


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
    Position pos0 = Position 
    {
        .qx = qx,
        .qy = qy,
        .qz = qz
    };
    StateVector vec0 = StateVector
    {
        .qx = qx,
        .qy = qy,
        .qz = qz,
        .vx = vx,
        .vy = vy,
        .vz = vz
    };

    // Tolerance for tests
    double tol_q = 1.0E-8;
    double tol_vec = 1.0E-8;

    // Calulate splined position and state vector
    Position pos = bv.interp_pos(mjd);
    StateVector vec = bv.interp_vec(mjd);

    // Report splined orbital elements
    print("\nSplined Sun state vectors at mjd {:8.2f} :\n", mjd);
    print("qx:         {:+9.6f}\n", vec.qx);
    print("qy:         {:+9.6f}\n", vec.qy);
    print("qz:         {:+9.6f}\n", vec.qz);
    print("vx:         {:+9.6f}\n", vec.vx);
    print("vy:         {:+9.6f}\n", vec.vy);
    print("vz:         {:+9.6f}\n", vec.vz);

    // Test that splined position matches expected results
    {
    bool is_ok = norm(pos0, pos) < tol_q;
    report_test("\nTest BodyVector::interp_pos() splined position matches database", is_ok);
    }

    // Test that splined state vectors match expected results
    {
    bool is_ok = norm(vec0, vec) < tol_vec;
    report_test("\nTest BodyVector::interp_vec() splined state vectors match database", is_ok);
    }
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
    double qx_ast = -2.485854955060206;
    double qy_ast = -0.6229239033462274;
    double qz_ast =  0.4383572489736212;
    double vx_ast =  0.0020617208493692073;
    double vy_ast = -0.010756337836425017;
    double vz_ast = -0.0007200628373353;
    // State vectors of the sun; copy / pasted from KS.GetStateVectors_Sun(58400, 58400, 1);
    double qx_sun = -0.0001095835967748;
    double qy_sun =  0.007235858951602957;
    double qz_sun = -0.0000736284237584;
    double vx_sun = -0.0000075730407095;
    double vy_sun =  0.0000026357733813;
    double vz_sun =  0.0000001892676823;
    // Expected results; heliocentric position of asteroid
    double qx = qx_ast - qx_sun;
    double qy = qy_ast - qy_sun;
    double qz = qz_ast - qz_sun;
    double vx = vx_ast - vx_sun;
    double vy = vy_ast - vy_sun;
    double vz = vz_ast - vz_sun;
    // Wrap expected results into Position and Vector objects
    Position pos0 = Position 
    {
        .qx = qx,
        .qy = qy,
        .qz = qz
    };
    StateVector vec0 = StateVector
    {
        .qx = qx,
        .qy = qy,
        .qz = qz,
        .vx = vx,
        .vy = vy,
        .vz = vz
    };

    // Tolerance for tests
    double tol_q = 1.0E-8;
    double tol_vec = 1.0E-8;
    // double tol_v = 1.0E-10;

    // Calulate splined position and state vector
    Position pos = ast_elt.interp_pos(asteroid_id, mjd);
    StateVector vec = ast_elt.interp_vec(asteroid_id, mjd);

    // Report splined orbital elements
    print("\nSplined Asteroid index {:d} at time index {:d} :\n", asteroid_idx, time_idx);
    print("AsteroidID: {:9d}\n", asteroid_id);
    print("mjd:        {:9.2f}\n", mjd);
    print("qx:         {:+9.6f}\n", vec.qx);
    print("qy:         {:+9.6f}\n", vec.qy);
    print("qz:         {:+9.6f}\n", vec.qz);
    print("vx:         {:+9.6f}\n", vec.vx);
    print("vy:         {:+9.6f}\n", vec.vy);
    print("vz:         {:+9.6f}\n", vec.vz);

    // Test that splined orbital elements match expected results (position only)
    {
    bool is_ok = norm(pos0, pos) < tol_q;
    report_test("\nTest AsteroidElement::interp_pos() splined state vectors match database", is_ok);
    }

    // Test that splined orbital elements match expected results
    {
    bool is_ok = norm(vec0, vec) < tol_vec;
    report_test("\nTest AsteroidElement::interp_vec() splined state vectors match database", is_ok);
    }
}    

// *****************************************************************************
void test_all()
{
    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors...\n");
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table.\n");

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Inputs to build DetectionCandidateTable, AsteroidSkypatchTable and AsteroidElement
    int d0 = 0;
    int d1 = 100;
    int n0 = 0;
    int n1 = 100;
    int mjd0 = 58000;
    int mjd1 = 61000;
    int time_step = 4;
    bool progbar = true;

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1, progbar);

    // Test DetectionTable
    test_detection_table(dt, d0);

    // Build AsteroidSkyPatch table
    print_newline();
    AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, n0, n1, progbar);

    // Rebuild BodyVectors from DB and save to disk
    // save_vectors("Sun");
    // save_vectors("Earth");

    // Initialize BodyVector for the Sun from disk
    BodyVector bv("Sun");

    // Test body vectors for Sun
    test_body_vector(bv);

    // Initialize AsteroidElement
    AsteroidElement elt(n0, n1, mjd0, mjd1, time_step);
    elt.load(conn, progbar);

    // Test asteroid elements - splining orbital elements
    test_asteroid_element(elt);

    // Test asteroid elements - state vectors from splined elements
    test_asteroid_element_vectors(elt);

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
}