/** @file   detection_near_elt_test.cpp
 *  @brief  Test harness for functions used in finding detections near candidate orbital elements.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-06
 * 
 * Example call:
 * ./detection_near_elt_test.x
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::Position;
    using ks::Velocity;
    using ks::print_orbital_element;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "PlanetElement.hpp"
    using ks::PlanetElement;    
#include "CandidateElement.hpp"
    using ks::CandidateElement;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "utils.hpp"
    using ks::norm;
    using ks::report_test;

// *****************************************************************************
// Declare functions defined in this module
int main(int argc, char* argv[]);
void print_detection(DetectionTable& dt);
void test_calc_traj(OrbitalElement& elt, DetectionTimeTable& dtt);
void test_massive_body();
void test_load_detection(db_conn_type& conn);
void test_all();

// *****************************************************************************
void print_detection(DetectionTable& dt)
{
    Detection d = dt[10];
    print("\nExample Detection:\n");
    print("detection_id = {:d}\n", d.detection_id);
    print("sky_patch_id = {:d}\n", d.sky_patch_id);
    print("time_id = {:d}\n", d.time_id);
    print("ux      = {:+8.6f}\n", d.ux);
    print("uy      = {:+8.6f}\n", d.uy);
    print("uz      = {:+8.6f}\n", d.uz);
    print("mag     = {:8.4f}\n", d.mag);
}

// *****************************************************************************
void test_calc_traj(OrbitalElement& elt, DetectionTimeTable& dtt)
{
    // Build CandidateElement for these elements
    CandidateElement ce(elt, dtt);

    // Choose an example date that is available in DB for testing
    double mjd_test = 58000.0;
    // Overwrite first slot of mjd array in CandidateElement object
    // Need to abuse const double pointer by casting it to non-const
    double* mjd = (double*) ce.get_mjd();
    mjd[0] = mjd_test;

    // Calculate trajectory of the asteroid in Kepler approximation
    ce.calc_trajectory();

    // Get the predicted position of the asteroid on the test date
    const double* q_ast = ce.get_q_ast();
    Position pos
    {
        .qx = q_ast[0],
        .qy = q_ast[1],
        .qz = q_ast[2]
    };
    // Get the predicted velocity of the asteroid on the test date
    const double* v_ast = ce.get_v_ast();
    Velocity vel
    {
        .vx = v_ast[0],
        .vy = v_ast[1],
        .vz = v_ast[2]
    };

    // Expected state vector components - location of Juno at this time.  
    // Copy / paste from KS.GetAsteroidVectors(3, 4, 58000, 58000);
    double qx =  1.0693547365201785;
    double qy = -2.684939245391761;
    double qz =  0.5674675777224312;
    double vx =  0.007764282851412018;
    double vy =  0.005217549084953882;
    double vz = -0.001498976011266847;

    // Wrap expected position and velocity objects
    Position pos0 {.qx = qx, .qy = qy, .qz = qz};
    Velocity vel0 {.vx = vx, .vy = vy, .vz = vz};

    // Calculate norm of position and velocity difference
    double dq = norm(pos0, pos);
    double dv = norm(vel0, vel);

    // Report results
    print("Distance between predicted state vectros in Kepler model and DB values for Juno @ {:9.4f}.\n", mjd_test);
    print("dq: {:8.2e} AU\n", dq);
    print("dv: {:8.2e} AU/day\n", dv);
}

// *****************************************************************************
void test_massive_body()
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Instantiate a MassiveBodyTable
    MassiveBodyTable mbt1(false);
    // Load it from DB
    mbt1.load(conn);
    // Save it
    mbt1.save();

    // Load from disk
    MassiveBodyTable mbt = MassiveBodyTable();

    // Print contents
    print("{:8s} : {:8s} : {:8s}\n", "BodyID", "M", "GM");
    for (int32_t body_id: mbt.get_body_id())
    {
        // Only print the planets; don't show the heavy asteroids
        if (body_id > 1000000) {continue;}
        // Access this body and print it
        MassiveBody mb = mbt[body_id];
        print("{:8d} : {:8.2e} : {:8.2e}\n", mb.body_id, mb.M, mb.GM);
    }

    // Close DB connection
    conn->close();
}

// *****************************************************************************
void test_planet_element(PlanetElement& pe)
{
    // Test state vectors of Earth @ 58000
    int32_t body_id=399;
    double mjd_test = 58000.0;
    // Calculate idx from body_id
    // int idx = pe.body_idx(body_id);
    // print("Earth has body_id={:d}, idx={:d}.\n", body_id, idx);

    // // Print the MJDs at the spline nodes
    // print("PlanetElement has {:d} MJDs forming spline nodes:\n", pe.N_t);
    // double* mjd = pe.get_mjd();
    // for (int i=0; i<pe.N_t; i++) {print("{:8.2f}, ", mjd[i]);}

    // // Print the element a at the spline nodes
    // print("\nValues of a at these nodes for Earth:\n");
    // double* a = pe.get_a(idx);
    // for (int i=0; i<pe.N_t; i++) {print("{:8.6f}, ", a[i]);}
    // print("\n");

    // // Manual spline of a
    // print("Manually splining a on gsl_spline.\n");
    // gsl_spline* gsl_interp_a = pe.elt_spline.a[idx];
    // gsl_interp_accel* acc = pe.acc;
    // double a_out = gsl_spline_eval(gsl_interp_a, mjd_test, acc);
    // print("a_out={:8.6f}.\n", a_out);

    // // Calculate interpolated elements
    // print("Splining orbital elements at {:8.4f}.\n", mjd_test);
    // OrbitalElement elt = pe.interp_elt(body_id, mjd_test);
    // print("a={:f}\n.", elt.a);

    // Calculate interpolated position and state vector of Earth
    Position pos = pe.interp_pos(body_id, mjd_test);
    StateVector vec = pe.interp_vec(body_id, mjd_test);
    // Wrap velocity of earth
    Velocity vel {.vx = vec.vx, .vy = vec.vy, .vz = vec.vz};

    // Expected state vector components - location of Juno at this time.  
    // Copy / paste from KS.GetStateVectors_Earth(58000, 58000, 1440);
    double qx =  0.9583774949482733;
    double qy = -0.31579917956842507;
    double qz = -0.0001266099206526;
    double vx =  0.005191609325627999;
    double vy =  0.016244467239388778;
    double vz = -0.0000000527623193;

    // Wrap expected position and velocity objects
    Position pos0 {.qx = qx, .qy = qy, .qz = qz};
    Velocity vel0 {.vx = vx, .vy = vy, .vz = vz};

    // Calculate norm of position and velocity difference
    double dq = norm(pos0, pos);
    double dv = norm(vel0, vel);

    // Report results
    print("Distance between interpolated state vectors and DB values for Earth @ {:9.4f}.\n", mjd_test);
    print("dq: {:8.2e} AU\n", dq);
    print("dv: {:8.2e} AU/day\n", dv);
}

// *****************************************************************************
void test_all()
{
    // Inputs used in testing
    int d0 = 0;
    int d1 = 1000000;
    int mjd0 = 58000;
    int mjd1 = 59000;
    int dt_min = 60;
    // int mjd0 = 57990;
    // int mjd1 = 58010;
    // int dt_min = 1440;

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Build PlanetElement
    PlanetElement pe(mjd0, mjd1, dt_min);
    pe.load(conn);
    pe.build_splines();
    print("\nBuilt PlanetElement object from mjd0 {:d} to mjd1 {:d} with time step {:d} minutes.\n", 
            mjd0, mjd1, dt_min);

    test_planet_element(pe);

    // Initialize DetectionTimeTable
    DetectionTimeTable dtt = DetectionTimeTable();
    print("Loaded DetectionTimeTable with {:d} detection times.\n", dtt.N());

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(d0, d1);
    dt.load();
    print("Loaded DetectionTable with detection_id in [{:d}, {:d}).\n", d0, d1);
    // print_detection(dt);

    // test_planet_element(pe);

    // Build orbital elements for asteroid 3 (Juno), which has 8 hits
    // Values copy / pasted from CALL KS.GetAsteroidElements(3, 4, 59000, 59000);
    OrbitalElement elt
    {
        .mjd    = 59000.0,
        .a      = 2.6682852999999986,
        .e      = 0.25693643000000027,
        .inc    = 0.22673642125828372,
        .Omega  = 2.9644675653853,
        .omega  =-1.9536135288017569,
        .f      = 46.517431678528794,
        .M      = 46.171557086433715
    };
    print("\nConstructed OrbitalElement for Juno @ {:f}\n", elt.mjd);
    print_orbital_element(elt);

    // Build CandidateElement for these elements
    CandidateElement ce(elt, dtt);
    print("\nConstructed CandidateElement from OrbitalElement.\n");

    // Calculate the trajectory of the elements matching Juno
    ce.calc_trajectory();
    print("Calculated trajectory of CandidateElement.\n");

    // Test the calculated trajectory
    test_calc_traj(elt, dtt);

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
    // test_massive_body();
}
