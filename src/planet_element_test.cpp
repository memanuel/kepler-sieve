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
#include "utils.hpp"
    using ks::norm;
    using ks::report_test;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::print_state_vector;
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

double mjd_test = 58000.0;

// *****************************************************************************
// Functions defined in this module
int main(int argc, char* argv[]);
void test_all();
void test_planet_vector(PlanetVector& pv);
void test_planet_element(PlanetElement& pe);

// *****************************************************************************
void test_planet_vector(PlanetVector& pv)
{
    // Test state vectors of Earth @ 58000
    int32_t body_id=body_id_earth;

    // Calculate interpolated position of Earth
    Position pos = pv.interp_pos(body_id, mjd_test);
    // Calculate and interpolated state vectors of Earth
    StateVector sv = pv.interp_vec(body_id, mjd_test);

    // Report splined state vectors
    print("\nSplined state vectors at {:8.4f}.\n", mjd_test);
    // print("{:9s} : {:9s} : {:9s} : {:9s} : {:9s} : {:9s}\n", 
    //     "qx", "qy", "qz", "vx", "vy", "vz");
    // print("{:+9.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f}\n", 
    //     sv.qx, sv.qx, sv.qz, sv.vx, sv.vx, sv.vz);
    print_state_vector(sv, true);
    print_state_vector(sv);

    // Wrap velocity of earth
    Velocity vel {.vx = sv.vx, .vy = sv.vy, .vz = sv.vz};

    // Expected state vector components
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
    // Relative difference
    // TODO

    // Report results
    print("\nDistance between interpolated state vectors and DB values for Earth @ {:9.4f}.\n", mjd_test);
    print("dq: {:8.2e} AU\n", dq);
    print("dv: {:8.2e} AU/day\n", dv);
}

// *****************************************************************************
void test_planet_element(PlanetElement& pe)
{
    // Test state vectors of Earth @ 58000
    int32_t body_id=399;   

    // Calculate idx from body_id
    // int idx = pe.body_idx(body_id);
    // print("Earth has body_id={:d}, idx={:d}.\n", body_id, idx);

    // Print the MJDs at the spline nodes
    // print("PlanetElement has {:d} MJDs forming spline nodes:\n", pe.N_t);
    // double* mjd = pe.get_mjd();
    // for (int i=0; i<pe.N_t; i++) {print("{:8.2f}, ", mjd[i]);}

    // Calculate and report interpolated elements on spline nodes
    OrbitalElement elt = pe.interp_elt(body_id, mjd_test);
    print("\nSplined orbital elements at {:8.4f}.\n", mjd_test);
    print("{:8s} : {:8s} : {:9s} : {:9s} : {:9s} : {:9s} : {:9s} \n", 
        "a", "e", "inc", "Omega", "omega", "f", "M");
    print("{:8.6f} : {:8.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f} : {:9.4f} : {:9.4f}\n", 
        elt.a, elt.e, elt.inc, elt.Omega, elt.omega, elt.f, elt.M);

    // double* mjd = pe.get_mjd();
    // for (int i=0; i<pe.N_t; i++) 
    // {
    //     OrbitalElement elt = pe.interp_elt(body_id, mjd[i]);
    //     print("{:10.4f} : {:8.6f} : {:8.6f}\n", mjd[i], elt.a, elt.e);
    // }

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
    print("\nDistance between interpolated state vectors and DB values for Earth @ {:9.4f}.\n", mjd_test);
    print("dq: {:8.2e} AU\n", dq);
    print("dv: {:8.2e} AU/day\n", dv);
}

// *****************************************************************************
void test_all()
{
    // Inputs used in testing
    int width = 100;
    int mjd0 = mjd_test - width;
    int mjd1 = mjd_test + width;
    int dt_min = 5;

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Timer object
    Timer t;

    // Build PlanetVector
    t.tick();
    PlanetVector pv(mjd0, mjd1, dt_min);
    pv.load();
    pv.build_splines();
    // PlanetVector pv = PlanetVector();
    print("\nBuilt PlanetVector object from mjd0 {:d} to mjd1 {:d} with time step {:d} minutes.\n", 
            pv.mjd0, pv.mjd1, pv.dt_min);
    t.tock_msg();

    // Test planet vectors
    test_planet_vector(pv);

    // // Build PlanetElement
    // t.tick();
    // PlanetElement pe(mjd0, mjd1, dt_min);
    // pe.load();
    // pe.build_splines();
    // // PlanetElement pe = PlanetElement();
    // print("\nBuilt PlanetElement object from mjd0 {:d} to mjd1 {:d} with time step {:d} minutes.\n", 
    //         pe.mjd0, pe.mjd1, pe.dt_min);
    // t.tock_msg();

    // // Test planet elements
    // test_planet_element(pe);

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
}
