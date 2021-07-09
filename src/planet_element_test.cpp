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
#include "PlanetElement.hpp"
    using ks::PlanetElement;    

// *****************************************************************************
// Functions defined in this module
int main(int argc, char* argv[]);
void test_all();
void test_planet_element(PlanetElement& pe);

// *****************************************************************************
void test_planet_element(PlanetElement& pe)
{
    // Test state vectors of Earth @ 58000
    int32_t body_id=399;
    double mjd_test = 58000.0;

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
    // int mjd0 = 57995;
    // int mjd1 = 58005;
    // int dt_min = 5;

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Timer object
    Timer t;

    // Build PlanetElement    
    t.tick();
    // PlanetElement pe(mjd0, mjd1, dt_min);
    // pe.load();
    // pe.build_splines();
    PlanetElement pe = PlanetElement();
    print("\nBuilt PlanetElement object from mjd0 {:d} to mjd1 {:d} with time step {:d} minutes.\n", 
            pe.mjd0, pe.mjd1, pe.dt_min);
    t.tock_msg();
    test_planet_element(pe);

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
}
