/** @file CandidateElement.cpp
 *  @brief Implmentation of CandidateElement class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-07
 */

// *****************************************************************************
// Local Dependencies
#include "CandidateElement.hpp"

// *****************************************************************************
// Local names used
using ks::CandidateElement;

// *****************************************************************************
// Constructor & destructor
// *****************************************************************************

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id, int N_t): 
    // Small data elements
    elt {elt},
    candidate_id {candidate_id},
    elt0 {elt},
    N_t {N_t},
    N_row {3*N_t},
    // Allocate arrays
    mjd   {new double[N_t]},
    q_obs {new double[N_row]},
    q_ast {new double[N_row]},
    v_ast {new double[N_row]},
    u_ast {new double[N_row]},
    r_ast {new double[N_t]},
    dq_ast {new double[N_row]},
    dv_ast {new double[N_row]}
    {}

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id, const double* mjd_, int N_t):
    // Delegate to main constructor using mjd array from detection time table
    CandidateElement(elt, candidate_id, N_t)
{
    // Initialize mjd array with input times
    init(mjd_);
}

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id):
    // Delegate to constructor taking an mjd array and size; use mjd array from detection time table
    CandidateElement(elt, candidate_id, dtt.mjd, dtt.N)
{}

// *****************************************************************************
// Delete manually allocated arrays
CandidateElement::~CandidateElement()
{
    delete [] mjd;
    delete [] q_obs;
    delete [] q_ast;
    delete [] v_ast;
    delete [] u_ast;
    delete [] r_ast;
    delete [] dq_ast;
    delete [] dv_ast;
}

// *****************************************************************************
void CandidateElement::init(const double* mjd_)
{
    // Copy from mjd_ input to mjd on the candidate element
    for (int i=0; i<N_t; i++) {mjd[i]=mjd_[i];}

    // Populate q_cal and v_cal from BodyVector of the Sun
    // This will be a pretty accurate first pass before running an optional numerical integration
    for (int i=0; i<N_t; i++)
    {
        // The time of this entry
        double t = mjd[i];
        // Interpolated position vector of the Earth at time i
        Position p_earth = bv_earth.interp_pos(t);
        // Interpolated state vector of the Sun at time i
        StateVector s_sun = bv_sun.interp_vec(t);
        // Index for arrays q_sun and v_sun
        const int j0 = 3*i;
        const int jx = j0+0;
        const int jy = j0+1;
        const int jz = j0+2;
        // Copy position components into q_obs
        q_obs[jx] = p_earth.qx;
        q_obs[jy] = p_earth.qy;
        q_obs[jz] = p_earth.qz;
        // Copy position components into dq_ast
        dq_ast[jx] = s_sun.qx;
        dq_ast[jy] = s_sun.qy;
        dq_ast[jz] = s_sun.qz;
        // Copy velocity components into dv_ast
        dv_ast[jx] = s_sun.vx;
        dv_ast[jy] = s_sun.vy;
        dv_ast[jz] = s_sun.vz;
    }   // for / i
}   // CandidateElement::init

// *****************************************************************************
// Calculate and calibrate trajectory in Kepler model
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calc_trajectory()
{
    // Calculate the mean motion, n
    double n = mean_motion(elt.a, mu_sun);
    // Unpack the five fields of the orbital element that don't change
    double a = elt.a;
    double e = elt.e;
    double inc = elt.inc;
    double Omega = elt.Omega;
    double omega = elt.omega;
    // Get the reference date and mean anomaly for the elements
    double mjd0 = elt.mjd;
    double M0 = elt.M;

    // Iterate through the mjds in the array
    for (int i=0; i<N_t; i++)
    {
        // The time shift from the reference time
        double dt = mjd[i] - mjd0;
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*dt;
        // Convert M to a true anomaly f
        double f = anomaly_M2f(M, e);
        // Convert orbital elements to a heliocentric position
        StateVector s = elt2vec(a, e, inc, Omega, omega, f, mu_sun);
        // Index for arrays q_ast and v_ast
        const int j = 3*i;
        const int jx = j+0;
        const int jy = j+1;
        const int jz = j+2;
        // Save into array q_ast
        q_ast[jx] = s.qx + dq_ast[jx];
        q_ast[jy] = s.qy + dq_ast[jy];
        q_ast[jz] = s.qz + dq_ast[jz];
        // Save into array v_ast
        v_ast[jx] = s.vx + dv_ast[jx];
        v_ast[jy] = s.vy + dv_ast[jy];
        v_ast[jz] = s.vz + dv_ast[jz];
    }
}

// *****************************************************************************
void CandidateElement::calibrate()
{
    // Get the reference time for the elements
    double mjd0 = elt.mjd;

    // Zero out old calibration
    for (int k=0; k<N_row; k++)
    {
        dq_ast[k] = 0.0;
        dv_ast[k] = 0.0;
    }

    // Calculate the trajectory without calibration
    calc_trajectory();

    // Build rebound simulation with planets at reference time
    Simulation sim = make_sim_planets(pv, mjd0);
    // Add asteroid to simulation with candidate elements
    sim.add_test_particle(elt, candidate_id);

    // Write the numerically integrated vectors from this simulation into q_cal and v_cal
    sim.write_vectors(dq_ast, dv_ast, mjd, N_t, candidate_id);

    // Iterate through all the rows; subtract the Kepler component from the numerical component
    // This is equivalent to writing
    // dq_ast[k] = numerical_trajectory[k] - kepler_trajectory[k]
    // It just uses less memory by using dq_ast to store the numerical trajectory and subtracting in place
    for (int k=0; k<N_row; k++)
    {
        // double q_num = dq_ast[k];
        dq_ast[k] -= q_ast[k];
        dv_ast[k] -= v_ast[k];
        // print("k={:d}, q_num = {:+10.6f}, q_kep = {:+10.6f}, dq = {:+10.6f}.\n", k, q_num, q_ast[k], dq_ast[k]);
    }
}

// *****************************************************************************
// Get predicted state vectors
// *****************************************************************************

// *****************************************************************************
const StateVector CandidateElement::state_vector(int i) const
{
    // The array index for the test date
    int j = 3*i;
    int jx = j+0;  
    int jy = j+1; 
    int jz = j+2;

    // Wrap the predicted position and velocity of the asteroid on the test date into a StateVector
    return StateVector 
    {
        .qx = q_ast[jx], 
        .qy = q_ast[jy],
        .qz = q_ast[jz],
        .vx = v_ast[jx],
        .vy = v_ast[jy],
        .vz = v_ast[jz]
    };
}

// *****************************************************************************
const Position CandidateElement::observer_pos(int i) const
{
    // The array index for the test date
    int j = 3*i;
    int jx = j+0;  
    int jy = j+1; 
    int jz = j+2;

    // Wrap the observer position on the test date into a Position
    return Position 
    {
        .qx = q_obs[jx], 
        .qy = q_obs[jy],
        .qz = q_obs[jz]
    };
}

// *****************************************************************************
// Calculate direction
// *****************************************************************************

// *****************************************************************************
// Calculate direction based on observer position q_obs and state vector q_ast, v_ast
// Caller is responsible to populate these before running calc_direction()
void CandidateElement::calc_direction()
{
    // Iterate over all N observation times
    for (int i=0; i<N_t; i++)
    {
        // // Wrap the observer position
        // Position q_obs {observer_pos(i)};
        // // Wrap the asteroid state vector
        // StateVector s_ast {state_vector(i)};
        // // Calculate direction and distance in Newtonian light model
        
    }
}

// *****************************************************************************
// Initialize the static member PlanetVector
// *****************************************************************************

BodyVector CandidateElement::bv_sun = BodyVector(SolarSystemBody_bv::sun);  

BodyVector CandidateElement::bv_earth = BodyVector(SolarSystemBody_bv::earth);

// Build DetectionTimeTable object as a variable first.
// Need to extract first and last time below.
// If the dtt is built immediately as a static member of CandidateElement, it is inaccessible here.
// The workaround is to (1) build it locally (2) extract first / last date (3) move it to the static member
// This sequence avoids building the object twice.
DetectionTimeTable dtt {DetectionTimeTable()};

// Inputs to build a PlanetVector table suited to CandidateElement
constexpr int pad = 32;
const int mjd_first = static_cast<int>(floor(dtt.mjd_first()));
const int mjd_last = static_cast<int>(ceil(dtt.mjd_last()));
const int mjd0 = std::min(mjd_first, 58000) - pad;
const int mjd1 = std::max(mjd_last,  60000) + pad;
constexpr int dt_min = 60;
bool load = true;
PlanetVector CandidateElement::pv {PlanetVector(mjd0, mjd1, dt_min, load)};
// PlanetVector CandidateElement::pv {PlanetVector(59000-32, 60000+32, 1440, true)};

// // Move the local copy of DetectionTable to avoid cost of building it twice
DetectionTimeTable CandidateElement::dtt = std::move(dtt);

// DEBUG - build a small Detection table quickly for testing
// DetectionTable CandidateElement::dt = DetectionTable()
int d1 = 85000000;
int d2 = 86000000;
DetectionTable CandidateElement::dt = DetectionTable(d1, d2, load);
