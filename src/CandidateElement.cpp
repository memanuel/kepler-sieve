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
    CandidateElement(elt, candidate_id, dtt.get_mjd(), dtt.N())
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
        // Interpolated position vector of the Earth at time i
        Position p = bv_earth.interp_pos(mjd[i]);
        // Interpolated state vector of the Sun at time i
        StateVector s = bv_sun.interp_vec(mjd[i]);
        // Index for arrays q_sun and v_sun
        const int j = 3*i;
        const int jx = j+0;
        const int jy = j+1;
        const int jz = j+2;
        // Copy position components into q_obs
        q_obs[jx] = p.qx;
        q_obs[jy] = p.qy;
        q_obs[jz] = p.qz;
        // Copy position components into q_cal
        dq_ast[jx] = s.qx;  
        dq_ast[jy] = s.qy;  
        dq_ast[jz] = s.qz;
        // Copy velocity components into v_cal
        dv_ast[jx] = s.vx;  
        dv_ast[jy] = s.vy;  
        dv_ast[jz] = s.vz;
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
        double t = mjd[i] - mjd0;
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*t;
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
void CandidateElement::calibrate(const PlanetVector& pv)
{
    // Get the reference time for the elements
    double mjd0 = elt.mjd;

    // Zero out old calibration
    for (int k=0; k<N_row; k++)
    {
        dq_ast[k]=0.0;
        dv_ast[k]=0.0;
    }

    // Calculate the trajectory without calibration
    calc_trajectory();
    // DEBUG
    print("CandidateElement::calibrate - ran calc_trajectory()\n");

    // Build rebound simulation with planets at reference time
    Simulation sim = make_sim_planets(pv, mjd0);
    // DEBUG
    print("CandidateElement::calibrate - created Simulation with sim = make_sim_planets(pv, mjd0)\n");
    // Add asteroid to simulation with candidate elements
    sim.add_test_particle(elt, candidate_id);
    // DEBUG
    print("CandidateElement::calibrate - added test particle for candidate elements.\n");
    // Write the numerically integrated vectors from this simulation into q_cal and v_cal
    sim.write_vectors(mjd, candidate_id, N_t, dq_ast, dv_ast);

    // Iterate through all the rows; subtract the Kepler component from the numerical component
    // This is equivalent to writing
    // dq_ast[k] = numerical_trajectory[k] - kepler_trajectory[k]
    // It just uses less memory by using dq_ast to store the numerical trajectory and subtracting in place
    for (int k=0; k<N_row; k++)
    {
        dq_ast[k] -= q_ast[k];
        dv_ast[k] -= v_ast[k];
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
// Calculate direction
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calc_direction()
{
    ;
}

