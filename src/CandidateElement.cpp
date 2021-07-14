/** @file CandidateElement.cpp
 *  @brief Implmentation of CandidateElement class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-07
 */

// *****************************************************************************
// Included files
#include "CandidateElement.hpp"

// *****************************************************************************
// Local names used
using ks::CandidateElement;

// *****************************************************************************
// Constructor & deestructor
// *****************************************************************************

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, DetectionTimeTable& dtt, int32_t candidate_id): 
    elt {elt},
    candidate_id {candidate_id},
    elt0 {elt},
    N_t {dtt.N()},
    mjd {new double[N_t]},
    q_obs {new double[3*N_t]},
    q_ast {new double[3*N_t]},
    v_ast {new double[3*N_t]},
    u_ast {new double[3*N_t]},
    q_cal {new double[3*N_t]},
    v_cal {new double[3*N_t]}
{
    // Populate mjd array with a copy taken from dtt
    size_t sz_mjd = N_t*sizeof(mjd[0]);
    memcpy((void*) dtt.get_mjd(), mjd, sz_mjd);

    // Populate q_obs with a copy taken from dtt
    size_t sz_q_obs = 3*N_t*sizeof(q_obs[0]);
    memcpy((void*) dtt.get_q_obs(), q_obs, sz_q_obs);

    // Intitialize calibration adjustments to zero
    for (int i=0; i<3*N_t; i++)
    {
        q_cal[i] = 0.0;
        v_cal[i] = 0.0;
    }
}

// *****************************************************************************
// Delete manually allocated arrays
CandidateElement::~CandidateElement()
{
    delete [] mjd;
    delete [] q_obs;
    delete [] q_ast;
    delete [] v_ast;
    delete [] u_ast;
    delete [] q_cal;
    delete [] v_cal;
}

// *****************************************************************************
// Calculate trajectory and direction
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calibrate(const PlanetVector& pv)
{
    // Get the reference time for the elements
    double mjd0 = elt.mjd;

    // Build rebound simulation with planets at reference time
    Simulation sim = make_sim_planets(pv, mjd0);
    // Add asteroid to simulation with candidate elements
    sim.add_test_particle(elt, candidate_id);
    // Write the numerically integrated vectors from this simulation into q_cal and v_cal
    // sim.write_vectors(mjd, N_t, q_cal, v_cal);

    // Unpack the five fields of the orbital element that don't change
    double a = elt.a;
    double e = elt.e;
    double inc = elt.inc;
    double Omega = elt.Omega;
    double omega = elt.omega;
    // Get the reference mean anomaly for the elements
    double M0 = elt.M;
    // Calculate the mean motion, n
    double n = mean_motion(elt.a, mu_sun);

    // Iterate through the times in the mjd array
    for (int i=0; i<N_t; i++)
    {
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*(mjd[i] - mjd0);
        // Convert M to a true anomaly f
        double f = anomaly_M2f(M, e);
        // Convert these elements to a position
        StateVector sv = elt2vec(a, e, inc, Omega, omega, f, mu_sun);
        // Index for arrays q_ast and v_ast
        int idx = 3*i;
        // Save into array q_ast with calibration adjustment
        q_ast[idx+0] = sv.qx + q_cal[idx+0];
        q_ast[idx+1] = sv.qy + q_cal[idx+1];
        q_ast[idx+2] = sv.qz + q_cal[idx+2];
        // Save into array v_ast
        v_ast[idx+0] = sv.vx + v_cal[idx+0];
        v_ast[idx+1] = sv.vy + v_cal[idx+1];
        v_ast[idx+2] = sv.vz + v_cal[idx+2];
    }
}

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
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*(mjd[i] - mjd0);
        // Convert M to a true anomaly f
        double f = anomaly_M2f(M, e);
        // Convert to a position
        StateVector sv = elt2vec(a, e, inc, Omega, omega, f, mu_sun);
        // Index for arrays q_ast and v_ast
        int idx = 3*i;
        // Save into array q_ast with calibration adjustment
        q_ast[idx+0] = sv.qx + q_cal[idx+0];
        q_ast[idx+1] = sv.qy + q_cal[idx+1];
        q_ast[idx+2] = sv.qz + q_cal[idx+2];
        // Save into array v_ast
        v_ast[idx+0] = sv.vx + v_cal[idx+0];
        v_ast[idx+1] = sv.vy + v_cal[idx+1];
        v_ast[idx+2] = sv.vz + v_cal[idx+2];
    }
}

// *****************************************************************************
void CandidateElement::calc_direction()
{
    ;
}

