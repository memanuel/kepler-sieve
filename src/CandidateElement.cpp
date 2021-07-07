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
CandidateElement::CandidateElement(OrbitalElement& elt, DetectionTimeTable& dtt): 
    elt(elt),
    N(dtt.N()),
    mjd(new double[N]),
    q_obs(new double[3*N]),
    q_ast(new double[3*N]),
    v_ast(new double[3*N]),
    u_ast(new double[3*N])
{
    // Populate mjds with a copy taken from dtt
    size_t sz_mjd = N*sizeof(mjd[0]);
    memcpy((void*) dtt.get_mjd(), mjd, sz_mjd);
    // Populate q_obs with a copy taken from dtt
    size_t sz_q_obs = 3*N*sizeof(q_obs[0]);
    memcpy((void*) dtt.get_q_obs(), q_obs, sz_q_obs);
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
}

// *****************************************************************************
// Calculate trajectory and direction
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calc_trajectory()
{
    // Calculate the mean motion, n
    double n = mean_motion(elt.a);
    // Unpack the five fields of orbital element that don't change
    double a = elt.a;
    double e = elt.e;
    double inc = elt.inc;
    double Omega = elt.Omega;
    double omega = elt.omega;
    // Get the reference date and mean anomaly for the elements
    double mjd0 = elt.mjd;
    double M0 = elt.M;

    // Iterate through the mjds in the array
    for (int i=0; i<N; i++)
    {
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*(mjd[i] - mjd0);
        // Convert M to a true anomaly f
        double f = anomaly_M2f(M, e);
        // Convert to a position
        StateVector sv = elt2vec(a, e, inc, Omega, omega, f);
        // Index for arrays q_ast and v_ast
        int idx = 3*i;
        // Save into array q_ast
        q_ast[idx+0] = sv.qx;
        q_ast[idx+1] = sv.qy;
        q_ast[idx+2] = sv.qz;
        // Save into array v_ast
        v_ast[idx+0] = sv.vx;
        v_ast[idx+1] = sv.vy;
        v_ast[idx+2] = sv.vz;
    }
}

// *****************************************************************************
void CandidateElement::calc_direction()
{
    ;
}

// *****************************************************************************
// Get read access to arrays
// *****************************************************************************

// *****************************************************************************
double* CandidateElement::get_mjd() const
{
    return mjd;
}

// *****************************************************************************
double* CandidateElement::get_q_obs() const
{
    return q_obs;
}

// *****************************************************************************
double* CandidateElement::get_q_ast() const
{
    return q_ast;
}

// *****************************************************************************
double* CandidateElement::get_v_ast() const
{
    return v_ast;
}

// *****************************************************************************
double* CandidateElement::get_u_ast() const
{
    return u_ast;
}
