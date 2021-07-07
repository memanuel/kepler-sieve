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
// Class CandidateElement
// *****************************************************************************

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement& elt, DetectionTimeTable& dtt): 
    elt(elt),
    N(dtt.N()),
    mjds(new double[N]),
    q_obs(new double[3*N]),
    q_ast(new double[3*N]),
    v_ast(new double[3*N]),
    u_ast(new double[3*N])
{
    // Populate mjds with a copy taken from dtt
    size_t sz_mjd = N*sizeof(mjds[0]);
    memcpy((void*) dtt.get_mjds(), mjds, sz_mjd);
    // Populate q_obs with a copy taken from dtt
    size_t sz_q_obs = 3*N*sizeof(q_obs[0]);
    memcpy((void*) dtt.get_q_obs(), q_obs, sz_q_obs);
}

// *****************************************************************************
// Delete manually allocated arrays
CandidateElement::~CandidateElement()
{
    delete [] mjds;
    delete [] q_obs;
    delete [] q_ast;
    delete [] v_ast;
    delete [] u_ast;
}

// *****************************************************************************
void CandidateElement::calc_trajectory()
{
    ;
}

// *****************************************************************************
void CandidateElement::calc_direction()
{
    ;
}
