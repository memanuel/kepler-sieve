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
CandidateElement::CandidateElement(OrbitalElement elt, DetectionTimeTable dtt): 
    N(dtt.N()),
    mjds(new double[N]),
    q_ast(new double[3*N]),
    v_ast(new double[3*N]),
    u_ast(new double[3*N])
{}

// *****************************************************************************
// Delete manually allocated arrays
CandidateElement::~CandidateElement()
{
    // DEBUG
    print("Entering destructor for CandidateElement.\n");
    ks::flush_console();
    delete [] mjds;
    print("Deleted mjd.\n");
    ks::flush_console();
    delete [] q_ast;
    delete [] v_ast;
    delete [] u_ast;
    print("Deleted q_ast, v_ast, u_ast.\n");
    ks::flush_console();
}
