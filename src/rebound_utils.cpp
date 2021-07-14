/** @file   rebound_utils.cpp
 *  @brief  Implementation of rebound_utils.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-10
 * 
 */

// *****************************************************************************
// Local dependencies
#include "rebound_utils.hpp"

// *****************************************************************************
namespace ks {
namespace reb {

// *****************************************************************************
// Utility functions
// *****************************************************************************

StateVector particle2state_vector(const Particle& p)
    {return StateVector{.qx=p.x, .qy=p.y, .qz=p.z, .vx=p.vx, .vy=p.vy, .vz=p.vz};}

Orbit particle2orbit(const Particle& target, const Particle& primary)
    {return reb_tools_particle_to_orbit(G, target, primary);}

// *****************************************************************************
} // Namespace reb
} // Namespace ks
