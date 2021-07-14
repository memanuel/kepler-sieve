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

// *****************************************************************************
Orbit particle2orbit(const Particle& target, const Particle& primary)
    {return reb_tools_particle_to_orbit(G, target, primary);}

// *****************************************************************************
Particle state_vector2particle(const StateVector& s, double m, int32_t body_id)
{
    // Convert body_id to a hash_id; type in rebound in unsigned rather than signed int32
    uint32_t hash = static_cast<uint32_t>(body_id);
    // Initialize a particle with this state vector and the correct mass
    return Particle 
    {
        .x  = s.qx, .y  = s.qy, .z  = s.qz,
        .vx = s.vx, .vy = s.vy, .vz = s.vz,
        .m  = m, 
        .hash = hash
    };
}

// *****************************************************************************
Particle elt2particle(const OrbitalElement& elt, double m, const Particle& primary)
    {return reb_tools_orbit_to_particle(G, primary, m, elt.a, elt.e, elt.inc, elt.Omega, elt.omega, elt.f);}

// *****************************************************************************
Particle orbit2particle(const Orbit& orb, double m, const Particle& primary)
    {return reb_tools_orbit_to_particle(G, primary, m, orb.a, orb.e, orb.inc, orb.Omega, orb.omega, orb.f);}

// *****************************************************************************
} // Namespace reb
} // Namespace ks
