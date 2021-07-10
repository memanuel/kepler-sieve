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


// *****************************************************************************
Simulation* make_sim()
{
    // Create an empty simulation
    Simulation* sim = reb_create_simulation();

    // Configure simulation
    // Unit system is (AU, Msun, day)

    // Set initial time step to 1.0 days
    sim->dt = 1.0;            // in days

    // Set gravitational constant in this unit system, i.e. in AU^3 / Msun / day^2
    sim->G  = G;

    // Use the IAS15 integrator
    sim->integrator = reb_simulation::REB_INTEGRATOR_IAS15;

    // Don't set a heartbeat for now
    // sim->heartbeat        = heartbeat;

    // Finish exactly at tmax in reb_integrate(). Default is already 1.
    sim->exact_finish_time = 1;

    // Return the assembled simulation object
    return sim;
}


// *****************************************************************************
Simulation* make_sim_planets(const PlanetVector& pv, double epoch)
{
    // Start with an empty simulation
    Simulation* sim = make_sim();

    // Set the time field in the simulation to match the epoch
    sim->t = epoch;

    // Load table of masses
    MassiveBodyTable mbt = MassiveBodyTable();

    // The body_ids from the PlanetVector object; useful in future for e.g. DE-435
    // const int32_t* body_ids = pv.get_body_id();

    // The number of bodies in the planets collection is known at compile time
    constexpr int N_body = 11;
    // The body_ids; add them to the simulation in Jacobi order rather than sorted by body_id
    constexpr int32_t body_ids[N_body] = {10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9};

    // Iterate through the bodies in the planet vector
    for (int i=0; i<N_body; i++)
    {
        // The body_id of this body
        int32_t body_id = body_ids[i];
    
        // Get the state vector of this particle at the epoch
        StateVector s = pv.interp_vec(body_id, epoch);

        // Initialize a particle with this state vector and the correct mass
        Particle p 
        {
            .x  = s.qx,
            .y  = s.qy,
            .z  = s.qz,
            .vx = s.vx,
            .vy = s.vy,
            .vz = s.vz,
            // Look up the mass on the MassiveBodyTable
            .m  = mbt.get_M(body_id)
        };

        // Add the particle to the simulation
        reb_add(sim, p); 
    }

    // Return the assmelbed simulation
    return sim;
}

// *****************************************************************************
} // Namespace ks
