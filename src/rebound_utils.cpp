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
// Class Simulation
// *****************************************************************************

// *****************************************************************************
Simulation::Simulation():
    // Create an empty simulation and save pointer to it
    prs {reb_create_simulation()}
{
    // Unit system is (AU, Msun, day)
    // Set gravitational constant in this unit system, i.e. in AU^3 / Msun / day^2
    prs->G  = G;

    // Set initial time step to 1.0 days
    prs->dt = 1.0;

    // Use the IAS15 integrator
    prs->integrator = reb_simulation::REB_INTEGRATOR_IAS15;

    // Finish exactly at tmax in reb_integrate(). Default is already 1.
    prs->exact_finish_time = 1;

    // Set number of active particles
    prs->N_active = 0;
}


// *****************************************************************************
Simulation::Simulation(const Simulation& s):
    // Create an empty simulation and save pointer to it
    prs {reb_create_simulation()}
{
    // Copy configuration from s
    prs->G = s.prs->G;
    prs->dt = s.prs->dt;
    prs->integrator = s.prs->integrator;
    prs->exact_finish_time = s.prs->exact_finish_time;
    prs->N_active = s.prs->N_active;

    // Copy particles from s
    for (int i=0; i<s.N(); i++)    
    {
        // The original particle
        // Particle p0 = s.particle(i);
        // The copy of the particle
        Particle p {s.particle(i)};
        // Add the copied particle to the new simulation
        reb_add(prs, p);
    }
}

// *****************************************************************************
Simulation::~Simulation()
{
    // Free memory in underlying rebound simulation
    reb_free_simulation(prs);
}

// *****************************************************************************
void Simulation::set_time(double t) {prs->t = t;}

// *****************************************************************************
void Simulation::add_particle(const StateVector& s, double m)
{
    // If this is a massive particle to be added, Check that no test particles present; 
    // must first add all massive particles, only then can test particle be added
    if ((m>0) && N_test()>0) 
    {
        throw(runtime_error("Simulation::add_particle(), "
        "can't add massive particle after any test particle added.\n"));
    }

    // Add the particle to the simulation
    add_particle_impl(s, m);
}

// *****************************************************************************
void Simulation::add_test_particle(const StateVector& s) 
{
    // Always safe to add a test particle; don't need to check    
    add_particle_impl(s, 0.0);
    // // Decrement active particles
    // prs->N_active--;
}

// *****************************************************************************
void Simulation::add_particle_impl(const StateVector& s, double m)
{
    // Initialize a particle with this state vector and the correct mass
    Particle p 
    {
        .x  = s.qx,
        .y  = s.qy,
        .z  = s.qz,
        .vx = s.vx,
        .vy = s.vy,
        .vz = s.vz,
        .m  = m
    };

    // Add the particle to the simulation
    reb_add(prs, p);

    // Increment active particles
    prs->N_active++;    
}

// *****************************************************************************
// Factor function make_sim_planets
// *****************************************************************************

// *****************************************************************************
Simulation make_sim_planets(const PlanetVector& pv, double epoch)
{
    // Start with an empty simulation
    Simulation sim;

    // Set the time field in the simulation to match the epoch
    sim.set_time(epoch);

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

        // Look up the mass on the MassiveBodyTable
        double m = mbt.get_M(body_id);

        // Add this particle
        sim.add_particle(s, m);
    }

    // Return the assembled simulation
    return sim;
}

// *****************************************************************************
} // Namespace reb
} // Namespace ks
