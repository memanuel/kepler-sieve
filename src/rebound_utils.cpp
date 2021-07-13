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
// Class Simulation
// *****************************************************************************

// *****************************************************************************
Simulation::Simulation():
    // Create an empty simulation and save pointer to it
    prs {reb_create_simulation()},
    // Initialize the vectors of body_ids (the particle and its primary)
    body_ids {},
    primary_body_ids {},
    // Initialize maps from body_id to body_idx (array position of in particles)
    body_map {},
    primary_body_map {}
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
    prs {reb_create_simulation()},
    // Copy the vectors of body_id and primary_body_id from s,
    body_ids {s.body_ids},
    primary_body_ids {s.body_ids},
    // Copy the body_map and primary_body_map from s
    body_map {s.body_map},
    primary_body_map {s.primary_body_map}
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
        // Look up body_id and primary_body_id
        const int32_t body_id = body_ids[i];
        const int32_t primary_body_id = primary_body_ids[i];
        // The particle
        Particle p {s.particle(i)};
        // Get state vector and mass of this particle
        StateVector s = particle2state_vector(p);
        double m = p.m;
        // Add this particle to the new simulation
        add_particle(s, m, body_id, primary_body_id);
    }
}

// *****************************************************************************
Simulation::~Simulation()
{
    // Free memory in underlying rebound simulation
    reb_free_simulation(prs);
}

// *****************************************************************************
void Simulation::add_particle(const StateVector& s, double m, int32_t body_id, int32_t primary_body_id)
{
    // If this is a massive particle to be added, Check that no test particles present; 
    // must first add all massive particles, only then can test particle be added
    if ((m>0) && N_test()>0) 
    {
        throw(runtime_error("Simulation::add_particle(), "
        "can't add massive particle after any test particle added.\n"));
    }

    // Add the particle to the simulation
    add_particle_impl(s, m, body_id, primary_body_id);
}

// *****************************************************************************
void Simulation::add_particle_impl(const StateVector& s, double m, int32_t body_id, int32_t primary_body_id)
{
    // Convert body_id to a hash_id; type in rebound in unsigned rather than signed int32
    uint32_t hash = static_cast<uint32_t>(body_id);
    // Initialize a particle with this state vector and the correct mass
    Particle p 
    {
        .x  = s.qx, .y  = s.qy, .z  = s.qz,
        .vx = s.vx, .vy = s.vy, .vz = s.vz,
        .m  = m, 
        .hash = hash
    };

    // Add the particle to the simulation
    reb_add(prs, p);
    // Increment active particles
    prs->N_active++;

    // Add the body_id and primary_body_id to vectors
    body_ids.push_back(body_id);
    primary_body_ids.push_back(primary_body_id);

    // Add the body_ids to body_map and primary_body_map
    int i = prs->N-1;
    body_map[body_id]=i;
    primary_body_map[body_id]=i;
}

// *****************************************************************************
void Simulation::add_test_particle(const StateVector& s, int32_t body_id, int32_t primary_body_id)
{
    // Always safe to add a test particle; don't need to check    
    add_particle_impl(s, 0.0, body_id, primary_body_id);
    // Decrement active particles
    prs->N_active--;
}

// *****************************************************************************
const StateVector Simulation::state_vector(int i) const
{
    // The particle
    Particle p {particle(i)};
    // Wrap this into a StateVector
    return StateVector {.qx=p.x, .qy=p.y, .qz=p.z, .vx=p.vx, .vy=p.vy, .vz=p.vz};
}

// *****************************************************************************
void Simulation::print_vectors() const
{
    // Display the particles
    print_newline();
    fmt::print("{:4s}: {:12s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s} \n",
          "   i", "body_name", " mass", "  qx", "  qy", "  qz", "  vx", "  vy", "  vz");
    for (int i=0; i<N(); i++)
    {
        Particle p = particle(i);
        const int32_t body_id = body_ids[i];
        const string body_name = get_body_name(body_id);
        fmt::print("{:4d}: {:12s}: {:10.3e}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f} \n", 
              i, body_name, p.m, p.x, p.y, p.z, p.vx, p.vy, p.vz);    
    }
    print_newline();
}

// *****************************************************************************
void Simulation::print_elements() const
{
    // Display the particles
    print_newline();
    fmt::print("{:4s}: {:12s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s} \n",
          "   i", "body_name", "  a", "e", "  inc", "  Omega", "  omega", "  f", "  M");
    for (int i=0; i<N(); i++)
    {
        // The body_id of the target and its primary
        const int32_t body_id = body_ids[i];
        const int32_t primary_body_id = primary_body_ids[i];
        const string body_name = get_body_name(body_id);
        // Particle of the target
        const Particle& target = particle(i);
        // Only print elements if this particle has a real primary, i.e. skip the Sun
        if (!primary_body_id) {continue;}
        // Particle of the primary
        int j = primary_body_map.at(primary_body_id);
        const Particle& primary = particle(j);

        // Print the orbit of this particle w.r.t. its primary
        Orbit orb = particle2orbit(target, primary);
        fmt::print("{:4d}: {:12s}: {:10.6f}: {:10.8f}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f} \n", 
              i, body_name, orb.a, orb.e, orb.inc, orb.Omega, orb.omega, orb.f, orb.M);
    }
    print_newline();
}

// *****************************************************************************
// Factory function make_sim_planets
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

    // Iterate through the bodies in this collection
    // Add a particle for each body with initial conditions based on interpolated state vector on PlanetVector
    for (int i=0; i<N_body_planets; i++)
    {
        // The body_id and primary_id of this body
        int32_t body_id = body_ids_planets[i];
        int32_t primary_body_id = get_primary_body_id(body_id);
        // Get the state vector of this particle at the epoch
        StateVector s = pv.interp_vec(body_id, epoch);
        // Look up the mass on the MassiveBodyTable
        double m = mbt.get_M(body_id);
        // Add this particle
        sim.add_particle(s, m, body_id, primary_body_id);
    }

    // Return the assembled simulation
    return sim;
}

// *****************************************************************************
Simulation make_sim_planets(const PlanetElement& pe, double epoch)
{
    // Start with an empty simulation
    Simulation sim;
    // Set the time field in the simulation to match the epoch
    sim.set_time(epoch);
    // Load table of masses
    MassiveBodyTable mbt = MassiveBodyTable();
    // The body_ids; add them to the simulation in Jacobi order rather than sorted by body_id
    constexpr int32_t body_ids[] = {10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9};

    // Iterate through the bodies in this collection
    // Add a particle for each body with initial conditions based on interpolated state vector on PlanetElement
    for (int i=0; i<N_body_planets; i++)
    {
        // The body_id and primary_id of this body
        int32_t body_id = body_ids[i];
        int32_t primary_body_id = get_primary_body_id(body_id);
        // Get the state vector of this particle at the epoch
        StateVector s = pe.interp_vec(body_id, epoch);
        // Look up the mass on the MassiveBodyTable
        double m = mbt.get_M(body_id);
        // Add this particle
        sim.add_particle(s, m, body_id, primary_body_id);
    }

    // Return the assembled simulation
    return sim;
}

// *****************************************************************************
} // Namespace reb
} // Namespace ks
