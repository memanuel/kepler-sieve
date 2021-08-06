/** @file   Simulation.cpp
 *  @brief  Implementation of Simulation class.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-10
 * 
 */

// *****************************************************************************
// Local dependencies
#include "Simulation.hpp"

// *****************************************************************************
namespace ks {
namespace reb {

// *****************************************************************************
// Constructor and Destructor
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
// Add Particles
// *****************************************************************************

// *****************************************************************************
void Simulation::add_particle(Particle& p, int32_t body_id, int32_t primary_body_id)
{
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
    // Convert the state vector to a particle
    Particle p = state_vector2particle(s, m, body_id);
    // Add the particle to the simulation
    add_particle(p, body_id, primary_body_id);
}

// *****************************************************************************
void Simulation::add_test_particle(const StateVector& s, int32_t body_id, int32_t primary_body_id)
{
    // Always safe to add a test particle; don't need to check. Just add it from its state vector.
    add_particle_impl(s, 0.0, body_id, primary_body_id);
    // Decrement active particles
    prs->N_active--;
}

// *****************************************************************************
void Simulation::add_test_particle(const OrbitalElement& elt, int32_t body_id, int32_t primary_body_id)
{
    // Look up the primary particle from its body_id
    Particle primary = particle_b(primary_body_id);
    // Build a test particle from this orbital element
    Particle p = elt2particle(elt, 0.0, primary);
    // Add the particle
    add_particle(p, body_id, primary_body_id);
    // Decrement active particles
    prs->N_active--;
}

// *****************************************************************************
void Simulation::add_test_particle(const Orbit& orb, int32_t body_id, int32_t primary_body_id)
{
    // Look up the primary particle from its body_id
    Particle primary = particle_b(primary_body_id);
    // Build a test particle from this orbit
    Particle p = orbit2particle(orb, 0.0, primary);
    // Add the particle
    add_particle(p, body_id, primary_body_id);
    // Decrement active particles
    prs->N_active--;
}

// *****************************************************************************
// Get state vectors, particles and orbits
// *****************************************************************************

// *****************************************************************************
const StateVector Simulation::state_vector(int i) const
{
    // The particle
    Particle p {particle(i)};
    // Wrap this into a StateVector
    return StateVector {.qx=p.x, .qy=p.y, .qz=p.z, .vx=p.vx, .vy=p.vy, .vz=p.vz};
}

// *****************************************************************************
const Orbit Simulation::orbit(int i) const
{
    // The body_id of the primary
    const int32_t primary_body_id = primary_body_ids[i];
    // Only return elements if this particle has a real primary, i.e. skip the Sun
    if (!primary_body_id) 
    {
        string msg = format("Simulation::orbit - bad particle i={:d}, body_id={:d} has no primary.");
        throw invalid_argument(msg);
    }
    // Particle of the target
    const Particle& target = particle(i);
    // Particle of the primary
    int j = primary_body_map.at(primary_body_id);
    const Particle& primary = particle(j);
    // Return the orbit of this particle w.r.t. its primary
    return particle2orbit(target, primary);
}

// *****************************************************************************
// Write arrays of vectors and orbital elements
// *****************************************************************************

// *****************************************************************************
// Wrap implementation functions into anonymous namespace
namespace {

// *****************************************************************************
// Implementation function - write out state vectors for one batch of bodies
void write_vector_row_batch(double* q, double* v, Simulation& sim, vector<int32_t> body_ids, int i, double t)
{
    // Number of bodies in output array
    const int N_body_out = body_ids.size();
    // Integrate to this time (either forward or backward)
    sim.integrate(t);
    // The array index for this page of output
    int j0 = 3*N_body_out*i;
    /// Iterate through the output particles; k is index into body_id
    for (int k=0; k<N_body_out; k++)
    {
        // The array index for this row of output
        int j = j0 + 3*k;
        // The body_id of this body
        int32_t body_id = body_ids[k];
        /// Get the StateVector for particle b
        StateVector s = sim.state_vector_b(body_id);
        /// Write the position array
        q[j+0] = s.qx;
        q[j+1] = s.qy;
        q[j+2] = s.qz;
        /// Write the velocity array
        v[j+0] = s.vx;
        v[j+1] = s.vy;
        v[j+2] = s.vz;
    }   // for / bodies
}   // function write_vector_row_batch
// *****************************************************************************
}   // anonymous namespace

// *****************************************************************************
/** Write an array of state vectors to q and v obtained by integrating this simulation.
 * \param[in] q - array of output positions written to; size is 3*N_t; layout (time, axis)
 * \param[in] v - array of output velocities written to; size is 3*N_t*N_body; layout (time, axis)
 * \param[in] body_id - integer ID of the body whose state vectors are being written 
 * \param[in] mjd - vector of dates when output is desired; size N_t
 * Caller is responsible to allocate arrays q and v of the correct size and delete them later. */
void Simulation::write_vectors(double* q, double* v, const double* mjd, const int N_t, int32_t body_id) const
{
    // Wrap body_id into a vector of length 1
    vector<int32_t> body_ids = {body_id};
    // Delegate to Simulation::write_vectors_batch
    write_vectors_batch(q, v, mjd, N_t, body_ids);
}

// *****************************************************************************
/** Write an array of state vectors to q and v obtained by integrating this simulation. 
 * \param[in] q - array of output positions; size is 3*N_t*N_body; layout (time, body, axis)
 * \param[in] v - array of output velocities; size is 3*N_t*N_body; layout (time, body, axis)
 * \param[in] mjd - vector of dates when output is desired; size N_t
 * \param[in] body_id - vector of integer body IDs whose position and velocity are requested, size N_body
 * Caller is responsible to allocate arrays q and v of the correct size and delete them later. */
void Simulation::write_vectors_batch(
        double* q, double* v, const double* mjd, const int N_t, const vector<int32_t>& body_ids) const
{
    // Alias the initial reference time
    const double t0 = t();
    // Number of times in the output
    // const int N_t = mjd.size();
    // Create a copy of the simulation for integrating forward
    Simulation sim_f = this->copy();
    // The last index that is prior to t; this will be integrated backwards down to 0
    int i0_back=-1;

    // *****************************************************************************
    // Loop forward for times t in mjd where t>= t0
    // *****************************************************************************

    // DEBUG
    print("Simulation::write_vectors_batch - entering forward loop...\n");
    print("t0={:.1f}, N_t={:d}.\n", t0, N_t);

    // Iterate forward through times in mjd on or after t0
    for (int i=0; i<N_t; i++)
    {
        // Time of this output
        double t = mjd[i];
        // If the output time is before the simulation time, update i0_back and keep going
        if (t < t0) 
        {
            i0_back = i;
            continue;
        }
        // DEBUG
        print("Simulation::write_vectors_batch - i={:d}, t={:f}.\n", i, mjd[i]);
        // If we get here, the output time is on or after the simulation time. 
        // Integrate forward to this time and save the output.
        write_vector_row_batch(q, v, sim_f, body_ids, i, t);
    } // for / forward over time (i)

    // DEBUG
    print("Simulation::write_vectors_batch - i0_back={:d}.\n", i0_back);
    return;

    // *****************************************************************************
    // Loop backward for times t in mjd where t< t0
    // *****************************************************************************

    // Create a copy of sim for the backwards integration
    Simulation sim_b = this->copy();
    // Set time step in sim_back to be negative for backwards integration
    // sim_b.prs->dt = -prs->dt;
    // Iterate backard through times in mjd prior to t0
    for (int i=i0_back; i>=0; i--)
    {
        // Time of this output
        double t = mjd[i];
        // Integrate backward to this time and save the output.
        write_vector_row_batch(q, v, sim_b, body_ids, i, t);
    }

    // Simulations sim_f and sim_b will be automatically cleaned up when they drop out of scope.
    // This is why it's better to initialize two copies than to overwrite sim twice.
}

// *****************************************************************************
// Print vectors and orbital elements
// *****************************************************************************

// *****************************************************************************
void Simulation::print_vectors() const
{
    // Display the particles
    print("\n{:4s}: {:12s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s} \n",
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
    fmt::print("\n{:4s}: {:12s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s}: {:10s} \n",
          "   i", "body_name", "  a", "e", "  inc", "  Omega", "  omega", "  f", "  M");

    // By convention, particle index 0 is the Sun and has no primary, so start the loop at i=1
    for (int i=1; i<N(); i++)
    {
        // The body_id of the target and its primary
        const int32_t body_id = body_ids[i];
        const int32_t primary_body_id = primary_body_ids[i];
        // Only print elements if this particle has a real primary, i.e. skip the Sun
        if (!primary_body_id) {continue;}
        // Calculate the orbit
        Orbit orb = orbit(i);
        // Get the body name
        const string body_name = get_body_name(body_id);
        fmt::print("{:4d}: {:12s}: {:10.6f}: {:10.8f}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f}: {:+10.6f} \n", 
              i, body_name, orb.a, orb.e, orb.inc, orb.Omega, orb.omega, orb.f, orb.M);
    }
    print_newline();
}

// *****************************************************************************
} // Namespace reb
} // Namespace ks
