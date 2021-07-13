/** @file   rebound_utils.hpp
 *  @brief  Utilities for working with rebound library for gravitational integrations.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-10
 * 
 */
/* ****************************************************************************/
#pragma once

// *****************************************************************************
// Library dependencies
#include <string>
    using std::string;
#include <stdexcept>
    using std::runtime_error;
#include <fmt/format.h>
    using fmt::print;

// C Library dependencies
extern "C" {
#include <rebound.hpp>
}

// *****************************************************************************
// Local depdendencies
#include "constants.hpp"
    using ks::cs::G;
    using ks::cs::body_id_sun;
    using ks::cs::body_id_earth;
    using ks::cs::body_id_moon;
#include "StateVector.hpp"
    using ks::StateVector;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "PlanetElement.hpp"
    using ks::PlanetElement;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;

// *****************************************************************************
namespace ks {
namespace reb {

// *****************************************************************************
// Alias commonly used Rebound types

/// A rebound particle in a simulation
using Particle = reb_particle;

/// A rebound orbit of a particle
using Orbit = reb_orbit;

// Function to integrate simulations
const auto integrate = reb_integrate;


// *****************************************************************************
// Utility functions

/// Convert a particle to a state vector
StateVector particle2state_vector(const Particle& p);

/// Calculate a rebound Orbit from a particle and a primary
Orbit particle2orbit(const Particle& target, const Particle& primary);

// *****************************************************************************
// Class Simulation
// *****************************************************************************

/// Class to encapsulate a rebound simulation structure
class Simulation
{
public:
    // *****************************************************************************
    // Constructor and destructor
    // *****************************************************************************

    /// Default constructor creates an empty rebound simulation with default configuration: units, G, integrator
    Simulation();

    /// Copy constructor clones all the particles
    Simulation(const Simulation& s);

    /// Destructor must call the reb_free function
    ~Simulation();

    // *****************************************************************************
    // Modify simulation: set time, add particles
    // *****************************************************************************

    /// Set the time parameter
    void set_time(double t) {prs->t = t;}

    /// Add one particle to a simulation given its state vector and mass
    void add_particle(const StateVector& s, double m, int32_t body_id=0);

    /// Add a test particle to a simulation given its state vector; MUST finish adding massive particles first
    void add_test_particle(const StateVector& s, int32_t body_id=0);

    /// Integrate to the given time
    void integrate(double t) {reb_integrate(prs, t);}

    /// Add mapping of particles to primary
    void set_primary(int32_t body_id, int32_t primary_body_id) {primary_body_map[body_id]=primary_body_id;}

    // *****************************************************************************
    // Get simulation properties: time, particle count, particles
    // *****************************************************************************

    /// Get the time
    const double t() const {return prs->t;}

    /// Get the number of particles
    const int N() const {return prs->N;}

    /// Get the number of active particles
    const int N_active() const {return prs->N_active;}

    /// Get the number of active particles
    const int N_test() const {return prs->N - prs->N_active;}

    /// Get the index number i of a particle from its body_id
    const int body_idx(int32_t body_id) {return body_map.at(body_id);};

    /// Get a particle by index in particles array
    const Particle& particle(int i) const {return prs->particles[i];}

    /// Get a state vector by index in particles array
    const StateVector state_vector(int i) const;

    /// Print the state vectors of all particles in the simulation
    void print_vectors() const;

    /// Print the orbital elements of all particles in the simulation
    void print_elements() const;

private:
    /// Pointer to the underlying rebound simulation structure; prs is "pointer to rebound simulation"
    reb_simulation* prs;

    /// Map with key=body_id, value=body_idx
    std::map<int32_t, int> body_map;

    /// Map with key=body_id, value=body_id_primary
    std::map<int32_t, int32_t> primary_body_map;

    /// Add one massive particle to a simulation given its state vector and mass
    void add_particle_impl(const StateVector& s, double m, int32_t body_id);
};

// *****************************************************************************

/// Create a rebound simulation with the planets; initialization from splined planet vectors at epoch.
Simulation make_sim_planets(const PlanetVector& pv, double epoch);

/// Create a rebound simulation with the planets; initialization from splined planet elements at epoch.
Simulation make_sim_planets(const PlanetElement& pe, double epoch);

// *****************************************************************************
} // Namespace reb
} // Namespace ks
