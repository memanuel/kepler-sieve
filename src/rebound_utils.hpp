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
#include "StateVector.hpp"
    using ks::StateVector;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "PlanetElement.hpp"
    using ks::PlanetElement;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "constants.hpp"
    using ks::cs::G;

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

    /// Add one rebound particle object to a simulation
    void add_particle(const Particle& p);

    /// Add one particle to a simulation given its state vector and mass
    void add_particle(const StateVector& s, double m);

    /// Add a test particle to a simulation given its state vector; MUST finish adding massive particles first
    void add_test_particle(const StateVector& s);

    /// Integrate to the given time
    void integrate(double t) {reb_integrate(prs, t);}

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

    /// Get a particle by index in particles array
    const Particle& particle(int i) const {return prs->particles[i];}

    /// Get a state vector by index in particles array
    const StateVector state_vector(int i) const;

    /// Print a simulation summary
    void print() const;

private:
    /// Pointer to the underlying rebound simulation structure; prs is "pointer to rebound simulation"
    reb_simulation* prs;

    /// Add one massive particle to a simulation given its state vector and mass
    void add_particle_impl(const StateVector& s, double m);
};

// *****************************************************************************

/// Create a rebound simulation with the planets; initialization from splined planet vectors at epoch.
Simulation make_sim_planets(const PlanetVector& pv, double epoch);

/// Create a rebound simulation with the planets; initialization from splined planet elements at epoch.
Simulation make_sim_planets(const PlanetElement& pe, double epoch);

// *****************************************************************************
} // Namespace reb
} // Namespace ks
