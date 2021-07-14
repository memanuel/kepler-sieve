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
#include <vector>
    using std::vector;
#include <map>
    using std::map;
#include <stdexcept>
    using std::runtime_error;
    using std::invalid_argument;
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
    using ks::cs::N_body_planets;
    using ks::cs::N_body_planets_ex_sun;
    using ks::cs::body_ids_planets;
    using ks::cs::body_ids_planets_ex_sun;
    using ks::cs::mjd0_hrzn;
    using ks::cs::mjd1_hrzn;
    using ks::cs::stride_hrzn_planets;
    using ks::cs::stride_hrzn_de435;
#include "db_utils.hpp"
    using ks::wrap_string;
#include "astro_utils.hpp"
    using ks::get_primary_body_id;
    using ks::get_body_name;
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

/// A rebound particle in a simulation. Fields include m, x, y, z, vx, vy, vz.
using Particle = reb_particle;

/// A rebound orbit of a particle.  Fields include a, e, inc, Omega, omega, f, M, d, v, h, p, n, l, theta, T. 
using Orbit = reb_orbit;

// Function to integrate simulations
const auto integrate = reb_integrate;

// *****************************************************************************
// Utility functions

/// Convert a particle to a state vector
StateVector particle2state_vector(const Particle& p);

/// Calculate a rebound Orbit from a particle and a primary
Orbit particle2orbit(const Particle& target, const Particle& primary);

/// Convert a rebound Orbit to a Particle
Particle orbit2particle(const Orbit& orb, double m, const Particle& primary);

/// Convert a rebound Orbit to a StateVector
// StateVector orbit2particle(const Orbit& orb, double m, const Particle& primary);

// *****************************************************************************
} // Namespace reb
} // Namespace ks
