/** @file   rebound_utils.hpp
 *  @brief  Utilities for working with rebound library.
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
// #include "PlanetElement.hpp"
//     using ks::PlanetElement;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "constants.hpp"
    using ks::cs::G;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Alias commonly used Rebound types

/// A rebound particle in a simulation
using Particle = reb_particle;

/// A rebound orbit of a particle
using Orbit = reb_orbit;

/// A rebound simulation
using Simulation = reb_simulation;

// *****************************************************************************

/// Create an empty rebound simulation with default configuration: units, G, integrator.
Simulation* make_sim();

/// Create a rebound simulation with the planets as of the specified date
Simulation* make_sim_planets(double epoch);

// *****************************************************************************
} // Namespace ks
