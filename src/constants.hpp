/** @file   constants.hpp
 *  @brief  Constants used more than once in Kepler Sieve project
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-09
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <cstdint>
    using std::int32_t;
#include <vector>
    using std::vector;
#include <numbers>

// *****************************************************************************
namespace ks {
namespace cs {

// *****************************************************************************
// Mathematical constants

/// Pi - Need we say more?
using std::numbers::pi;

/// Tau is an alias for 2pi
constexpr double tau = 2*pi;

// *****************************************************************************
// Time conversions between minutes and days

/// Number of minutes in one day
constexpr int mpd = 1440;

/// Number of days in one minute
constexpr double dpm = 1.0 / mpd;

/// Number of seconds in one day
constexpr int spd = mpd * 60;

// *****************************************************************************
// Physical constants

/// The gravitational constant in unit system (days, AU, Msun), i.e. AU^3 / MSun / day^2
constexpr double G {2.959122082855910945e-04};
// see rebound documentation for exact value        
// sim = make_sim_planets(epoch=59000)
// G_ = sim.G_

/// gravitational field strength: mu = G * (m0 + m1) in AU^3 / day^2
constexpr double m_sun {1.0};
constexpr double mu_sun {G * m_sun};
// here m0 = 1.0 (Sun) and we assume m1 is light, i.e. m1 = 0.0; works well for asteroids / test bodies.
// do NOT use this for orbital elements of planets!

// gravitational field strength of the Earth; for elements of the moon
constexpr double m_earth_moon {0.000003040432648};
constexpr double mu_earth_moon {G * m_earth_moon};

/// speed of light in meters per second - exact by definition (!)
/// see https://en.wikipedia.org/wiki/Speed_of_light
constexpr double c_m_s {299792458.0};

/// one astronomical unit (AU) in meters - exact by definition (!)
/// see https://en.wikipedia.org/wiki/Astronomical_unit
constexpr double au_m {149597870700};

/// speed of light in AU / day - exact by definition
/// quoted value (Google) 173.145
constexpr double c {c_m_s * static_cast<double>(spd) / au_m};

/// inverse speed of light in day / AU - exact by definition
/// quoted value (Google) 0.00577552
constexpr double c_inv {au_m / (c_m_s * static_cast<double>(spd))};

// sanity check on Google done by searching for "speed of light in AU per day"

// *****************************************************************************
// Time span for solar system integrations saved to database

/// Start date for planets integration saved to database
constexpr int mjd0_db = 48000;
/// End date for planets integration saved to database
constexpr int mjd1_db = 63000;

/// Stride in minutes for planets integration in database
constexpr int stride_db_min = 5;
/// Number of times in database integration of planets
constexpr int N_t_db = (mjd1_db-mjd0_db)*mpd/stride_db_min+1;

// *****************************************************************************
// Time span for solar system state vectors from horizons

/// Start date for planets integration saved to database
constexpr int mjd0_hrzn = 48000;
/// End date for planets integration saved to database
constexpr int mjd1_hrzn = 63000;
/// Stride in days for planets integration
constexpr int stride_hrzn_planets = 1;
/// Stride in days for DE435 integration
constexpr int stride_hrzn_de435 = 5;

// *****************************************************************************
// Body IDs for bodies with special handling in planets collection or used in examples

/// The body_id of the Sun
constexpr int32_t body_id_sun = 10;
/// The body_id of the Earth
constexpr int32_t body_id_earth = 399;
/// The body_id of the Moon
constexpr int32_t body_id_moon = 301;

// the body_ids of the planet barycenters
constexpr int32_t body_id_mercury_bc = 1;       // Mercury barycenter; same as Mercury
constexpr int32_t body_id_venus_bc = 2;         // Venus barycenter; Same as Venus
constexpr int32_t body_id_earth_bc = 3;         // Earth barycenter; skipped in favor of Earth + Moon
constexpr int32_t body_id_mars_bc = 4;          // Mars barycenter
constexpr int32_t body_id_jupiter_bc = 5;       // Jupiter barycenter
constexpr int32_t body_id_saturn_bc = 6;        // Saturn barycenter
constexpr int32_t body_id_uranus_bc = 7;        // Uranus barycenter
constexpr int32_t body_id_neptune_bc = 8;       // Neptune barycenter
constexpr int32_t body_id_pluto_bc = 9;         // Pluto barycenter
// placeholder for when a body_id is not applicable, e.g. primary of the Sun
constexpr int32_t body_id_null = 0;

/// The number of bodies in the planets collection
const int N_body_planets = 11;
/// The body_ids in the planets collection in Jacobi order: Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
const vector<int32_t> body_ids_planets = 
{body_id_sun, body_id_mercury_bc, body_id_venus_bc, body_id_earth, body_id_moon, body_id_mars_bc, 
body_id_jupiter_bc, body_id_saturn_bc, body_id_uranus_bc, body_id_neptune_bc, body_id_pluto_bc};

/// The number of bodies in planets ex sun collection (entities in planets collection that have orbital elements)
constexpr int N_body_planets_ex_sun = N_body_planets-1;
const vector<int32_t> body_ids_planets_ex_sun =
{body_id_mercury_bc, body_id_venus_bc, body_id_earth, body_id_moon, body_id_mars_bc, 
body_id_jupiter_bc, body_id_saturn_bc, body_id_uranus_bc, body_id_neptune_bc, body_id_pluto_bc};

/// The asteroid_id of Ceres
constexpr int32_t asteroid_id_ceres = 1;

// *****************************************************************************
} // namespace cs
} // namespace ks
