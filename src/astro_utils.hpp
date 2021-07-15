/** @file astro_utils.cpp
 *  @brief Utilities relating to astronomy.
 *
 * jd_to_mjd(jd) converts a julian date to a modified julian date
 * mjd_to_jd(mjd) converts a modified julian date to a julian date
 * dist2rad(s) convers a Cartesian distance on the unit sphere to radians
 * rad2dist(s_rad) convers a distance on the unit in radians sphere to Cartesian distance
 * dist2deg(s) convers a Cartesian distance on the unit sphere to degrees
 * deg2dist(s_deg) convers a distance on the unit in degrees sphere to Cartesian distance
 * dist2sec(s) convers a Cartesian distance on the unit sphere to arc seoncds
 * sec2dist(s_sec) convers a distance on the unit in arc seconds sphere to Cartesian distance
 * 
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */
// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <cmath>
#include <numbers>
#include <stdexcept>
    using std::invalid_argument;
#include <string>
    using std::string;
#include <cstring>
#include <fmt/format.h>
    using fmt::print;
    using fmt::format;

// *****************************************************************************
// Local dependencies
#include "constants.hpp"
    using ks::cs::body_id_sun;
    using ks::cs::body_id_mercury_bc;
    using ks::cs::body_id_venus_bc;
    using ks::cs::body_id_earth;
    using ks::cs::body_id_moon;
    using ks::cs::body_id_mars_bc;
    using ks::cs::body_id_jupiter_bc;
    using ks::cs::body_id_saturn_bc;
    using ks::cs::body_id_uranus_bc;
    using ks::cs::body_id_neptune_bc;
    using ks::cs::body_id_pluto_bc;
    using ks::cs::body_id_null;

// *****************************************************************************
// Put all functions into namespace ks
namespace ks {

// *****************************************************************************
// Set the grid size for the sky patch at compile time
namespace astro_utils{
    /// Pi
    using std::numbers::pi;
    /// Number of degrees in one radian
    constexpr double deg_per_rad = 180.0 / pi;
    /// Number of radians in one degree
    constexpr double rad_per_deg= pi / 180.0;

    /// The modified Julian Date is 2400000.5 less than the Julian Base Number
    constexpr double modified_julian_offset = 2400000.5;
    /// Number of seconds in one day
    constexpr double day2sec= 24.0 * 3600.0;
    /// Number of days in one second
    constexpr double sec2day = 1.0 / day2sec;
}

// Use these names
using astro_utils::pi;
using astro_utils::deg_per_rad;
using astro_utils::rad_per_deg;
using astro_utils::modified_julian_offset;

// *****************************************************************************
// Dates and times
// *****************************************************************************

/// Convert julian date (jd) to a modified julian date (mjd)
double jd_to_mjd(double jd);

/// Convert a modified julian date (mjd) to a julian date (jd)
double mjd_to_jd(double mjd);

// *****************************************************************************
// Angles and distances
// *****************************************************************************

// *****************************************************************************
/// Convert a distance on the unit sphere in [0, 2] to radians in [0, pi]
double dist2rad(double s);

/// Convert a distance on the unit sphere in [0, 2] to radians in [0, pi]
double rad2dist(double s_rad);

/// Convert a cartesian distance on unit sphere in [0, 2] to degrees in [0, 180]
double dist2deg(double s);

/// Convert a distance on unit sphere from degrees in [0, 180] to cartesian distance in [0, 2]
double deg2dist(double s_deg);

/// Convert a cartesian distance on unit sphere in [0, 2] to arc seconds in [0, 180*3600]
double dist2sec(double s);

/// Convert a distance on unit sphere from degrees in [0, 180] to cartesian distance in [0, 2]
double sec2dist(double s_sec);

// *****************************************************************************
// Solar system bodies
// *****************************************************************************

/// Enumeration of solar system bodies keyed by body_id; limited to planets collection
enum class SolarSystemBody : int32_t
{
    sun             = body_id_sun,
    mercury_bc      = body_id_mercury_bc,
    venus_bc        = body_id_venus_bc,
    earth           = body_id_earth,
    moon            = body_id_moon,
    mars_bc         = body_id_mars_bc,
    jupiter_bc      = body_id_jupiter_bc,
    saturn_bc       = body_id_saturn_bc,
    uranus_bc       = body_id_uranus_bc,
    neptune_bc      = body_id_neptune_bc,
    pluto_bc        = body_id_pluto_bc
};

/// Enumeration of solar system bodies keyed by body_id with cached BodyVector inputs
enum class SolarSystemBody_bv : int32_t
{
    sun             = body_id_sun,
    earth           = body_id_earth,
};

/// Get the primary body_id of a body in the solar system
const int32_t get_primary_body_id(int32_t body_id);

/// Get body index (row number) of a body in the planets collection given its SolarSystemBody enum
const int get_body_idx(SolarSystemBody body);

/// Get body index (row number) of a body in the planets collection given its body_id
const int get_body_idx(int32_t body_id);

/// Get the name of a solar system body from its SolarSystemBody enum
const char* get_body_name(SolarSystemBody body);

/// Get the name of a solar system body from its SolarSystemBody_bv enum
const char* get_body_name(SolarSystemBody_bv body);

/// Get the name of a solar system body from its body_id; limited to planets collection
const string get_body_name(int32_t body_id);

/// Get the ID of a solar system body from its name; limited to planets collection
const int32_t get_body_id(const char* body_name);

// *****************************************************************************
} // Namespace ks
