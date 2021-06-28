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

#pragma once

// *****************************************************************************
// Included libraries
#include <cmath>
#include <string>

// *****************************************************************************
// Names used

// *****************************************************************************
// Put all functions into namespace ks
namespace ks {

// *****************************************************************************
// Set the grid size for the sky patch at compile time
namespace astro_utils{
    /// Pi
    constexpr double pi = M_PI;
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
/// Convert julian date (jd) to a modified julian date (mjd)
double jd_to_mjd(double jd);

/// Convert a modified julian date (mjd) to a julian date (jd)
double mjd_to_jd(double mjd);

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
} // Namespace ks
