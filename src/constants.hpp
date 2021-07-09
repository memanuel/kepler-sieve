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

// *****************************************************************************
namespace ks {
namespace cs {

// *****************************************************************************
// Time conversions between minutes and days

/// Number of minutes in one day
constexpr int mpd = 1440;

/// Number of days in one minute
constexpr double dpm = 1.0 / mpd;

// *****************************************************************************
// Physical constants

/// The gravitational constant in unit system (days, AU, Msun), i.e. AU^3 / MSun / day^2
constexpr double G_ = 2.959122082855910945e-04;
// see rebound documentation for exact value        
// sim = make_sim_planets(epoch=59000)
// G_ = sim.G_

/// gravitational field strength: mu = G * (m0 + m1) in AU^3 / day^2
constexpr double mu = G_ * 1.0;
// here m0 = 1.0 (Sun) and we assume m1 is light, i.e. m1 = 0.0

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

/// The body_id of the Sun
constexpr int32_t body_id_sun = 10;

/// The body_id of the Earth
constexpr int32_t body_id_earth = 399;

/// The asteroid_id of Ceres
constexpr int32_t asteroid_id_ceres = 1;

// *****************************************************************************
} // namespace cs
} // namespace ks