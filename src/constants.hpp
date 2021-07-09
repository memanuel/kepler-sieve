/** @file   constants.hpp
 *  @brief  Constants used more than once in Kepler Sieve project
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-09
 */

// *****************************************************************************
// Library dependencies
#include <cstdint>
    using std::int32_t;

// *****************************************************************************
namespace ks {
namespace cs {

/// Number of minutes in one day
constexpr int mpd = 1440;

/// Number of days in one minute
constexpr double dpm = 1.0 / mpd;

/// Start date for planets integration saved to database
constexpr int mjd0_db = 48000;

/// End date for planets integration saved to database
constexpr int mjd1_db = 63000;

/// Stride in minutes for planets integration in database
constexpr int stride_db_min = 5;

/// Number of times in database integration of planets
constexpr int N_t_db = (mjd1_db-mjd0_db)*mpd/stride_db_min+1;

// The body_id of the sun
constexpr int32_t body_id_sun = 10;


// *****************************************************************************
} // namespace cs
} // namespace ks
