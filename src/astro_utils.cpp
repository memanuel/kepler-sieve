/** @file astro_utils.cpp
 *  @brief Implmentation of astronomy utilities.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

/*****************************************************************************
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/

// *****************************************************************************
// Included files
#include "astro_utils.hpp"

// *****************************************************************************
namespace ks {

// *****************************************************************************
double jd_to_mjd(double jd)
    {return jd - modified_julian_offset;}

// *****************************************************************************
double mjd_to_jd(double mjd)
    {return mjd + modified_julian_offset;}

// *****************************************************************************
double dist2rad(double s)
    {return asin(0.5*s) * 2.0;}

// *****************************************************************************
double rad2dist(double s_rad)
    {return sin(0.5*s_rad) * 2.0;}

// *****************************************************************************
double dist2deg(double s)
{
    double s_rad = dist2rad(s);
    return s_rad*deg_per_rad;
}

// *****************************************************************************
double deg2dist(double s_deg)
{
    double s_rad = s_deg*rad_per_deg;
    return rad2dist(s_rad);
}

// *****************************************************************************
double dist2sec(double s)
{
    double s_deg = dist2deg(s);
    return s_deg*3600.0;
}

// *****************************************************************************
double sec2dist(double s_sec)
{
    double s_deg = s_sec / 3600.0;
    return deg2dist(s_deg);
}

// *****************************************************************************
}; // namespace
