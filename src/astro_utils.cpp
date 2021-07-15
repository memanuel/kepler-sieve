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
// Local dependencies
#include "astro_utils.hpp"

// *****************************************************************************
// Names used
using ks::SolarSystemBody;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Dates and times
// *****************************************************************************

// *****************************************************************************
double jd_to_mjd(double jd)
    {return jd - modified_julian_offset;}

// *****************************************************************************
double mjd_to_jd(double mjd)
    {return mjd + modified_julian_offset;}

// *****************************************************************************
// Angles and distances
// *****************************************************************************

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
// Solar system bodies
// *****************************************************************************

// *****************************************************************************
const int32_t get_primary_body_id(int32_t body_id)
{
    if (body_id == body_id_sun)
        {return body_id_null;}
    else if (body_id == body_id_moon)
        {return body_id_earth;}
    else 
        {return body_id_sun;}
}

// *****************************************************************************
const char* get_body_name(SolarSystemBody body)
{
    switch (body)
    {
        // Most common bodies first
        case SolarSystemBody::sun           : return "Sun";
        case SolarSystemBody::earth         : return "Earth";
        case SolarSystemBody::moon          : return "Moon";
        // Remaining bodies in Jacobi order
        case SolarSystemBody::mercury_bc    : return "Mercury_bc";
        case SolarSystemBody::venus_bc      : return "Venus_bc";
        case SolarSystemBody::mars_bc       : return "Mars_bc";
        case SolarSystemBody::jupiter_bc    : return "Jupiter_bc";
        case SolarSystemBody::saturn_bc     : return "Saturn_bc";
        case SolarSystemBody::uranus_bc     : return "Uranus_bc";
        case SolarSystemBody::neptune_bc    : return "Neptune_bc";
        case SolarSystemBody::pluto_bc      : return "Pluto_bc";
    }
    // Impossible to get here b/c all enum cases covered above.
    // Put something to suppress compiler warning
    throw invalid_argument("get_body_name - bad body enum!");
}

// *****************************************************************************
const int get_body_idx(SolarSystemBody body)
{
    // Bodies are laid out in Jacobi order 0, 1, 2, 301, 399, 4, 5, 6, 7, 8, 9
    // This corresponds to Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto.
    switch (body)
    {
        case SolarSystemBody::sun           : return 0;
        case SolarSystemBody::mercury_bc    : return 1;
        case SolarSystemBody::venus_bc      : return 2;
        case SolarSystemBody::earth         : return 3;
        case SolarSystemBody::moon          : return 4;
        case SolarSystemBody::mars_bc       : return 5;
        case SolarSystemBody::jupiter_bc    : return 6;
        case SolarSystemBody::saturn_bc     : return 7;
        case SolarSystemBody::uranus_bc     : return 8;
        case SolarSystemBody::neptune_bc    : return 9;
        case SolarSystemBody::pluto_bc      : return 10;
    }
    // Impossible to get here b/c all enum cases covered above.
    // Put something to suppress compiler warning
    throw invalid_argument("get_body_name - bad body enum!");
}

// *****************************************************************************
const int get_body_idx(int32_t body_id)
    {return get_body_idx(static_cast<SolarSystemBody>(body_id));}

// *****************************************************************************
const string get_body_name(int32_t body_id)
{
    switch (body_id)
    {
        // Put the most common bodies first
        case body_id_sun:           return "Sun";
        case body_id_earth:         return "Earth";
        case body_id_moon:          return "Moon";
        // Put the remaining bodies in Jacobi order
        case body_id_mercury_bc:    return "Mercury_bc";
        case body_id_venus_bc:      return "Venus_bc";
        case body_id_mars_bc:       return "Mars_bc";
        case body_id_jupiter_bc:    return "Jupiter_bc";
        case body_id_saturn_bc:     return "Saturn_bc";
        case body_id_uranus_bc:     return "Uranus_bc";
        case body_id_neptune_bc:    return "Neptune_bc";
        case body_id_pluto_bc:      return "Pluto_bc";
        // If an unknown body_id passed, 
        default:                    return "Unknown";        
    }   // switch
}

// *****************************************************************************
const char* get_body_name(SolarSystemBody_bv body)
    {return get_body_name(static_cast<SolarSystemBody>(body) );}

// *****************************************************************************
const int32_t get_body_id(const char* body_name)
{
    // Convert the body_name 
    // Put the most common bodies first
    if (strcmp(body_name, "Sun")==0)            {return body_id_sun;}
    if (strcmp(body_name, "Earth")==0)          {return body_id_earth;}
    if (strcmp(body_name, "Moon")==0)           {return body_id_moon;}
    // Put the remaining bodies in Jacobi order
    if (strcmp(body_name, "Mercury bc")==0)     {return body_id_mercury_bc;}
    if (strcmp(body_name, "Venus bc")==0)       {return body_id_venus_bc;}
    if (strcmp(body_name, "Marx bc")==0)        {return body_id_mars_bc;}
    if (strcmp(body_name, "Jupiter bc")==0)     {return body_id_jupiter_bc;}
    if (strcmp(body_name, "Saturn bc")==0)      {return body_id_saturn_bc;}
    if (strcmp(body_name, "Uranus bc")==0)      {return body_id_uranus_bc;}
    if (strcmp(body_name, "Neptune bc")==0)     {return body_id_neptune_bc;}
    if (strcmp(body_name, "Pluto bc")==0)       {return body_id_pluto_bc;}
    // If we get here, it was an error
    print("get_body_id - bad body_name! got {:s}, must be in planets collection.\n", string(body_name));
    throw invalid_argument("get_body_id - bad body_name!");
}

// *****************************************************************************
}; // namespace ks
