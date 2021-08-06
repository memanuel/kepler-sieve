/** @file   Direction.hpp
 *  @brief  Class to store direction from an observer to a target in the solar system.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-08-06
 * 
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
// #include <cmath>
#include <fmt/format.h>
    using fmt::print;
#include <gsl/gsl_spline.h>
    // gsl_spline

// Local dependencies
#include "StateVector.hpp"
    using ks::Position;
    using ks::StateVector;
#include "constants.hpp"
    using ks::cs::c_inv;    
#include "utils.hpp"
    using ks::sign;
    using ks::sqr;
#include "astro_utils.hpp"
    using ks::dist2rad;
    using ks::dist2deg;
    using ks::dist2sec;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Data types for Direction
// *****************************************************************************

// *****************************************************************************
/// Direction from an observer to a target in the solar system
struct Direction
{
    /// x coordinate on the unit sphere in the barcycentric mean ecliptic (BME) frame
    double ux;
    /// y coordinate on the unit sphere in the barcycentric mean ecliptic (BME) frame
    double uy;
    /// z coordinate on the unit sphere in the barcycentric mean ecliptic (BME) frame
    double uz;
};

/// The results of one observation - a direction and a distance
struct ObservationResult
{
    /// The direction from observer to target in the BME frame
    Direction u;
    /// The distance from observer to target in AU
    double r;
};

// *****************************************************************************
/// Encapsulate three GSL interpolators, one for each component, into one structure
struct DirectionSpline
{
    gsl_spline* const ux;
    gsl_spline* const uy;
    gsl_spline* const uz;
};

// *****************************************************************************
// Add and subtract direction vectors
// *****************************************************************************

/// Add two direction vectors; result no longer on unit sphere
Direction operator+ (const Direction& u1, const Direction& u2);

/// Subtract two direction vectors; result no longer on unit sphere
Direction operator- (const Direction& u1, const Direction& u2);

// *****************************************************************************
// Norm and distance of direction vectors
// *****************************************************************************

/// Return the norm of a direction vector
double norm(const Direction& u);

/// Return the distance between two direction vectors
double dist(const Direction& u1, const Direction& u2);

/// Return the distance between two direction vectors in radians
inline double dist_rad(const Direction& u1, const Direction& u2) {return dist2rad(dist(u1, u2));}

/// Return the distance between two direction vectors in degrees
inline double dist_deg(const Direction& u1, const Direction& u2) {return dist2deg(dist(u1, u2));}

/// Return the distance between two direction vectors in arc seconds
inline double dist_sec(const Direction& u1, const Direction& u2) {return dist2sec(dist(u1, u2));}

// *****************************************************************************
// Check if two directions are close
// *****************************************************************************

/// Test if two position vectors are close within the given tolerance
bool is_close(const Direction& u1, const Direction& u2, double tol);

// *****************************************************************************
// Calculate a direction from position of observer and state vector of target
// *****************************************************************************

/// Create a direction by normalizing a position (difference) vector
ObservationResult normalize_pos(const Position& pos);

/// Calculate observed direction and distance between two objects with Newtonian light model
ObservationResult observe(const Position& q_obs, const StateVector& s_tgt);

/// Calculate only the direction between two objects with Newtonian light model
Direction observe_dir(const Position& q_obs, const StateVector& s_tgt);

// *****************************************************************************
// Print a direction
// *****************************************************************************

// *****************************************************************************
/// Print a column headers for one line description of a direction
void print_direction_headers(const string prefix="");

/// Print a one line description of a direction in scientific notation
void print_direction(const Direction& u, const string prefix="");

/// Print a one line description of a state vector in scientific notation
void print_direction_sci(const Direction& u, const string prefix="");

// *****************************************************************************
} // namespace ks
