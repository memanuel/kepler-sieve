/** @file   StateVector.hpp
 *  @brief  Class to store one position, velocity or state vector (position and velocity)
 *          of a single object in the solar system.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-02
 * 
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <cmath>
#include <vector>
    using std::vector;
#include <fmt/format.h>
    using fmt::print;
#include <gsl/gsl_spline.h>
    // gsl_spline

// Local dependencies
#include "utils.hpp"
    using ks::sign;
    using ks::sqr;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Data types for Position, Velocity, StateVector
// *****************************************************************************

// *****************************************************************************
/// Position of an object in the solar system
struct Position
{
    /// Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame
    double qx;
    /// Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame
    double qy;
    /// Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame
    double qz;
};

// *****************************************************************************
/// Velocity of an object in the solar system
struct Velocity
{
    /// Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame
    double vx;
    /// Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame
    double vy;
    /// Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame
    double vz;
};

// *****************************************************************************
/// State Vector of an object in the solar system
struct StateVector
{
    /// Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame
    double qx;
    /// Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame
    double qy;
    /// Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame
    double qz;
    /// Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame
    double vx;
    /// Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame
    double vy;
    /// Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame
    double vz;
};

// *****************************************************************************
/// Encapsulate six GSL interpolators, one for each component, into one structure
struct StateVectorSpline
{
    gsl_spline* qx;
    gsl_spline* qy;
    gsl_spline* qz;
    gsl_spline* vx;
    gsl_spline* vy;
    gsl_spline* vz;
};

// *****************************************************************************
/** Encapsulate six vectors of GSL interpolators into one structure for code legibility
 *  One vector for each component qx, qy, qz, vx, vy, vz.
 *  Each body has one one entry in each vector. */
struct StateVectorSplines
{
    vector<gsl_spline*> qx;
    vector<gsl_spline*> qy;
    vector<gsl_spline*> qz;
    vector<gsl_spline*> vx;
    vector<gsl_spline*> vy;
    vector<gsl_spline*> vz;
};


// *****************************************************************************
} // namespace ks
