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
// Functions for manupulating positions, velocities and and state vectors
// *****************************************************************************

// *****************************************************************************
/// Add two position vectors
Position operator+ (const Position& p1, const Position& p2);

/// Add two velocity vectors
Velocity operator+ (const Velocity& p1, const Velocity& p2);

/// Add two state vectors
StateVector operator+ (const StateVector& p1, const StateVector& p2);

/// Multiply a veocity vector by a scalar
Velocity operator* (const Velocity& v, const double alpha);
Velocity operator* (const double alpha, const Velocity& v);

/// Extract the posiiton from a state vector
Position sv2pos(StateVector& sv);

/// Extract the velocity from a state vector
Velocity sv2vel(StateVector& sv);

// *****************************************************************************
// Utility function - norm of two positions
// *****************************************************************************

/// Return the distance between two position vectors
double norm(Position p1, Position p2);

/// Return the distance between two position vectors
double norm(Velocity v1, Velocity v2);

/// Return a norm between two state vectors
double norm(StateVector v1, StateVector v2);

// *****************************************************************************
// Print positions and state vectlors
// *****************************************************************************

// *****************************************************************************
/// Print a one line description of an orbital element
void print_state_vector(StateVector& sv, bool header=false);

/// Print a multi-line description of an orbital element
void print_state_vector_long(StateVector& sv);

// *****************************************************************************
} // namespace ks
