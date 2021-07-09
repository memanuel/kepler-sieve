/** @file OrbitalElement.hpp
 *  @brief Class to store one Keplerian orbital element, with 7 entries a, e, inc, Omega, omega, f, M.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-02
 * 
 * See DB table KS.AsteroidElements and stored procedure KS.GetAsteroidElements.
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <cmath>
#include <numbers>
    using std::numbers::pi;    
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
    using ks::cube;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Data types for OrbitalElement and a collection of orbital element splines
// *****************************************************************************

// *****************************************************************************
/**Data contained in one orbital element.
 * See DB table KS.Detection and stored procedure KS.GetDetectionObs.*/
struct OrbitalElement
{
    /// The Modified Julian Date in the TDB (barycentric dynamical time) frame
    double mjd;
    /// The semimajor axis in AU
    double a;
    /// The eccentricity; dimensionless
    double e;
    /// The inclination in radians
    double inc;
    /// The longitude of the ascending node in radians
    double Omega;
    /// The argument of periapsis in radians
    double omega;
    /// The true anomaly in radians
    double f;
    /// The mean anomaly in radians
    double M;
};

// *****************************************************************************
/** Encapsulate all seven vectors of GSL interpolators into one structure for code legibility
 *  One vector for each of seven orbital elements a, e, inc, Omega, omega, f, M.
 *  Each asteroid has one one entry in each vector. */
struct ElementSplines
{
    vector<gsl_spline*> a;
    vector<gsl_spline*> e;
    vector<gsl_spline*> inc;
    vector<gsl_spline*> Omega;
    vector<gsl_spline*> omega;
    vector<gsl_spline*> f;
    vector<gsl_spline*> M;
};

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
// Functions for converting one type of anomaly to another
// *****************************************************************************

// *****************************************************************************
/// Convert the eccentric anomaly E to the true anomaly f; e is the eccentricity
double anomaly_E2f(double E, double e);

// *****************************************************************************
/// Convert the true anomaly f to the eccentric anomaly E; e is the eccentricity
double anomaly_f2E(double f, double e);

// *****************************************************************************
/// Convert the true anomaly f to the eccentric anomaly E; e is the eccentricity
double anomaly_E2M(double E, double e);

// *****************************************************************************
/// Convert the true anomaly f to the mean anomaly M; e is the eccentricity
double anomaly_f2M(double f, double e);

// *****************************************************************************
/// Convert the mean anomaly M to the eccentric anomaly E using Danby iterations
double anomaly_M2E(double M, double e);

// *****************************************************************************
/// Convert the mean anomaly M to the true anomaly f using Danby iterations
double anomaly_M2f(double M, double e);

// *****************************************************************************
// Functions for calculating additional elements e.g period T, mean motion n
// *****************************************************************************
/// Calculate the orbital period, T, from the semimajor axis; uses mu for Sun
double period(double a);

/// Calculate the mean motion, n, from the semimajor axis; uses mu for Sun
double mean_motion(double a);

// *****************************************************************************
// Functions for converting between orbital elements and state vectors
// *****************************************************************************

/// Convert from orbital elements (six doubles) to a position vector. See SSD page 51, equation 2.122.
Position elt2pos(double a, double e, double inc, double Omega, double omega, double f);

/// Convert from orbital elements to a position vector. See SSD page 51, equation 2.122.
Position elt2pos(OrbitalElement& elt);

/// Convert from orbital elements (six doubles) to a state vector. See SSD page 51, equation 2.122.
StateVector elt2vec(double a, double e, double inc, double Omega, double omega, double f);

/// Convert from orbital elements to a state vector. See SSD page 51, equation 2.122.
StateVector elt2vec(OrbitalElement& elt);

// *****************************************************************************
/// Print a description of an orbital element
void print_orbital_element(OrbitalElement& elt);

// *****************************************************************************
// Utility function - norm of two positions
// *****************************************************************************
// Return the distance between two position vectors
double norm(Position p1, Position p2);

// Return the distance between two position vectors
double norm(Velocity v1, Velocity v2);

// Return a norm between two state vectors
double norm(StateVector v1, StateVector v2);

// *****************************************************************************
} // namespace ks
