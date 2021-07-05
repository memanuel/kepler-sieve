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

// Local dependencies
#include "utils.hpp"
    using ks::sign;
    using ks::sqr;

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Data types for OrbitalElement, Position, Velocity, StateVector
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
// Functions for converting between orbital elements and state vectors
// *****************************************************************************

// Convert from orbital elements to a position vector
Position elt2pos(OrbitalElement& elt);

// Convert from orbital elements to a state vector
StateVector elt2vec(OrbitalElement& elt);

// *****************************************************************************
} // namespace ks