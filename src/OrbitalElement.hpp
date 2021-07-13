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
#include "constants.hpp"
#include "utils.hpp"
    using ks::sign;
    using ks::sqr;
    using ks::cube;
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::StateVector;
    
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
/// Orbital element with each angles replaced by its cos and sin; this form is splinable.
struct OrbitalElementSplinable
{
    /// The Modified Julian Date in the TDB (barycentric dynamical time) frame
    double mjd;
    /// The semimajor axis in AU
    double a;
    /// The eccentricity; dimensionless
    double e;
    /// Cosine of the inclination in radians
    double cos_inc;
    /// Sine of the inclination in radians
    double sin_inc;
    /// Cosine of the longitude of the ascending node in radians
    double cos_Omega;
    /// Sine of the longitude of the ascending node in radians
    double sin_Omega;
    /// Cosine of the argument of periapsis in radians
    double cos_omega;
    /// Sine of the argument of periapsis in radians
    double sin_omega;
    /// Cosine of the true anomaly in radians
    double cos_f;
    /// Sine of the true anomaly in radians
    double sin_f;
    /// Cosine of the mean anomaly in radians
    double cos_M;
    /// Sine of the mean anomaly in radians
    double sin_M;
};


// *****************************************************************************
/** Encapsulate vectors of GSL interpolators into one structure for code legibility
 *  One vector for orbital elements a, e.
 *  Angle elements inc, Omega, omega, f, M each get two splines, for cosine and sine respectively.
 *  Each asteroid has one one entry in each vector. */
struct ElementSplines
{
    vector<gsl_spline*> a;
    vector<gsl_spline*> e;
    vector<gsl_spline*> cos_inc;
    vector<gsl_spline*> sin_inc;
    vector<gsl_spline*> cos_Omega;
    vector<gsl_spline*> sin_Omega;
    vector<gsl_spline*> cos_omega;
    vector<gsl_spline*> sin_omega;
    vector<gsl_spline*> cos_f;
    vector<gsl_spline*> sin_f;
    vector<gsl_spline*> cos_M;
    vector<gsl_spline*> sin_M;
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

/// Calculate the orbital period, T, from the semimajor axis; requires gravitational field mu
double period(double a);

/// Calculate the mean motion, n, from the semimajor axis; requires gravitational field mu
double mean_motion(double a, double mu);

// *****************************************************************************
// Functions for converting between orbital elements and state vectors
// *****************************************************************************

/// Convert from orbital elements (six doubles) to a position vector. See SSD page 51, equation 2.122.
Position elt2pos(double a, double e, double inc, double Omega, double omega, double f);

/// Convert from orbital elements to a position vector. See SSD page 51, equation 2.122.
Position elt2pos(const OrbitalElement& elt);

/// Convert from orbital elements (six doubles) to a position vector. This version has slightly more steps but exactly matches elt2vec.
Position elt2pos_vec(double a, double e, double inc, double Omega, double omega, double f);

/// Convert from orbital elements to a position vector. This version has slightly more steps but exactly matches elt2vec.
Position elt2pos_vec(const OrbitalElement& elt);

/// Convert from orbital elements (six doubles) to a state vector. See SSD page 51, equation 2.122.
StateVector elt2vec(double a, double e, double inc, double Omega, double omega, double f, double mu);

/// Convert from orbital elements to a state vector. See SSD page 51, equation 2.122.
StateVector elt2vec(const OrbitalElement& elt, double mu);

// *****************************************************************************
// Add a perturbation to an orbital element
// *****************************************************************************

/// Add the components of two orbital elements (typically a real element and a small shift)
OrbitalElement operator+ (const OrbitalElement& e1, const OrbitalElement& e2);

// *****************************************************************************
// Functions for converting between standard and splinable orbital elements
// *****************************************************************************

/// Convert a standard OrbitalElement elt to a SplinableOrbitalElement elts
OrbitalElementSplinable elt2elts(const OrbitalElement& elt);

/// Convert a SplinableOrbitalElement selt to a standard OrbitalElement elt
OrbitalElement elts2elt(const OrbitalElementSplinable& elts);

// *****************************************************************************
// Print description of orbital elements
// *****************************************************************************

/// Print a one line description of an orbital element
void print_orbital_element(const OrbitalElement& elt, bool header=false);

/// Print a multi-line description of an orbital element
void print_orbital_element_long(const OrbitalElement& elt);

// *****************************************************************************
} // namespace ks
