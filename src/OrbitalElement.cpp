/** @file OrbitalElement.hpp
 *  @brief Class to store one Keplerian orbital element, with 7 entries a, e, inc, Omega, omega, f, M.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-02
 * 
 * See DB table KS.AsteroidElements and stored procedure KS.GetAsteroidElements.
 */

// *****************************************************************************
// Local dependencies
#include "OrbitalElement.hpp"

// Constants
// tau = 2pi
// constexpr double tau = 2.0 * pi;

// *****************************************************************************
namespace ks {

// The gravitational field strength: mu = G * (m0 + m1)
using ks::cs::mu;
using ks::cs::tau;

// *****************************************************************************
// Functions for converting one type of anomaly to another
// *****************************************************************************

// *****************************************************************************
/// Convert the eccentric anomaly E to the true anomaly f; e is the eccentricity
double anomaly_E2f(double E, double e)
{
    // SSD equation 2.46; solve for f in terms of E
    double tan_half_f = sqrt((1.0 + e) / (1.0 - e)) * tan(E/2.0);
    return atan(tan_half_f) * 2.0;
}

// *****************************************************************************
/// Convert the true anomaly f to the eccentric anomaly E; e is the eccentricity
double anomaly_f2E(double f, double e)
{
    // SSD equation 2.46; solve for E in terms of f
    double tan_half_f = tan(f * 0.5);
    double tan_half_E = tan_half_f / sqrt((1.0 + e) / (1.0 - e));
    return atan(tan_half_E) * 2.0;
}

// *****************************************************************************
/// Convert the true anomaly f to the eccentric anomaly E; e is the eccentricity
double anomaly_E2M(double E, double e)
{
    // SSD equation 2.52; this is Kepler's Equation
    return E - e * sin(E);
}

// *****************************************************************************
/// Convert the true anomaly f to the mean anomaly M; e is the eccentricity
double anomaly_f2M(double f, double e)
{
    // SSD equation 2.52; this is Kepler's Equation
    // Delegate to anomaly_f2E
    double E = anomaly_f2E(f, e);
    // Delegate to anomaly_E2M
    return anomaly_E2M(E, e);
}

// *****************************************************************************
/** Initial guess for E0 for iterative calculation of E from M
    See SSD equation 2.64 on page 36. */
double danby_guess(double M, double e)
{
    double k = 0.85;
    return M + sign(sin(M)) * k * e;
}

// *****************************************************************************
/** Perform one iteration of the Danby algorithm for computing E fomr M.
 *  See SSD equation 2.62 on page 36. */
double danby_iteration(double M, double e, double E)
{
    // The objective function that is set to zero using Newton-Raphson is
    // This is just the error term if we take Kepler's Equation and subtract M from both sides
    // This function will be zero at the correct value of E
    // f(E) = E - e Sin E - M
    
    // Save two intermediate arrays that are reused
    double eSinE = e * sin(E);
    double eCosE = e * cos(E);
    
    // Save the value of the function and its first three derivatives
    double f0 = E - eSinE - M;
    double f1 = 1.0 - eCosE;
    double f2 = eSinE;
    double f3 = eCosE;

    // The three delta adjustments; see SSD Equation 2.62
    double d1 = -f0 / f1;
    double d2 = -f0 / (f1 + 0.5*d1*f2);
    double d3 = -f0 / (f1 + 0.5*d2*f2 + (1.0/6.0)*sqr(d2)*f3);

    // Return E_next
    return E + d3;
}

// *****************************************************************************
/// Convert mean anomaly M to eccentric anomaly E using Danby iterations
double anomaly_M2E_danby(double M, double e, int n)
{
    // The initial guess
    double E = danby_guess(M=M, e=e);
    // n iterations of Danby algorithm
    for (int i=0; i<n; i++)
    {
        E = danby_iteration(M, e, E);
    }
    return E;
}

// *****************************************************************************
/// Delegate to anomaly_M2E with 3 iterations
double anomaly_M2E(double M, double e)
{
    return anomaly_M2E_danby(M, e, 3);
}

// *****************************************************************************
double anomaly_M2f(double M, double e)
{
    // Delegate to anomaly_M2E
    double E = anomaly_M2E(M, e);

    // Delegate to anomaly_E2f
    return anomaly_E2f(E, e);
}

// *****************************************************************************
// Functions for calculating additional elements e.g period T, mean motion n
// *****************************************************************************

// *****************************************************************************
double period(double a)
{
    // T^2 = 4pi^2 a^3 / mu (SSD 2.22)
    // Solve for T and obtain
    // T = 2pi * Sqrt(a^3 / mu)
    return tau * sqrt(cube(a) / mu);
}

// *****************************************************************************
double mean_motion(double a)
{
    // mu = n^2 * a^3 (SSD 2.26)
    // Solve this for n and obtain
    // n = sqrt(mu / a^3)
    return sqrt(mu / cube(a));
}

// *****************************************************************************
// Functions for converting between orbital elements and state vectors
// *****************************************************************************

// *****************************************************************************
Position elt2pos(double a, double e, double inc, double Omega, double omega, double f)
{
    // Distance from the center, r; SSD equation 2.20
    double r = a * (1.0 - sqr(e) ) / (1.0 + e * cos(f));

    // Intermediate resluts used for angular rotations
    // The angle in the elliptic plane, measured from the reference direction
    double theta = omega + f;

    // Trigonometric functions of the angles
    double cos_inc = cos(inc);
    double sin_inc = sin(inc);
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    // The Cartesian position coordinates; see SSD equation 2.122
    return Position 
    {
        .qx = r * (cos_Omega*cos_theta - sin_Omega*sin_theta*cos_inc),
        .qy = r * (sin_Omega*cos_theta + cos_Omega*sin_theta*cos_inc),
        .qz = r * (sin_theta*sin_inc)
    };
}

// *****************************************************************************
Position elt2pos(const OrbitalElement& elt)
    {return elt2pos(elt.a, elt.e, elt.inc, elt.Omega, elt.omega, elt.f);}

// *****************************************************************************
Position elt2pos_vec(double a, double e, double inc, double Omega, double omega, double f)
{
    // These equations taken directly from rebound library. 
    // File tools.c; function reb_tools_orbit_to_particle_err. Line 868.

    // The distance to the center
	double r = a*(1.0-e*e)/(1.0 + e*cos(f));

    // Components of the rotation matrix
	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
    // Initialize empty Position
    Position p;
    
	// Murray & Dermott Eq 2.122
	p.qx = r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.qy = r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.qz = r*(so*cf+co*sf)*si;

    return p;
}

// *****************************************************************************
Position elt2pos_vec(const OrbitalElement& elt)
    {return elt2pos_vec(elt.a, elt.e, elt.inc, elt.Omega, elt.omega, elt.f);}

// *****************************************************************************
StateVector elt2vec(double a, double e, double inc, double Omega, double omega, double f)
{
    // The equations used here were taken from the rebound library
    // The position calculations are equivalent to the ones above from SSD.
    // The velocity calculations are a bit more involved, and I did not see them with explicit equations in SSD.

    // These equations adapted directly from rebound library. 
    // File tools.c; function reb_tools_orbit_to_particle_err. Line 868.

    // The distance to the center
	double r = a*(1.0-e*e)/(1.0 + e*cos(f));
    // The speed of the body
	double v0 = sqrt(mu/a/(1.0-e*e)); // in this form it works for elliptical and hyperbolic orbits

    // Components of the rotation matrix
	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
    // Initialize empty StateVector
    StateVector s;
    
	// Murray & Dermott Eq 2.122
	s.qx = r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	s.qy = r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	s.qz = r*(so*cf+co*sf)*si;

	// Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
	s.vx = v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	s.vy = v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	s.vz = v0*((e+cf)*co*si - sf*si*so);
  
    // Return the assembled StateVector
    return s;
}

// *****************************************************************************
StateVector elt2vec(const OrbitalElement& elt)
{
    return elt2vec(elt.a, elt.e, elt.inc, elt.Omega, elt.omega, elt.f);
}

// *****************************************************************************
// Add a perturbation to an orbital element
// *****************************************************************************

// *****************************************************************************
OrbitalElement operator+ (const OrbitalElement& e1, const OrbitalElement& e2)
{
    // Add the two positions componentwise
    return OrbitalElement
    {
        .a      = e1.a      + e2.a,
        .e      = e1.e      + e2.e,
        .inc    = e1.inc    + e2.inc,
        .Omega  = e1.Omega  + e2.Omega,
        .omega  = e1.omega  + e2.omega,
        .f      = e1.f      + e2.f,
        .M      = e1.M      + e2.M
    };
}

// *****************************************************************************
// Print description of orbital elements
// *****************************************************************************

// *****************************************************************************
void print_orbital_element(const OrbitalElement& elt, bool header)
{
    if (header)
    {print("{:8s} : {:8s} : {:9s} : {:9s} : {:9s} : {:9s} : {:9s} \n", 
        "a", "e", "inc", "Omega", "omega", "f", "M");}
    else
    {print("{:8.6f} : {:8.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f} : {:9.4f} : {:9.4f}\n", 
        elt.a, elt.e, elt.inc, elt.Omega, elt.omega, elt.f, elt.M);}
}

// *****************************************************************************
void print_orbital_element_long(const OrbitalElement& elt)
{
    print("mjd      = {:9.3f}\n", elt.mjd);
    print("a        = {:9.6f}\n", elt.a);
    print("e        = {:9.6f}\n", elt.e);
    print("inc      = {:+8.6f}\n", elt.inc);
    print("Omega    = {:+8.6f}\n", elt.Omega);
    print("omega    = {:+8.6f}\n", elt.omega);
    print("f        = {:+8.6f}\n", elt.f);
    print("M        = {:+8.6f}\n", elt.M);
}

// *****************************************************************************
} // namespace ks
