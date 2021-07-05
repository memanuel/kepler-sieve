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

// *****************************************************************************
namespace ks {

// The gravitational constant in unit system (days, AU, Msun)
// see rebound documentation for exact value        
// sim = make_sim_planets(epoch=59000)
// G_ = sim.G_
constexpr double G_ = 2.959122082855910945e-04;

// mu is the gravitational field strength: mu = G * (m0 + m1)
// here m0 = 1.0 (Sun) and we assume m1 is light, i.e. m0 = 0.0
constexpr double mu = G_ * 1.0;

// *****************************************************************************
// Functions for working with orbital elements
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
Position elt2pos(OrbitalElement& elt)
{
    // Distance from the center, r; SSD equation 2.20
    double r = elt.a * (1.0 - sqr(elt.e) ) / (1.0 + elt.e * cos(elt.f));

    // Intermediate resluts used for angular rotations
    // The angle in the elliptic plane, measured from the reference direction
    double theta = elt.omega + elt.f;

    // Trigonometric functions of the angles
    double cos_inc = cos(elt.inc);
    double sin_inc = sin(elt.inc);
    double cos_Omega = cos(elt.Omega);
    double sin_Omega = sin(elt.Omega);
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
StateVector elt2vec(OrbitalElement& elt)
{
    // The equations used here were taken from the rebound library
    // The position calculations are equivalent to the ones above from SSD.
    // The velocity calculations are a bit more involved, and I did not see them with explicit equations in SSD.

    // sine and cosine of the angles inc, Omega, omega, and f
    double ci = cos(elt.inc);
    double si = sin(elt.inc);
    double cO = cos(elt.Omega);
    double sO = sin(elt.Omega);
    double co = cos(elt.omega);
    double so = sin(elt.omega);
    double cf = cos(elt.f);
    double sf = sin(elt.f);

    // Distance from center
    double one_minus_e2 = 1.0 - sqr(elt.e);
    double one_plus_e_cos_f = 1.0 + elt.e*cos(elt.f);
    double r = elt.a * one_minus_e2 / one_plus_e_cos_f;

    // Current speed
    double v0 = sqrt(mu / elt.a / one_minus_e2);

    // Position
    // qx = r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    // qy = r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    // qz = r*(so*cf+co*sf)*si
    // the term cos_omega*cos_f - sin_omega*sin_f appears 2 times
    // the term sin_omega*cos_f + cos_omega*sin_f appears 3 times
    double cocf_sosf = co*cf - so*sf;
    double socf_cosf = so*cf + co*sf;
    double qx = r*(cO*cocf_sosf - sO*socf_cosf*ci);
    double qy = r*(sO*cocf_sosf + cO*socf_cosf*ci);
    double qz = r*(socf_cosf*si);
    
    // Velocity
    // vx = v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    // vy = v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    // vz = v0*((e+cf)*co*si - sf*si*so)
    // The term e+cf appears three times
    double epcf = elt.e + cf;
    // The term cocO appears twice
    double cocO = co*cO;
    // The term cosO appears twice
    double cosO = co*sO;
    // The term so*sO appears twice
    double sosO = so*sO;
    // The terms socO appears twice
    double socO = so*cO;

    // Simplified expression for velocity with substitutions
    double vx = v0*(epcf*(-ci*cosO - socO) - sf*(cocO - ci*sosO));
    double vy = v0*(epcf*(ci*cocO - sosO)  - sf*(cosO + ci*socO));
    double vz = v0*(epcf*co*si - sf*si*so);
    
    // Wrap the components into one StateVector object
    return StateVector
    {
        .qx = qx,
        .qy = qy,
        .qz = qz,
        .vx = vx,
        .vy = vy,
        .vz = vz
    };
}

// *****************************************************************************
double norm(Position p1, Position p2)
{
    double dq2 = sqr(p1.qx-p2.qx) + sqr(p1.qy-p2.qy) + sqr(p1.qz-p2.qz);
    return sqrt(dq2);
}

// *****************************************************************************
double norm(Velocity v1, Velocity v2)
{
    double dv2 = sqr(v1.vx-v2.vx) + sqr(v1.vy-v2.vy) + sqr(v1.vz-v2.vz);
    return sqrt(dv2);
}

// *****************************************************************************
double norm(StateVector v1, StateVector v2)
{
    double dq2 = sqr(v1.qx-v2.qx) + sqr(v1.qy-v2.qy) + sqr(v1.qz-v2.qz);
    double dv2 = sqr(v1.vx-v2.vx) + sqr(v1.vy-v2.vy) + sqr(v1.vz-v2.vz);
    return sqrt(dq2 + dv2);
}

// *****************************************************************************
} // namespace ks
