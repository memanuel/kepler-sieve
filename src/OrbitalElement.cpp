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
} // namespace ks
