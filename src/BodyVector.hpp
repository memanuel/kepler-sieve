/** @file BodyVector.hpp
 *  @brief Class to state vectors of one astronomical body in memory and perform interpolation on them.
 *  See DB tables KS.StateVectors_Sun, KS.StateVectors_Earth.
 *  See stored procedures Ks.GetStateVectors_Sun and KS.GetStateVectors_Earth.
 * 
 *  @author Michael S. Emanuel
 *  @date 2021-07-05
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Included libraries
#include <string>
    using std::string;
    using std::to_string;
#include <vector>
    using std::vector;
#include <numeric>
    using std::lcm;
#include <stdexcept>
    using std::domain_error;
#include <gsl/gsl_spline.h>

// Local dependencies
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::Position;
    // using ks::Velocity;
    using ks::StateVector;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

// *****************************************************************************
namespace ks {

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
class BodyVector
{
public:
    // ********************************************************************************************
    // Constructor and destructor
    // ********************************************************************************************

    /// Constructor takes time range, time step and allocates memory
    BodyVector(int mjd0, int mjd1, int dt_min);

    /// Constructor - take a DB connection and a body_name which must be one of "Sun", "Earth"
    BodyVector(db_conn_type& conn, string body_name);

    /// Destructor - delete manually created arrays
    ~BodyVector();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************

    // ********************************************************************************************
    // Public Methods
    // ********************************************************************************************

    /// Get the array of times; this is shared by all the asteroids
    double* get_mjd() const;

    /// Calculate the interpolated position of the given asteroid at time mjd
    Position interp_pos(double mjd);

    /// Calculate the interpolated state vector of the given asteroid at time mjd
    StateVector interp_vec(double mjd);

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

    /// The number of times
    const int N_t;
    /// First date loaded (inclusive); an integer divisible by time_step
    int mjd0;
    /// Last date loaded (inclusive); an integer divisible by time_step
    int mjd1;
    /// Time step in minutes
    int dt_min;

    /// One shared array for the times as of which orbital elements apply (every 4th day)
    double* mjd;

    // One array for state vector component; array size is N_t
    double* qx;
    double* qy;
    double* qz;
    double* vx;
    double* vy;
    double* vz;

    // GSL spline interpolators for splined orbital elements
    StateVectorSpline vec_spline;
    // Get a GSL cubic spline accelerator for lookups on orbital element splines
    gsl_interp_accel* acc;
    
    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    /// Load data from the database and construct interpolating splines
    void load(db_conn_type& conn, const string sp_name);
    // Function to return the row index given a time_id
    int32_t row(int32_t time_id) const;

    // Build GSL splines
    void build_splines();
    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
