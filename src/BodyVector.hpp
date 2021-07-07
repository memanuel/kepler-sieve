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
#include <fstream>
    using std::ifstream;
    using std::ofstream;
#include <vector>
    using std::vector;
#include <numeric>
    using std::lcm;
#include <algorithm>
    using std::count;
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
    BodyVector(int mjd0, int mjd1, int dt_min, string body_name);

    /// Constructor - load data from file on disk
    BodyVector(string body_name);

    /// Constructor - take a DB connection and a body_name which must be one of "Sun", "Earth"
    BodyVector(db_conn_type& conn, string body_name);

    /// Destructor - delete manually created arrays
    ~BodyVector();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************
    const string body_name;

    // ********************************************************************************************
    // Public Methods
    // ********************************************************************************************

    /// Get the array of times used to build the spline (times with known state vectors)
    double* get_mjd() const;

    /// Calculate the interpolated position of the body at time mjd
    Position interp_pos(double mjd) const;

    /// Calculate the interpolated state vector of the body at time mjd
    StateVector interp_vec(double mjd) const;

    /// Load data from the database and construct interpolating splines
    void load(db_conn_type& conn);
    /// Load data from disk and construct interpolating splines
    void load();
    /// Save this object to disk
    void save() const;

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

    /// The number of times
    const int N_t;
    /// First date loaded (inclusive); an integer
    int mjd0;
    /// Last date loaded (inclusive); an integer
    int mjd1;
    /// Time step in minutes
    int dt_min;

    /// One shared array for the times when splined elements are available; size N_t
    double* mjd;

    // One array for state vector component; size N_t
    double* qx;
    double* qy;
    double* qz;
    double* vx;
    double* vy;
    double* vz;

    // GSL spline interpolators for splined orbital elements
    StateVectorSpline vec_spline;
    // Get a GSL cubic spline accelerator for lookups on orbital element splines; shared by all six splines
    gsl_interp_accel* acc;
    
    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    // Function to return the row index given a time_id
    int32_t row(int32_t time_id) const;
    /// Name of database stored procedure to fetch state vectors for this body
    const string sp_name_from_body() const;
    /// Name of filename with serialized data for this body
    const string file_name_from_body() const;

    /// Process a batch of rows from the databse SP
    void process_rows(db_conn_type& conn, int mjd0, int mjd1);
    // Build GSL splines
    void build_splines();
    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
/// Helper function - build and save vectors for Sun and Earth
void save_vectors(string body_name);

// *****************************************************************************
} // namespace ks
