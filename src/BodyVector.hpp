/** @file   BodyVector.hpp
 *  @brief  Class to state vectors of one astronomical body in memory and perform interpolation on them.
 *  See DB tables KS.StateVectors_Sun, KS.StateVectors_Earth.
 *  See stored procedures Ks.GetStateVectors_Sun and KS.GetStateVectors_Earth.
 * 
 *  @author Michael S. Emanuel
 *  @date   2021-07-05
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Included libraries
#include <string>
    using std::string;
    using std::string_view;
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
    using std::invalid_argument;
    using std::runtime_error;
#include <gsl/gsl_spline.h>

// Local dependencies
#include "constants.hpp"
    using ks::cs::mpd;
    using ks::cs::dpm;
    using ks::cs::mjd0_db;
    using ks::cs::mjd1_db;
    using ks::cs::stride_db_min;
    using ks::cs::body_id_sun;
    using ks::cs::body_id_earth;
#include "astro_utils.hpp"
    using ks::SolarSystemBody_bv;
    using ks::get_body_name;
#include "StateVector.hpp"
    using ks::Position;
    using ks::StateVector;
    using ks::StateVectorSpline;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

// *****************************************************************************
namespace ks {

// *****************************************************************************

// *****************************************************************************
class BodyVector
{
public:
    // ********************************************************************************************
    // Constructor and destructor
    // ********************************************************************************************

    /// Constructor takes time range, time step and allocates memory
    BodyVector(SolarSystemBody_bv ssb, int mjd0, int mjd1, int dt_min);

    /// Constructor gives caller option to load data from disk
    BodyVector(SolarSystemBody_bv ssb, int mjd0, int mjd1, int dt_min, bool load_);

    /// Constructor - load data from file on disk
    BodyVector(SolarSystemBody_bv ssb);

    /// Constructor - take a DB connection and a body_name which must be one of "Sun", "Earth"
    BodyVector(SolarSystemBody_bv ssb, db_conn_type& conn);

    /// Destructor - delete manually created arrays
    ~BodyVector();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************

    /// The body whose vectors are saved as a SolarSystemBody enum
    const SolarSystemBody_bv ssb;
    /// The name of this body as a const char*
    const char* body_name;

    // ********************************************************************************************
    // Public Methods
    // ********************************************************************************************

    /// Get the array of times used to build the spline (times with known state vectors)
    const double* get_mjd() const;

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

    // ********************************************************************************************
    // Data Members
    // ********************************************************************************************
public:
    /// The number of times
    const int N_t;
    /// First date loaded (inclusive); an integer
    const int mjd0;
    /// Last date loaded (inclusive); an integer
    const int mjd1;
    /// Time step in minutes
    const int dt_min;

private:
    /// One shared array for the times when splined elements are available; size N_t
    double* const mjd;

    // One array for each state vector component; size N_t
    double* const qx;
    double* const qy;
    double* const qz;
    double* const vx;
    double* const vy;
    double* const vz;

    // GSL spline interpolators for splined orbital elements
    StateVectorSpline vec_spline;
    // Get a GSL cubic spline accelerator for lookups on orbital element splines; shared by all six splines
    gsl_interp_accel* const acc;
    
    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    // Function to return the row index given a time_id
    int32_t row(int32_t time_id) const;
    /// Name of database stored procedure to fetch state vectors for this body
    const char* sp_name_() const;
    /// Name of filename with serialized data for this body
    const char* file_name_() const;

    /// Process a batch of rows from the databse SP
    void process_rows(db_conn_type& conn, int mjd0, int mjd1);
    // Build GSL splines
    void build_splines();
    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
