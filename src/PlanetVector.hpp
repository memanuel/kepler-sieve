/** @file   PlanetVector.hpp
 *  @brief  Class to load state vectors for all the planets. Supports interpolation of state vectors.
 *  See DB table KS.StateVectors_Planet and stored procedure KS.GetStateVectors_Planets.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-08
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <numeric>
    using std::lcm;
#include <string>
    using std::string;
    using std::to_string;
#include <fstream>
    using std::ifstream;
    using std::ofstream;
#include <vector>
    using std::vector;
#include <stdexcept>
    using std::domain_error;
#include <gsl/gsl_spline.h>
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "constants.hpp"
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::StateVector;
    using ks::StateVectorSpline;
#include "utils.hpp"
    using ks::flush_console;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "Timer.hpp"
    using ks::Timer;

// *****************************************************************************
namespace ks {

// *****************************************************************************
class PlanetVector
{
public:
    // ********************************************************************************************
    // Constructor and destructor
    // ********************************************************************************************

    /// Constructor takes time range, time step and allocates memory
    PlanetVector(int mjd0, int mjd1, int dt_min);

    /// Default Constructor - load data from file on disk
    PlanetVector();

    /// Constructor - take a DB connection
    PlanetVector(db_conn_type& conn);

    /// Destructor - delete manually created arrays
    ~PlanetVector();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************

    /// The number of bodies
    const int N_body;
    /// The number of times
    const int N_t;
    /// First date loaded (inclusive); an integer
    int mjd0;
    /// Last date loaded (inclusive); an integer
    int mjd1;
    /// Time step in minutes
    int dt_min;

    // ********************************************************************************************
    // Public Methods - Load and Save; build splines
    // ********************************************************************************************

    /// Load data from the database and construct interpolating splines
    void load(db_conn_type &conn);
    /// Load data from disk and construct interpolating splines
    void load();
    /// Save this object to disk
    void save() const;

    // Build GSL splines
    void build_splines();

    // ********************************************************************************************
    // Public Methods - Get Data
    // ********************************************************************************************

    /// Get the array of body IDs whose elements are in this table
    int32_t* get_body_id() const;

    /// Get the array of times; this is shared by all the bodies
    double* get_mjd() const;

    /// Calculate the interpolated position of the given body at time mjd in the BME frame
    Position interp_pos(int32_t body_id, double mjd) const;

    /// Calculate the interpolated state vector of the given body at time mjd in the BME frame
    StateVector interp_vec(int32_t body_id, double mjd) const;

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

    /// Number of rows of data
    int N_row; 
    /// The time_id of the first row corresponding to mjd0
    const int time_id0;
    /// The time_id of the last row corresponding to mjd1
    const int time_id1;

    // One shared array for the body_id; always the same here, the sun, 9 planets, and the moon
    int32_t* body_id;
    /// One shared array for the times as of which orbital elements apply
    double* mjd;

    // One array for each state vector component; array size is N_body * N_t
    // Array is laid out first by body, then by time (same order that SP returns data).
    // This is the required layout to spline each body vs. time.
    double* qx;
    double* qy;
    double* qz;
    double* vx;
    double* vy;
    double* vz;

    /// GSL spline interpolators for splined orbital elements
    StateVectorSplines vec_spline;
    /// Get a GSL cubic spline accelerator for lookups on orbital element splines
    gsl_interp_accel* acc;
    
    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    // Process a batch of rows
    void process_rows(db_conn_type& conn, int t0, int t1);
    /// Function to return the body_idx given a body_id; this is the row number of the body on a sorted list of bodies
    const int body_idx(int32_t body_id) const;
    /// Function to return the row index of a body_id; this is the index into the data arrays
    const int body_row(int32_t body_id) const;
    /// The row index of a time_id in a block for the same body
    const int time_row(int32_t time_id) const;
    /// The row index into the data array as a function of both body_id and time_id
    inline const int row_id(int32_t body_id, int32_t time_id) {return body_row(body_id) + time_row(time_id);};

    // Get an array (pointer to double) of state vector components given a body index number
    // Note the argument is body_idx (ranging from 0 to 9), NOT body_id!
    double* get_qx(int idx) const;
    double* get_qy(int idx) const;
    double* get_qz(int idx) const;
    double* get_vx(int idx) const;
    double* get_vy(int idx) const;
    double* get_vz(int idx) const;

    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
