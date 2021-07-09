/** @file AsteroidElement.hpp
 *  @brief Class to load a batch of orbital elements for a block of asteroids into memory.
 *  Supports interpolation of orbital elements and state vectors.
 *  See DB table KS.AsteroidElements and stored procedure KS.GetAsteroidElements.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-02
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
#include <gsl/gsl_spline.h>

// Local dependencies
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::StateVector;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;    
    using ks::ElementSplines;
    using ks::elt2pos;
    using ks::elt2vec;
#include "BodyVector.hpp"
    using ks::BodyVector;
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
class AsteroidElement
{
public:
    // ********************************************************************************************
    // Constructor and destructor
    // ********************************************************************************************

    /// Constructor takes a range of asteroids and dates
    AsteroidElement(int n0, int n1, int mjd0, int mjd1, int dt);

    /// Destructor - delete manually created arrays
    ~AsteroidElement();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************

    /// The number of asteroids
    const int N_ast;
    /// The number of times
    const int N_t;

    // ********************************************************************************************
    // Public Methods
    // ********************************************************************************************

    /// Load data from the database and construct interpolating splines
    void load(db_conn_type &conn, bool progbar);

    /// Get the array of asteroid IDs whose elements are in this table
    int32_t* get_asteroid_id() const;

    /// Get the array of times; this is shared by all the asteroids
    double* get_mjd() const;

    /// Calculate the interpolated orbital elements of the given asteroid at time mjd
    OrbitalElement interp_elt(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated position of the given asteroid at time mjd in the heliocentric frame
    Position interp_pos_hel(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated state vector of the given asteroid at time mjd in the heliocentric frame
    StateVector interp_vec_hel(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated position of the given asteroid at time mjd in the BME frame
    Position interp_pos(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated state vector of the given asteroid at time mjd in the BME frame
    StateVector interp_vec(int32_t asteroid_id, double mjd) const;

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

    /// First asteroid ID loaded (inclusive)
    int n0;
    /// Last asteroid ID loaded (exclusive)
    int n1;
    /// First date loaded (inclusive); an integer divisible by time_step
    int mjd0;
    /// Last date loaded (inclusive); an integer divisible by time_step
    int mjd1;
    /// Time step in days
    int dt;

    // One shared array for the distinct asteroid IDs (typically a sequence, possibly with some holes)
    int32_t* asteroid_id;
    /// One shared array for the times as of which orbital elements apply (every 4th day)
    double* mjd;

    // One array for each orbital element; array size is N_ast * N_t
    // Array is laid out first by asteroid, then by time (same order that SP returns data).
    // This is the required layout to spline each asteroid vs. time.
    double* elt_a;
    double* elt_e;
    double* elt_inc;
    double* elt_Omega;
    double* elt_omega;
    double* elt_f;
    double* elt_M;

    /// GSL spline interpolators for splined orbital elements
    ElementSplines elt_spline;
    /// Get a GSL cubic spline accelerator for lookups on orbital element splines
    gsl_interp_accel* acc;
    /// Interpolated state vectors of the Sun; used to calculate state vectors in the BME frame
    const BodyVector bv_sun;
    
    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    // Process a batch of rows
    void process_rows(db_conn_type& conn, int i0, int i1);
    // Function to return the asteroid index given an asteroid_id
    const int asteroid_idx(int32_t asteroid_id) const;
    // Function to return the row index given an asteroid_id
    const int asteroid_row(int32_t asteroid_id) const;

    // Get an array (pointer to double) of orbital elements given an asteroid_id
    double* get_a(int32_t asteroid_id) const;
    double* get_e(int32_t asteroid_id) const;
    double* get_inc(int32_t asteroid_id) const;
    double* get_Omega(int32_t asteroid_id) const;
    double* get_omega(int32_t asteroid_id) const;
    double* get_f(int32_t asteroid_id) const;
    double* get_M(int32_t asteroid_id) const;

    // Build GSL splines
    void build_splines();
    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
