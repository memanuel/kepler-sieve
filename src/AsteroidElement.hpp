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
#include "constants.hpp"
    using ks::cs::mu_sun;
    using ks::cs::mpd;
#include "utils.hpp"
    using ks::flush_console;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "Timer.hpp"
    using ks::Timer;
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::StateVector;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::OrbitalAngle;
    using ks::OrbitalElementSplines;
    using ks::OrbitalAngleSplines;
    using ks::elt2pos;
    using ks::elt2vec;
    using ks::elt2pos_vec;
#include "BodyVector.hpp"
    using ks::BodyVector;

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
    const int32_t* get_asteroid_id() const {return asteroid_id;}

    /// Get the array of times; this is shared by all the asteroids
    const double* get_mjd() const {return mjd;}

    /// Calculate the interpolated orbital elements of the given asteroid at time mjd
    const OrbitalElement interp_elt(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated position of the given asteroid at time mjd in the heliocentric frame
    const Position interp_pos_hel(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated state vector of the given asteroid at time mjd in the heliocentric frame
    const StateVector interp_vec_hel(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated position of the given asteroid at time mjd in the BME frame
    const Position interp_pos(int32_t asteroid_id, double mjd) const;

    /// Calculate the interpolated state vector of the given asteroid at time mjd in the BME frame
    const StateVector interp_vec(int32_t asteroid_id, double mjd) const;

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

    /// The number of rows of data
    const int N_row;

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

    // One array for each orbital element; array size is N_body * N_t
    // Array is laid out first by body, then by time (same order that SP returns data).
    // This is the required layout to spline each asteroid vs. time.

    // Seven standard orbital elements
    double* elt_a;
    double* elt_e;
    double* elt_inc;
    double* elt_Omega;
    double* elt_omega;
    double* elt_f;
    double* elt_M;

    // Ten orbital angles
    double* elt_cos_inc;
    double* elt_sin_inc;
    double* elt_cos_Omega;
    double* elt_sin_Omega;
    double* elt_cos_omega;
    double* elt_sin_omega;
    double* elt_cos_f;
    double* elt_sin_f;
    double* elt_cos_M;
    double* elt_sin_M;

    /// GSL spline interpolators for orbital elements
    OrbitalElementSplines elt_spline;
    /// GSL spline interpolators for orbital angles
    OrbitalAngleSplines ang_spline;
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
    const int asteroid_idx(int32_t asteroid_id) const {return asteroid_id - n0;}
    // Function to return the row index given an asteroid_id
    const int asteroid_row(int32_t asteroid_id) const {return asteroid_idx(asteroid_id)*N_t;}

    // Seven traditional orbital elements
    const double* get_a(        int idx) const {return elt_a         + N_t*idx;}
    const double* get_e(        int idx) const {return elt_e         + N_t*idx;}
    const double* get_inc(      int idx) const {return elt_inc       + N_t*idx;}
    const double* get_Omega(    int idx) const {return elt_Omega     + N_t*idx;}
    const double* get_omega(    int idx) const {return elt_omega     + N_t*idx;}
    const double* get_f(        int idx) const {return elt_f         + N_t*idx;}
    const double* get_M(        int idx) const {return elt_M         + N_t*idx;}

    // Five pairs of cosine / sine of angle orbital elements
    const double* get_cos_inc(  int idx) const {return elt_cos_inc   + N_t*idx;}
    const double* get_sin_inc(  int idx) const {return elt_sin_inc   + N_t*idx;}
    const double* get_cos_Omega(int idx) const {return elt_cos_Omega + N_t*idx;}
    const double* get_sin_Omega(int idx) const {return elt_sin_Omega + N_t*idx;}
    const double* get_cos_omega(int idx) const {return elt_cos_omega + N_t*idx;}
    const double* get_sin_omega(int idx) const {return elt_sin_omega + N_t*idx;}
    const double* get_cos_f(    int idx) const {return elt_cos_f     + N_t*idx;}
    const double* get_sin_f(    int idx) const {return elt_sin_f     + N_t*idx;}
    const double* get_cos_M(    int idx) const {return elt_cos_M     + N_t*idx;}
    const double* get_sin_M(    int idx) const {return elt_sin_M     + N_t*idx;}

    // Build GSL splines
    void build_splines();
    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
