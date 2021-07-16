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
#include "astro_utils.hpp"
    using ks::SolarSystemBody_bv;
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
    AsteroidElement(int32_t n0, int32_t n1, int mjd0, int mjd1, int dt);

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
    void load(db_conn_type& conn, bool progbar);

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
    int32_t n0;
    /// Last asteroid ID loaded (exclusive)
    int32_t n1;

    /// First date loaded (inclusive); an integer divisible by time_step
    const int mjd0;
    /// Last date loaded (inclusive); an integer divisible by time_step
    const int mjd1;
    /// Time step in days
    const int dt;

    // One shared array for the distinct asteroid IDs (typically a sequence, possibly with some holes)
    int32_t* const asteroid_id;
    /// One shared array for the times as of which orbital elements apply (every 4th day)
    double* const mjd;

    // One array for each orbital element; array size is N_body * N_t
    // Array is laid out first by body, then by time (same order that SP returns data).
    // This is the required layout to spline each asteroid vs. time.

    // Seven standard orbital elements
    double* const elt_a;
    double* const elt_e;
    double* const elt_inc;
    double* const elt_Omega;
    double* const elt_omega;
    double* const elt_f;
    double* const elt_M;

    // Orbital angles splined via cosine / sine
    double* const elt_cos_inc;
    double* const elt_sin_inc;
    double* const elt_cos_Omega;
    double* const elt_sin_Omega;
    double* const elt_cos_omega;
    double* const elt_sin_omega;

    /// GSL spline interpolators for orbital elements
    OrbitalElementSplines elt_spline;
    /// Get a GSL cubic spline accelerator for lookups on orbital element splines
    gsl_interp_accel* const acc;
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

    // Orbital elements splined directly
    const double* get_a(        int idx) const {return elt_a         + N_t*idx;}
    const double* get_e(        int idx) const {return elt_e         + N_t*idx;}
    const double* get_M(        int idx) const {return elt_M         + N_t*idx;}

    // Orbital angles cosine / sine
    const double* get_cos_inc(  int idx) const {return elt_cos_inc   + N_t*idx;}
    const double* get_sin_inc(  int idx) const {return elt_sin_inc   + N_t*idx;}
    const double* get_cos_Omega(int idx) const {return elt_cos_Omega + N_t*idx;}
    const double* get_sin_Omega(int idx) const {return elt_sin_Omega + N_t*idx;}
    const double* get_cos_omega(int idx) const {return elt_cos_omega + N_t*idx;}
    const double* get_sin_omega(int idx) const {return elt_sin_omega + N_t*idx;}

    // Build GSL splines
    void build_splines();
    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
