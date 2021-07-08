/** @file   PlanetElement.hpp
 *  @brief  Class to load orbital elements for all the planets.
 *  Supports interpolation of orbital elements and state vectors.
 *  See DB table KS.AsteroidElements and stored procedure KS.GetAsteroidElements.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-08
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <string>
    using std::string;
    using std::to_string;
#include <vector>
    using std::vector;
#include <stdexcept>
    using std::domain_error;
#include <gsl/gsl_spline.h>

// Local dependencies
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::ElementSpline;
    using ks::Position;
    using ks::Velocity;
    using ks::StateVector;
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
class PlanetElement
{
public:
    // ********************************************************************************************
    // Constructor and destructor
    // ********************************************************************************************

    /// Constructor
    PlanetElement(int mjd0, int mjd1, int dt_min);

    /// Destructor - delete manually created arrays
    ~PlanetElement();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************

    /// The number of bodies
    const int N_body;
    /// The number of times
    const int N_t;

    // ********************************************************************************************
    // Public Methods - Load and Save
    // ********************************************************************************************

    /// Load data from the database and construct interpolating splines
    void load(db_conn_type &conn);
    /// Load data from disk and construct interpolating splines
    void load();
    /// Save this object to disk
    void save() const;

    // ********************************************************************************************
    // Public Methods - Get Data
    // ********************************************************************************************

    /// Get the array of body IDs whose elements are in this table
    int32_t* get_body_id() const;

    /// Get the array of times; this is shared by all the bodies
    double* get_mjd() const;

    /// Calculate the interpolated orbital elements of the given body at time mjd
    OrbitalElement interp_elt(int32_t body_id, double mjd) const;

    /// Calculate the interpolated position of the given body at time mjd in the heliocentric frame
    Position interp_pos_hel(int32_t body_id, double mjd) const;

    /// Calculate the interpolated state vector of the given body at time mjd in the heliocentric frame
    StateVector interp_vec_hel(int32_t body_id, double mjd) const;

    /// Calculate the interpolated position of the given body at time mjd in the BME frame
    Position interp_pos(int32_t body_id, double mjd) const;

    /// Calculate the interpolated state vector of the given body at time mjd in the BME frame
    StateVector interp_vec(int32_t body_id, double mjd) const;

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

    /// First date loaded (inclusive); an integer divisible by time_step
    int mjd0;
    /// Last date loaded (inclusive); an integer divisible by time_step
    int mjd1;
    /// Time step in minutes
    int dt_min;

    // One shared array for the body_id; always the same here, the sun, 9 planets, and the moon
    int32_t* body_id;
    /// One shared array for the times as of which orbital elements apply
    double* mjd;

    // One array for each orbital element; array size is N_body * N_t
    // Array is laid out first by body, then by time (same order that SP returns data).
    // This is the required layout to spline each asteroid vs. time.
    double* elt_a;
    double* elt_e;
    double* elt_inc;
    double* elt_Omega;
    double* elt_omega;
    double* elt_f;
    double* elt_M;

    /// GSL spline interpolators for splined orbital elements
    ElementSpline elt_spline;
    /// Get a GSL cubic spline accelerator for lookups on orbital element splines
    gsl_interp_accel* acc;
    /// Interpolated state vectors of the Sun; used to calculate state vectors in the BME frame
    const BodyVector bv;
    
    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    // Process a batch of rows
    void process_rows(db_conn_type& conn, int t0, int t1);
    // Function to return the body_idx given a body_id; this is the row number of the body on a sorted list of bodies
    const int body_idx(int32_t body_id) const;
    // Function to return the row index of a body_id; this is the index into the data arrays
    const int body_row(int32_t body_id) const;

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
