/** @file   PlanetElement.hpp
 *  @brief  Class to load orbital elements for all the planets.
 *  Supports interpolation of orbital elements and state vectors.
 *  See DB table KS.OrbitalElements_Planets and stored procedure KS.GetElements_Planets.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-08
 */
// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <cmath>
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
    using std::invalid_argument;
    using std::runtime_error;
#include <gsl/gsl_spline.h>
#include <fmt/format.h>
    using fmt::print;

// *****************************************************************************
// Local dependencies
#include "constants.hpp"
    using ks::cs::tau;
    using ks::cs::mpd;
    using ks::cs::dpm;
    using ks::cs::body_id_sun;
    using ks::cs::body_id_earth;
    using ks::cs::body_id_moon;
    using ks::cs::N_body_planets_ex_sun;
    using ks::cs::body_ids_planets_ex_sun;
    using ks::cs::m_sun;
    using ks::cs::G;
    using ks::cs::mjd0_db;
    using ks::cs::mjd1_db;
    using ks::cs::stride_db_min;
    using ks::cs::N_t_db;
#include "astro_utils.hpp"
    using ks::get_primary_body_id;
    using ks::get_body_idx;
#include "StateVector.hpp"
    using ks::Position;
    using ks::Velocity;
    using ks::StateVector;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::OrbitalAngle;
    using ks::OrbitalElementSplines;
    using ks::OrbitalAngleSplines;
    using ks::eltsp2elt;
    using ks::elt2pos;
    using ks::elt2vec;
    using ks::elt2pos_vec;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
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

    /// Constructor takes time range, time step and allocates memory
    PlanetElement(int mjd0, int mjd1, int dt_min);

    /// Default Constructor - load data from file on disk
    PlanetElement();

    /// Constructor - take a DB connection
    PlanetElement(db_conn_type& conn);

    /// Destructor - delete manually created arrays
    ~PlanetElement();

    // ********************************************************************************************
    // Public Data Members
    // ********************************************************************************************

    /// The number of bodies
    const int N_body;
    /// The number of times
    const int N_t;
    /// First date loaded (inclusive); an integer
    const int mjd0;
    /// Last date loaded (inclusive); an integer
    const int mjd1;
    /// Time step in minutes
    const int dt_min;

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
    const int32_t* get_body_id() const {return body_id;}

    /// Get the array of gravitational field strength mu
    const double* get_mu() const {return mu;}

    /// Get the array of times; this is shared by all the bodies
    const double* get_mjd() const {return mjd;}

    /// Calculate the interpolated orbital elements of the given body at time mjd
    const OrbitalElement interp_elt(int32_t body_id, double mjd) const;

    /// Calculate the interpolated position of the given body at time mjd relative to its primary
    const Position interp_pos_rel(int32_t body_id, double mjd) const;

    /// Calculate the interpolated state vector of the given body at time mjd relative to its primary
    const StateVector interp_vec_rel(int32_t body_id, double mjd) const;

    /// Calculate the interpolated position of the given body at time mjd in the BME frame
    const Position interp_pos(int32_t body_id, double mjd) const;

    /// Calculate the interpolated state vector of the given body at time mjd in the BME frame
    const StateVector interp_vec(int32_t body_id, double mjd) const;

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
    /// One shared array for mu of each interaction
    double* mu;

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
    // double* elt_cos_f;
    // double* elt_sin_f;
    // double* elt_cos_M;
    // double* elt_sin_M;

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
    void process_rows(db_conn_type& conn, int t0, int t1);
    /// Function to return the body_idx given a body_id; this is the row number of the body on a sorted list of bodies
    const int body_idx(int32_t body_id) const;
    /// Function to return the row index of a body_id; this is the index into the data arrays
    const int body_row(int32_t body_id) const;
    /// The row index of a time_id in a block for the same body
    const int time_row(int32_t time_id) const;
    /// The row index into the data array as a function of both body_id and time_id
    inline const int row_id(int32_t body_id, int32_t time_id) {return body_row(body_id) + time_row(time_id);};

    // Get an array (pointer to double) of orbital elements given a body index number
    // Note the argument is body_idx (ranging from 0 to 9), NOT body_id!

    // Seven traditional orbital elements
    const double* get_a(        int idx) const {return elt_a         + N_t*idx;}
    const double* get_e(        int idx) const {return elt_e         + N_t*idx;}
    // const double* get_inc(      int idx) const {return elt_inc       + N_t*idx;}
    // const double* get_Omega(    int idx) const {return elt_Omega     + N_t*idx;}
    // const double* get_omega(    int idx) const {return elt_omega     + N_t*idx;}
    // const double* get_f(        int idx) const {return elt_f         + N_t*idx;}
    const double* get_M(        int idx) const {return elt_M         + N_t*idx;}

    // Five pairs of cosine / sine of angle orbital elements
    const double* get_cos_inc(  int idx) const {return elt_cos_inc   + N_t*idx;}
    const double* get_sin_inc(  int idx) const {return elt_sin_inc   + N_t*idx;}
    const double* get_cos_Omega(int idx) const {return elt_cos_Omega + N_t*idx;}
    const double* get_sin_Omega(int idx) const {return elt_sin_Omega + N_t*idx;}
    const double* get_cos_omega(int idx) const {return elt_cos_omega + N_t*idx;}
    const double* get_sin_omega(int idx) const {return elt_sin_omega + N_t*idx;}
    // const double* get_cos_f(    int idx) const {return elt_cos_f     + N_t*idx;}
    // const double* get_sin_f(    int idx) const {return elt_sin_f     + N_t*idx;}
    // const double* get_cos_M(    int idx) const {return elt_cos_M     + N_t*idx;}
    // const double* get_sin_M(    int idx) const {return elt_sin_M     + N_t*idx;}

    /// Calculate the interpolated orbital elements of the given body at time mjd; idx is the body index
    const OrbitalElement interp_elt_by_idx(int idx, double mjd) const;

    // Free up GSL resources
    void gsl_free();
};

// *****************************************************************************
} // namespace ks
