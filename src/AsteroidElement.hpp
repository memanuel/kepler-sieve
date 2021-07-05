/** @file AsteroidSkyPatch.hpp
 *  @brief Class to load a batch of AsteroidSkyPatch entries in memory and look them up by SkyPatchID.
 *  See DB table AsteroidSkyPatch and stored procedure GetDetectionDire.
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
/** Encapsulate all seven vectors of GSL interpolators into one structure for code legibility
 *  One vector for each of seven orbital elements a, e, inc, Omega, omega, f, M.
 *  Each asteroid has one one entry in each vector. */
struct ElementSpline
{
    vector<gsl_spline*> a;
    vector<gsl_spline*> e;
    vector<gsl_spline*> inc;
    vector<gsl_spline*> Omega;
    vector<gsl_spline*> omega;
    vector<gsl_spline*> f;
    vector<gsl_spline*> M;
};

// *****************************************************************************
class AsteroidElement
{
public:
    // ********************************************************************************************
    // Constructor and destructor
    // ********************************************************************************************

    /// Constructor in terms of ranges delegate to native constructor for memory allocation
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

    /// First asteroid ID loaded (inclusive)
    int n0;
    /// Last asteroid ID loaded (exclusive)
    int n1;
    /// First date loaded (inclusive); an integer divisible by time_step
    int mjd0;
    /// Last date loaded (inclusive); an integer divisible by time_step
    int mjd1;
    /// Time step
    int dt;

    // ********************************************************************************************
    // Public Methods
    // ********************************************************************************************

    /// Load data from the database
    void load(db_conn_type &conn, bool progbar);

    /// Get the array of asteroid IDs whose elements are in this table
    int32_t* get_asteroid_id() const;
    /// Get the array of times; this is shared by all the asteroids
    double* get_mjd() const;

    // Get an array (pointer to double) of orbital elements given an asteroid_id
    double* get_a(int32_t asteroid_id) const;
    double* get_e(int32_t asteroid_id) const;
    double* get_inc(int32_t asteroid_id) const;
    double* get_Omega(int32_t asteroid_id) const;
    double* get_omega(int32_t asteroid_id) const;
    double* get_f(int32_t asteroid_id) const;
    double* get_M(int32_t asteroid_id) const; 

private:
    // ********************************************************************************************
    // Private Data Members
    // ********************************************************************************************

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

    // GSL objects for splined orbital elements
    gsl_interp_accel* acc;
    ElementSpline elt_spline;

    // ********************************************************************************************
    // Private Methods
    // ********************************************************************************************

    // Process a batch of rows
    void process_rows(db_conn_type& conn, int i0, int i1);
    // Function to return the asteroid index given an asteroid_id
    int32_t asteroid_idx(int32_t asteroid_id) const;
    // Function to return the row index given an asteroid_id
    int32_t asteroid_row(int32_t asteroid_id) const;

    // Free up GSL resources
    void gsl_free();

};

// // *****************************************************************************
// /// Encapsulate all seven 2D arrays into one structure for code legibility
// struct ElementArray
// {
//     double* a;
//     double* e;
//     double* inc;
//     double* Omega;
//     double* omega;
//     double* f;
//     double* M;
// };

// *****************************************************************************
} // namespace ks
