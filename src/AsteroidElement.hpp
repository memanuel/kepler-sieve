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
/// Data contained in one AsteroidElement entry from output of SP KS.GetAsteroidElements.
struct AsteroidElementEntry
{
    /// Integer ID for the timestamp of these state vectors; FK to KS.IntegrationTime
    int32_t time_id;
    /// The asteroid whose orbital elements are described; FK to KS.Asteroid
    int32_t asteroid_id;
    /// The Modified Julian Date in the TDB (barycentric dynamical time) frame
    double mjd;
    /// The semimajor axis in AU
    double a;
    /// The eccentricity; dimensionless
    double e;
    /// The inclination in radians
    double inc;
    /// The longitude of the ascending node in radians
    double Omega;
    /// The argument of periapsis in radians
    double omega;
    /// The true anomaly in radians
    double f;
    /// The mean anomaly in radians
    double M;
};

// *****************************************************************************
/// Encapsulate all seven 2D arrays into one structure for code legibility
struct ElementArrays
{
    double* a;
    double* e;
    double* inc;
    double* Omega;
    double* omega;
    double* f;
    double* M;
};

// *****************************************************************************
class AsteroidElement
{
public:
    /// Native constructor - just allocate memory
    AsteroidElement(int n0, int N_ast, int mjd0, int N_t);

    /// Database constructor - delegate to native constructor for memory allocation, then write from DB
    AsteroidElement(db_conn_type &conn, int n0, int n1, int mjd0, int mjd1, bool progbar);

    /// Destructor - delete manually created arrays
    ~AsteroidElement();

    // Data
    /// First asteroid ID loaded (inclusive)
    const int n0;
    /// First date loaded (inclusive); an integer and divisible by 4
    const int mjd0;
    /// The number of asteroids
    const int N_ast;
    /// The number of times
    const int N_t;

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

    // Process a batch of rows
    void process_rows(db_conn_type& conn, int i0, int i1);
    // Function to return the asteroid index given an asteroid_id
    int32_t asteroid_idx(int32_t asteroid_id);
    // Function to return the row index given an asteroid_id
    int32_t row_idx(int32_t asteroid_id);

};

// *****************************************************************************
} // namespace ks
