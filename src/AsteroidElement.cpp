/** @file AsteroidElement.cpp
 *  @brief Implmentation of AsteroidElement class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-02
 */

// *****************************************************************************
// Included files
#include "AsteroidElement.hpp"

// *****************************************************************************
// Constants used in this module
namespace {
// Set batch size; this is the number of ASTEROIDS, not rows!
constexpr int batch_size = 1000;
}

// *****************************************************************************
// Local names used
using ks::AsteroidElement;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
AsteroidElement::AsteroidElement(int32_t n0_, int32_t n1_, int mjd0_, int mjd1_, int dt_):
    // Data size
    N_ast {n1_ - n0_},
    N_t {(mjd1_ - mjd0_)/dt_ +1},
    N_row {N_ast*N_t},
    // Asteroid range
    n0 {n0_},
    n1 {n1_},
    // Date range
    mjd0 {mjd0_},
    mjd1 {mjd1_},
    dt {dt_},
    // Allocate the one dimensional arrays for asteroid_id and mjd
    asteroid_id {new int32_t[N_ast]},
    mjd {new double[N_t]},
    // Allocate a two dimensional array for seven traditional orbital elements
    elt_a          {new double[N_row]},
    elt_e          {new double[N_row]},
    elt_inc        {new double[N_row]},
    elt_Omega      {new double[N_row]},
    elt_omega      {new double[N_row]},
    elt_f          {new double[N_row]},
    elt_M          {new double[N_row]},
    // Allocate a two dimensional array for orbital angles splined via cosine/sine
    elt_cos_inc    {new double[N_row]},
    elt_sin_inc    {new double[N_row]},
    elt_cos_Omega  {new double[N_row]},
    elt_sin_Omega  {new double[N_row]},
    elt_cos_omega  {new double[N_row]},
    elt_sin_omega  {new double[N_row]},
    // Initialize GSL objects
    // elt_spline has one member for every orbital quantity that is splined
    elt_spline 
    {
        .a          =    {vector<gsl_spline*>(0)},
        .e          =    {vector<gsl_spline*>(0)},
        .M          =    {vector<gsl_spline*>(0)},
        .cos_inc    =    {vector<gsl_spline*>(0)},
        .sin_inc    =    {vector<gsl_spline*>(0)},
        .cos_Omega  =    {vector<gsl_spline*>(0)},
        .sin_Omega  =    {vector<gsl_spline*>(0)},
        .cos_omega  =    {vector<gsl_spline*>(0)},
        .sin_omega  =    {vector<gsl_spline*>(0)}
    },
    // The GSL accelerator
    acc {gsl_interp_accel_alloc() }
{
    // Populate asteroid_id
    for (int i=0; i<N_ast; i++) {asteroid_id[i] = n0+i;}
    // Populate mjd
    for (int i=0; i<N_t; i++) {mjd[i] = mjd0 + i*dt;}

    // Reserve space in vector of splines
    elt_spline.a.reserve(         N_ast);
    elt_spline.e.reserve(         N_ast);
    elt_spline.M.reserve(         N_ast);
    elt_spline.cos_inc.reserve(   N_ast);
    elt_spline.sin_inc.reserve(   N_ast);
    elt_spline.cos_Omega.reserve( N_ast);
    elt_spline.sin_Omega.reserve( N_ast);
    elt_spline.cos_omega.reserve( N_ast);
    elt_spline.sin_omega.reserve( N_ast);

    // Initialize the splines for the elements of each asteroid; one spline for each asteroid and element
    for (int i=0; i<N_ast; i++)
    {
        elt_spline.a.push_back(         gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.e.push_back(         gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.M.push_back(         gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.cos_inc.push_back(   gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.sin_inc.push_back(   gsl_spline_alloc(gsl_interp_cspline, N_t));        
        elt_spline.cos_Omega.push_back( gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.sin_Omega.push_back( gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.cos_omega.push_back( gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.sin_omega.push_back( gsl_spline_alloc(gsl_interp_cspline, N_t));
    }   // for / asteroid number in batch

    // DEBUG
    print("N_ast={:d}, N_t={:d}, N_row={:d}.\n", N_ast, N_t, N_row);
    print("n0={:d}, n1={:d}.\n", n0, n1);
    print("mjd0={:d}, mjd1={:d}, dt={:d}.\n", mjd0, mjd1, dt);

} // end function

// *****************************************************************************
/// Need a non-trivial destructor to release all manually allocated memory
AsteroidElement::~AsteroidElement()
{
    // Delete two 1D arrays
    delete [] asteroid_id;
    delete [] mjd;

    // Delete seven 2D arrays, one for each orbital element
    delete [] elt_a;
    delete [] elt_e;
    delete [] elt_inc;
    delete [] elt_Omega;
    delete [] elt_omega;
    delete [] elt_f;
    delete [] elt_M;

    // Delete 2D arrays, with cosine / sine of the angle orbital elements
    delete [] elt_cos_inc;
    delete [] elt_sin_inc;
    delete [] elt_cos_Omega;
    delete [] elt_sin_Omega;
    delete [] elt_cos_omega;
    delete [] elt_sin_omega;

    // Free GSL resources
    gsl_free();   
}

// *****************************************************************************
// GSL resources that need to be freed before exit are
// (1) One spline for each asteroid / element pair
// (2) The interpolation accelerator
void AsteroidElement::gsl_free()
{
    // One interpolator for each asteroid and each element
    for (int i=0; i<N_ast; i++)
    {
        gsl_spline_free(elt_spline.a[i]);
        gsl_spline_free(elt_spline.e[i]);
        gsl_spline_free(elt_spline.M[i]);
        gsl_spline_free(elt_spline.cos_inc[i]);
        gsl_spline_free(elt_spline.sin_inc[i]);
        gsl_spline_free(elt_spline.cos_Omega[i]);
        gsl_spline_free(elt_spline.sin_Omega[i]);
        gsl_spline_free(elt_spline.cos_omega[i]);
        gsl_spline_free(elt_spline.sin_omega[i]);
    }

    // Just one accelerator object
    gsl_interp_accel_free(acc);
}

// *****************************************************************************
void AsteroidElement::load(db_conn_type &conn, bool progbar)
{
    // Status update
    if (progbar) 
    {
        int batch_count = std::max(N_ast/batch_size, 1);
        print("Processing {:d} asteroids of AsteroidElement data from {:d} to {:d} in {:d} batches of size {:d}...\n",
                N_ast, n0, n1, batch_count, batch_size);
    }

    // Timer for processing from DB
	Timer t;
    t.tick();

    // Iterate over the batches; i0 is the first asteroid number of the batch being processed
    for (int i0=n0; i0<n1; i0+=batch_size)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+batch_size, n1);
        // Process SQL data in this batch
        process_rows(conn, i0, i1);
        // Progress bar
        if (progbar) 
        {
            print(".");
            flush_console();            
        }
    } // for / i (batches of asteroids)
    if (progbar) 
    {
        print("\nLoaded AsteroidSkyPatch table.\n");
        t.tock_msg();
    }

    // Delegate to build_splines method to initialize the gsl_spline* objects
    build_splines();
}   // end function

// *****************************************************************************
void AsteroidElement::process_rows(db_conn_type& conn, int i0, int i1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetAsteroidElements";
    vector<string> params = {to_string(i0), to_string(i1), to_string(mjd0), to_string(mjd1)};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // The time_id corresponding to mjd0
    const int time_id0 = mpd * mjd0;
    // The space in minutes between entries in mjd array.
    // Used to calculate the index in the time array from the mjd field
    const int dt_min = mpd * dt;

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset; 10 total fields
        // Two integer IDs
        int32_t time_id     = rs->getInt("TimeID");
        int32_t asteroid_id = rs->getInt("AsteroidID");
        // The time
        // double mjd = rs->getDouble("mjd");
        // Seven orbital elements
        // MariaDB SQL Connector appears to match strings using case insensitive collation.
        // This brutally disregards the database collation (case sensitive) and is also different from Pandas.
        // This causes field "omega" to receive the value of "Omega".
        // As a workaround, use the integer column number instead of the column name.
        double a        = rs->getDouble("a");
        double e        = rs->getDouble("e");
        double inc      = rs->getDouble("inc");
        double Omega    = rs->getDouble("Omega");
        // double omega = rs->getDouble("omega");
        double omega    = rs->getDouble(8);
        double f        = rs->getDouble("f");
        double M        = rs->getDouble("M");

        // The index (row number) has two terms.
        // The base is the asteroid_row, which is asteroid index times the number of times each asteroid.
        // The offset is the time_row, which is the number of rows written for earlier times on this asteroid.
        int j = asteroid_row(asteroid_id) + (time_id-time_id0)/dt_min;

        // Save seven standard orbital elements
        elt_a[j]      = a;
        elt_e[j]      = e;
        elt_inc[j]    = inc;
        elt_Omega[j]  = Omega;
        elt_omega[j]  = omega;
        elt_f[j]      = f;
        elt_M[j]      = M;

        // Save pairs of cosine / sine of angle orbital elements
        elt_cos_inc[j]    = cos(inc);
        elt_sin_inc[j]    = sin(inc);
        elt_cos_Omega[j]  = cos(Omega);
        elt_sin_Omega[j]  = sin(Omega);
        elt_cos_omega[j]  = cos(omega);
        elt_sin_omega[j]  = sin(omega);
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
void AsteroidElement::build_splines()
{
    // One interpolator for each asteroid and each element
    for (int idx=0; idx<N_ast; idx++)
    {
        // Initialize an interpolation spline for seven traditional orbital elements of this body
        gsl_spline_init(elt_spline.a[idx]           , mjd, get_a(idx),         N_t);
        gsl_spline_init(elt_spline.e[idx]           , mjd, get_e(idx),         N_t);
        gsl_spline_init(elt_spline.M[idx]           , mjd, get_M(idx),         N_t);

        // Initialize an interpolation spline for ten orbital angles of this body
        gsl_spline_init(elt_spline.cos_inc[idx]     , mjd, get_cos_inc(idx),   N_t);
        gsl_spline_init(elt_spline.sin_inc[idx]     , mjd, get_sin_inc(idx),   N_t);
        gsl_spline_init(elt_spline.cos_Omega[idx]   , mjd, get_cos_Omega(idx), N_t);
        gsl_spline_init(elt_spline.sin_Omega[idx]   , mjd, get_sin_Omega(idx), N_t);
        gsl_spline_init(elt_spline.cos_omega[idx]   , mjd, get_cos_omega(idx), N_t);
        gsl_spline_init(elt_spline.sin_omega[idx]   , mjd, get_sin_omega(idx), N_t);
    }
}

// *****************************************************************************
// The asteroids are laid out in contiguous order.
// The asteroid_id determines the offset into the array.
// If there is a hole in the sequence, those memory locations are left in whatever
// random state they were initialized with by operator new.
// This allows for fastest memory access in O(1) steps, bypassing the need
// to do a binary search when accessing the elements of an asteroid identified
// by its asteroid_id field.  
// There are very few gaps in the asteroid_id sequence except between 545135 
// (last numbered asteroid) and 1000001 (first unnumbered asteroid).

// *****************************************************************************
// Calculate splined orbital elements
// *****************************************************************************

// *****************************************************************************
const OrbitalElement AsteroidElement::interp_elt(int32_t asteroid_id, double mjd) const
{
    // Get the index number of this asteroid_id
    int idx = asteroid_idx(asteroid_id);

    // Seven traditional elements
    gsl_spline* interp_a         = elt_spline.a[idx];
    gsl_spline* interp_e         = elt_spline.e[idx];
    gsl_spline* interp_M         = elt_spline.M[idx];

    // Ten orbital angles
    gsl_spline* interp_cos_inc   = elt_spline.cos_inc[idx];
    gsl_spline* interp_sin_inc   = elt_spline.sin_inc[idx];
    gsl_spline* interp_cos_Omega = elt_spline.cos_Omega[idx];
    gsl_spline* interp_sin_Omega = elt_spline.sin_Omega[idx];
    gsl_spline* interp_cos_omega = elt_spline.cos_omega[idx];
    gsl_spline* interp_sin_omega = elt_spline.sin_omega[idx];

    // Spline a and e at selected time
    double a = gsl_spline_eval(interp_a, mjd, acc);
    double e = gsl_spline_eval(interp_e, mjd, acc);

    // Spline cosine and sine of inc, Omega, omega
    double cos_inc    = gsl_spline_eval(interp_cos_inc,   mjd, acc);
    double sin_inc    = gsl_spline_eval(interp_sin_inc,   mjd, acc);
    double cos_Omega  = gsl_spline_eval(interp_cos_Omega, mjd, acc);
    double sin_Omega  = gsl_spline_eval(interp_sin_Omega, mjd, acc);
    double cos_omega  = gsl_spline_eval(interp_cos_omega, mjd, acc);
    double sin_omega  = gsl_spline_eval(interp_sin_omega, mjd, acc);

    // Calculate inc, Omega, omega from its components
    double inc      = atan2(sin_inc,    cos_inc);
    double Omega    = atan2(sin_Omega,  cos_Omega);
    double omega    = atan2(sin_omega,  cos_omega);

    // Directly spline M - this is safe b/c DB source includes winding for f and M
    double M = gsl_spline_eval(interp_M, mjd, acc);

    // Calculate f from M and e
    double f = anomaly_M2f(M, e);

    // Return the an assembled orbital element
    return OrbitalElement {.a=a, .e=e, .inc=inc, .Omega=Omega, .omega=omega, .f=f, .M=M};
}

// *****************************************************************************
// This calculation returns the position in the HELIOCENTRIC frame, NOT the BME!
// The orbital element describes the relative position of the asteroid vs. the primary, which is the Sun here.
const Position AsteroidElement::interp_pos_hel(int32_t asteroid_id, double mjd) const
{
    // Delegate to interp_elt to spline the elements
    OrbitalElement elt = interp_elt(asteroid_id, mjd);
    // Call elt2pos to calculate a position
    return elt2pos(elt);
}

// *****************************************************************************
// This calculation returns the state vector in the HELIOCENTRIC frame, NOT the BME!
// The orbital element describes the relative vectors of the asteroid vs. the primary, which is the Sun here.
const StateVector AsteroidElement::interp_vec_hel(int32_t asteroid_id, double mjd) const
{
    // Delegate to interp_elt to spline the elements
    OrbitalElement elt = interp_elt(asteroid_id, mjd);
    // Call elt2pos to calculate a position
    return elt2vec(elt, mu_sun);   
}

// *****************************************************************************
// Add the heliocentric position from splined orbital elements to the Sun's state vectors
// to get the position in the BME frame.
const Position AsteroidElement::interp_pos(int32_t asteroid_id, double mjd) const
{
    // Delegate to interp_pos_hel for the relative position of asteroid vs. the Sun
    Position ast = interp_pos_hel(asteroid_id, mjd);
    // Delegate to bv_sun to get interpolated position of Sun in BME
    Position sun = bv_sun.interp_pos(mjd);
    // Add the two components
    return ast + sun;
}

// *****************************************************************************
// Add the heliocentric state vectors from splined orbital elements to the Sun's state vectors
// to get the state vectors in the BME frame.
const StateVector AsteroidElement::interp_vec(int32_t asteroid_id, double mjd) const
{
    // Delegate to interp_pos_hel for the relative position of asteroid vs. the Sun
    StateVector ast = interp_vec_hel(asteroid_id, mjd);
    // Delegate to bv_sun to get interpolated state vector of Sun in BME
    StateVector sun = bv_sun.interp_vec(mjd);
    // Add the two components
    return ast + sun;
}
