/** @file AsteroidElement.cpp
 *  @brief Implmentation of AsteroidElement class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-02
 */

// *****************************************************************************
// Included files
#include "AsteroidElement.hpp"

// Set batch size; this is the number of ASTEROIDS, not rows!
constexpr int batch_size = 100;
// Number of minutes in one day
constexpr int mpd = 1440;

// *****************************************************************************
// Local names used
using ks::AsteroidElement;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
AsteroidElement::AsteroidElement(int n0, int n1, int mjd0, int mjd1, int dt):
    // Data size
    N_ast(n1-n0),
    N_t((mjd1-mjd0)/dt+1),
    // Asteroid range
    n0(n0),
    n1(n1),
    // Date range
    mjd0((mjd0/dt)*dt),
    mjd1((mjd1/dt)*dt),
    dt(dt),
    // Allocate the one dimensional arrays for asteroid_id and mjd
    asteroid_id(new int32_t[N_ast]),
    mjd(new double[N_t]),
    // Allocate a two dimensional array for each orbital element
    elt_a(new double[N_t*N_ast]),
    elt_e(new double[N_t*N_ast]),
    elt_inc(new double[N_t*N_ast]),
    elt_Omega(new double[N_t*N_ast]),
    elt_omega(new double[N_t*N_ast]),
    elt_f(new double[N_t*N_ast]),
    elt_M(new double[N_t*N_ast]),
    // Initialize GSL objects
    acc(gsl_interp_accel_alloc() ),
    // BodyVector object for interpolated Sun state vectors
    bv_sun(BodyVector("Sun"))
{
    // Populate asteroid_id
    for (int i=0; i<N_ast; i++) {asteroid_id[i] = n0+i;}
    // Populate mjd
    for (int i=0; i<N_t; i++) {mjd[i] = mjd0 + i*dt;}

    // elt_spline is a structure with one member for each of seven elements
    // Reserve space in vector of splines for each element    
    elt_spline.a.reserve(N_ast);
    elt_spline.e.reserve(N_ast);
    elt_spline.inc.reserve(N_ast);
    elt_spline.Omega.reserve(N_ast);
    elt_spline.omega.reserve(N_ast);
    elt_spline.f.reserve(N_ast);
    elt_spline.M.reserve(N_ast);

    // Initialize the splines for the elements of each asteroid; one spline for each asteroid and element
    for (int i=0; i<N_ast; i++)
    {
        elt_spline.a.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.e.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.inc.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.Omega.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.omega.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.f.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.M.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
    }

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
        gsl_spline_free(elt_spline.inc[i]);
        gsl_spline_free(elt_spline.Omega[i]);
        gsl_spline_free(elt_spline.omega[i]);
        gsl_spline_free(elt_spline.f[i]);
        gsl_spline_free(elt_spline.M[i]);
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
    int mjd1 = mjd0 + N_t*4;
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
        int idx = asteroid_row(asteroid_id) + (time_id-time_id0)/dt_min;

        // Save the data fields to the member arrays
        elt_a[idx]      = a;
        elt_e[idx]      = e;
        elt_inc[idx]    = inc;
        elt_Omega[idx]  = Omega;
        elt_omega[idx]  = omega;
        elt_f[idx]      = f;
        elt_M[idx]      = M;

    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
void AsteroidElement::build_splines()
{
    // One interpolator for each asteroid and each element
    for (int i=0; i<N_ast; i++)
    {
        gsl_spline_init(elt_spline.a[i],     mjd, get_a(i), N_t);
        gsl_spline_init(elt_spline.e[i],     mjd, get_e(i), N_t);
        gsl_spline_init(elt_spline.inc[i],   mjd, get_inc(i), N_t);
        gsl_spline_init(elt_spline.Omega[i], mjd, get_Omega(i), N_t);
        gsl_spline_init(elt_spline.omega[i], mjd, get_omega(i), N_t);
        gsl_spline_init(elt_spline.f[i],     mjd, get_f(i), N_t);
        gsl_spline_init(elt_spline.M[i],     mjd, get_M(i), N_t);
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
const int AsteroidElement::asteroid_idx(int32_t asteroid_id) const 
{
    return asteroid_id - n0;
}

// *****************************************************************************
const int AsteroidElement::asteroid_row(int32_t asteroid_id) const
{
    return asteroid_idx(asteroid_id)*N_t;
}

// *****************************************************************************
// Get 1D arrays of asteroid_id and times
// *****************************************************************************

// *****************************************************************************
const int32_t* AsteroidElement::get_asteroid_id() const 
    {return asteroid_id;}

// *****************************************************************************
const double* AsteroidElement::get_mjd() const 
    {return mjd;}

// *****************************************************************************
// Get 1D array of each orbital element given an asteroid_id
// *****************************************************************************

// *****************************************************************************
const double* AsteroidElement::get_a(int32_t asteroid_id) const
    {return elt_a + asteroid_row(asteroid_id);}

// *****************************************************************************
const double* AsteroidElement::get_e(int32_t asteroid_id) const
    {return elt_e + asteroid_row(asteroid_id);}

// *****************************************************************************
const double* AsteroidElement::get_inc(int32_t asteroid_id) const
    {return elt_inc + asteroid_row(asteroid_id);}

// *****************************************************************************
const double* AsteroidElement::get_Omega(int32_t asteroid_id) const
    {return elt_Omega + asteroid_row(asteroid_id);}

// *****************************************************************************
const double* AsteroidElement::get_omega(int32_t asteroid_id) const
    {return elt_omega + asteroid_row(asteroid_id);}

// *****************************************************************************
const double* AsteroidElement::get_f(int32_t asteroid_id) const
    {return elt_f + asteroid_row(asteroid_id);}

// *****************************************************************************
const double* AsteroidElement::get_M(int32_t asteroid_id) const
    {return elt_M + asteroid_row(asteroid_id);}

// *****************************************************************************
const OrbitalElement AsteroidElement::interp_elt(int32_t asteroid_id, double mjd) const
{
    // The interpolators for each orbital element for this asteroid
    gsl_spline* gsl_interp_a = elt_spline.a[asteroid_id];
    gsl_spline* gsl_interp_e = elt_spline.e[asteroid_id];
    gsl_spline* gsl_interp_inc = elt_spline.inc[asteroid_id];
    gsl_spline* gsl_interp_Omega = elt_spline.Omega[asteroid_id];
    gsl_spline* gsl_interp_omega = elt_spline.omega[asteroid_id];
    gsl_spline* gsl_interp_f = elt_spline.f[asteroid_id];
    gsl_spline* gsl_interp_M = elt_spline.M[asteroid_id];

    // Evaluate the splines at the selected time
    OrbitalElement elt 
    {
        .a      = gsl_spline_eval(gsl_interp_a, mjd, acc),
        .e      = gsl_spline_eval(gsl_interp_e, mjd, acc),
        .inc    = gsl_spline_eval(gsl_interp_inc, mjd, acc),
        .Omega  = gsl_spline_eval(gsl_interp_Omega, mjd, acc),
        .omega  = gsl_spline_eval(gsl_interp_omega, mjd, acc),
        .f      = gsl_spline_eval(gsl_interp_f, mjd, acc),
        .M      = gsl_spline_eval(gsl_interp_M, mjd, acc)
    };
    // Replace splined f value with conversion from M (M splines better than f)
    elt.f = anomaly_M2f(elt.M, elt.e);   
    return elt;
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
    return elt2vec(elt);   
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
