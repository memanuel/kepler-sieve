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
// constexpr int batch_size = 100;
// DEBUG
constexpr int batch_size = 10;
// Number of minutes in one day
constexpr int mpd = 1440;

// *****************************************************************************
// Local names used
using ks::AsteroidElement;
// using ks::AsteroidElementEntry;
// using ks::ElementArrays;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
AsteroidElement::AsteroidElement(
    int n0, int n1, int mjd0, int mjd1, int dt):
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
    elt_M(new double[N_t*N_ast])
{
    // Populate asteroid_id
    for (int i=0; i<N_ast; i++) {asteroid_id[i] = n0+i;}
    // Populate mjd
    for (int i=0; i<N_t; i++) {mjd[i] = mjd0 + i*dt;}
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
        int32_t time_id = rs->getInt("TimeID");
        int32_t asteroid_id = rs->getInt("AsteroidID");
        // The time
        // double mjd = rs->getDouble("mjd");
        // Seven orbital elements
        double a = rs->getDouble("a");
        double e = rs->getDouble("e");
        double inc = rs->getDouble("inc");
        double Omega = rs->getDouble("Omega");
        double omega = rs->getDouble("omega");
        double f = rs->getDouble("f");
        double M = rs->getDouble("M");

        // The index (row number) has two terms.
        // The base is the asteroid_row, which is asteroid index times the number of times each asteroid.
        // The offset is the time_row, which is the number of rows written for earlier times on this asteroid.
        int idx = asteroid_row(asteroid_id) + (time_id-time_id0)/dt_min;

        // Save the data fields to the member arrays
        elt_a[idx] = a;
        elt_e[idx] = e;
        elt_inc[idx] = inc;
        elt_Omega[idx] = Omega;
        elt_omega[idx] = omega;
        elt_f[idx] = f;
        elt_M[idx] = M;

    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;    
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
int32_t AsteroidElement::asteroid_idx(int32_t asteroid_id) const
{
    return asteroid_id - n0;
}

// *****************************************************************************
int32_t AsteroidElement::asteroid_row(int32_t asteroid_id) const
{
    return asteroid_idx(asteroid_id)*N_t;
}

// *****************************************************************************
// Get 1D arrays of asteroid_id and times
// *****************************************************************************

// *****************************************************************************
int32_t* AsteroidElement::get_asteroid_id() const
{
    return asteroid_id;
}

// *****************************************************************************
double* AsteroidElement::get_mjd() const
{
    return mjd;
}

// *****************************************************************************
// Get 1D array of each orbital element given an asteroid_id
// *****************************************************************************

// *****************************************************************************
double* AsteroidElement::get_a(int32_t asteroid_id) const
{
    return elt_a + asteroid_row(asteroid_id);
}

// *****************************************************************************
double* AsteroidElement::get_e(int32_t asteroid_id) const
{
    return elt_e + asteroid_row(asteroid_id);
}

// *****************************************************************************
double* AsteroidElement::get_inc(int32_t asteroid_id) const
{
    return elt_inc + asteroid_row(asteroid_id);
}

// *****************************************************************************
double* AsteroidElement::get_Omega(int32_t asteroid_id) const
{
    return elt_Omega + asteroid_row(asteroid_id);
}

// *****************************************************************************
double* AsteroidElement::get_omega(int32_t asteroid_id) const
{
    return elt_omega + asteroid_row(asteroid_id);
}

// *****************************************************************************
double* AsteroidElement::get_f(int32_t asteroid_id) const
{
    return elt_f + asteroid_row(asteroid_id);
}

// *****************************************************************************
double* AsteroidElement::get_M(int32_t asteroid_id) const
{
    return elt_M + asteroid_row(asteroid_id);
}
