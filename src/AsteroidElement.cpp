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

// *****************************************************************************
// Local names used
using ks::AsteroidElement;
using ks::AsteroidElementEntry;
using ks::ElementArrays;

// *****************************************************************************
// Native constructor just allocates memory
AsteroidElement::AsteroidElement(int n0, int N_ast, int mjd0, int N_t):
    // Initialize the integer members from the arguments
    n0(n0),
    mjd0(mjd0),
    N_ast(N_ast),
    N_t(N_t),
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
{}

// *****************************************************************************
/**Helper function: Process a batch of rows, 
 * writing data from stored preocedure output to vector of detections. */
void AsteroidElement::process_rows(db_conn_type& conn, int i0, int i1)
{
    ;
}

// *****************************************************************************
AsteroidElement::AsteroidElement(db_conn_type &conn, int n0, int n1, int mjd0, int mjd1, bool progbar):
    // Delegate to the native constructor
    AsteroidElement(n0, n1-n0, (mjd0/4)*4, (mjd1/4)-(mjd0/4)) 
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

    // Wrap the element arrays into the stucture elt
    ElementArrays elt {
        .a = elt_a,
        .e = elt_e,
        .inc = elt_inc,
        .Omega = elt_Omega,
        .omega = elt_omega,
        .f = elt_f,
        .M = elt_M
    };

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
    }
    if (progbar) 
    {
        print("\nLoaded AsteroidSkyPatch table.\n");
        t.tock_msg();
    }

}

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