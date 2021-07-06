/** @file BodyVector.cpp
 *  @brief Implmentation of BodyVector class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-05
 */

// *****************************************************************************
// Included files
#include "BodyVector.hpp"

// Number of minutes in one day
constexpr int mpd = 1440;
// One day as a floating point number of minutes
constexpr double dpm = 1.0 / mpd;
// Stride in minutes for database vectors for Sun and Earth
constexpr int stride_db_min = 5;

// *****************************************************************************
// Local names used
using ks::BodyVector;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
BodyVector::BodyVector(int mjd0, int mjd1, int dt_min):
    // number of times (data size)
    N_t((mjd1-mjd0)*mpd/dt_min+1),
    // Date range
    mjd0(mjd0),
    mjd1(mjd1),
    // Stride in minutes between consecutive entries
    // database resolution is every 5 minutes, so apply LCM function to the input to ensure validity.
    dt_min(lcm(dt_min, stride_db_min)),
    // Allocate array for mjd
    mjd(new double[N_t]),
    // Allocate arrays for each component of the state vector
    qx(new double[N_t]),
    qy(new double[N_t]),
    qz(new double[N_t]),
    vx(new double[N_t]),
    vy(new double[N_t]),
    vz(new double[N_t]),

    // Initialize GSL accelerator
    acc(gsl_interp_accel_alloc() )
{
    // Populate mjd
    for (int i=0; i<N_t; i++) {mjd[i] = mjd0 + i*dt_min*dpm;}

    // Initialize the splines for the elements of each asteroid; one spline for each asteroid and element
    vec_spline.qx = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.qy = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.qz = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.vx = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.vy = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.vz = gsl_spline_alloc(gsl_interp_cspline, N_t);

} // end function

// *****************************************************************************
/// Need a non-trivial destructor to release all manually allocated memory
BodyVector::~BodyVector()
{
    // Delete arrays for mjd and six vector components
    delete [] mjd;
    delete [] qx;
    delete [] qy;
    delete [] qz;
    delete [] vx;
    delete [] vy;
    delete [] vz;

    // Free GSL resources
    gsl_free();
}

// *****************************************************************************
// GSL resources that need to be freed before exit are
// (1) One spline for each state vector component
// (2) The interpolation accelerator
void BodyVector::gsl_free()
{
    // One interpolator for each component of state vector
    gsl_spline_free(vec_spline.qx);
    gsl_spline_free(vec_spline.qy);
    gsl_spline_free(vec_spline.qz);
    gsl_spline_free(vec_spline.vx);
    gsl_spline_free(vec_spline.vy);
    gsl_spline_free(vec_spline.vz);

    // Just one accelerator object
    gsl_interp_accel_free(acc);
}

// *****************************************************************************
void BodyVector::load(db_conn_type &conn, bool progbar)
{
    // Status update
    if (progbar) 
    {
        print("Processing body vectors data from mjd {:d} to {:d}...\n", mjd0, mjd1);
    }

    // Timer for processing from DB
	// Timer t;
    // t.tick();

    // Process SQL data
    process_rows(conn);

    if (progbar) 
    {
        print("\nLoaded BodyVector data.\n");
        // t.tock_msg();
    }
    // Delegate to build_splines method to initialize the gsl_spline* objects
    build_splines();
}   // end function

// *****************************************************************************
void BodyVector::process_rows(db_conn_type& conn)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetStateVectors_Sun";
    vector<string> params = {to_string(mjd0), to_string(mjd1), to_string(dt_min)};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // The time_id corresponding to mjd0
    const int time_id0 = mpd * mjd0;

    // Loop through resultset
    while (rs->next()) 
    {
        // The TimeID
        int32_t time_id = rs->getInt("TimeID");
        // The row where this data is written
        int idx = (time_id-time_id0)/dt_min;
        // Write to mjd array
        mjd[idx] = rs->getDouble("mjd");
        // Write state vector components to arrays
        qx[idx] = rs->getDouble("qx");
        qy[idx] = rs->getDouble("qy");
        qz[idx] = rs->getDouble("qz");
        vx[idx] = rs->getDouble("vx");
        vy[idx] = rs->getDouble("vy");
        vz[idx] = rs->getDouble("vz");

    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
void BodyVector::build_splines()
{
    // One interpolator for each element component
    gsl_spline_init(vec_spline.qx, mjd, qx, N_t);
    gsl_spline_init(vec_spline.qy, mjd, qy, N_t);
    gsl_spline_init(vec_spline.qz, mjd, qz, N_t);
    gsl_spline_init(vec_spline.vx, mjd, vx, N_t);
    gsl_spline_init(vec_spline.vy, mjd, vy, N_t);
    gsl_spline_init(vec_spline.vz, mjd, vz, N_t);
}

// *****************************************************************************
// Get interpolated position and state vectors
// *****************************************************************************

// *****************************************************************************
double* BodyVector::get_mjd() const
{
    return mjd;
}

// *****************************************************************************
Position BodyVector::interp_pos(double mjd)
{
    return Position
    {
        .qx = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .qy = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .qz = gsl_spline_eval(vec_spline.qx, mjd, acc)
    };
}

// *****************************************************************************
StateVector BodyVector::interp_vec(double mjd)
{
    return StateVector
    {
        .qx = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .qy = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .qz = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .vx = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .vy = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .vz = gsl_spline_eval(vec_spline.qx, mjd, acc)
    };
}
