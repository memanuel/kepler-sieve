/** @file   PlanetElement.cpp
 *  @brief  Implmentation of PlanetElement class.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-08
 */

// *****************************************************************************
// Included files
#include "PlanetElement.hpp"

// *****************************************************************************
// Constants used in this module

// Number of minutes in one day
constexpr int mpd = 1440;
// Number of days in one minute
constexpr double dpm = 1.0 / mpd;

// Number of massive bodies in Planets collection: Sun + 9 planets + moon = 11 bodies
// But the Sun is the primary so it has no elements.  There are *10* bodies with orbital elements.
constexpr int N_body_planets = 10;
// Array of body_id in this problem
constexpr int32_t body_id_planets[10] = {1, 2, 4, 5, 6, 7, 8, 9, 301, 399};
// The body_id of the sun
constexpr int32_t body_id_sun = 10;
// The position and state vector of the sun in the heliocentric frame are zero by definition
constexpr Position pos_sun_hel {.qx=0.0, .qy=0.0, .qz=0.0};
constexpr StateVector vec_sun_hel {.qx=0.0, .qy=0.0, .qz=0.0, .vx=0.0, .vy=0.0, .vz=0.0};

// Stride in minutes for database vectors for Sun and planets
constexpr int stride_db_min = 5;
// Batch size for loading dates from DB
constexpr int batch_size = 1000;

// *****************************************************************************
// Local names used
using ks::PlanetElement;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
PlanetElement::PlanetElement(int mjd0, int mjd1, int dt_min):
    // Data size
    N_body(N_body_planets),
    // number of times (data size)
    N_t((mjd1-mjd0)*mpd/lcm(dt_min, stride_db_min)+1),
    // Date range
    mjd0(mjd0),
    mjd1(mjd1),
    // Stride in minutes between consecutive entries
    // database resolution is every 5 minutes, so apply LCM function to the input to ensure validity.
    dt_min(lcm(dt_min, stride_db_min)),
    // Allocate the one dimensional arrays for body_id and mjd
    body_id(new int32_t[N_body]),
    mjd(new double[N_t]),
    // Allocate a two dimensional array for each orbital element
    elt_a(      new double[N_t*N_body]),
    elt_e(      new double[N_t*N_body]),
    elt_inc(    new double[N_t*N_body]),
    elt_Omega(  new double[N_t*N_body]),
    elt_omega(  new double[N_t*N_body]),
    elt_f(      new double[N_t*N_body]),
    elt_M(      new double[N_t*N_body]),
    // Initialize GSL objects
    acc(gsl_interp_accel_alloc() ),
    // BodyVector object for interpolated Sun state vectors
    bv_sun(BodyVector("Sun"))
{
    // Populate body_id - this is a copy from body_id_planets
    for (int i=0; i<N_body; i++) {body_id[i]=body_id_planets[i];}
    // Populate mjd; these are on a fixed schedule spaced dt_min minutes apart from mjd0 to mjd1
    const double dt = dt_min * dpm;
    for (int i=0; i<N_t; i++) {mjd[i] = mjd0 + i*dt;}

    // elt_spline is a structure with one member for each of seven elements
    // Reserve space in vector of splines for each element    
    elt_spline.a.reserve(     N_body);
    elt_spline.e.reserve(     N_body);
    elt_spline.inc.reserve(   N_body);
    elt_spline.Omega.reserve( N_body);
    elt_spline.omega.reserve( N_body);
    elt_spline.f.reserve(     N_body);
    elt_spline.M.reserve(     N_body);

    // Initialize the splines for the elements of each body; one spline for each asteroid and element
    for (int i=0; i<N_body; i++)
    {
        elt_spline.a.push_back(     gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.e.push_back(     gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.inc.push_back(   gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.Omega.push_back( gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.omega.push_back( gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.f.push_back(     gsl_spline_alloc(gsl_interp_cspline, N_t));
        elt_spline.M.push_back(     gsl_spline_alloc(gsl_interp_cspline, N_t));
    }   // for / i (loop over massive bodies)
} // end function

// *****************************************************************************
/// Need a non-trivial destructor to release all manually allocated memory
PlanetElement::~PlanetElement()
{
    // Delete two 1D arrays
    delete [] body_id;
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
void PlanetElement::gsl_free()
{
    // One interpolator for each body and each element
    for (int i=0; i<N_body; i++)
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
void PlanetElement::load(db_conn_type &conn)
{
    // Iterate over the batches
    for (int t0=mjd0; t0<mjd1; t0+=batch_size)
    {
        // Upper limit for this batch
        int t1 = std::min(t0+batch_size, mjd1);
        // Process SQL data in this batch
        process_rows(conn, t0, t1);
        // Status
        print("Loaded rows for mjd in [{:d}, {:d}].\n", t0, t1);
    }   // for / t0
}   // end function

// *****************************************************************************
void PlanetElement::process_rows(db_conn_type& conn, int t0, int t1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetElements_Planets";
    // Run the stored procedure to get planet elements in [t0, t1]
    vector<string> params = {to_string(t0), to_string(t1), to_string(dt_min)};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // The time_id corresponding to mjd0
    const int time_id0 = mpd * mjd0;

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset; 10 total fields
        // Two integer IDs
        int32_t time_id = rs->getInt("TimeID");
        int32_t body_id = rs->getInt("BodyID");
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
        int idx = body_row(body_id) + (time_id-time_id0)/dt_min;

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
void PlanetElement::build_splines()
{
    // One interpolator for each asteroid and each element
    // idx is the body index (i.e. counter from 0 to 9), NOT the body_id!
    for (int idx=0; idx<N_body; idx++)
    {
        // Evaluate the interpolation spline for each orbital element of this body
        gsl_spline_init(elt_spline.a[idx],     mjd, get_a(idx),     N_t);
        gsl_spline_init(elt_spline.e[idx],     mjd, get_e(idx),     N_t);
        gsl_spline_init(elt_spline.inc[idx],   mjd, get_inc(idx),   N_t);
        gsl_spline_init(elt_spline.Omega[idx], mjd, get_Omega(idx), N_t);
        gsl_spline_init(elt_spline.omega[idx], mjd, get_omega(idx), N_t);
        gsl_spline_init(elt_spline.f[idx],     mjd, get_f(idx),     N_t);
        gsl_spline_init(elt_spline.M[idx],     mjd, get_M(idx),     N_t);
    }
}

// *****************************************************************************
const int PlanetElement::body_idx(int32_t body_id) const
{
    // Bodies are laid out in the order 1, 2, 4, 5, 6, 7, 8, 9, 301, 399
    if ((0 < body_id) && (body_id < 3)) {return body_id-1;}
    else if ((3 < body_id) && (body_id < 10)) {return body_id-2;}
    else if (body_id==301) {return 8;}
    else if (body_id==399) {return 9;}
    else {throw domain_error(
        "Bad body_id! Must be one of 1, 2, 4, 5, 6, 7, 8, 9, 301, 399 (planet barycenters, earth and moon.\n");}
}

// *****************************************************************************
const int PlanetElement::body_row(int32_t body_id) const
{
    return body_idx(body_id)*N_t;
}

// *****************************************************************************
// Get 1D arrays of body_id and times
// *****************************************************************************

// *****************************************************************************
int32_t* PlanetElement::get_body_id() const {return body_id;}

// *****************************************************************************
double* PlanetElement::get_mjd() const {return mjd;}

// *****************************************************************************
// Get 1D array of each orbital element given a body_id
// *****************************************************************************

// *****************************************************************************
double* PlanetElement::get_a(int idx) const {return elt_a + N_t*idx;}

// *****************************************************************************
double* PlanetElement::get_e(int idx) const {return elt_e + N_t*idx;}

// *****************************************************************************
double* PlanetElement::get_inc(int idx) const {return elt_inc + N_t*idx;}

// *****************************************************************************
double* PlanetElement::get_Omega(int idx) const {return elt_Omega + N_t*idx;}

// *****************************************************************************
double* PlanetElement::get_omega(int idx) const {return elt_omega + N_t*idx;}

// *****************************************************************************
double* PlanetElement::get_f(int idx) const {return elt_f + N_t*idx;}

// *****************************************************************************
double* PlanetElement::get_M(int idx) const {return elt_M + N_t*idx;}

// *****************************************************************************
OrbitalElement PlanetElement::interp_elt(int32_t body_id, double mjd) const
{
    // Get the index number of this body
    int idx = body_idx(body_id);

    // The interpolators for each orbital element for this body
    // The array index into the elt_spline members is idx, NOT body_id;
    // body_id will go over the end, e.g. for Earth body_id=399 and idx=9.
    gsl_spline* gsl_interp_a     = elt_spline.a[idx];
    gsl_spline* gsl_interp_e     = elt_spline.e[idx];
    gsl_spline* gsl_interp_inc   = elt_spline.inc[idx];
    gsl_spline* gsl_interp_Omega = elt_spline.Omega[idx];
    gsl_spline* gsl_interp_omega = elt_spline.omega[idx];
    gsl_spline* gsl_interp_f     = elt_spline.f[idx];
    gsl_spline* gsl_interp_M     = elt_spline.M[idx];

    // Evaluate the splines at the selected time
    OrbitalElement elt 
    {
        .a      = gsl_spline_eval(gsl_interp_a,     mjd, acc),
        .e      = gsl_spline_eval(gsl_interp_e,     mjd, acc),
        .inc    = gsl_spline_eval(gsl_interp_inc,   mjd, acc),
        .Omega  = gsl_spline_eval(gsl_interp_Omega, mjd, acc),
        .omega  = gsl_spline_eval(gsl_interp_omega, mjd, acc),
        .f      = gsl_spline_eval(gsl_interp_f,     mjd, acc),
        .M      = gsl_spline_eval(gsl_interp_M,     mjd, acc)
    };
    return elt;
}

// *****************************************************************************
// This calculation returns the position in the HELIOCENTRIC frame, NOT the BME!
// The orbital element describes the relative position of the body vs. the primary, which is the Sun here.
Position PlanetElement::interp_pos_hel(int32_t body_id, double mjd) const
{
    // Peel off special case that the body is the sun; then the position is zero
    if (body_id==body_id_sun){return pos_sun_hel;}
    // Delegate to interp_elt to spline the elements
    OrbitalElement elt = interp_elt(body_id, mjd);
    // Call elt2pos to calculate a position
    return elt2pos(elt);
}

// *****************************************************************************
// This calculation returns the state vector in the HELIOCENTRIC frame, NOT the BME!
// The orbital element describes the relative vectors of the body vs. the primary, which is the Sun here.
StateVector PlanetElement::interp_vec_hel(int32_t body_id, double mjd) const
{
    // Peel off special case that the body is the sun; then the state vector is zero
    if (body_id==body_id_sun){return vec_sun_hel;}
    // Delegate to interp_elt to spline the elements
    OrbitalElement elt = interp_elt(body_id, mjd);
    // Call elt2pos to calculate a position
    return elt2vec(elt);   
}

// *****************************************************************************
// Add the heliocentric position from splined orbital elements to the Sun's state vectors
// to get the position in the BME frame.
Position PlanetElement::interp_pos(int32_t body_id, double mjd) const
{
    // Delegate to interp_pos_hel for the relative position of asteroid vs. the Sun
    Position ast = interp_pos_hel(body_id, mjd);
    // Delegate to bv_sun to get interpolated position of Sun in BME
    Position sun = bv_sun.interp_pos(mjd);
    // Add the two components
    return Position
    {
        .qx = ast.qx + sun.qx,
        .qy = ast.qy + sun.qy,
        .qz = ast.qz + sun.qz
    };
}

// *****************************************************************************
// Add the heliocentric state vectors from splined orbital elements to the Sun's state vectors
// to get the state vectors in the BME frame.
StateVector PlanetElement::interp_vec(int32_t body_id, double mjd) const
{
    // Delegate to interp_pos_hel for the relative position of asteroid vs. the Sun
    StateVector ast = interp_vec_hel(body_id, mjd);
    // Delegate to bv_sun to get interpolated state vectors of Sun in BME
    StateVector sun = bv_sun.interp_vec(mjd);
    // Add the two components
    return StateVector
    {
        .qx = ast.qx + sun.qx,
        .qy = ast.qy + sun.qy,
        .qz = ast.qz + sun.qz,
        .vx = ast.vx + sun.vx,
        .vy = ast.vy + sun.vy,
        .vz = ast.vz + sun.vz
    };
}
