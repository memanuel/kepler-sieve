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
using ks::cs::mpd;
// Number of days in one minute
using ks::cs::dpm;

// Number of massive bodies in Planets collection: Sun + 9 planets + moon = 11 bodies
// But the Sun is the primary so it has no elements.  There are *10* bodies with orbital elements.
constexpr int N_body_planets = 10;
// Array of body_id in this problem
constexpr int32_t body_id_planets[N_body_planets] = {1, 2, 4, 5, 6, 7, 8, 9, 301, 399};
// The body_id of the Sun, Earth and Moon
using ks::cs::body_id_sun;
using ks::cs::body_id_earth;
using ks::cs::body_id_moon;
// Mass of the Sun and gravitational constant
using ks::cs::m_sun;
using ks::cs::G;
// The position and state vector of the sun in the heliocentric frame are zero by definition
// constexpr Position pos_sun_hel {.qx=0.0, .qy=0.0, .qz=0.0};
// constexpr StateVector vec_sun_hel {.qx=0.0, .qy=0.0, .qz=0.0, .vx=0.0, .vy=0.0, .vz=0.0};

// Start and end date for database
using ks::cs::mjd0_db;
using ks::cs::mjd1_db;
// Stride in minutes for database vectors for Sun and planets
using ks::cs::stride_db_min;
// Number of times expected in database
using ks::cs::N_t_db;

// Batch size for loading dates from DB; this is a block of dates e.g. [48000, 49000]
constexpr int batch_size = 1000;
// Location of file with serialized data
const string file_name = "data/cache/PlanetElement.bin";

// *****************************************************************************
// Local names used
using ks::PlanetElement;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
PlanetElement::PlanetElement(int mjd0, int mjd1, int dt_min):
    // Data size: number of bodies and number of times
    N_body(N_body_planets),
    N_t((mjd1-mjd0)*mpd/lcm(dt_min, stride_db_min)+1),
    // Date range
    mjd0(mjd0),
    mjd1(mjd1),
    // Stride in minutes between consecutive entries
    // database resolution is every 5 minutes, so apply LCM function to the input to ensure validity.
    dt_min(lcm(dt_min, stride_db_min)),
    // Number of rows- used for array initialization
    N_row(N_body*N_t),
    // time_id range - for fast loading from file
    time_id0(mjd0*mpd),
    time_id1(mjd1*mpd),
    // Allocate the one dimensional arrays for body_id and mjd
    body_id(new int32_t[N_body]),
    mjd(new double[N_t]),
    // Gravitational interaction strength mu for each pair of primary / target
    mu(new double[N_body]),
    // Allocate a two dimensional array for each orbital element
    elt_a(      new double[N_row]),
    elt_e(      new double[N_row]),
    elt_inc(    new double[N_row]),
    elt_Omega(  new double[N_row]),
    elt_omega(  new double[N_row]),
    elt_f(      new double[N_row]),
    elt_M(      new double[N_row]),
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

    // Populate the gravitational field strength mu for each pair of primary / target
    // This array is aligned with body_id
    MassiveBodyTable mbt;
    for (int i=0; i<N_body; i++)
    {
        // The target body
        int32_t body_id = this->body_id[i];
        // The primary body is always the sun, EXCEPT for moon orbiting around the Earth
        int32_t body_id_primary = (body_id == body_id_moon) ? body_id_earth : body_id_sun;
        // Look up mass of primary and target body on the massive body table
        double m_pri = mbt.get_M(body_id_primary);
        double m_tgt = mbt.get_M(body_id);
        // The gravitational strength of this interaction
        mu[i] = G * (m_pri + m_tgt);
    }

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
PlanetElement::PlanetElement() :
    // Delegate to memory allocating constructor with database inputs 
    PlanetElement(mjd0_db, mjd1_db, stride_db_min)
{
    // Load data from disk
    load();

    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
PlanetElement::PlanetElement(db_conn_type& conn) :
    // Delegate to memory allocating constructor with database inputs 
    PlanetElement(mjd0_db, mjd1_db, stride_db_min)
{
    // Load data from database
    load(conn);

    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
/// Need a non-trivial destructor to release all manually allocated memory
PlanetElement::~PlanetElement()
{
    // Delete three 1D arrays; body_id and mu are aligned
    delete [] body_id;
    delete [] mu;
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
    // Run the stored procedure to get orbital elements of the planets
    string sp_name = "KS.GetElements_Planets";
    // Run the stored procedure to get planet elements in [t0, t1]
    vector<string> params = {to_string(t0), to_string(t1), to_string(dt_min)};
    ResultSet* rs = sp_run(conn, sp_name, params);

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

        // The row number j into the data arrays has two terms.
        // The base is the body_row, which is body index times the number of times each body.
        // The offset is the time_row, which is the number of rows written for earlier times on this body.
        // int j = body_row(body_id) + time_row(time_id);
        int j = row_id(body_id, time_id);

        // Save the data fields to the member arrays
        elt_a[j]      = a;
        elt_e[j]      = e;
        elt_inc[j]    = inc;
        elt_Omega[j]  = Omega;
        elt_omega[j]  = omega;
        elt_f[j]      = f;
        elt_M[j]      = M;
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
void PlanetElement::build_splines()
{
    // One interpolator for each body and each element
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
    // Bodies are laid out in the strategic order 1, 2, 301, 399, 4, 5, 6, 7, 8, 9
    // This corresponds to Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
    if ((0 < body_id) && (body_id < 3)) {return body_id-1;}
    else if ((3 < body_id) && (body_id < 10)) {return body_id;}
    else if (body_id==301) {return 8;}
    else if (body_id==399) {return 9;}
    else {throw domain_error(
        "Bad body_id! Must be one of 1, 2, 301, 399, 4, 5, 6, 7, 8, 9 (planet barycenters, Earth and Moon.\n");}
}

// *****************************************************************************
const int PlanetElement::body_row(int32_t body_id) const
{
    return body_idx(body_id)*N_t;
}

// *****************************************************************************
const int PlanetElement::time_row(int32_t time_id) const
    {return (time_id-time_id0)/dt_min;}

// *****************************************************************************
// Get 1D arrays of body_id and times
// *****************************************************************************

// *****************************************************************************
const int32_t* PlanetElement::get_body_id() const 
    {return body_id;}

// *****************************************************************************
const double* PlanetElement::get_mjd() const 
    {return mjd;}

// *****************************************************************************
// Get 1D array of each orbital element given a body_id
// *****************************************************************************

// *****************************************************************************
const double* PlanetElement::get_a(int idx) const 
    {return elt_a + N_t*idx;}

// *****************************************************************************
const double* PlanetElement::get_e(int idx) const 
    {return elt_e + N_t*idx;}

// *****************************************************************************
const double* PlanetElement::get_inc(int idx) const 
    {return elt_inc + N_t*idx;}

// *****************************************************************************
const double* PlanetElement::get_Omega(int idx) const 
    {return elt_Omega + N_t*idx;}

// *****************************************************************************
const double* PlanetElement::get_omega(int idx) const 
    {return elt_omega + N_t*idx;}

// *****************************************************************************
const double* PlanetElement::get_f(int idx) const 
    {return elt_f + N_t*idx;}

// *****************************************************************************
const double* PlanetElement::get_M(int idx) const 
    {return elt_M + N_t*idx;}

// *****************************************************************************
const OrbitalElement PlanetElement::interp_elt_by_idx(int idx, double mjd) const
{
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
    // Replace splined f value with conversion from M (M splines better than f)
    elt.f = anomaly_M2f(elt.M, elt.e);
    // Return the splined orbital element
    return elt;
}

// *****************************************************************************
const OrbitalElement PlanetElement::interp_elt(int32_t body_id, double mjd) const
{
    // Get the index number of this body
    int idx = body_idx(body_id);

    // Delegate to method keyed by body index
    return interp_elt_by_idx(idx, mjd);
}

// *****************************************************************************
// This calculation returns the position RELATIVE to the primary, NOT in the BME!
// The orbital element describes the relative position of the body vs. the primary, which is the Sun here.
const Position PlanetElement::interp_pos_rel(int32_t body_id, double mjd) const
{
    // Peel off special case that the body is the sun; then the position is zero
    // if (body_id==body_id_sun){return pos_sun_hel;}

    // Delegate to interp_elt to spline the elements
    OrbitalElement elt = interp_elt(body_id, mjd);

    // Call elt2pos to calculate a position. 
    // Use the version elt2pos_vec to guarantee identical positions when using either method.
    return elt2pos_vec(elt);
}

// *****************************************************************************
// This calculation returns the state vector RELATIVE to the primary, NOT in the BME!
// The orbital element describes the relative vectors of the body vs. the primary, which is the Sun here.
const StateVector PlanetElement::interp_vec_rel(int32_t body_id, double mjd) const
{
    // Peel off special case that the body is the sun; then the state vector is zero
    // if (body_id==body_id_sun){return vec_sun_hel;}

    // Get the index number of this body
    int idx = body_idx(body_id);

    // Delegate to interp_elt to spline the elements
    OrbitalElement elt = interp_elt_by_idx(idx, mjd);

    // Convert to state vectors using the gravitational field strength mu for this interaction
    return elt2vec(elt, mu[idx]);
}

// *****************************************************************************
// Add the heliocentric position from splined orbital elements to the Sun's state vectors
// to get the position in the BME frame.
const Position PlanetElement::interp_pos(int32_t body_id, double mjd) const
{
    // Delegate to bv_sun to get interpolated position of Sun in BME
    Position sun = bv_sun.interp_pos(mjd);
    // Delegate to interp_pos_rel for the relative position of body vs. the primary
    Position tgt = interp_pos_rel(body_id, mjd);
    // Peel off the special case that the body is the moon; then the primary is Earth, not the Sun
    if (body_id==body_id_moon)
    {
        Position earth = sun + interp_pos_rel(body_id_earth, mjd);
        return earth + tgt;
    }
    // Add the position of the primary (sun) and the relative position
    return sun + tgt;
}

// *****************************************************************************
// Add the heliocentric state vectors from splined orbital elements to the Sun's state vectors
// to get the state vectors in the BME frame.
const StateVector PlanetElement::interp_vec(int32_t body_id, double mjd) const
{
    // Delegate to bv_sun to get interpolated position of Sun in BME
    StateVector sun = bv_sun.interp_vec(mjd);
    // Delegate to interp_vec_rel for the relative state vector of body vs. the primary
    StateVector tgt = interp_vec_rel(body_id, mjd);
    // Peel off the special case that the body is the moon; then the primary is Earth, not the Sun
    if (body_id==body_id_moon)
    {
        StateVector earth = sun + interp_vec_rel(body_id_earth, mjd);
        return earth + tgt;
    }
    // Add the position of the primary (sun) and the relative position
    return sun + tgt;
}

// *****************************************************************************
/// Structure for one row of serializing the PlanetElement data
struct PlanetElementEntry
{
    int32_t body_id;
    int32_t time_id;
    double mjd;
    double a;
    double e;
    double inc;
    double Omega;
    double omega;
    double f;
    double M;
};

namespace 
{
/// Size of PlanetElementEntry; used in save() and load() methods.
int entry_sz = sizeof(PlanetElementEntry);
// Wrap in anonymous namespace to avoid naming collisions with other modules.
}

// *****************************************************************************
void PlanetElement::save() const
{
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of bodies and times in binary as int
    fs.write((char*) &N_body, sizeof(N_body));
    fs.write((char*) &N_t, sizeof(N_t));
    // Write the date range; used to support loading in a time interval
    fs.write((char*) &mjd0, sizeof(mjd0));
    fs.write((char*) &mjd1, sizeof(mjd1));
    fs.write((char*) &dt_min, sizeof(dt_min));

    // Write the rows to the file in binary
    // Iterate over the N_body bodies in this collection
    for (int idx=0; idx<N_body; idx++)
    {
        // The body_id for this idx
        int32_t body_id = (this->body_id)[idx];
        // The row offset for this body
        int body_row = N_t * idx;
        // Iterate over the N_t times
        for (int i=0; i<N_t; i++)
        {
            // Calculate the time_id
            int32_t time_id = time_id0 + i*dt_min;
            // The array address for the 2D data elements
            int j = body_row + i;
            // Initialize a PlanetElementEntry
            PlanetElementEntry pe
            {
                .body_id    =body_id,
                .time_id    =time_id,
                .mjd        =mjd[i],
                .a          =elt_a[j],
                .e          =elt_e[j],
                .inc        =elt_inc[j],
                .Omega      =elt_Omega[j],
                .omega      =elt_omega[j],
                .f          =elt_f[j],
                .M          =elt_M[j]
            };
            // Write the PlanetElementEntry object
            fs.write((char*) &pe, entry_sz);
        } // for / i (loop over times)
    } // for / idx (loop over bodies)

    // Close output filestream
    fs.close();
} // function save()

// *****************************************************************************
void PlanetElement::load()
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of bodies and times
    auto N_body_file=N_body;
    auto N_t_file=N_t;
    fs.read( (char*) &N_body_file, sizeof(N_body));
    fs.read( (char*) &N_t_file, sizeof(N_t));

    // Read the date range
    auto mjd0_file=mjd0;
    auto mjd1_file=mjd1;
    auto dt_min_file=dt_min;
    fs.read( (char*) &mjd0_file, sizeof(mjd0));
    fs.read( (char*) &mjd1_file, sizeof(mjd1));
    fs.read( (char*) &dt_min_file, sizeof(dt_min));

    // Check that the number of rows agrees with the constexpr specification in this file
    if (N_body_file != N_body_planets)
    {throw domain_error(
        format("Bad data file! N_body={:d} does not match specification {:d}.\n", N_body_file, N_body_planets));}
    if (N_t_file != N_t_db)
    {throw domain_error(
        format("Bad data file! N_t={:d} does not match specification {:d}.\n", N_t_file, N_t_db));}

    // Number of rows in the file
    int N_row_file = N_body_file * N_t_file;
    // PlanetElementEntry object used to read rows
    PlanetElementEntry pe;

    // Iterate through the rows in the file
    for (int i=0; i<N_row_file; i++)
    {
        // Read the PlanetElementEntry
        fs.read( (char*) &pe, entry_sz);
        // Only load dates in the selected time_id range
        bool in_range = (time_id0 <= pe.time_id) && (pe.time_id <= time_id1);
        bool in_stride = (pe.time_id % dt_min == 0);
        bool process_row = in_range && in_stride;
        if (!process_row ) {continue;}
        // Calculate the row index from the ID fields in the entry
        int j = body_row(pe.body_id) + (pe.time_id-time_id0)/dt_min;

        // Save the data fields to the member arrays with the calculated row index
        elt_a[j]      = pe.a;
        elt_e[j]      = pe.e;
        elt_inc[j]    = pe.inc;
        elt_Omega[j]  = pe.Omega;
        elt_omega[j]  = pe.omega;
        elt_f[j]      = pe.f;
        elt_M[j]      = pe.M;
    } // for / i (loop over rows in the file)

    // Close input filestream
    fs.close();
}
