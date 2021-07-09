/** @file   PlanetVector.cpp
 *  @brief  Implmentation of PlanetVector class.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-09
 */

// *****************************************************************************
// Local dependencies
#include "PlanetVector.hpp"

// *****************************************************************************
// Constants used in this module

// Number of minutes in one day
using ks::cs::mpd;
// Number of days in one minute
using ks::cs::dpm;

// Number of massive bodies in Planets collection: Sun + 9 planets + moon = 11 bodies
constexpr int N_body_planets = 11;
// Array of body_id in this problem
constexpr int32_t body_id_planets[N_body_planets] = {1, 2, 4, 5, 6, 7, 8, 9, 10, 301, 399};

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
const string file_name = "data/cache/PlanetVector.bin";

// *****************************************************************************
// Local names used
using ks::PlanetVector;

// *****************************************************************************
// The constructor just allocates memory.  
// It does not load data from database, that is done with the load() method.
PlanetVector::PlanetVector(int mjd0, int mjd1, int dt_min):
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
    // Allocate a two dimensional array for each state vector component
    qx(new double[N_row]),
    qy(new double[N_row]),
    qz(new double[N_row]),
    vx(new double[N_row]),
    vy(new double[N_row]),
    vz(new double[N_row]),
    // Initialize GSL objects
    acc(gsl_interp_accel_alloc() )
{
    // Populate body_id - this is a copy from body_id_planets
    for (int i=0; i<N_body; i++) {body_id[i]=body_id_planets[i];}
    // Populate mjd; these are on a fixed schedule spaced dt_min minutes apart from mjd0 to mjd1
    const double dt = dt_min * dpm;
    for (int i=0; i<N_t; i++) {mjd[i] = mjd0 + i*dt;}

    // vec_spline is a structure with one member for each state vector component
    // Reserve space in vector of splines for each element    
    vec_spline.qx.reserve(N_body);
    vec_spline.qy.reserve(N_body);
    vec_spline.qz.reserve(N_body);
    vec_spline.vx.reserve(N_body);
    vec_spline.vy.reserve(N_body);
    vec_spline.vz.reserve(N_body);

    // Initialize the splines for the elements of each body; one spline for each body and component
    for (int i=0; i<N_body; i++)
    {
        vec_spline.qx.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        vec_spline.qy.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        vec_spline.qz.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        vec_spline.vx.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        vec_spline.vy.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
        vec_spline.vz.push_back(gsl_spline_alloc(gsl_interp_cspline, N_t));
    }   // for / i (loop over massive bodies)
} // end function

// *****************************************************************************
PlanetVector::PlanetVector() :
    // Delegate to memory allocating constructor with database inputs 
    PlanetVector(mjd0_db, mjd1_db, stride_db_min)
{
    // Load data from disk
    load();

    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
PlanetVector::PlanetVector(db_conn_type& conn) :
    // Delegate to memory allocating constructor with database inputs 
    PlanetVector(mjd0_db, mjd1_db, stride_db_min)
{
    // Load data from database
    load(conn);

    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
/// Need a non-trivial destructor to release all manually allocated memory
PlanetVector::~PlanetVector()
{
    // Delete two 1D arrays
    delete [] body_id;
    delete [] mjd;

    // Delete six 2D arrays, one for each state vector component
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
// (1) One spline for each asteroid / element pair
// (2) The interpolation accelerator
void PlanetVector::gsl_free()
{
    // One interpolator for each body and each element
    for (int i=0; i<N_body; i++)
    {
        gsl_spline_free(vec_spline.qx[i]);
        gsl_spline_free(vec_spline.qy[i]);
        gsl_spline_free(vec_spline.qz[i]);
        gsl_spline_free(vec_spline.vx[i]);
        gsl_spline_free(vec_spline.vy[i]);
        gsl_spline_free(vec_spline.vz[i]);
    }

    // Just one accelerator object
    gsl_interp_accel_free(acc);
}

// *****************************************************************************
void PlanetVector::load(db_conn_type &conn)
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
void PlanetVector::process_rows(db_conn_type& conn, int t0, int t1)
{
    // Run the stored procedure to get state vectors of the planets
    string sp_name = "KS.GetStateVectors_Planets";
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
        // Six state vector components
        double qx = rs->getDouble("qx");
        double qy = rs->getDouble("qy");
        double qz = rs->getDouble("qz");
        double vx = rs->getDouble("vx");
        double vy = rs->getDouble("vy");
        double vz = rs->getDouble("vz");

        // The row number j into the data arrays has two terms.
        // The base is the body_row, which is body index times the number of times each body.
        // The offset is the time_row, which is the number of rows written for earlier times on this body.
        // int j = body_row(body_id) + time_row(time_id);
        int j = row_id(body_id, time_id);

        // Save the data fields to the member arrays
        this->qx[j] = qx;
        this->qy[j] = qy;
        this->qz[j] = qz;
        this->vx[j] = vx;
        this->vy[j] = vy;
        this->vz[j] = vz;
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
void PlanetVector::build_splines()
{
    // One interpolator for each body and each element
    // idx is the body index (i.e. counter from 0 to 10), NOT the body_id!
    for (int idx=0; idx<N_body; idx++)
    {
        // Evaluate the interpolation spline for each state vector component of this body
        gsl_spline_init(vec_spline.qx[idx],  mjd, get_qx(idx), N_t);
        gsl_spline_init(vec_spline.qy[idx],  mjd, get_qy(idx), N_t);
        gsl_spline_init(vec_spline.qz[idx],  mjd, get_qz(idx), N_t);
        gsl_spline_init(vec_spline.vx[idx],  mjd, get_vx(idx), N_t);
        gsl_spline_init(vec_spline.vy[idx],  mjd, get_vy(idx), N_t);
        gsl_spline_init(vec_spline.vz[idx],  mjd, get_vz(idx), N_t);
    }
}

// *****************************************************************************
const int PlanetVector::body_idx(int32_t body_id) const
{
    // Bodies are laid out in the order 1, 2, 4, 5, 6, 7, 8, 9, 10, 301, 399
    if ((0 < body_id) && (body_id < 3)) {return body_id-1;}
    else if ((3 < body_id) && (body_id < 11)) {return body_id-2;}
    else if (body_id==301) {return 9;}
    else if (body_id==399) {return 10;}
    else {throw domain_error(
        "Bad body_id! Must be one of 1, 2, 4, 5, 6, 7, 8, 9, 10, 301, 399 "
        "(sun, planet barycenters, earth and moon.\n");}
}

// *****************************************************************************
const int PlanetVector::body_row(int32_t body_id) const
{
    return body_idx(body_id)*N_t;
}

// *****************************************************************************
const int PlanetVector::time_row(int32_t time_id) const
{
    return (time_id-time_id0)/dt_min;
}

// *****************************************************************************
// Get 1D arrays of body_id and times
// *****************************************************************************

// *****************************************************************************
int32_t* PlanetVector::get_body_id() const {return body_id;}

// *****************************************************************************
double* PlanetVector::get_mjd() const {return mjd;}

// *****************************************************************************
// Get 1D array of each orbital element given a body_id
// *****************************************************************************

// *****************************************************************************
double* PlanetVector::get_qx(int idx) const {return qx + N_t*idx;}

// *****************************************************************************
double* PlanetVector::get_qy(int idx) const {return qy + N_t*idx;}

// *****************************************************************************
double* PlanetVector::get_qz(int idx) const {return qz + N_t*idx;}

// *****************************************************************************
double* PlanetVector::get_vx(int idx) const {return vx + N_t*idx;}

// *****************************************************************************
double* PlanetVector::get_vy(int idx) const {return vy + N_t*idx;}

// *****************************************************************************
double* PlanetVector::get_vz(int idx) const {return vz + N_t*idx;}

// *****************************************************************************
// Return the interpolated position in the BME frame of this body at the input time.
Position PlanetVector::interp_pos(int32_t body_id, double mjd) const
{
    // Get the index number of this body
    int idx = body_idx(body_id);

    // The interpolators for each orbital element for this body
    // The array index into the elt_spline members is idx, NOT body_id;
    // body_id will go over the end, e.g. for Earth body_id=399 and idx=9.
    gsl_spline* gsl_interp_qx = vec_spline.qx[idx];
    gsl_spline* gsl_interp_qy = vec_spline.qy[idx];
    gsl_spline* gsl_interp_qz = vec_spline.qz[idx];

    // Evaluate the splines at the selected time
    return Position 
    {
        .qx= gsl_spline_eval(gsl_interp_qx, mjd, acc), 
        .qy= gsl_spline_eval(gsl_interp_qy, mjd, acc), 
        .qz= gsl_spline_eval(gsl_interp_qz, mjd, acc)
    };
}

// *****************************************************************************
// Return the interpolated state vector in the BME frame of this body at the input time.
StateVector PlanetVector::interp_vec(int32_t body_id, double mjd) const
{
    // Get the index number of this body
    int idx = body_idx(body_id);

    // The interpolators for each orbital element for this body
    // The array index into the elt_spline members is idx, NOT body_id;
    // body_id will go over the end, e.g. for Earth body_id=399 and idx=9.
    gsl_spline* gsl_interp_qx = vec_spline.qx[idx];
    gsl_spline* gsl_interp_qy = vec_spline.qy[idx];
    gsl_spline* gsl_interp_qz = vec_spline.qz[idx];
    gsl_spline* gsl_interp_vx = vec_spline.vx[idx];
    gsl_spline* gsl_interp_vy = vec_spline.vy[idx];
    gsl_spline* gsl_interp_vz = vec_spline.vz[idx];

    // Evaluate the splines at the selected time
    return StateVector 
    {
        .qx= gsl_spline_eval(gsl_interp_qx, mjd, acc), 
        .qy= gsl_spline_eval(gsl_interp_qy, mjd, acc), 
        .qz= gsl_spline_eval(gsl_interp_qz, mjd, acc),
        .vx= gsl_spline_eval(gsl_interp_vx, mjd, acc), 
        .vy= gsl_spline_eval(gsl_interp_vy, mjd, acc), 
        .vz= gsl_spline_eval(gsl_interp_vz, mjd, acc)
    };
}

// *****************************************************************************
/// Structure for one row of serializing the PlanetElement data
struct PlanetVectorEntry
{
    int32_t body_id;
    int32_t time_id;
    double mjd;
    double qx;
    double qy;
    double qz;
    double vx;
    double vy;
    double vz;
};

// Size of this
int entry_sz = sizeof(PlanetVectorEntry);

// *****************************************************************************
void PlanetVector::save() const
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
            PlanetVectorEntry pe
            {
                .body_id = body_id,
                .time_id = time_id,
                .mjd     = mjd[i],
                .qx      = qx[j],
                .qy      = qy[j],
                .qz      = qz[j],
                .vx      = vx[j],
                .vy      = vy[j],
                .vz      = vz[j]
            };
            // Write the PlanetElementEntry object
            fs.write((char*) &pe, entry_sz);
        } // for / i (loop over times)
    } // for / idx (loop over bodies)

    // Close output filestream
    fs.close();
} // function save()

// *****************************************************************************
void PlanetVector::load()
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
    // PlanetVectorEntry object used to read rows
    PlanetVectorEntry pe;

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
        qx[j] = pe.qx;
        qy[j] = pe.qy;
        qz[j] = pe.qz;
        vx[j] = pe.vx;
        vy[j] = pe.vy;
        vz[j] = pe.vz;
    } // for / i (loop over rows in the file)

    // Close input filestream
    fs.close();
}
