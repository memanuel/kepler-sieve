/** @file BodyVector.cpp
 *  @brief Implmentation of BodyVector class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-05
 */

// *****************************************************************************
// Included files
#include "BodyVector.hpp"

// *****************************************************************************
// Constants used in this module
// Batch size for loading dates from DB
constexpr int batch_size = 1000;

// List of supported body names
constexpr std::array<const char*, 2> body_names_bv = {"Sun", "Earth"};

// Location of file with serialized data
constexpr const char* file_name_sun   = "data/cache/BodyVector_Sun.bin";
constexpr const char* file_name_earth = "data/cache/BodyVector_Earth.bin";
// constexpr std::array<const char*, 2> file_names_bv = {file_name_sun, file_name_earth};

// Name of stored procedure as format string to be evaluated with body_name
constexpr const char* sp_name_sun   = "KS.GetStateVectors_Sun";
constexpr const char* sp_name_earth = "KS.GetStateVectors_Earth";

// *****************************************************************************
// Local names used in this module
using ks::BodyVector;

// *****************************************************************************
// The constructor just allocates memory and initializes GSL splines.  
// It does not load data from database or file, that is done with the load() method.
BodyVector::BodyVector(SolarSystemBody_bv ssb_, int mjd0_, int mjd1_, int dt_min_):
    // The body and its name
    ssb {ssb_},
    body_name {get_body_name(ssb_)},
    // number of times (data size)
    N_t {(mjd1_-mjd0_)*mpd/dt_min_+1},
    // Date range
    mjd0 {mjd0_},
    mjd1 {mjd1_},
    // Stride in minutes between consecutive entries
    // database resolution is every 5 minutes, so apply LCM function to the input to ensure validity.
    dt_min {lcm(dt_min_, stride_db_min)},
    // Allocate array for mjd
    mjd {new double[N_t]},
    // Allocate arrays for each component of the state vector
    qx {new double[N_t]},
    qy {new double[N_t]},
    qz {new double[N_t]},
    vx {new double[N_t]},
    vy {new double[N_t]},
    vz {new double[N_t]},
    // Initialize the splines for the six state vector components
    vec_spline 
    {
        .qx = gsl_spline_alloc(gsl_interp_cspline, N_t),
        .qy = gsl_spline_alloc(gsl_interp_cspline, N_t),
        .qz = gsl_spline_alloc(gsl_interp_cspline, N_t),
        .vx = gsl_spline_alloc(gsl_interp_cspline, N_t),
        .vy = gsl_spline_alloc(gsl_interp_cspline, N_t),
        .vz = gsl_spline_alloc(gsl_interp_cspline, N_t),
    },
    // Initialize GSL accelerator
    acc {gsl_interp_accel_alloc()}
    {} // end function

// *****************************************************************************
BodyVector::BodyVector(SolarSystemBody_bv ssb, int mjd0, int mjd1, int dt_min, bool load_) :
    // Delegate to memory allocating constructor with these inputs 
    BodyVector(ssb, mjd0, mjd1, dt_min)
{
    if (load_)
    {
        // Load data from disk
        load();
        // Initialize the gsl_spline objects
        build_splines();
    }
}

// *****************************************************************************
BodyVector::BodyVector(SolarSystemBody_bv ssb) :
    // Delegate to memory allocating constructor with database inputs and load_ = true
    BodyVector(ssb, mjd0_db, mjd1_db, stride_db_min, true)
    {}

// *****************************************************************************
BodyVector::BodyVector(SolarSystemBody_bv ssb, db_conn_type& conn) :
    // Delegate to memory allocating constructor with database inputs 
    BodyVector(ssb, mjd0_db, mjd1_db, stride_db_min)
{
    // Load data from database
    load(conn);
    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
const char* BodyVector::sp_name_() const
{
    switch (ssb)
    {
        case SolarSystemBody_bv::sun:      return sp_name_sun;
        case SolarSystemBody_bv::earth:    return sp_name_earth;
    }   // switch (ssb)
    // Impossible to get here b/c all enum cases covered above.
    // Put something to suppress compiler warning
    throw invalid_argument("get_body_name - bad body enum!");
}

// *****************************************************************************
const char* BodyVector::file_name_() const
{
    switch (ssb)
    {
        case SolarSystemBody_bv::sun:      return file_name_sun;
        case SolarSystemBody_bv::earth:    return file_name_earth;
    }   // switch (ssb)
    // Impossible to get here b/c all enum cases covered above.
    // Put something to suppress compiler warning
    throw invalid_argument("get_body_name - bad body enum!");
}

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
void BodyVector::load(db_conn_type &conn)
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
void BodyVector::process_rows(db_conn_type& conn, int t0, int t1)
{
    // Choose correct SP name from body_name
    const string sp_name = string(sp_name_());

    // Run the stored procedure to get state vectors in [t0, t1]
    vector<string> params = {to_string(t0), to_string(t1), to_string(dt_min)};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // The time_id corresponding to mjd0
    const int time_id0 = mpd * mjd0;

    // Loop through resultset
    while (rs->next()) 
    {   // The TimeID
        int32_t time_id = rs->getInt("TimeID");
        // The row where this data is written
        int i = (time_id-time_id0)/dt_min;
        // Write to mjd array
        mjd[i] = rs->getDouble("mjd");
        // Write state vector components to arrays
        qx[i] = rs->getDouble("qx");
        qy[i] = rs->getDouble("qy");
        qz[i] = rs->getDouble("qz");
        vx[i] = rs->getDouble("vx");
        vy[i] = rs->getDouble("vy");
        vz[i] = rs->getDouble("vz");
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;   
}   // end function

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
const double* BodyVector::get_mjd() const
    {return mjd;}

// *****************************************************************************
Position BodyVector::interp_pos(double mjd) const
{
    return Position
    {
        .qx = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .qy = gsl_spline_eval(vec_spline.qy, mjd, acc),
        .qz = gsl_spline_eval(vec_spline.qz, mjd, acc)
    };
}

// *****************************************************************************
StateVector BodyVector::interp_vec(double mjd) const
{
    return StateVector
    {
        .qx = gsl_spline_eval(vec_spline.qx, mjd, acc),
        .qy = gsl_spline_eval(vec_spline.qy, mjd, acc),
        .qz = gsl_spline_eval(vec_spline.qz, mjd, acc),
        .vx = gsl_spline_eval(vec_spline.vx, mjd, acc),
        .vy = gsl_spline_eval(vec_spline.vy, mjd, acc),
        .vz = gsl_spline_eval(vec_spline.vz, mjd, acc)
    };
}

// *****************************************************************************
// Save to and load from disk (much faster than DB)
// *****************************************************************************

// *****************************************************************************
void BodyVector::save() const
{
    // Build file_name from file_name_base and body_name
    const char* file_name = file_name_();
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of rows in binary as long int
    fs.write((char*) &N_t, sizeof(N_t));

    // The time_id corresponding to mjd0
    // const int time_id0 = mpd * mjd0;
    // All the data elements except time_id are doubles
    int data_sz = sizeof(double);
    // Write the rows to the file in binary
    for (int i=0; i<N_t; i++)
    {
        // Write the time as a time_id
        // int32_t time_id = time_id0 + i*dt_min;
        // fs.write((char*) &time_id, sizeof(time_id));

        // Write the time as an mjd
        fs.write((char*) (mjd+i), data_sz);
        // Write the six components of the state vector
        fs.write((char*) (qx+i), data_sz);
        fs.write((char*) (qy+i), data_sz);
        fs.write((char*) (qz+i), data_sz);
        fs.write((char*) (vx+i), data_sz);
        fs.write((char*) (vy+i), data_sz);
        fs.write((char*) (vz+i), data_sz);
    }

    // Status
    // print("Wrote {:d} rows of data to {:s}\n", N_t, file_name);

    // Close output filestream
    fs.close();
}

// *****************************************************************************
void BodyVector::load()
{
    // Build file_name from file_name_base and body_name
    const char* file_name = file_name_();
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows according to the file
    int N_t_file=-1;
    fs.read( (char*) &N_t_file, sizeof(N_t_file));
    // Check that the number of rows agrees with the constexpr specification in this file
    if (N_t_file != N_t)
    {
        string msg = format(
            "BodyVector::load - Bad data file! N_t={:d} does not match specification {:d}.\n", N_t_file, N_t);
        // DEBUG
        print("file_name");
        throw runtime_error(msg);
    }

    // All the data elements except time_id are doubles
    int data_sz = sizeof(double);

    // The time_id corresponding to mjd0
    // const int time_id0 = mpd * mjd0;
    // The time_id on this loop
    // int32_t time_id

    // Loop through the file
    for (int i=0; i<N_t; i++)
    {
        // Read the time_id
        // fs.read( (char*) &time_id, sizeof(time_id));
        // Write to mjd array
        fs.read( (char*) (mjd+i), data_sz);
        // Write state vector components to arrays
        fs.read( (char*) (qx+i), data_sz);
        fs.read( (char*) (qy+i), data_sz);
        fs.read( (char*) (qz+i), data_sz);
        fs.read( (char*) (vx+i), data_sz);
        fs.read( (char*) (vy+i), data_sz);
        fs.read( (char*) (vz+i), data_sz);
    }

    // Close input filestream
    fs.close();
}   
