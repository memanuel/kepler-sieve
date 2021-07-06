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
// Number of minutes in one day
constexpr int mpd = 1440;
// One day as a floating point number of minutes
constexpr double dpm = 1.0 / mpd;
// Start and end date for database
constexpr int mjd0_db = 48000;
constexpr int mjd1_db = 63000;
// Stride in minutes for database vectors for Sun and Earth
constexpr int stride_db_min = 5;
// Batch size for loading dates from DB
constexpr int batch_size = 1000;
// Location of file with serialized data
const string file_name_fmt = "data/cache/BodyVector_{:s}.bin";

// *****************************************************************************
// Local names used in this module
using ks::BodyVector;

// *****************************************************************************
// The constructor just allocates memory and initializes GSL splines.  
// It does not load data from database or file, that is done with the load() or load_db() methods.
BodyVector::BodyVector(int mjd0, int mjd1, int dt_min, string body_name):
    // Name of the body
    body_name(body_name),
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
    // Validate the body_name field

    // Initialize the splines for the six state vector components
    vec_spline.qx = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.qy = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.qz = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.vx = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.vy = gsl_spline_alloc(gsl_interp_cspline, N_t);
    vec_spline.vz = gsl_spline_alloc(gsl_interp_cspline, N_t);

} // end function

// *****************************************************************************
BodyVector::BodyVector(string body_name) :
    // Delegate to memory allocating constructor with database inputs 
    BodyVector(mjd0_db, mjd1_db, stride_db_min, body_name)
{
    // Load data from disk
    load();

    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
BodyVector::BodyVector(db_conn_type& conn, string body_name) :
    // Delegate to memory allocating constructor with database inputs 
    BodyVector(mjd0_db, mjd1_db, stride_db_min, body_name)
{
    // Load data from database
    load_db(conn);

    // Initialize the gsl_spline objects
    build_splines();
}

// *****************************************************************************
const string BodyVector::sp_name_from_body()
{
    // List of supported body names
    std::array<std::string, 2> body_names= {"Sun", "Earth"};
    // Count number of matches
    int match_count = std::count(std::begin(body_names), std::end(body_names), body_name);
    // If we get a match, return the SP name
    if (match_count > 0)
        {return format("KS.GetStateVectors_{:s}", body_name);}
    // If it's not a match, throw a domain error
    else
        {throw domain_error("Bad body_name! Must be one of \"Sun\", \"Earth\".\n");}
}

// *****************************************************************************
const string BodyVector::file_name_from_body()
{
    // List of supported body names
    std::array<std::string, 2> body_names= {"Sun", "Earth"};
    // Count number of matches
    int match_count = std::count(std::begin(body_names), std::end(body_names), body_name);
    // If we get a match, return the file_name
    if (match_count > 0)
        {return format(file_name_fmt, body_name);}
    // If it's not a match, throw a domain error
    else
        {throw domain_error("Bad body_name! Must be one of \"Sun\", \"Earth\".\n");}
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
void BodyVector::load_db(db_conn_type &conn)
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
    string sp_name = sp_name_from_body();

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
        .qy = gsl_spline_eval(vec_spline.qy, mjd, acc),
        .qz = gsl_spline_eval(vec_spline.qz, mjd, acc)
    };
}

// *****************************************************************************
StateVector BodyVector::interp_vec(double mjd)
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
void BodyVector::save()
{
    // Build file_name from file_name_base and body_name
    string file_name = file_name_from_body();
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
    string file_name = file_name_from_body();
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows according to the file
    int N_t_file=-1;
    fs.read( (char*) &N_t_file, sizeof(N_t_file));
    // Check that the number of rows agrees with the constexpr specification in this file
    if (N_t_file != N_t)
    {throw domain_error("Bad data file! N_t does not match specification.\n");}

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

// *****************************************************************************
// Additional functions for working with BodyVector objects
// *****************************************************************************

/// Helper function - build and save vectors
void ks::save_vectors(string body_name)
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Initialize BodyVector for this body using DataBase
    print("Loading BodyVector for {:s} from DB...\n", body_name);
    BodyVector bv(conn, body_name);
    // Save to disk
    bv.save();
    print("Saved BodyVector to disk for {:s}.\n", body_name);

    // Close DB connection
    conn->close();
}

// *****************************************************************************
/// Factory function - load vectors from disk
// BodyVector ks::load_vectors(string body_name)
// {
//     // Build file_name from file_name_base and body_name
//     string file_name = format("KS.GetStateVectors_{:s}", body_name);
//     // Open input filestream in binary mode
//     std::ifstream fs;
//     std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
//     fs.open(file_name, file_mode);

//     // Read the number of rows
//     int N_t_file=-1;
//     fs.read( (char*) &N_t_file, sizeof(N_t_file));

// }
