/** @file DetectionTime.cpp
 *  @brief Implmentation of DetectionTime class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-07
 */

// *****************************************************************************
// Included files
#include "DetectionTime.hpp"

// *****************************************************************************
// Local names used
using ks::DetectionTime;
using ks::DetectionTimeTable;
// Declare helper functions
int file_length();

// Location of file with serialized data
namespace {
const string file_name = "data/cache/DetectionTimeTable.bin";
}

// *****************************************************************************
// The "real" constructor that initializes all the members
DetectionTimeTable::DetectionTimeTable(int N):
    N {N},
    // dtv is a vector of DetectionTimes
    // size is N+1 because we index by detection_time_id, which starts at 1.
    dtv {vector<DetectionTime>(N+1)},
    // dtm is an empty map keyed by TimeID; it will be populated later
    dtm {map<int32_t, vector<int32_t> >{}},
    // mjds is an array of doubles of size N+1 (size of dtv)
    mjd_ {new double[N+1]},
    // q_obs_ is an array of doubles with 3(N+1) entries for the HELIOCENTRIC observatory position
    q_obs_ {new double[3*(N+1)]},
    // mjd is a read-only copy of mjd_ that is publicly accessible
    mjd {const_cast<double*>(mjd_)},
    // q_obs is a read-only copy of q_obs_ that is publicly accessible
    q_obs {const_cast<double*>(q_obs_)}
    {}

// *****************************************************************************
// The default constructor loads from the disk; it needs to know row count before delegating.
DetectionTimeTable::DetectionTimeTable():
    // Delegate to size constructor, looking up the file size
    DetectionTimeTable {file_length()}
{
    // Load data from disk
    load();
    // Write out arrays mjd_ and q_obs_
    write_arrays();
}

// *****************************************************************************
DetectionTimeTable::DetectionTimeTable(db_conn_type& conn):
    DetectionTimeTable(sp_run_int(conn, "KS.GetDetectionTimeMaxID"))
    {
        // Now load data from the database
        load(conn);
        // Write out arrays mjd_ and q_obs_
        write_arrays();
    }

// *****************************************************************************
// Delete manually allocated arrays
DetectionTimeTable::~DetectionTimeTable()
{
    delete [] mjd_;
    delete [] q_obs_;
}

// *****************************************************************************
void DetectionTimeTable::load(db_conn_type& conn)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetDetectionTimes";
    vector<string> params = {};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next())
    {
        // Get the DetectionTimeID and TimeID
        int32_t detection_time_id = rs->getInt("DetectionTimeID");
        int32_t time_id = rs->getInt("TimeID");
        // Get the doubles
        double mjd = rs->getDouble("mjd");
        double q_obs_x = rs->getDouble("qObs_x");
        double q_obs_y = rs->getDouble("qObs_y");
        double q_obs_z = rs->getDouble("qObs_z");
        double q_sun_x = rs->getDouble("qSun_x");
        double q_sun_y = rs->getDouble("qSun_y");
        double q_sun_z = rs->getDouble("qSun_z");
        // Wrap this entry into one DetectionTime
        DetectionTime dt = DetectionTime
        {
            .detection_time_id = detection_time_id,
            .time_id = time_id,
            .data_source_id = rs->getByte("DataSourceID"),
            .observatory_id = rs->getByte("ObservatoryID"),
            .mjd        = mjd,
            .q_obs_x    = q_obs_x,
            .q_obs_y    = q_obs_y,
            .q_obs_z    = q_obs_z,
            .q_sun_x    = q_sun_x,
            .q_sun_y    = q_sun_y,
            .q_sun_z    = q_sun_z
        };

        // Write this DetectionTime into the vector dtv
        dtv[detection_time_id] = dt;
        // Add this detection time to the map keyed by TimeID
        (dtm[time_id]).push_back(detection_time_id);
    }
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
const DetectionTime& DetectionTimeTable::operator[](int32_t id) const
    {return dtv[id];}

// *****************************************************************************
const vector<int32_t> DetectionTimeTable::get_time(int32_t time_id) const
    {return dtm.at(time_id);}

// *****************************************************************************
const vector<DetectionTime> DetectionTimeTable::detection_times() const
    {return dtv;}

// *****************************************************************************
void DetectionTimeTable::save()
{
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of rows in binary as long int
    fs.write((char*) &N, sizeof(N));

    // Write the rows to the file in binary
    for (DetectionTime dt : dtv)
    {
        fs.write((char*) &dt, sizeof(dt));
    }

    // Close output filestream
    fs.close();
}

// *****************************************************************************
const int DetectionTimeTable::file_length() const
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);
    // Read the number of rows in the file
    int N_file=-1;
    fs.read( (char*) &N_file, sizeof(N_file));
    // Close input filestream
    fs.close();
    return N_file;
}

// *****************************************************************************
void DetectionTimeTable::load()
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows in the file
    int N_file=-1;
    fs.read( (char*) &N_file, sizeof(N_file));
    // Need to read one extra row because we index from detection_id which starts at 1
    int sz = N_file +1;
    // Status
    // print("DetectionTimeTable::load() - Opened file {:s} with {:d} rows of DetectionTime data (N={:d}).\n", 
    //     file_name, sz, N_file);

    // Read the rows from the file in binary
    DetectionTime dt;
    for (int i=0; i<sz; i++)
    {
        // Read this row into dt
        fs.read( (char*) &dt, sizeof(dt));
        // Get the detection_time_id
        int32_t detection_time_id = dt.detection_time_id;
        // Write this DetectionTime into the vector dtv
        dtv[detection_time_id] = dt;
        // Add this detection time to the map keyed by TimeID
        (dtm[dt.time_id]).push_back(detection_time_id);
    }

    // Write out the arrays mjd_ and q_obs_
    write_arrays();

    // Close input filestream
    fs.close();
}

// *****************************************************************************
void DetectionTimeTable::write_arrays()
{
    // Iterate through the detection times in the vector 
    for (DetectionTime dt : dtv)
    {
        // Save the mjd to the array mjd_
        mjd_[dt.detection_time_id] = dt.mjd;

        // Save the observatory coordinates to the array q_obs_
        int idx = dt.detection_time_id*3;
        q_obs_[idx+0] = dt.q_obs_x;
        q_obs_[idx+1] = dt.q_obs_y;
        q_obs_[idx+2] = dt.q_obs_z;
    }
}
