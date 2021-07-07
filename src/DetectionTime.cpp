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

// Location of file with serialized data
const string file_name = "data/cache/DetectionTimeTable.bin";

// *****************************************************************************
DetectionTimeTable::DetectionTimeTable():
    // Delegate to size constructor with a placeholder size of zero
    DetectionTimeTable(0)
{
    // Load data from disk; this resizes dtv to match size on disk
    load();
}

// *****************************************************************************
DetectionTimeTable::DetectionTimeTable(int max_id):
    // dtv is a vector of DetectionTimes; size is max_id+1 because we index by detection_time_id
    dtv(vector<DetectionTime>(max_id+1)),
    // dtm is an empty map keyed by TimeID; it will be populated later
    dtm(map<int32_t, vector<int32_t> >{}),
    // mjds is an array of doubles of size N (size of dtv)
    mjds(new double[max_id+1]),
    // q_obs is an array of doubles with 3N entries for the HELIOCENTRIC observatory position
    q_obs(new double[max_id+1])
    {}

// *****************************************************************************
DetectionTimeTable::DetectionTimeTable(db_conn_type& conn):
    DetectionTimeTable(sp_run_int(conn, "KS.GetDetectionTimeMaxID"))
    {
        // Now load data from the database
        load(conn);
    }

// *****************************************************************************
// Delete manually allocated arrays
DetectionTimeTable::~DetectionTimeTable()
{
    delete [] mjds;
    delete [] q_obs;
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
        dtv.push_back(dt);
        // Add this detection time to the map keyed by TimeID
        (dtm[time_id]).push_back(detection_time_id);
        // Save mjd to the array
        mjds[detection_time_id] = mjd;
    }
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
DetectionTime DetectionTimeTable::operator[](int32_t id) const
{
    return dtv[id];
}

// *****************************************************************************
vector<int32_t> DetectionTimeTable::get_time(int32_t time_id)
{
    return dtm[time_id];
}

// *****************************************************************************
const int DetectionTimeTable::N() const
{
    return dtv.size();
}

// *****************************************************************************
const vector<DetectionTime> DetectionTimeTable::detection_times() const
{
    return dtv;
}

// *****************************************************************************
void DetectionTimeTable::save()
{
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of rows in binary as long int
    long sz = dtv.size();
    fs.write((char*) &sz, sizeof(sz));

    // Write the rows to the file in binary
    for (DetectionTime dt : dtv)
    {
        fs.write((char*) &dt, sizeof(dt));
    }

    // Close output filestream
    fs.close();
}

// *****************************************************************************
void DetectionTimeTable::load()
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows in the file
    long sz=-1;
    fs.read( (char*) &sz, sizeof(sz));
    // Status
    // print("Opened file {:s} with {:d} rows of DetectionTime data.\n", file_name, sz);

    // Resize the detection time vector
    dtv.resize(sz+1);
    // Need to free memory from the placeholder arrays of mjds and q_obs to avoid a memory leak
    delete [] mjds;
    delete [] q_obs;
    // Allocate new storage for arrays mjds and q_obs now that we know their size
    mjds = new double[sz+1];
    q_obs = new double[3*(sz+1)];

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
        // Save the mjd to the array mjds
        mjds[detection_time_id] = dt.mjd;
    }

    // Close input filestream
    fs.close();

    // Calculate observatory position in heliocentric frame
    calc_q_obs();
}

// *****************************************************************************
void DetectionTimeTable::calc_q_obs()
{
    // Iterate through the detection times in the vector 
    for (DetectionTime dt : dtv)
    {
        // The observatory position in the heliocentric frame
        double qx = dt.q_obs_x - dt.q_sun_x;
        double qy = dt.q_obs_y - dt.q_sun_y;
        double qz = dt.q_obs_z - dt.q_sun_z;
        // Save these entries to the q_obs array
        int idx = dt.detection_time_id*3;
        q_obs[idx+0] = qx;
        q_obs[idx+1] = qy;
        q_obs[idx+2] = qz;
    }
}
