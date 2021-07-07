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
    dtv(vector<DetectionTime>(max_id+1)),
    dtm(map<int32_t, vector<int32_t> >{})
    {}

// *****************************************************************************
DetectionTimeTable::DetectionTimeTable(db_conn_type& conn):
    DetectionTimeTable(sp_run_int(conn, "KS.GetDetectionTimeMaxID"))
    {
        // Now load data from the database
        load(conn);
    }

// Default destructor is OK
// *****************************************************************************
DetectionTimeTable::~DetectionTimeTable()
{}

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
const int DetectionTimeTable::size() const
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

    // Reserve space for the detection time vector
    dtv.reserve(sz);

    // // Read the rows from the file in binary
    DetectionTime dt;
    for (int i=0; i<sz; i++)
    {
        // Read this row into d
        fs.read( (char*) &dt, sizeof(dt));
        // Write this DetectionTime into the vector dtv
        dtv.push_back(dt);
        // Add this detection time to the map keyed by TimeID
        (dtm[dt.time_id]).push_back(dt.detection_time_id);
    }

    // Close input filestream
    fs.close();
}
