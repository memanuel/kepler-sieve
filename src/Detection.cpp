/** @file Detection.cpp
 *  @brief Implmentation of Detection class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-28
 */

// *****************************************************************************
// Included files
#include "Detection.hpp"

// *****************************************************************************
// Local names used
using ks::Detection;
using ks::DetectionTable;
using ks::DetectionTime;
using ks::DetectionTimeTable;

// Set batch size to one million
constexpr int batch_size = 1000000;

// Location of file with serialized data
const string file_name = "data/cache/DetectionTable.bin";

// *****************************************************************************
// Class DetectionTable
// *****************************************************************************

// *****************************************************************************
DetectionTable::DetectionTable(): 
    d0(0),
    d1(0),
    dt(vector<Detection>(0)),
    dtsp(vector<vector<int32_t>>(N_sp))
    {}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type& conn, int d0, int d1, bool progbar): 
    d0(d0),
    d1(d1),
    // Initialize dt to a vector with sz entries, one for each possible detection in the interval
    dt(vector<Detection>(d1-d0)),
    // Initialize dtsp to a vector with N_sp entries, one for each SkyPatch (whether occupied or not)
    dtsp(vector<vector<int32_t>>(N_sp))
{
    // Write 0 into DetectionID field so we can later ignore any holes (e.g. DetectionID=0)
    int sz = d1-d0;
    for (int i=0; i<sz; i++)
    {
        dt[i].detection_id=0;
    }

    // Calculate number of batches
    int batch_count = std::max((sz+1)/batch_size, 1);
    // Status update
    if (progbar) 
    {
        print("Processing {:d} detections from {:d} to {:d} in {:d} batches of {:d}...\n",
                sz, d0, d1, batch_count, batch_size);
    }

    // Timer for processing detections from DB
    Timer t;
    t.tick();

    // Iterate over the batches
    for (int i0=d0; i0<d1; i0+=batch_size)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+batch_size, d1);
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
        print("\nLoaded detection candidates table with {:d} entries.\n", dt.size());
        t.tock_msg();
    }
}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, bool progbar): 
    // Delegate to range constructor, using DB to compute d1
    DetectionTable(conn, 0, sp_run_int(conn, "KS.GetMaxDetectionID"), progbar) {}

// *****************************************************************************
/**Helper function: Process a batch of rows, 
 * writing data from stored procedure output to vector of detections. */
void DetectionTable::process_rows(db_conn_type& conn, int i0, int i1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetDetectionsObs";
    vector<string> params = {to_string(i0), to_string(i1)};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset; 10 total fields
        // 3 key fields
        int32_t detection_id = rs->getInt("DetectionID");
        int32_t sky_patch_id = rs->getInt("SkyPatchID");
        int32_t time_id = rs->getInt("TimeID");
        // 1 time field
        double mjd = rs->getDouble("mjd");
        // 3 direction components
        double ux = rs->getDouble("ux");
        double uy = rs->getDouble("uy");
        double uz = rs->getDouble("uz");
        // 3 observatory position components
        double q_obs_x = rs->getDouble("qObs_x");
        double q_obs_y = rs->getDouble("qObs_y");
        double q_obs_z = rs->getDouble("qObs_z");

        // Initialize the Detection at this location in dt
        int idx = detection_id - d0;
        dt[idx] = {
            .detection_id = detection_id,
            .sky_patch_id = sky_patch_id,
            .time_id = time_id,
            .mjd = mjd,
            .ux = ux,
            .uy = uy,
            .uz = uz,
            .q_obs_x = q_obs_x,
            .q_obs_y = q_obs_y,
            .q_obs_z = q_obs_z,
        };

        // Write this DetectionID to the vector keyed by this SkyPatchID
        (dtsp[sky_patch_id]).push_back(detection_id);
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
///Default destructor is OK here
DetectionTable::~DetectionTable() {}

// *****************************************************************************
const int DetectionTable::size() const
{
    return (d1-d0);
}

// *****************************************************************************
Detection DetectionTable::operator[](int32_t id) const
{
    return dt[id];
}

// *****************************************************************************
vector<int32_t> DetectionTable::get_skypatch(int32_t spid) const
{
    return dtsp[spid];
} // function

// *****************************************************************************
void DetectionTable::save()
{
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of rows in binary as long int
    long sz = dt.size();
    fs.write((char*) &sz, sizeof(sz));

    // Write the rows to the file in binary
    for (Detection d : dt)
    {
        fs.write((char*) &d, sizeof(d));
    }

    // Status
    // print("Wrote data to {:s}\n", file_name);

    // Close output filestream
    fs.close();
}

// *****************************************************************************
void DetectionTable::load()
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows
    long sz=-1;
    fs.read( (char*) &sz, sizeof(sz));

    // Reserve space for the detection table
    dt.reserve(sz);

    // Read the rows from the file in binary
    Detection d;
    for (int i=0; i<sz; i++)
    {
        // Read this row into d
        fs.read( (char*) &d, sizeof(d));
        // Save this element to dt
        dt.push_back(d);
        // Write this DetectionID to the vector keyed by this SkyPatchID
        (dtsp[d.sky_patch_id]).push_back(d.detection_id);
    }

    // Update paramater d1
    d1 = dt.size();

    // Close input filestream
    fs.close();
}

// *****************************************************************************
// Class DetectionTimeTable
// *****************************************************************************

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
