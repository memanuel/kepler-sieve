/** @file AsteroidSkyPatch.cpp
 *  @brief Implmentation of AsteroidSkyPatch class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

// *****************************************************************************
// Included files
#include "AsteroidSkyPatch.hpp"

// *****************************************************************************
// Local names used
using ks::AsteroidSkyPatch;
using ks::AsteroidSkyPatchTable;

// Set batch size
constexpr int bs = 10000;

// *****************************************************************************
/** Helper function: Process a batch of rows, 
 * writing data from stored preocedure output to vector of detections. */
void process_rows(db_conn_type &conn, vector<AsteroidSkyPatch> &asp, int n0, int n1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetAsteroidSkyPatch";
    vector<string> params = {to_string(n0), to_string(n1)};
    ResultSet *rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset; 10 total fields
        int32_t asteroid_id = rs->getInt("AsteroidID");
        int32_t segment = rs->getInt("Segment");
        int32_t sky_patch_id = rs->getInt("SkyPatchID");
        int32_t time_id_0 = rs->getInt("TimeID_0");
        int32_t time_id_1 = rs->getInt("TimeID_1");

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
    }   // while rs
    // Close the resultset
    rs->close();
}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, int d0, int d1, bool progbar): 
    d0(d0),
    d1(d1),
    dt(vector<Detection>(d1-d0))
{
    // Write -1 into DetectionID field so we can later ignore any holes (e.g. DetectionID=0)
    int sz = d1-d0;
    for (int i=0; i<sz; i++)
    {
        dt[i].detection_id=-1;
    }

    // Set batch size
    int batch_count = sz / bs;

    // Status update
    if (progbar) 
    {
        print("Processing {:d} detections from {:d} to {:d} in {:d} batches of {:d}...\n",
                sz, d0, d1, batch_count, bs);
    }

    // Timer for processing detections from DB
	Timer t;
    t.tick();

    // Iterate over the batches
    for (int i0=0; i0<d1; i0+=bs)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+bs, d1);
        // Process SQL data in this batch
        process_rows(conn, dt, i0, i1);
        // Progress bar
        if (progbar) {print("."); }
    }
    if (progbar) 
    {
        print("\nLoaded detection table.\n");
        t.tock_msg();
    }
}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, bool progbar): 
    // Delegate to range constructor, using DB to compute d1
    DetectionTable(conn, 0, sp_run_int(conn, "KS.GetMaxDetectionID"), progbar) {}

// *****************************************************************************
///Default destructor is OK here
DetectionTable::~DetectionTable() {}

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
