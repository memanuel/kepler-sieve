/*****************************************************************************
 * Michael S. Emanuel
 * 2021-06-28
 * ****************************************************************************/

// *****************************************************************************
// Included files
#include "Detection.hpp"

// *****************************************************************************
// Local names used
using ks::Detection;
using ks::DetectionTable;

// *****************************************************************************
/** Helper function: Process a batch of rows, 
 * writing data from stored preocedure output to vector of detections. */
void process_rows(db_conn_type &conn, vector<Detection> &dt, int d0, int d1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetDetectionsObs";
    vector<string> params = {to_string(d0), to_string(d1)};
    ResultSet *rs = sp_run(conn, sp_name, params);

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
    }   // while rs
    // Close the resultset
    rs->close();
}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, int d0, int d1): 
    d0(d0),
    d1(d1),
    dt(vector<Detection>(0))
{
    // Allocate a vector of detections with the appropriate size
    int sz = d1-d0;
    dt.resize(sz);

    // Write -1 into DetectionID field so we can later ignore any holes (e.g. DetectionID=0)
    for (int i=0; i<sz; i++)
    {
        dt[i].detection_id=-1;
    }

    // Write vector of detections using SQL data
    process_rows(conn, dt, d0, d1);
}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, bool progbar): 
    // Provisionally populate d0 and d1
    d0(0),
    d1(1),
    dt(vector<Detection>(0))
{

    // Calculate d1 from the database
    // TODO: replace hard coded values with real DB Calls!
    d1 = 1000000;
    // d1 = 159000000;

    // Allocate a vector of detections with the appropriate size
    int sz = d1-d0;
    dt.resize(sz);

    // Write -1 into DetectionID field so we can later ignore any holes (e.g. DetectionID=0)
    for (int i=0; i<sz; i++)
    {
        dt[i].detection_id=-1;
    }

    // Set batch size
    int b = 10000;
    int n = (d1-d0);
    int batch_count = n / b;
    // Status update
    if (progbar) 
    {
        print("Processing detections from {:d} to {:d} in {:d} batches of {:d}...\n", d0, d1, batch_count, b);
    }

    // Iterate over the batches
    for (int i0=0; i0<d1; i0+=b)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+b, d1);
        // Process SQL data in this batch
        process_rows(conn, dt, i0, i1);
        // Progress bar
        if (progbar) {print("."); }
    }
    if (progbar) {print("Loaded detections.\n");}

}
// *****************************************************************************
//*Default destructor is fine
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
