/*****************************************************************************
 * Michael S. Emanuel
 * 2021-06-28
 * ****************************************************************************/

// *****************************************************************************
// Included files
#include "Detection.h"

// *****************************************************************************
// Local names used
using ks::Detection;
using ks::DetectionTable;

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, int d0, int d1): 
    d0(d0),
    dt(vector<Detection>(0))
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetDetectionsObs";
    vector<string> params = {to_string(d0), to_string(d1)};
    ResultSet *rs = sp_run(conn, sp_name, params);

    // Allocate a vector with of detections with the appropriate size
    int sz = d1-d0;
    dt.resize(sz);
    // Write -1 into DetectionID field so we can later ignore any holes (e.g. DetectionID=0)
    for (int i=0; i<sz; i++)
    {
        dt[i].detection_id=-1;
    }

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
    }
};

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
