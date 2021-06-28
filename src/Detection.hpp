/*****************************************************************************
 * Class to load a batch of detections in memory and look them up by SkyPatchID.
 * See DB table Detection and stored procedure GetDetectionDire.
 * 
 * Michael S. Emanuel
 * 2021-06-28
 * ****************************************************************************/
#pragma once

// *****************************************************************************
// Included libraries
#include <string>

// Local dependencies
#include "db_utils.hpp"

// *****************************************************************************
// Standard library class names used
using std::string;
using std::to_string;
using std::vector;

// *****************************************************************************
// Local names used
using ks::db_conn_type;
using ks::get_db_conn;
using ks::sp_run;

// *****************************************************************************
namespace ks {

// *****************************************************************************
namespace detection{
} // namespace ks::detection

// *****************************************************************************
/**Data contained in one detection.
 * See DB table KS.Detection and stored procedure KS.GetDetectionObs.*/
struct Detection
{
    //*Integer ID of this detection
    int32_t detection_id;
    //*Integer ID of the SkyPatch where this detection is located in the sky
    int32_t sky_patch_id;
    //*Integer ID of the detection time; FK to KS.HiResTime; the mjd converted to number of minutes.
    int32_t time_id;
    //*Time of this detection
    double mjd;
    //*Direction of this position; x component
    double ux;
    //*Direction of this position; y component
    double uy;
    //*Direction of this position; z component
    double uz;
    //*Position of observatory; x component
    double q_obs_x;
    //*Position of observatory; x component
    double q_obs_y;
    //*Position of observatory; x component
    double q_obs_z;
};

// *****************************************************************************
class DetectionTable
{
public:
    //*Initialize a DetectionTable object with all available detections
    DetectionTable(db_conn_type &conn, bool progbar);
    //*Initialize a DetectionTable object with detections in the given range
    DetectionTable(db_conn_type &conn, int d0, int d1);
    //*Destructor for DetectionTable.
    ~DetectionTable();

    //*Get a detection given its ID
    Detection operator[](int32_t id) const;

    //*Get vector of DetectionIDs matching a given SkyPatchID
    vector<int32_t> get_skypatch(int32_t spid) const;

private:
    // Data
    //*First detection ID loaded; base for indexing into arrays
    const int d0;
    //*Last detection ID loaded
    int d1;
    //*Vector of detections; dt stands for "Detection Table"
    vector<Detection> dt;
    //*Vector of detection ID vectors keyed by SkyPatchID; dtsp stands for "Detection Table [keyed by] SkyPatchID"
    vector< vector<int32_t> > dtsp;
};

// *****************************************************************************
} // namespace ks
