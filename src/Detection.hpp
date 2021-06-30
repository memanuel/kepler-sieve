/** @file Detection.hpp
 *  @brief Class to load a batch of detections in memory and look them up by SkyPatchID.
 *  See DB table Detection and stored procedure GetDetectionDire.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-06-28
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Included libraries
#include <string>
    using std::string;
    using std::to_string;
#include <vector>
    using std::vector;

// Local dependencies
#include "utils.hpp"
    using ks::flush_console;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "SkyPatch.hpp"
    using ks::sky_patch::N_sp;
#include "Timer.hpp"
    using ks::Timer;

// *****************************************************************************
namespace ks {

// *****************************************************************************
/**Data contained in one detection.
 * See DB table KS.Detection and stored procedure KS.GetDetectionObs.*/
struct Detection
{
    /// Integer ID of this detection
    int32_t detection_id;
    /// Integer ID of the SkyPatch where this detection is located in the sky
    int32_t sky_patch_id;
    /// Integer ID of the detection time; FK to KS.HiResTime; the mjd converted to number of minutes.
    int32_t time_id;
    /// Time of this detection
    double mjd;
    /// Direction of this position; x component
    double ux;
    /// Direction of this position; y component
    double uy;
    /// Direction of this position; z component
    double uz;
    /// Position of observatory; x component
    double q_obs_x;
    /// Position of observatory; x component
    double q_obs_y;
    /// Position of observatory; x component
    double q_obs_z;
};

// *****************************************************************************
class DetectionTable
{
public:
    /// Initialize a DetectionTable object with detections in the given range
    DetectionTable(db_conn_type &conn, int d0, int d1, bool progbar);
    /// Initialize a DetectionTable object with all available detections
    DetectionTable(db_conn_type &conn, bool progbar);
    /// Destructor for DetectionTable.
    ~DetectionTable();

    /// Get a detection given its ID
    Detection operator[](int32_t id) const;

    /// Get vector of DetectionIDs matching a given SkyPatchID
    vector<int32_t> get_skypatch(int32_t spid) const;

    /// First detection ID loaded; base for indexing into arrays
    const int d0;
    /// Last detection ID loaded
    const int d1;
    /// The size
    const int size() const;

private:
    /// Vector of detections; dt stands for "Detection Table"
    vector<Detection> dt;
    /// Vector of detection ID vectors keyed by SkyPatchID; dtsp stands for "Detection Table [keyed by] SkyPatchID"
    vector< vector<int32_t> > dtsp;
    /// Helper function to process a batch of rows
    void process_rows(db_conn_type& conn, int i0, int i1);
};

// *****************************************************************************
} // namespace ks
