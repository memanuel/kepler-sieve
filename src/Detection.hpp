/** @file Detection.hpp
 *  @brief Class to load a batch of detections in memory and look them up by SkyPatchID.
 *  See DB table Detection and stored procedures KS.GetDetections and KS.GetDetectionsObs.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-06-28
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <string>
    using std::string;
    using std::to_string;
#include <fstream>
    using std::ifstream;
    using std::ofstream;
#include <vector>
    using std::vector;
#include <algorithm>
    using std::min;

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
    /// Integer ID of the detection time; FK to KS.DetectionTime; unique ID for times
    int32_t detection_time_id;
    /// Time of this detection
    double mjd;
    /// Direction of this detection; x component
    double ux;
    /// Direction of this detection; y component
    double uy;
    /// Direction of this detection; z component
    double uz;
    /// Magnitude of this detection
    double mag;
};

// *****************************************************************************
class DetectionTable
{
public:
    /// Default constructor builds an empty table, then loads from disk
    DetectionTable();
    /// Initialize an empty DetectionTable object with detections in the given range
    DetectionTable(int d0, int d1);
    /// Initialize a DetectionTable object with all available detections
    DetectionTable(db_conn_type& conn, bool progbar);
    /// Destructor for DetectionTable.
    ~DetectionTable();

    /// Get a detection given its ID
    Detection operator[](int32_t id) const;
    /// Get vector of DetectionIDs matching a given SkyPatchID
    vector<int32_t> get_skypatch(int32_t spid) const;

    /// First detection ID loaded; base for indexing into arrays
    const int d0;
    /// Last detection ID loaded
    int d1;
    /// The size
    const int size() const;

    /// Load this object with from database
    void load(db_conn_type& conn, bool progbar);
    /// Save this object to disk
    void save();
    /// Load this object from disk
    void load();

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
