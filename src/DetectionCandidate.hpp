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
#include <fstream>
    using std::ifstream;
    using std::ofstream;
#include <vector>
    using std::vector;
#include <algorithm>
    using std::min;
#include <fmt/format.h>
    using fmt::print;

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
/**Subset of data contained in one detection needed to search for candidate detections near an asteroid.
 * See DB table KS.Detection and stored procedure KS.GetDetectionObs.*/
struct DetectionCandidate
{
    /// Integer ID of this detection
    int32_t detection_id;
    /// Integer ID of the SkyPatch where this detection is located in the sky
    int32_t sky_patch_id;
    /// Integer ID of the detection time; FK to KS.HiResTime; the mjd converted to number of minutes.
    int32_t time_id;
};

// *****************************************************************************
class DetectionCandidateTable
{
public:
    /// Default constructor builds an empty table, then loads from disk
    DetectionCandidateTable();
    /// Initialize an empty DetectionCandidateTable object with detections in the given range
    DetectionCandidateTable(int d0, int d1);
    /// Initialize a DetectionCandidateTable object with all available detections
    DetectionCandidateTable(db_conn_type& conn, bool progbar);
    /// Destructor for DetectionCandidateTable.
    ~DetectionCandidateTable();

    /// Get a detection candidate given its ID
    const DetectionCandidate operator[](int32_t id) const;

    /// Get vector of DetectionIDs matching a given SkyPatchID
    const vector<int32_t> get_skypatch(int32_t spid) const;

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
    vector<DetectionCandidate> dt;
    /// Vector of detection ID vectors keyed by SkyPatchID; dtsp stands for "Detection Table [keyed by] SkyPatchID"
    vector< vector<int32_t> > dtsp;
    /// Process a batch of rows, writing data from stored procedure output to vector of detections.
    void process_rows(db_conn_type& conn, int i0, int i1);
};

// *****************************************************************************
} // namespace ks
