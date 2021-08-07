/** @file DetectionTime.hpp
 *  @brief Class to load a batch of detection times in memory and look them up by TimeID.
 *  See DB table Detection and stored procedure GetDetectionDire.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-07
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
#include <map>
    using std::map;

// Local dependencies
#include "utils.hpp"
    using ks::flush_console;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "Timer.hpp"
    using ks::Timer;

// *****************************************************************************
namespace ks {

// *****************************************************************************
/**Data contained in one detection time.
 * See DB table KS.DetectionTime and stored procedure KS.GetDetectionTimes.*/
struct DetectionTime
{
    /// Integer ID of this detection time
    int32_t detection_time_id;
    /// Integer ID of the detection time; FK to KS.HiResTime; the mjd converted to number of minutes.
    int32_t time_id;
    /// Data source of this detection
    int8_t data_source_id;
    /// Observatory where this detection was made
    int8_t observatory_id;
    /// Time of this detection
    double mjd;
    /// Position of observatory; x component
    double q_obs_x;
    /// Position of observatory; y component
    double q_obs_y;
    /// Position of observatory; z component
    double q_obs_z;
    /// Position of Sun; x component
    double q_sun_x;
    /// Position of Sun; y component
    double q_sun_y;
    /// Position of Sun; z component
    double q_sun_z;
};

// *****************************************************************************
class DetectionTimeTable
{

// Constructor and destructor
public:
    /// Default constructor builds a table and populates it from disk
    DetectionTimeTable();
    /// Constructor builds an empty table with the requested size
    DetectionTimeTable(int N);
    /// This constructor builds an empty table, then loads it using the given database connection
    DetectionTimeTable(db_conn_type& conn);
    /// Destructor for DetectionTimeTable.
    ~DetectionTimeTable() ;

// Data elements
public:
    // Size of this table
    const int N;
private:
    /// Vector of detection times; dtv stands for "DetectionTime Vector"
    vector<DetectionTime> dtv;
    /// Map of detection ID vectors keyed by TimeID; dtm stands for "DetectionTime map"
    /// Given a detection_time_id (the key), returns a vector of all the detection_id's made at that time (the value)
    map<int32_t, vector<int32_t> > dtm;
    /// Array of mjds when detections taken; size N+1
    double* const mjd_;
    /// Array of observatory positions in HELIOCENTRIC frame; size 3(N+1)
    double* const q_obs_;
public:
    /// Read-only copy of mjd
    const double* const mjd;
    /// Read-only copy of q_obs
    const double* const q_obs;

// Member functions
public:
    /// Get a detection time given its ID
    const DetectionTime operator[](int32_t id) const;
    /// Get vector of DetectionIDs matching a given TimeID
    const vector<int32_t> get_time(int32_t time_id) const;
    /// Vector of all detection time objects
    const vector<DetectionTime> detection_times() const;
    /// Get date of first detection; this is in array slot 1 because indexing matches detection_id
    const double mjd_first() const {return mjd_[1];}
    /// Get date of last detection; this is in array slot N because indexing matches detection_id
    const double mjd_last() const {return mjd_[N];}
   
    /// Load all available detections using the DB connection
    void load(db_conn_type& conn);
    /// Save this object to disk
    void save();
    /// Load this object from disk
    void load();

    // Write out the arrays mjd_ and q_obs from the vector of detection time structs
    void write_arrays();
    /// Number of rows in data file
    const int file_length() const;
};

// *****************************************************************************
} // namespace ks
