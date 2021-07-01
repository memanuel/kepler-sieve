/** @file Detection.cpp
 *  @brief Implmentation of Detection class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-28
 */

// *****************************************************************************
// Included files
#include "DetectionCandidate.hpp"

// *****************************************************************************
// Local names used
using ks::DetectionCandidate;
using ks::DetectionCandidateTable;

// Set batch size to one million
constexpr int batch_size = 1000000;

// *****************************************************************************
void DetectionCandidateTable::process_rows(db_conn_type& conn, int i0, int i1)
{
    // Run the stored procedure to get detections; does not need the observatory positions
    string sp_name = "KS.GetDetectionCandidates";
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

        // Initialize the DetectionCandidate at this location in dt
        int idx = detection_id - d0;
        dt[idx] = {
            .detection_id = detection_id,
            .sky_patch_id = sky_patch_id,
            .time_id = time_id,
        };

        // Write this DetectionID to the vector keyed by this SkyPatchID
        (dtsp[sky_patch_id]).push_back(detection_id);
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    // delete rs;
}

// *****************************************************************************
DetectionCandidateTable::DetectionCandidateTable(): 
    d0(0),
    d1(0),
    dt(vector<DetectionCandidate>(0)),
    dtsp(vector<vector<int32_t>>(0))
    {}

// *****************************************************************************
DetectionCandidateTable::DetectionCandidateTable(db_conn_type &conn, int d0, int d1, bool progbar): 
    d0(d0),
    d1(d1),
    // Initialize dt to a vector with sz entries, one for each possible detection in the interval
    dt(vector<DetectionCandidate>(d1-d0)),
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
        print("Processing {:d} detection candidates from {:d} to {:d} in {:d} batches of {:d}...\n",
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
DetectionCandidateTable::DetectionCandidateTable(db_conn_type &conn, bool progbar): 
    // Delegate to range constructor, using DB to compute d1
    DetectionCandidateTable(conn, 0, sp_run_int(conn, "KS.GetMaxDetectionID"), progbar) {}

// *****************************************************************************
///Default destructor is OK here
DetectionCandidateTable::~DetectionCandidateTable() {}

// *****************************************************************************
const int DetectionCandidateTable::size() const
{
    return (d1-d0);
}

// *****************************************************************************
DetectionCandidate DetectionCandidateTable::operator[](int32_t id) const
{
    // Need to subtract d0 from id to get the index into the vector dt
    return dt[id-d0];
}

// *****************************************************************************
vector<int32_t> DetectionCandidateTable::get_skypatch(int32_t spid) const
{
    return dtsp[spid];
}

// *****************************************************************************
void serialize()
{
    // Name of file where this is saved
    string file_path = __FILE__;
    print("file_path={:s}\n", file_path);
}