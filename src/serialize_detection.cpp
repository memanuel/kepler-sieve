/** @file serialize_detection.cpp
 *  @brief Batch program to save output of DB queries of detections to disk.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-01
 * 
 * Example call:
 * ./serialize_detection.x
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "utils.hpp"
    using ks::report_test;

#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;

#include "DetectionCandidate.hpp"
    using ks::DetectionCandidate;
    using ks::DetectionCandidateTable;

// *****************************************************************************
// Declare functions
void serialize_Detection(db_conn_type& conn);
void serialize_DetectionCandidate(db_conn_type& conn);
void test(db_conn_type& conn);

// *****************************************************************************
int main()
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run test on DetectionCandidate
    test(conn);

    // Close DB connection
    conn->close();

    // Normal program exit
    return 0;
}

// *****************************************************************************
void serialize_Detection(db_conn_type conn)
{
    // Load the detection candidate table
    bool progbar = true;
    DetectionTable dt = DetectionTable(conn, progbar);

    // Save detection table to disk
    dt.save();
}

// *****************************************************************************
void serialize_DetectionCandidate(db_conn_type conn)
{
    // Load the detection candidate table
    bool progbar = true;
    DetectionCandidateTable dt = DetectionCandidateTable(conn, progbar);

    // Save detection table to disk
    dt.save();
}

// *****************************************************************************
void test(db_conn_type& conn)
{
    // Load the detection candidate table
    int d0 = 0;
    int d1 = 1000000;
    bool progbar = true;
    DetectionCandidateTable dt1 = DetectionCandidateTable(conn, d0, d1, progbar);

    // Example row
    int i=10;

    // Version loaded from DB
    DetectionCandidate dc1 = dt1[i];
    print("\nRow {:d} of dt1 (loaded from database):\n", i);
    print("dc.detection_id = {:d}\n", dc1.detection_id);
    print("dc.sky_patch_id = {:d}\n", dc1.sky_patch_id);
    print("dc.time_id = {:d}\n", dc1.time_id);

    // Load detection table from disk
    DetectionCandidateTable dt2;
    dt2.load();

    // Version loaded from DB
    DetectionCandidate dc2 = dt2[i];
    print("\nRow {:d} of dt2 (loaded from disk:\n", i);
    print("dc.detection_id = {:d}\n", dc2.detection_id);
    print("dc.sky_patch_id = {:d}\n", dc2.sky_patch_id);
    print("dc.time_id = {:d}\n", dc2.time_id);

    // Test if the two records are identical
    bool is_ok =    (dc1.detection_id == dc2.detection_id) && 
                    (dc1.sky_patch_id == dc2.sky_patch_id) &&
                    (dc1.time_id      == dc2.time_id);
    report_test("Load DetectionCandidate", is_ok);
}
