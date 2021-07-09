/** @file detection_test.cpp
 *  @brief Test harness for classes DetectionTime, DetectionCandidate and Detection.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-01
 * 
 * Example calls:
 * ./detection_test.x --DetectionTime
 * ./save_detection.x --DetectionCandidate
 * ./save_detection.x --Detection
 * ./save_detection.x --test
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <filesystem>
    using std::filesystem::exists;
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "utils.hpp"
    using ks::print_stars;
    using ks::report_test;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "DetectionCandidate.hpp"
    using ks::DetectionCandidate;
    using ks::DetectionCandidateTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;

// *****************************************************************************
// Declare functions
bool test_DetectionTime(db_conn_type& conn);
bool test_DetectionCandidate(db_conn_type& conn);
bool test_Detection(db_conn_type& conn);
bool test_all(db_conn_type& conn);

// *****************************************************************************
int main(int argc, char* argv[])
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run test on DetectionCandidate if requested
    test_all(conn);

    // Close DB connection
    conn->close();

    // Normal program exit
    return 0;
}

// *****************************************************************************
// Test that loading tables from disk matches DB
// *****************************************************************************

// *****************************************************************************
bool test_all(db_conn_type& conn)
{
    // Flag with overall test result
    bool is_ok = true;

    // Test DetectionTime
    print_stars();
    is_ok = is_ok && test_DetectionTime(conn);
    // Test DetectionCandidate
    print_stars();    
    is_ok = is_ok && test_DetectionCandidate(conn);
    // Test Detection
    print_stars();
    is_ok = is_ok && test_Detection(conn);

    // Report overall results
    print_stars();
    report_test("\nOverall DetectionTime, Detection, and DetectionCandidate", is_ok);
    return is_ok;
}

// *****************************************************************************
bool test_DetectionTime(db_conn_type& conn)
{
    // Example row to test
    int i=10;

    // Load DetectionTimeTable from database
    DetectionTimeTable dtt1 = DetectionTimeTable(conn);
    print("\nLoaded DetectionTimeTable from DB with {:d} detection times.\n", dtt1.N());

    // Version loaded from DB
    DetectionTime dt1 = dtt1[i];
    print("Row {:d} of dtt1 (loaded from database):\n", i);
    print("detection_time_id = {:d}\n", dt1.detection_time_id);
    print("time_id           = {:d}\n", dt1.time_id);
    print("mjd               = {:9.4f}\n", dt1.mjd);

    // Load DetectionTimeTable from disk
    DetectionTimeTable dtt2 = DetectionTimeTable();
    print("\nLoaded DetectionTimeTable from disk with {:d} detection times.\n", dtt2.N());

    // Version loaded from disk
    DetectionTime dt2 = dtt2[i];
    print("Row {:d} of dtt2 (loaded from disk):\n", i);
    print("detection_time_id = {:d}\n", dt2.detection_time_id);
    print("time_id           = {:d}\n", dt2.time_id);
    print("mjd               = {:9.4f}\n", dt2.mjd);

    // Test if the two records are identical
    bool is_ok =    (dt1.detection_time_id == dt2.detection_time_id) && 
                    (dt1.time_id           == dt2.time_id) &&
                    (dt1.mjd               == dt2.mjd);

    // Report test results and return status
    report_test("\nLoad DetectionTime", is_ok);
    return is_ok;
}

// *****************************************************************************
bool test_DetectionCandidate(db_conn_type& conn)
{
    // Range of detections to test
    int d0 = 0;
    int d1 = 20;
    bool progbar = true;
    // Example row to test
    int i=10;

    // Load the detection candidate table from the database
    DetectionCandidateTable dt1 = DetectionCandidateTable(d0, d1);
    dt1.load(conn, progbar);

    // Version loaded from DB
    DetectionCandidate dc1 = dt1[i];
    print("\nRow {:d} of dt1 (loaded from database):\n", i);
    print("detection_id = {:d}\n", dc1.detection_id);
    print("sky_patch_id = {:d}\n", dc1.sky_patch_id);
    print("time_id = {:d}\n", dc1.time_id);

    // Load detection table from disk
    DetectionCandidateTable dt2 = DetectionCandidateTable(d0, d1);
    dt2.load();
    // DetectionCandidateTable dt2 = DetectionCandidateTable();

    // Version loaded from disk
    DetectionCandidate dc2 = dt2[i];
    print("\nRow {:d} of dt2 (loaded from disk):\n", i);
    print("detection_id = {:d}\n", dc2.detection_id);
    print("sky_patch_id = {:d}\n", dc2.sky_patch_id);
    print("time_id = {:d}\n", dc2.time_id);

    // Test if the two records are identical
    bool is_ok =    (dc1.detection_id == dc2.detection_id) && 
                    (dc1.sky_patch_id == dc2.sky_patch_id) &&
                    (dc1.time_id      == dc2.time_id);

    // Report test results and return status
    report_test("\nLoad DetectionCandidate", is_ok);
    return is_ok;
}

// *****************************************************************************
bool test_Detection(db_conn_type& conn)
{
    // Range of detections to test
    int d0 = 0;
    int d1 = 20;
    bool progbar = true;
    // Example row to test
    int i=10;

    // Load the detection table from the database
    DetectionTable dt1 = DetectionTable(d0, d1);
    dt1.load(conn, progbar);

    // Version loaded from DB
    Detection det1 = dt1[i];
    print("\nRow {:d} of dt1 (loaded from database):\n", i);
    print("d.detection_id = {:d}\n", det1.detection_id);
    print("d.sky_patch_id = {:d}\n", det1.sky_patch_id);
    print("d.time_id = {:d}\n", det1.time_id);
    print("d.ux      = {:8.6f}\n", det1.ux);
    print("d.uy      = {:8.6f}\n", det1.uy);
    print("d.uz      = {:8.6f}\n", det1.uz);

    // Load detection table from disk
    DetectionTable dt2 = DetectionTable(d0, d1);
    dt2.load();

    // Version loaded from disk
    Detection det2 = dt2[i];
    print("\nRow {:d} of dt2 (loaded from disk):\n", i);
    print("d.detection_id = {:d}\n", det2.detection_id);
    print("d.sky_patch_id = {:d}\n", det2.sky_patch_id);
    print("d.time_id = {:d}\n", det2.time_id);
    print("d.ux      = {:8.6f}\n", det2.ux);
    print("d.uy      = {:8.6f}\n", det2.uy);
    print("d.uz      = {:8.6f}\n", det2.uz);

    // Test if the two records are identical
    bool is_ok =    (det1.detection_id == det2.detection_id) && 
                    (det1.sky_patch_id == det2.sky_patch_id) &&
                    (det1.time_id      == det2.time_id) &&
                    (det1.ux           == det2.ux) &&
                    (det1.uy           == det2.uy) &&
                    (det1.uz           == det2.uz);

    // Report test results and return status
    report_test("\nLoad Detection", is_ok);
    return is_ok;
}
