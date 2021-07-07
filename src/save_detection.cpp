/** @file serialize_detection.cpp
 *  @brief Batch program to save output of DB queries of detections to disk.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-01
 * 
 * Example call:
 * ./save_detection.x
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <filesystem>
    using std::filesystem::exists;
#include <fmt/format.h>
    using fmt::print;

#include <boost/program_options.hpp>
    namespace po = boost::program_options;

// Local dependencies
#include "utils.hpp"
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
void save_DetectionTime(db_conn_type& conn);
void save_DetectionCandidate(db_conn_type& conn);
void save_Detection(db_conn_type& conn);

void test_DetectionCandidate(db_conn_type& conn);
void test_Detection(db_conn_type& conn);
void test_all(db_conn_type& conn);

// *****************************************************************************
int main(int argc, char* argv[])
{
    // *****************************************************************************
    // Process commandline arguments.
    // *****************************************************************************

    // Flags from commandline arguments
    bool run_test = false;
    bool run_save_DetectionTime = false;
    bool run_save_Detection = false;
    bool run_save_DetectionCandidate = false;

    // Set up parser for named commandline arguments ("options")
    po::options_description desc("Find detections near asteroids");
    desc.add_options()
        ("help,h", "Produce help message")
        ("test", po::bool_switch(&run_test), "Run test")
        ("DetectionTime", po::bool_switch(&run_save_DetectionTime), 
            "Save contents of stored procedure KS.GetDetectionTimes")
        ("Detection", po::bool_switch(&run_save_Detection), 
            "Save contents of stored procedure KS.GetDetections on full range")
        ("DetectionCandidate", po::bool_switch(&run_save_DetectionCandidate), 
            "Save contents of stored procedure KS.GetCandidateDetections on full range")
    ;
    po::variables_map vm;

    // Parse commandline arguments including positional arguments
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // If program called with no arguments, assume user wanted to run tests
    bool run_one_file = (run_save_DetectionTime || run_save_Detection || run_save_DetectionCandidate);
    run_test = run_test || (!run_one_file);

    // *****************************************************************************
    // Program body
    // *****************************************************************************

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run test on DetectionCandidate if requested
    if (run_test) {test_all(conn);}

    // Save DetectionTime if requested
    if (run_save_DetectionTime) {save_DetectionTime(conn);}

    // Save Detection if requested
    if (run_save_Detection) {save_Detection(conn);}

    // Save DetectionCandidate if requested
    if (run_save_DetectionCandidate) {save_DetectionCandidate(conn);}

    // Close DB connection
    conn->close();

    // Normal program exit
    return 0;
}

// *****************************************************************************
// Save a file to disk with cached data from database.
// *****************************************************************************

// *****************************************************************************
void save_DetectionTime(db_conn_type& conn)
{
    // Load the detection time table from database
    DetectionTimeTable dtt = DetectionTimeTable(conn);
    // Save detection time table to disk
    dtt.save();
    print("Saved DetectionTime to disk.");
}

// *****************************************************************************
void save_DetectionCandidate(db_conn_type& conn)
{
    // Load the detection candidate table from database
    bool progbar = true;
    DetectionCandidateTable dt = DetectionCandidateTable(conn, progbar);
    // Save detection candidate table to disk
    dt.save();
}

// *****************************************************************************
void save_Detection(db_conn_type& conn)
{
    // Load the detection candidate table from database
    bool progbar = true;
    DetectionTable dt = DetectionTable(conn, progbar);
    // Save detection table to disk
    dt.save();
}

// *****************************************************************************
// Test that loading tables from disk matches DB
// *****************************************************************************

// *****************************************************************************
void test_all(db_conn_type& conn)
{
    test_DetectionCandidate(conn);
    test_Detection(conn);
}

// *****************************************************************************
void test_DetectionCandidate(db_conn_type& conn)
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
    print("dc.detection_id = {:d}\n", dc1.detection_id);
    print("dc.sky_patch_id = {:d}\n", dc1.sky_patch_id);
    print("dc.time_id = {:d}\n", dc1.time_id);

    // Load detection table from disk
    DetectionCandidateTable dt2 = DetectionCandidateTable(d0, d1);
    dt2.load();
    // DetectionCandidateTable dt2 = DetectionCandidateTable();

    // Version loaded from disk
    DetectionCandidate dc2 = dt2[i];
    print("\nRow {:d} of dt2 (loaded from disk):\n", i);
    print("dc.detection_id = {:d}\n", dc2.detection_id);
    print("dc.sky_patch_id = {:d}\n", dc2.sky_patch_id);
    print("dc.time_id = {:d}\n", dc2.time_id);

    // Test if the two records are identical
    bool is_ok =    (dc1.detection_id == dc2.detection_id) && 
                    (dc1.sky_patch_id == dc2.sky_patch_id) &&
                    (dc1.time_id      == dc2.time_id);

    // Report test results
    report_test("\nLoad DetectionCandidate", is_ok);
}

// *****************************************************************************
void test_Detection(db_conn_type& conn)
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

    // Report test results
    report_test("\nLoad Detection", is_ok);
}
