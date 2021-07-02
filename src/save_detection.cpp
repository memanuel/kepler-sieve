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

#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;

#include "DetectionCandidate.hpp"
    using ks::DetectionCandidate;
    using ks::DetectionCandidateTable;

// *****************************************************************************
// Declare functions
void save_Detection(db_conn_type& conn);
void save_DetectionCandidate(db_conn_type& conn);
void test(db_conn_type& conn);

// *****************************************************************************
int main(int argc, char* argv[])
{

    // *****************************************************************************
    // Process commandline arguments.
    // *****************************************************************************

    // Flags from commandline arguments
    bool run_test = false;
    bool run_save_Detection = false;
    bool run_save_DetectionCandidate = false;

    // Set up parser for named commandline arguments ("options")
    po::options_description desc("Find detections near asteroids");
    desc.add_options()
        ("help,h", "Produce help message")
        ("test", po::bool_switch(&run_test), "Run test")
        ("Detection", po::bool_switch(&run_save_Detection), 
            "Serialize contents of stored procedure KS.GetDetections on full range")
        ("DetectionCandidate", po::bool_switch(&run_save_DetectionCandidate), 
            "Serialize contents of stored procedure KS.GetCandidateDetections on full range")
    ;
    po::variables_map vm;

    // Parse commandline arguments including positional arguments
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // If program called with no arguments, assume user wanted to run tests
    // if !(run_save_Detection | run_save_DetectionCandidate) {run_test = true;}
    run_test =  run_test | (!(run_save_Detection | run_save_DetectionCandidate));

    // *****************************************************************************
    // Program body
    // *****************************************************************************

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run test on DetectionCandidate if requested
    if (run_test) {test(conn);}

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
void save_Detection(db_conn_type& conn)
{
    // Load the detection candidate table
    bool progbar = true;
    DetectionTable dt = DetectionTable(conn, progbar);

    // Save detection table to disk
    dt.save();
}

// *****************************************************************************
void save_DetectionCandidate(db_conn_type& conn)
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

    // Does the file with saved data exist?
    const string file_name = "data/cache/DetectionCandidateTable.bin";
    bool file_exists = std::filesystem::exists(file_name);
    // If the file does not exist, save it
    if (!file_exists)
    {
        dt1.save();
    }

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

    // Report test results
    report_test("\nLoad DetectionCandidate", is_ok);
}
