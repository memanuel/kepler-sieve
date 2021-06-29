/*****************************************************************************
 * Batch program to calculate candidate pairs (AsteroidID, DetectionID)
 * where an asteroid is in a neighboring skypatch of the detection.
 * Results further refined downstream in detection_near_ast.py.
 * 
 * Commandline arguments:
 * jn: the job number
 * sz: the job size; defaults to 100000
 * This will process asteroids with asteroid numbers in [n0, n1), where
 * n0 = jn*sz
 * n1 = n0+sz
 * 
 * Example calls:
 * ./detection_near_ast.x --jn 0 --sz 1000
 * ./detection_near_ast.x 0 1000
 * both of these will process asteroids in [0, 1000)
 * 
 * Michael S. Emanuel
 * 2021-06-25
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;
#include <boost/program_options.hpp>
    namespace po = boost::program_options;

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::sql_stmt_type;
    using ks::sql_prepared_stmt_type;
    using ks::get_db_conn;
    using ks::sp_run;
    using ks::sp_run_int;

#include "utils.hpp"
    using ks::print_stars;
    using ks::print_newline;
    using ks::time2hms;

#include "SkyPatchNeighbor.hpp"
    using ks::SkyPatch;
    using ks::SkyPatchNeighbor;
    using ks::sky_patch::N_sp;

#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;

#include "AsteroidSkyPatch.hpp"
    using ks::AsteroidSkyPatch;
    using ks::AsteroidSkyPatchTable;

// *****************************************************************************
// Data type to describe one detection / asteroid candidate match
struct DetectionNearAsteroidCandidate
{
    int32_t detection_id;
    int32_t asteroid_id;
};

// Batch size used for processing asteroids
constexpr int batch_size = 10;

// *****************************************************************************
// Declare functions
void search_asteroid_detection(
    DetectionTable& dt, AsteroidSkyPatchTable& aspt, 
    SkyPatchNeighbor& spn, vector<DetectionNearAsteroidCandidate>& cv);
void write_candidates_db(db_conn_type& conn, const vector<DetectionNearAsteroidCandidate>& cv, int k0, int k1);
void test_detection_table(DetectionTable& dt, int detection_id);
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id);
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt);
void test_all();
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);

// *****************************************************************************
/** Perform a search over detections in the detection table dt and
 *  AsteroidSkyPatch entries in the AsteroidSkyPatchTable aspt.
 * \param[in] dt - table of detection data to search over
 * \param[in] aspt - table of asteroid sky patch segments to search over
 * \param[in] spn- table with the 9 neighbors of each sky patch
 * \param[in] cv - vector of candidate detections near an asteroid
 * */
void search_asteroid_detection(
    DetectionTable& dt, AsteroidSkyPatchTable& aspt, 
    SkyPatchNeighbor& spn, vector<DetectionNearAsteroidCandidate>& cv)
{
    // Counter for the number of matches found and pairs examined
    // int pairs = 0;
    // Reusable pointer to the 9 neighbors of a SkyPatch
    int32_t* ns;

    // Iterate through all the AsteroidSkyPatch entries
    for (int i=0; i<aspt.size(); i++)
    {
        // The AsteroidSkyPatch
        AsteroidSkyPatch asp = aspt[i];
        // The AsteroidID of this segment
        int32_t asteroid_id = asp.asteroid_id;
        // Get the neighbors of the SkyPatch occupied by the asteroid in this segment
        int32_t sky_patch_ast = asp.sky_patch_id;
        ns = spn[sky_patch_ast];
        // Get the start and end time_id
        int32_t time_id_0 = asp.time_id_0;
        int32_t time_id_1 = asp.time_id_1;

        // Iterate through all 9 neighbors
        for (int j=0; j<9; j++)
        {
            // The neighboring SkyPatch
            int32_t spid = ns[j];
            // The vector of Detections in the neighboring SkyPatch
            vector<int32_t> ds = dt.get_skypatch(spid);
            for (int32_t detection_id : ds)
            {
                // The time of this detection
                int32_t time_id = dt[detection_id].time_id;
                // Check each detection to see if its time falls in the relevant window
                if ((time_id_0 <= time_id) && (time_id <= time_id_1))
                {
                    // Assemble the candidate
                    DetectionNearAsteroidCandidate c = {
                        .detection_id = detection_id,
                        .asteroid_id = asteroid_id,
                    };
                    // Add the candidate to the table
                    cv.push_back(c);
                } // if / time_id
                // Increment the pairs counter
                // pairs++;
            } // for / detection_id            
        } // for / j (neighbors of asteroid sky patch i)        
    } // for / i (starting asteroid  sky patch)
}

// *****************************************************************************
/// Write one batch of detection candidates to the databaes.
void write_candidates_db(
    db_conn_type& conn, const vector<DetectionNearAsteroidCandidate>& cv,
    int k0, int k1)
{
     // SQL to insert one record into DetectionNearAsteroid
    string sql = R"(
    REPLACE INTO KS.DetectionNearAsteroidCandidate 
    (DetectionID, AsteroidID) 
    VALUES
    (?, ?);
    )";
    // print("\nSQL string:\n");
    // print(sql);
    // Create a SQL PreparedStatement for the insert
    sql_prepared_stmt_type stmt(conn->prepareStatement(sql));

    // Iterate through the selected block of candidates in this batch
    for (int k=k0; k<k1; k++)
    {
        // The candidate
        DetectionNearAsteroidCandidate c = cv[k];
        // Bind two integer parameters to the prepared statement    
        stmt->setInt(1, c.detection_id);
        stmt->setInt(2, c.asteroid_id);
        // Execute the SQL statement
        stmt->executeQuery();
    } // for / k 
}

// *****************************************************************************
/// Parse the commandline arguments; report them and wrap into a map.
std::map<string, int> parse_args(int argc, char* argv[])
{
    // Set up parser for commandline arguments
    po::options_description desc(
        "./detection_near_ast.x jn sz.\n"
        "Process asteroids from n0=jn*sz to n1=(jn+1)*sz.\n");
    desc.add_options()
        ("jn", po::value<int>()->default_value(0), "Job number")
        ("sz", po::value<int>()->default_value(1000), "Batch size for jobs");
    po::variables_map vm;

    po::positional_options_description p;
    p.add("jn", 1);
    p.add("sz", 2);

    // Unpack commandline arguments
    // po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    notify(vm);
    po::notify(vm);

    int jn = vm["jn"].as<int>();
    int sz = vm["sz"].as<int>();
    int n0 = jn*sz;
    int n1 = n0+sz;

    // Report commandline arguments
    print("Commandline arguments\n");
    print("Job Number jn       : {:d}\n", jn);
    print("Batch Size sz       : {:d}\n", sz);
    print("Start AsteroidID n0 : {:d}\n", n0);
    print("End AsteroidID   n1 : {:d}\n", n1);
    print_newline();

    // Wrap into standard map
    std::map<string, int> args = {{"n0", n0},{"n1", n1}};
    return args;
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    // Parse command line arguments
    std::map<string, int> args = parse_args(argc, argv);
    int n0 = args["n0"];
    int n1 = args["n1"];

    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors...\n");
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table.\n\n");

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Get last detection in database
    int d_max = sp_run_int(conn, "KS.GetMaxDetectionID");
    // print("Max DetectionID: {:d}.\n", d_max);

    // Inputs to build DetectionTable and AsteroidSkypatchTable
    int d0 = 0;
    int d1 = d_max;
    bool progbar = true;

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1, progbar);
    // DetectionTable dt = DetectionTable(conn, progbar);

    // Initialize an empty vector of candidates, cv (for "candidate vector")
    vector<DetectionNearAsteroidCandidate> cv;

    // Status message for main loop over asteroid blocks
    int n_ast = (n1-n0);
    int n_batch = std::max((n_ast+1)/batch_size, 1);
    print("\nProcessing {:d} asteroids with AsteroidID in [{:d}, {:d}) in {:d} batches of size {:d}.\n", 
        n_ast, n0, n1, n_batch, batch_size);

    // Counter for the number of candidates found; used for inserting batches of records into DB
    int n_cand = 0;

    // Timer for progress bar
    Timer t;
    t.tick();

    // Column headers for progress updates
    print("{:17s} : {:10s} : {:17s} :\n", "Processed", "Average", "Remaining");
    print("{:6s} : {:8s} : {:10s} : {:6s} : {:8s} : {:7s}\n", "N_ast", "time", "sec/batch", "N_ast", "time", "Matches" );

    // Iterate over batches of asteroids
    for (int i0=n0; i0<n1; i0+=batch_size)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+batch_size, n1);
        // Build AsteroidSkyPatch table for this asteroid block
        bool progbar=false;
        AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, i0, i1, progbar);
        // Search for detection matching this asteroid block
        search_asteroid_detection(dt, aspt, spn, cv);
        // Range of records to insert
        int k0 = n_cand;
        int k1 = cv.size();
        // Insert this batch
        write_candidates_db(conn, cv, k0, k1);
        // Update the candidate counter
        n_cand = cv.size();

        // Progress indicator
        // print(".");
        // flush_console();

        // Estimate remaining time
        double t_proc = t.tock();
        int n_proc = i1 - n0;
        double t_per_ast = t_proc / n_proc;
        int n_left = n1 - i1;
        double t_left = n_left * t_per_ast;
        double t_per_batch = t_per_ast * batch_size;
        print("{:6d} : {:8s} : {:10.2f} : {:6d} : {:8s} : {:7d}\n", 
            n_proc, time2hms(t_proc), t_per_batch, n_left, time2hms(t_left), n_cand);
    } // for / i0

    // Report matches
    int matches = cv.size();
    print("\nFound {:d} matches.\n", matches);
    if (matches) {print("{:12s} {:12s}\n", "DetectionID", "AsteroidID");}
    for (DetectionNearAsteroidCandidate c: cv)
    {
        print("{:11d} {:11d}\n", c.detection_id, c.asteroid_id);
    }

    // Close DB connection
    conn->close();

    // Normal program exit
    return 0;
}

// *****************************************************************************
// Tests: Detection Table, AsteroidSkyPatch, an
// *****************************************************************************

// *****************************************************************************
void test_detection_table(DetectionTable& dt, int detection_id)
{
    // Print first 10 detections
    int i0=detection_id;
    int i1=std::min(dt.d1, i0+10);
    print("\nSample data: detections with IDs in [{:d}, {:d}) :\n", i0, i1);
    print("{:12s} {:10s} {:8s}  {:9s}  {:9s}  {:9s}\n", "DetectionID", "SkyPatchID", "mjd", "ux", "uy", "uz");
    for (int i=i0; i<i1;i++)
    {
        // The ith detection on the table
        Detection d = dt[i];
        // Holes in the detection vector indicated by d.detectio_id=0; skip these
        if (!d.detection_id) {continue;}
        // Print this detection
        print("{:11d} {:10d} {:8.3f} {:+9.6f} {:+9.6f} {:+9.6f}.\n", 
            d.detection_id, d.sky_patch_id, d.mjd, d.ux, d.uy, d.uz);
    }
}

// *****************************************************************************
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id)
{
    // Demonstrate searching by SkyPatchID
    // Detection d = dt[sky_patch_id];
    // int sky_patch_id = d.sky_patch_id;
    print("\nSearch detections with SkyPatchID={:d}:\n", sky_patch_id);
    for (int detection_id: dt.get_skypatch(sky_patch_id)) {print("{:d},", detection_id);}
    print_newline();
}

// *****************************************************************************
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt)
{
    // Display size
    print("AsteroidSkyPatchTable contains {:d} entries.\n", aspt.size());
    // Print first 10 entries
    int i0 = 0;
    int i1 = std::min(aspt.size(), 10);

    print("\nSample data: first {:d} AsteroidSkyPatch entries:\n", i1-i0);
    print("{:12s} {:12s} {:10s} {:10s}\n", "AsteroidID", "SkyPatchID", "TimeID_0", "TimeID_1");
    for (int i=i0; i<i1;i++)
    {
        AsteroidSkyPatch asp = aspt[i];
        print("{:10d} {:12d} {:10d} {:10d}\n", 
            asp.asteroid_id, asp.sky_patch_id, asp.time_id_0, asp.time_id_1);
    } // for / i
}

// *****************************************************************************
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn)
{
    print("\nRunning search function on {:d} asteroid segments and {:d} detections...\n", aspt.size(), dt.size());
    vector<DetectionNearAsteroidCandidate> cv;
    search_asteroid_detection(dt, aspt, spn, cv);
    int matches = cv.size();
    long pairs = dt.size() * aspt.size() * 9;
    print("Found {:d} matches in {:d} million possible detection / asteroid pairs.\n", matches, pairs/1000000);
    if (matches) {print("{:12s} {:12s}\n", "DetectionID", "AsteroidID");}
    for (DetectionNearAsteroidCandidate c: cv)
    {
        print("{:11d} {:11d}\n", c.detection_id, c.asteroid_id);
    }
}

// *****************************************************************************
void test_all()
{
    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors...\n");
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table.\n");

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Get last detection in database
    int d_max = sp_run_int(conn, "KS.GetMaxDetectionID");
    print("Max DetectionID: {:d}.\n", d_max);

    // Inputs to build DetectionTable and AsteroidSkypatchTable
    int d0 = 0;
    int d1 = 100;
    int n0 = 51000;
    int n1 = 51100;
    bool progbar = true;

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1, progbar);

    // Build AsteroidSkyPatch table
    print_newline();
    AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, n0, n1, progbar);

    // Test values of detection_id and skypatch_id
    int test_detection_id = 10;
    int test_skypatch_id = dt[test_detection_id].sky_patch_id;

    // Test detection table and search of detections by SkyPatchID
    test_detection_table(dt, test_detection_id);
    test_detection_table_by_sp(dt, test_skypatch_id);

    // Test AsteroidSkyPatch table
    test_asteroid_skypatch(aspt);

    // Test search function for asteroid / detection matches
    test_search(dt, aspt, spn);

    // Close DB connection
    conn->close();
}
