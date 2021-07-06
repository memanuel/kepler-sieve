/** @file detection_near_ast.cpp
 *  @brief Batch program to search for detections that are near a known asteroid's trajectory.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-25
 * 
 * Batch program to calculate candidate pairs (AsteroidID, DetectionID)
 * where an asteroid is in a neighboring skypatch of the detection.
 * Results further refined downstream in detection_near_ast.py.
 * 
 * Commandline arguments:
 * jn: the job number
 * sz: the job size; defaults to 100000
 * n0: the first asteroid to process
 * n1: the last asteroid to process
 * This will process asteroids with asteroid numbers in [n0, n1).
 * When asteroid range is specified using jn and sz, the relationship is
 * n0 = jn*sz
 * n1 = n0+sz
 * 
 * Example calls:
 * ./detection_near_ast.x --jn 0 --sz 1000
 * ./detection_near_ast.x 0 1000
 * both of these will process asteroids in [0, 1000)
 * ****************************************************************************/

// Detailed help message for this program
#include <string>
    using std::string;
std::string help_message = 
R"(
detection_near_ast - 
Batch program to interactions where a known asteroid is within a close threshold
(e.g. 10 arc seconds) of a detection.
Data written to DB table KS.DetectionNearAsteroid, with one record taking the form
(DetectionID, AsteroidID, s, LightTime)
The two IDs identify the detection and asteroid.
s is the cartesian distance between the two directions and light time is measured in minutes.
The direction from asteroid to observatory calculated with splined orbital elements.

Commandline arguments:
jn: the job number
sz: the job size; defaults to 100000
n0: the first asteroid to process
n1: the last asteroid to process
This will process asteroids with asteroid numbers in [n0, n1).
When asteroid range is specified using jn and sz, the relationship is
n0 = jn*sz
n1 = n0+sz

Example calls:
./detection_near_ast.x --jn 0 --sz 1000
./detection_near_ast.x 0 1000
both of these will process asteroids in [0, 1000)
)";

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;
#include <boost/program_options.hpp>
    namespace po = boost::program_options;
#include <gsl/gsl_spline.h>

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::sql_stmt_type;
    using ks::sql_prepared_stmt_type;
    using ks::get_db_conn;
    using ks::sp_run;
    using ks::sp_run_int;

#include "utils.hpp"
    using ks::is_close_abs;
    using ks::print_stars;
    using ks::print_newline;
    using ks::report_test;
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

#include "BodyVector.hpp"
    using ks::BodyVector;

#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::norm;

#include "AsteroidElement.hpp"
    using ks::AsteroidElement;

// *****************************************************************************
// Data type to describe one detection near one known asteroid
struct DetectionNearAsteroid
{
    /// The integer ID of the detection
    int32_t detection_id;
    /// The integer ID of the asteroid
    int32_t asteroid_id;
    /// The cartesian distance on the unit sphere between directions of the asteroid and detection
    double s;
    /// The time in minutes that light traveled before reaching the detector
    double light_time;
};

// Batch size used for processing asteroids
constexpr int batch_size = 1000;

// *****************************************************************************
// Functions defined in this module
void search_candidates(
    DetectionTable& dt, AsteroidSkyPatchTable& aspt, 
    SkyPatchNeighbor& spn, vector<DetectionNearAsteroid>& cv);
void calculate_directions(
    DetectionTable& dt, vector<DetectionNearAsteroid>& dv, int k0, int k1);
void write_detections_db(db_conn_type& conn, const vector<DetectionNearAsteroid>& cv, int k0, int k1);

// Test functions defined in detection_near_ast_test.cpp
void test_detection_table(DetectionTable& dt, int detection_id);
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id);
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt);
void test_all();
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);

// *****************************************************************************
// Main
// *****************************************************************************
int main(int argc, char* argv[])
{
    // *****************************************************************************
    // Process commandline arguments.
    // Result of this stage:
    // Set variables n0, n1 with asteroid range.
    // Set flag 
    // *****************************************************************************

    // The job number and batch size from commandline arguments 
    int jn = -1;
    int sz = -1;
    // The asteroid range from commandline arguments
    int n0 = -1;
    int n1 = -1;
    // The detection range from commandline arguments; usually omitted to process all detections
    int d0 = -1;
    int d1 = -1;
    // Flags from commandline arguments
    bool is_dry_run = false;
    bool print_matches = false;

    // Set up parser for named commandline arguments ("options")
    po::options_description desc("Find detections near asteroids");
    desc.add_options()
        ("help,h", "Produce help message")
        ("jn,j", po::value<int>(&jn)->default_value(0), "Job number")
        ("sz,s", po::value<int>(&sz)->default_value(100000), "Batch size for jobs")
        ("n0", po::value<int>(&n0)->default_value(-1), "First asteroid to process; omit to set using jn and sz")
        ("n1", po::value<int>(&n1)->default_value(-1), "Last asteroid to process; omit to set using jn and sz")
        ("d0", po::value<int>(&d0)->default_value(-1), "First detection to process; omit to process all detections")
        ("d1", po::value<int>(&d1)->default_value(-1), "Last detection to process; omit to process all detections")
        ("dryrun", po::bool_switch(&is_dry_run), "Dry run - report commandline arguments and quit early")
        ("print_matches", po::bool_switch(&print_matches), "Print matches to console at end")
    ;
    po::variables_map vm;

    // Add positional arguments for jn and sz in first two positions
    po::positional_options_description pos;
    pos.add("jn", 1);
    pos.add("sz", 2);

    // Parse commandline arguments including positional arguments
    // po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
    po::notify(vm);

    // If help requested, print description and quit early
    extern std::string help_message;
    if (vm.count("help")) 
    {
        std::cout << desc;
        print(help_message);
        return 0;
    }
    
    // Set flag indicating whether asteroid range was set manually
    bool is_manual_ast_range = (n0>=0) && (n1>=0);
    // Set flag indicating whether detection range was set
    bool is_manual_detection_range = (d0>=0) && (d1>=0);

    // Calculate asteroid range from jn and sz if not set directly
    if (n0<0) {n0 = jn*sz;}
    if (n1<0) {n1 = n0+sz;}
    // Calculate batch size from input asteroid range if applicable
    if (is_manual_ast_range) {sz = n1-n0;}
    // Report commandline arguments and the resulting asteroid range
    print("Commandline arguments\n");
    if (!is_manual_ast_range) 
    {
        print("Job Number jn       : {:d}\n", jn);
        print("Batch Size sz       : {:d}\n", sz);
    }
    print("Start AsteroidID n0 : {:d}\n", n0);
    print("End AsteroidID   n1 : {:d}\n", n1);
    if (is_manual_ast_range)
    {
        print("Batch Size sz       : {:d}\n", sz);
    }
    if (is_manual_detection_range)
    {
        print("First DetectionID   : {:d}\n", d0);
        print("Last DetectionID    : {:d}\n", d1);        
    }
    else 
    {
        print("Detection range     : ALL\n");
    }
    print_newline();

    // Quit early if it was a dry run 
    if (is_dry_run) 
    {
        print("This was a dry run. Bye!\n");
        return 0;
    }

    // *****************************************************************************
    // Program body.
    // Result of this stage:
    // Build SkyPatch neighbor table spn
    // Build DetectionCandidateTable dt
    // Loop over blocks of asteroids in range [n0, n1).
    // For each block build an AsteroidSkyPatchTable table aspt and then
    // perform a search over dt and aspt.
    // During the search, write output to both the candidate vector cv in memort
    // as well as directly to the database.
    // *****************************************************************************

    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors...\n");
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table.\n\n");

    // Initialize DetectionCandidateTable
    // bool progbar = true;
    // Load the detection candidate table; process either the selected manual range or all detections (default)
    // DetectionTable dt = is_manual_detection_range ? 
    //     DetectionTable(conn, d0, d1, progbar) :
    //     DetectionTable(conn, progbar);
    DetectionTable dt;
    dt.load();
    print("Loaded Detection table from disk with {:d} detections.\n", dt.size());
    
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Initialize an empty vector of candidates, cv (for "candidate vector")
    vector<DetectionNearAsteroid> cv;

    // Status message for main loop over asteroid blocks
    int n_ast = (n1-n0);
    int n_batch = std::max(n_ast/batch_size, 1);
    print("\nProcessing {:d} asteroids with AsteroidID in [{:d}, {:d}) in {:d} batches of size {:d}.\n", 
        n_ast, n0, n1, n_batch, batch_size);

    // Counter for the number of candidates found; used for inserting batches of records into DB
    int n_cand = 0;
    // Timer for progress bar
    Timer t;
    t.tick();
    // Column headers for progress updates
    print("{:17s} : {:10s} : {:17s} : {:7s}\n", 
          "Processed", "Average", "Remaining", "Matches");
    print("{:6s} : {:8s} : {:10s} : {:6s} : {:8s} : {:7s}\n", 
          "N_ast", "time", "sec/batch", "N_ast", "time", "Found" );

    // Iterate over batches of asteroids; i0 is first asteroid ID of the current batch
    for (int i0=n0; i0<n1; i0+=batch_size)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+batch_size, n1);
        // Build AsteroidSkyPatch table for this asteroid block
        bool progbar=false;
        AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, i0, i1, progbar);
        // Search for camdidate detections matching this asteroid block
        search_candidates(dt, aspt, spn, cv);
        // Range of records to insert
        int k0 = n_cand;
        int k1 = cv.size();
        // Insert this batch of results to the database
        write_detections_db(conn, cv, k0, k1);
        // Update the candidate counter
        n_cand = cv.size();

        // Estimate remaining time and print status
        double t_proc = t.tock();
        int n_proc = i1 - n0;
        double t_per_ast = t_proc / n_proc;
        int n_left = n1 - i1;
        double t_left = n_left * t_per_ast;
        double t_per_batch = t_per_ast * batch_size;
        print("{:6d} : {:8s} : {:10.2f} : {:6d} : {:8s} : {:7d}\n", 
            n_proc, time2hms(t_proc), t_per_batch, n_left, time2hms(t_left), n_cand);
    } // for / i0

    // Report number of matches
    int matches = cv.size();
    print("\nFound {:d} matches.\n", matches);
    // Report the matches themselves if requested
    if (print_matches && matches) 
    {
        // Maximum number of matches to report; dont' want to spam the console
        constexpr int max_report_size = 100;
        int j_max = std::min(matches, max_report_size);
        print("{:12s} {:12s}\n", "DetectionID", "AsteroidID");
        for (int j=0; j<j_max; j++)
        {
            DetectionNearAsteroid c = cv[j];
            print("{:11d} {:11d}\n", c.detection_id, c.asteroid_id);
        } // for over matches
    } // if reporting matches

    // Close DB connection
    conn->close();

    // Normal program exit
    return 0;

}

// *****************************************************************************
// Functions to perform search and write results to database.
// *****************************************************************************

// *****************************************************************************
/** Perform a search over detections in the detection table dt and
 *  AsteroidSkyPatch entries in the AsteroidSkyPatchTable aspt.
 *  \param[in] dt - table of detection data to search over; "detection table"
 *  \param[in] aspt - table of asteroid sky patch segments to search over; "asteroid sky patch table"
 *  \param[in] spn - table with the 9 neighbors of each sky patch; "sky patch neighbor"
 *  \param[in] cv - vector of candidate detections near an asteroid; "candidate vector"
 * */
void search_candidates(
    DetectionTable& dt, AsteroidSkyPatchTable& aspt, 
    SkyPatchNeighbor& spn, vector<DetectionNearAsteroid>& cv)
{
    // Counter for the number of matches found and pairs examined
    // int pairs = 0;
    // Reusable pointer to the 9 neighbors of a SkyPatch
    int32_t* ns;

    // Iterate through all the AsteroidSkyPatch entries in the table
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

        // Iterate through all 9 neighbors of this sky patch
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
                    DetectionNearAsteroid c = {
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
void write_detections_db(
    db_conn_type& conn, const vector<DetectionNearAsteroid>& dv,
    int k0, int k1)
{
     // SQL to insert one record into DetectionNearAsteroid
    string sql = R"(
    REPLACE INTO KS.DetectionNearAsteroid
    (DetectionID, AsteroidID, s, LightTime) 
    VALUES
    (?, ?, ?, ?);
    )";
    // Create a SQL PreparedStatement for the insert
    sql_prepared_stmt_type stmt(conn->prepareStatement(sql));

    // Iterate through the selected block of candidates in this batch
    for (int k=k0; k<k1; k++)
    {
        // The candidate c is row k of the candidate vector cv
        DetectionNearAsteroid d = dv[k];
        // Bind two integer parameters to the prepared statement    
        stmt->setInt(1, d.detection_id);
        stmt->setInt(2, d.asteroid_id);
        stmt->setDouble(3, d.s);
        stmt->setDouble(4, d.light_time);
        // Execute the SQL statement
        stmt->executeQuery();
    } // for / k 
}
