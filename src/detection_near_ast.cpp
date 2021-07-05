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
// Declare functions
void search_candidates(
    DetectionTable& dt, AsteroidSkyPatchTable& aspt, 
    SkyPatchNeighbor& spn, vector<DetectionNearAsteroid>& cv);
void calculate_directions(
    DetectionTable& dt, vector<DetectionNearAsteroid>& dv, int k0, int k1);
void write_detections_db(db_conn_type& conn, const vector<DetectionNearAsteroid>& cv, int k0, int k1);
void test_detection_table(DetectionTable& dt, int detection_id);
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id);
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt);
void test_all();
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);

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
    bool is_test = false;
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
        ("test", po::bool_switch(&is_test), "Test - run tests and quit early")
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

    // Run tests then quit early if in test mode
    if (is_test) 
    {
        print("Running in test mode:\n\n");
        test_all();
        return 0;
    }

    // *****************************************************************************
    // Main body of program
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

// *****************************************************************************
// Tests: DetectionTable and AsteroidElement
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
/// Test that AsteroidElement instance splines orbital elements that match copy / pasted values for Ceres @ 58400.
void test_asteroid_element(AsteroidElement& ast_elt)
{
    // This test is designed to check Ceres (asteroid_id=1) at mjd 58400 (time_idx=100)
    int asteroid_idx = 1;
    int time_idx = 100;

    // Expected results, copy / pasted from database from KS.GetAsteroidElements(1, 2, 58400, 58400);
    double a0     = 2.7673528257126296;
    double e0     = 0.07561068735641437;
    double inc0   = 0.18489327222145555;
    double Omega0 = 1.4016429441591765;
    double omega0 = 1.2779179332926616;
    double f0     = 38.4030882294433;
    double M0     = 38.30932258771573;
    
    // Tolerance for tests
    double tol_a = 1.0E-14;
    double tol_e = 1.0E-14;
    double tol_angle = 1.0E-12;

    // Read off some asteroid elements
    print("\nAsteroidElement properties:\n");
    print("N_ast    : {:d}\n", ast_elt.N_ast);
    print("N_t      : {:d}\n", ast_elt.N_t);

    // Read the two 1D arrays
    int32_t* asteroid_ids = ast_elt.get_asteroid_id();
    double* mjds = ast_elt.get_mjd();
    // The selected asteroid_id and time
    int32_t asteroid_id = asteroid_ids[asteroid_idx];
    double mjd = mjds[time_idx];

    // Calulate splined orbital elements; these should match
    OrbitalElement elt = ast_elt.interp_elt(asteroid_id, mjd);

    // Report splined orbital elements
    print("\nSplined Asteroid index {:d} at time index {:d} :\n", asteroid_idx, time_idx);
    print("AsteroidID: {:9d}\n", asteroid_id);
    print("mjd:        {:9.2f}\n", mjd);
    print("a:          {:9.6f}\n", elt.a);
    print("e:          {:9.6f}\n", elt.e);
    print("inc:        {:9.6f}\n", elt.inc);
    print("Omega:      {:9.6f}\n", elt.Omega);
    print("omega:      {:9.6f}\n", elt.omega);
    print("f:          {:9.6f}\n", elt.f);
    print("M:          {:9.6f}\n", elt.M);

    // Test that splined orbital elements match expected results
    bool is_ok = true;
    is_ok = is_ok && is_close_abs(a0,     elt.a, tol_a);
    is_ok = is_ok && is_close_abs(e0,     elt.e, tol_e);
    is_ok = is_ok && is_close_abs(inc0,   elt.inc, tol_angle);
    is_ok = is_ok && is_close_abs(Omega0, elt.Omega, tol_angle);
    is_ok = is_ok && is_close_abs(omega0, elt.omega, tol_angle);
    is_ok = is_ok && is_close_abs(f0,     elt.f, tol_angle);
    is_ok = is_ok && is_close_abs(M0,     elt.M, tol_angle);
    report_test("\nTest AsteroidElement::interp_elt() splined elements match database", is_ok);
}

// *****************************************************************************
/// Test that AsteroidElement instance splines state vectors that match copy / pasted values for Ceres @ 58400.
void test_asteroid_element_vectors(AsteroidElement& ast_elt)
{
    // This test is designed to check Ceres (asteroid_id=1) at mjd 58400 (time_idx=100)
    int asteroid_idx = 1;
    int time_idx = 100;

    // Read the two 1D arrays
    int32_t* asteroid_ids = ast_elt.get_asteroid_id();
    double* mjds = ast_elt.get_mjd();
    // The selected asteroid_id and time
    int32_t asteroid_id = asteroid_ids[asteroid_idx];
    double mjd = mjds[time_idx];
    
    // Asteroid position; copy / pasted from database from KS.GetAsteroidElements(1, 2, 58400, 58400);
    double qx_ast = -2.485854955060206;
    double qy_ast = -0.6229239033462274;
    double qz_ast =  0.4383572489736212;
    double vx_ast =  0.0020617208493692073;
    double vy_ast = -0.010756337836425017;
    double vz_ast = -0.0007200628373353;
    // State vectors of the sun; copy / pasted from KS.GetStateVectors_Sun(58400, 58400, 1);
    double qx_sun = -0.0001095835967748;
    double qy_sun =  0.007235858951602957;
    double qz_sun = -0.0000736284237584;
    double vx_sun = -0.0000075730407095;
    double vy_sun =  0.0000026357733813;
    double vz_sun =  0.0000001892676823;
    // Expected results; heliocentric position of asteroid
    double qx = qx_ast - qx_sun;
    double qy = qy_ast - qy_sun;
    double qz = qz_ast - qz_sun;
    double vx = vx_ast - vx_sun;
    double vy = vy_ast - vy_sun;
    double vz = vz_ast - vz_sun;
    // Wrap expected results into Position and Vector objects
    Position pos0 = Position 
    {
        .qx = qx,
        .qy = qy,
        .qz = qz
    };
    StateVector vec0 = StateVector
    {
        .qx = qx,
        .qy = qy,
        .qz = qz,
        .vx = vx,
        .vy = vy,
        .vz = vz
    };

    // Tolerance for tests
    double tol_q = 1.0E-8;
    double tol_vec = 1.0E-8;
    // double tol_v = 1.0E-10;

    // Calulate splined position and state vector
    Position pos = ast_elt.interp_pos(asteroid_id, mjd);
    StateVector vec = ast_elt.interp_vec(asteroid_id, mjd);

    // Report splined orbital elements
    print("\nSplined Asteroid index {:d} at time index {:d} :\n", asteroid_idx, time_idx);
    print("AsteroidID: {:9d}\n", asteroid_id);
    print("mjd:        {:9.2f}\n", mjd);
    print("qx:         {:+9.6f}\n", vec.qx);
    print("qy:         {:+9.6f}\n", vec.qy);
    print("qz:         {:+9.6f}\n", vec.qz);
    print("vx:         {:+9.6f}\n", vec.vx);
    print("vy:         {:+9.6f}\n", vec.vy);
    print("vz:         {:+9.6f}\n", vec.vz);

    // Test that splined orbital elements match expected results (position only)
    {
    bool is_ok = norm(pos0, pos) < tol_q;
    report_test("\nTest AsteroidElement::interp_pos() splined state vectors match database", is_ok);
    }

    // Test that splined orbital elements match expected results
    {
    bool is_ok = norm(vec0, vec) < tol_vec;
    report_test("\nTest AsteroidElement::interp_vec() splined state vectors match database", is_ok);
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

    // Inputs to build DetectionCandidateTable, AsteroidSkypatchTable and AsteroidElement
    int d0 = 0;
    int d1 = 100;
    int n0 = 0;
    int n1 = 100;
    int mjd0 = 58000;
    int mjd1 = 61000;
    int time_step = 4;
    bool progbar = true;

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1, progbar);

    // Test DetectionTable
    test_detection_table(dt, d0);

    // Build AsteroidSkyPatch table
    print_newline();
    AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, n0, n1, progbar);

    // Initialize AsteroidElement
    AsteroidElement elt = AsteroidElement(n0, n1, mjd0, mjd1, time_step);
    elt.load(conn, progbar);

    // Test asteroid elements - splining orbital elements
    test_asteroid_element(elt);

    // Test asteroid elements - state vectors from splined elements
    test_asteroid_element_vectors(elt);

    // Close DB connection
    conn->close();
}
