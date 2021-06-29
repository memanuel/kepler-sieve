/*****************************************************************************
 * Batch program to calculate candidate pairs (AsteroidID, DetectionID)
 * where an asteroid is in a neighboring skypatch of the detection.
 * Results further refined downstream in detection_near_ast.py.
 * 
 * Michael S. Emanuel
 * 2021-06-25
 * ****************************************************************************/

// *****************************************************************************
// Included libraries
// #include <utility>
//     using std::pair;

#include <fmt/format.h>
    using fmt::print;
    using fmt::format;

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

// *****************************************************************************
// Declare functions
vector<DetectionNearAsteroidCandidate>
    search_asteroid_detection(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);
void test_detection_table(DetectionTable& dt, int detection_id);
void test_detection_table_by_sp(DetectionTable& dt, int sky_patch_id);
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt);
void test_all();
void test_search(DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn);
int main();

// *****************************************************************************
/** Perform a search over detections in the detection table dt and
 *  AsteroidSkyPatch entries in the AsteroidSkyPatchTable aspt
 * */
vector<DetectionNearAsteroidCandidate> search_asteroid_detection(
    DetectionTable& dt, AsteroidSkyPatchTable& aspt, SkyPatchNeighbor& spn)
{
    // Counter for the number of matches found and pairs examined
    // int pairs = 0;
    // Reusable pointer to the 9 neighbors of a SkyPatch
    int32_t* ns;

    // Initialize an empty vector of candidates, cv (for "candidate vector")
    vector<DetectionNearAsteroidCandidate> cv;    

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

        // Iterate through all the neighbors
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
    // Return the candidates vector
    return cv;
}

// *****************************************************************************
int main()
{
    /*
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
    // DetectionTable dt = DetectionTable(conn, progbar);

    // Build AsteroidSkyPatch table
    print_newline();
    AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, n0, n1, progbar);

    // Test values of detection_id and skypatch_id
    // int test_detection_id = 348750;
    // int test_sky_patch_id = 14499016;
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
    */

    test_all();
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
    vector<DetectionNearAsteroidCandidate> cv = search_asteroid_detection(dt, aspt, spn);
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
