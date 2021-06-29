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
// #include <cmath>
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
void test_detection_table(DetectionTable& dt)
{
    // Print first 10 detections
    int i0=dt.d0;
    int i1=std::min(dt.d1, 10);
    print("\nSample data: first {:d} detections:\n", i1-i0);
    print("{:12s} {:10s} {:8s}  {:9s}  {:9s}  {:9s}\n", "DetectionID", "SkyPatchID", "mjd", "ux", "uy", "uz");
    for (int i=i0; i<i1;i++)
    {
        Detection d = dt[i];
        if (d.detection_id < 0) {continue;}
        print("{:11d} {:10d} {:8.3f} {:+9.6f} {:+9.6f} {:+9.6f}.\n", 
            d.detection_id, d.sky_patch_id, d.mjd, d.ux, d.uy, d.uz);
    }

}

// *****************************************************************************
void test_detection_table_by_sp(DetectionTable& dt)
{
    // Demonstrate searching by SkyPatchID
    Detection d = dt[1];
    int sky_patch_id = d.sky_patch_id;
    print("\nSearch detections with SkyPatchID={:d}:\n", sky_patch_id);
    for (int detection_id: dt.get_skypatch(sky_patch_id)) {print("{:d},", detection_id);}
    print_newline();
}

// *****************************************************************************
void test_asteroid_skypatch(AsteroidSkyPatchTable& aspt)
{
    // Display size
    print("AsteroidSkyPatchTable aspt contains {:d} entries.\n", aspt.size());
    // Print first 10 entries
    int i0 = 0;
    int i1 = std::min(aspt.size(), 10);

    print("\nSample data: first {:d} AsteroidSkyPatch entries:\n", i1-i0);
    print("{:12s} {:12s} {:10s}  {:10s}\n", "AsteroidID", "SkyPatchID", "TimeID_0", "TimeID_1");
    for (int i=i0; i<i1;i++)
    {
        AsteroidSkyPatch asp = aspt[i];
        print("{:12d} {:12d} {:10d} {:10d}\n", 
            asp.asteroid_id, asp.sky_patch_id, asp.time_id_0, asp.time_id_1);
    }

}


// *****************************************************************************
/** Perform a search over detections in the detection table dt and
 *  AsteroidSkyPatch entries in the AsteroidSkyPatchTable aspt
 * */
int search(DetectionTable dt, AsteroidSkyPatchTable aspt)
{
    return 0;
}

// *****************************************************************************
int main()
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

    // Inputs for DetectionTable constructor
    int d0 = 0;
    int d1 = 10000;
    bool progbar = true;
    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1, progbar);
    // print("Loaded detection table with detections {:d} to {:d}.\n", d0, d1);
    // DetectionTable dt = DetectionTable(conn, progbar);

    // Build AsteroidSkyPatch table
    int n0=0;
    int n1=10;
    print_newline();
    AsteroidSkyPatchTable aspt = AsteroidSkyPatchTable(conn, n0, n1, progbar);

    // Test detection table and search by SkyPatchID
    test_detection_table(dt);
    test_detection_table_by_sp(dt);

    // Test AsteroidSkyPatch table
    test_asteroid_skypatch(aspt);

    // Close Connection
    conn->close();

    // Normal program exit
    return 0;
}
