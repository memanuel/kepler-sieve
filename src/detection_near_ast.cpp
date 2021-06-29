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

// *****************************************************************************
void test_detection_table(DetectionTable& dt)
{
    // Print first 10 detections
    int i0=dt.d0;
    int i1=std::min(dt.d1, 10);
    print("\nSample data: first {:d} detections:\n", i1);
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
    
    // Demonstrate searching by SkyPatchID
    Detection d = dt[1];
    int sky_patch_id = d.sky_patch_id;
    print("\nSearch detections with SkyPatchID={:d}:\n", sky_patch_id);
    for (int detection_id: dt.get_skypatch(sky_patch_id)) {print("{:d},", detection_id);}
    print_newline();

    // Close Connection
    conn->close();

    // Normal program exit
    return 0;
}
