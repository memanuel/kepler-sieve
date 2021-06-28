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
#include <cmath>
#include <fmt/format.h>

// Local dependencies
#include "db_utils.h"
#include "utils.h"
#include "SkyPatchNeighbor.h"
#include "Detection.h"

// *****************************************************************************
// Names used
using fmt::print;
using fmt::format;

// DB utilities
using ks::db_conn_type;
using ks::sql_stmt_type;
using ks::sql_prepared_stmt_type;
using ks::get_db_conn;
using ks::sp_run;
using ks::result_set_size;

// SkyPatch
using ks::SkyPatch;
using ks::SkyPatchNeighbor;
using ks::Detection;
using ks::DetectionTable;
using ks::sky_patch::N_sp;
// using ks::sky_patch::N_spn;

// Utilities
using ks::print_stars;

// *****************************************************************************
int main()
{
    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors...\n");
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table spn.\n");

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Range of detections
    int d0 = 0;
    int d1 = 10;
    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1);

    // Get detections
    print("{:12s} {:8s}  {:9s}  {:9s}  {:9s}\n", "DetectionID", "mjd", "ux", "uy", "uz");
    for (int i=d0; i<d1;i++)
    {
        Detection d = dt[i];
        if (d.detection_id < 0) {continue;}
        print("{:11d} {:8.3f} {:+9.6f} {:+9.6f} {:+9.6f}.\n", d.detection_id, d.mjd, d.ux, d.uy, d.uz);
    }

    // Close Connection
    conn->close();

    // Normal program exit
    return 0;
}
