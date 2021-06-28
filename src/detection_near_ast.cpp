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
// #include "sky_patch.h"
#include "SkyPatchNeighbor.h"

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
using ks::CubeFace;
using ks::SkyPatch;
using ks::sky_patch::N;
using ks::sky_patch::N_sp;
using ks::sky_patch::N_spn;
using ks::SkyPatch_from_id;
// using ks::spn_type;
// using ks::make_sky_patch_neighbor_table;
// Utilities
using ks::print_stars;

// *****************************************************************************
int main()
{
    /*
    // Build the SkyPatchNeighbor table
    print("Building SkyPatch neighbors for N = {:d}...\n", N);
    spn_type spn = make_sky_patch_neighbor_table();
    */

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Run the stored procedure
    string sp_name = "KS.GetDetections";
    vector<string> params = {"0", "10"};
    ResultSet *rs = sp_run(conn, sp_name, params);
    // ResultSet *rs = sp_run(conn);

    // Get size of resultset
    int rows = result_set_size(rs);
    print("rows={:d}.\n", rows);

    // Loop through resultset
    print("{:12s} {:8s}  {:8s}  {:8s}  {:8s}\n", "DetectionID", "mjd", "ux", "uy", "uz");
    while (rs->next()) 
    {
        int32_t detection_id = rs->getInt("DetectionID");
        double mjd = rs->getDouble("mjd");
        double ux = rs->getDouble("ux");
        double uy = rs->getDouble("uy");
        double uz = rs->getDouble("uz");
        print("{:12d} {:8.3f} {:8.6f} {:8.6f} {:8.6f}.\n", detection_id, mjd, ux, uy, uz);
    }
    // Close Connection
    conn->close();

    /*
    // Initialize a starting SkyPatch
    int32_t spid0 = 0;
    SkyPatch sp0 = SkyPatch_from_id(spid0);

    // Read off neighbors of first row
    print("Starting SkyPatch:\n{:s}", sp0.str());
    print("Neighbors of this SkyPatch:\n");
    // Offset into table for sp0
    int32_t idx0 = spid0*9;
    for (int j=0; j<9; j++)
    {
        // The jth neighbor
        int32_t spid1 = spn[idx0+j];
        // Only process *real* neighbors with non-negative spids
        if (spid1 >= 0)
        {
            SkyPatch sp1 = SkyPatch_from_id(spid1);
            print(sp1.str());
        }
    }
    */
    // Normal program exit
    return 0;
}
