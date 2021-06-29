/** @file AsteroidSkyPatch.cpp
 *  @brief Implmentation of AsteroidSkyPatch class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

// *****************************************************************************
// Included files
#include "AsteroidSkyPatch.hpp"

// *****************************************************************************
// Local names used
using ks::AsteroidSkyPatch;
using ks::AsteroidSkyPatchTable;

// Set batch size; this is in terms of the number of asteroids, not rows
constexpr int bs = 100;

// *****************************************************************************
/** Helper function: Process a batch of rows, 
 * writing data from stored preocedure output to vector of detections. */
void process_rows(db_conn_type &conn, vector<AsteroidSkyPatch> &aspt, int n0, int n1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetAsteroidSkyPatch";
    vector<string> params = {to_string(n0), to_string(n1)};
    ResultSet *rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset; 10 total fields
        int32_t asteroid_id = rs->getInt("AsteroidID");
        int32_t segment = rs->getInt("Segment");
        int32_t sky_patch_id = rs->getInt("SkyPatchID");
        int32_t time_id_0 = rs->getInt("TimeID_0");
        int32_t time_id_1 = rs->getInt("TimeID_1");

        // Build the AsteroidSkyPatch at this location
        AsteroidSkyPatch asp = {
            .asteroid_id = asteroid_id,
            .segment = segment,
            .sky_patch_id = sky_patch_id,
            .time_id_0 = time_id_0,
            .time_id_1 = time_id_1,
        };
        // Save it to the asteroid sky patch table
        aspt.push_back(asp);
    }   // while rs
    // Close the resultset
    rs->close();
}

// *****************************************************************************
AsteroidSkyPatchTable::AsteroidSkyPatchTable(db_conn_type &conn, int n0, int n1, bool progbar): 
    n0(n0),
    n1(n1),
    aspt(vector<AsteroidSkyPatch>())
{
    // Number of asteroids to be processed
    int ast_count = n1-n0;

    // Status update
    if (progbar) 
    {
        print("Processing {:d} asteroids of AsteroidSkyPatch data from {:d} to {:d} in batches of {:d}...\n",
                ast_count, n0, n1, bs);
    }

    // Timer for processing AsteroidSkyPatch from DB
	Timer t;
    t.tick();

    // Iterate over the batches
    for (int i0=0; i0<n1; i0+=bs)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+bs, n1);
        // Process SQL data in this batch
        process_rows(conn, aspt, i0, i1);
        // Progress bar
        if (progbar) {print("."); }
    }
    if (progbar) 
    {
        print("\nLoaded detection table.\n");
        t.tock_msg();
    }
}

// *****************************************************************************
///Default destructor is OK here
AsteroidSkyPatchTable::~AsteroidSkyPatchTable() {}

// *****************************************************************************
AsteroidSkyPatch AsteroidSkyPatchTable::operator[](int32_t i) const
{
    return aspt[i];
}
