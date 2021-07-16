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

// Set batch size; this is the number of ASTEROIDS, not rows!
constexpr int batch_size = 100;

// *****************************************************************************
/**Helper function: Process a batch of rows, 
 * writing data from stored preocedure output to vector of detections. */
void process_rows(db_conn_type& conn, vector<AsteroidSkyPatch>& aspt, int n0, int n1)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetAsteroidSkyPatch";
    vector<string> params = {to_string(n0), to_string(n1)};
    ResultSet* rs = sp_run(conn, sp_name, params);

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
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
AsteroidSkyPatchTable::AsteroidSkyPatchTable(db_conn_type &conn, int n0, int n1, bool progbar): 
    n0 {n0},
    n1 {n1},
    aspt {vector<AsteroidSkyPatch>()}
{
    // Status update
    if (progbar) 
    {
        int ast_count = n1-n0;
        int batch_count = std::max(ast_count/batch_size, 1);
        print("Processing {:d} asteroids of AsteroidSkyPatch data from {:d} to {:d} in {:d} batches of size {:d}...\n",
                ast_count, n0, n1, batch_count, batch_size);
    }

    // Timer for processing AsteroidSkyPatch from DB
	Timer t;
    t.tick();

    // Iterate over the batches; i0 is the first asteroid number of the batch being processed
    for (int i0=n0; i0<n1; i0+=batch_size)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+batch_size, n1);
        // Process SQL data in this batch
        process_rows(conn, aspt, i0, i1);
        // Progress bar
        if (progbar) 
        {
            print(".");
            flush_console();            
        }
    }
    if (progbar) 
    {
        print("\nLoaded AsteroidSkyPatch table.\n");
        t.tock_msg();
    }
}

// *****************************************************************************
///Default destructor is OK here
AsteroidSkyPatchTable::~AsteroidSkyPatchTable() {}

// *****************************************************************************
const AsteroidSkyPatch AsteroidSkyPatchTable::operator[](int32_t i) const
    {return aspt.at(i);}

// *****************************************************************************
const int AsteroidSkyPatchTable::size() const
    {return aspt.size();}
