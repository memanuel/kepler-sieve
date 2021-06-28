/** @file AsteroidSkyPatch.hpp
 *  @brief Class to load a batch of AsteroidSkyPatch entries in memory and look them up by SkyPatchID.
 *  See DB table AsteroidSkyPatch and stored procedure GetDetectionDire.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Included libraries
#include <string>
    using std::string;
    using std::to_string;
#include <vector>
    using std::vector;

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "Timer.hpp"
    using ks::Timer;

// *****************************************************************************
namespace ks {

// *****************************************************************************
/**Data contained in one AsteroidSkyPatch entry.
 * See DB table KS.AsteroidSkyPatch and stored procedure KS.AsteroidSkyPatch.*/
struct AsteroidSkyPatch
{
    /// Integer ID of the asteroid whose path is described
    int32_t asteroid_id;
    /// The segment number of this asteroid's path in the sky
    int32_t segment;
    /// Integer ID of the SkyPatch where this detection is located in the sky on this segment
    int32_t sky_patch_id;
    /// Integer ID of the time when this segment starts; nearest multiple of 15 minutes
    int32_t time_id_0;
    /// Integer ID of the time when this segment ends; nearest multiple of 15 minutes
    int32_t time_id_1;
};

// *****************************************************************************
class AsteroidSkyPatchTable
{
public:
    /// Initialize an AsteroidSkyPatchTable object for asteroids in the given range
    AsteroidSkyPatchTable(db_conn_type &conn, int n0, int n1, bool progbar);

    /// Destructor for AsteroidSkyPatchTable.
    ~AsteroidSkyPatchTable();

    /// Get an AsteroidSkyPatch entry given its row number i
    AsteroidSkyPatch operator[](int32_t i) const;

private:
    // Data
    //*First asteroid ID loaded
    const int n0;
    //*Last asteroid ID loaded
    const int n1;
    //*Vector of AsteroidSkyPatch entries
    vector<AsteroidSkyPatch> asp;
};

// *****************************************************************************
} // namespace ks
