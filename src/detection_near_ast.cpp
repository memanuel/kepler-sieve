// *****************************************************************************
// Included libraries
#include <cmath>
#include <iostream>
#include <boost/format.hpp>

// Local dependencies
#include "utils.h"
#include "sky_patch.h"

// *****************************************************************************
// Names used
using std::cout;
using boost::format;
using ks::CubeFace;
using ks::SkyPatch;
using ks::sky_patch::N;
using ks::sky_patch::N_sp;
using ks::sky_patch::N_spn;
using ks::SkyPatch_from_id;
// using ks::spn_type;
using ks::make_sky_patch_neighbor_table;
using ks::print_stars;


// *****************************************************************************
int main()
{
    // Build the SkyPatchNeighbor table
    cout << format("Building SkyPatch neighbors for N = %d...\n") % N;
    /*
    spn_type spn = make_sky_patch_neighbor_table();

    // Initialize a starting SkyPatch
    int32_t spid0 = 0;
    SkyPatch sp0 = SkyPatch_from_id(spid0);

    // Read off neighbors of first row
    cout << format("Starting SkyPatch:\n%s") % sp0.str();
    cout << format("Neighbors of this SkyPatch:\n");
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
            // cout << format("spid: %d\n") % sp1.id();            
            cout << sp1.str();
        }
    }
    */
    // Normal program exit
    return 0;
}