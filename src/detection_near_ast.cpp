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
using ks::SkyPatch_from_id;
using ks::sky_patch_count;
using ks::write_sky_patch_neighbor_table;
using ks::print_stars;

// The grid size
constexpr int N = ks::N_sky_patch;

// *****************************************************************************
int main()
{
    // Total number of sky patches is known at compile time
    int32_t N_spc = sky_patch_count();
    // Allocate an array of size 9*N_spc to hold the 9 neighbors of each patch
    int32_t *spn = new int32_t [N_spc*9];

    // Build the SkyPatchNeighbor table
    cout << format("Building SkyPatch neighbors for N = %d...\n") % N;
    write_sky_patch_neighbor_table(spn);

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

    // Normal program exit
    return 0;
}