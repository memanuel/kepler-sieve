/** @file   SkyPatchNeighbor.cpp
 *  @brief  Implmentation of SkyPatchNeighbor class.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-06-28
 */

// *****************************************************************************
// Local dependencies
#include "SkyPatchNeighbor.hpp"
    using ks::SkyPatch;
    using ks::SkyPatchNeighbor;
    using ks::make_SkyPatch_table;
    using ks::sky_patch::N_spn;

// *****************************************************************************
SkyPatchNeighbor::SkyPatchNeighbor(): 
    // Allocate an array of size N_spn (number of sky patch neighbors) to hold the 9 neighbors of each patch
    spn {new int32_t[N_spn]},
    // Allocate array to hold the neighbor distances
    spnd {new double[N_spn]}
{
    // Populate the neighbor table
    // Loop through the starting sky patch, with ID sky_patch_id_1
    for (int8_t f=0; f<6; f++)
    {
        // The offset for this grid face is M2*f
        int idx_f = M2*f;
        for (int16_t i0=0; i0<M; i0++)
        {
            for (int16_t j0=0; j0<M; j0++)
            {
                // The starting SkyPatch and its ID
                SkyPatch sp = SkyPatch(f, i0, j0);
                int32_t spid0 = sp.id();

                // The starting index into the neighbors table
                int idx = spid0*9;

                // Get grid coordinates of 9 candidate neighbors.
                // First iterate over the shift in the i index, di               
                for (int16_t di=-1; di<=1; di++)
                {
                    // The grid entry i1 for the three candidate neighbors in this row
                    int16_t i1 = i0+di;
                    // The additive term to the index into the sky patch table for this row
                    int idx_i1 = M*i1;
                    // Now iterate over the shift in the j index, dj
                    for (int16_t dj=-1; dj<=1; dj++)
                    {
                        // The grid entry j1 for this candidate neighbor
                        int16_t j1 = j0+dj;
                        // Storage for the new SkyPatchID; initialize to signify missing neighbors (e.g. for corners)
                        int32_t spid1=-1;

                        // Determine if we're wrapping in the i and j directions
                        bool is_on_grid_i = (0 <= i1) && (i1 < M);
                        bool is_on_grid_j = (0 <= j1) && (j1 < M);
                        // Is this the simple case where we are on the same cube face?
                        bool is_on_grid = is_on_grid_i & is_on_grid_j;

                        // Simple case; we're on the grid, use fast calculation
                        if (is_on_grid) 
                        {
                            spid1 = idx_f+idx_i1+j1;
                        }
                        // Handle shifts that wrap to another face, including diagonals
                        else 
                        {
                            // Is this the special case where we are trying to double wrap around a corner?
                            bool is_corner = (!is_on_grid_i) & (!is_on_grid_j);
                            // Exclude case of corners; want to write -1 here, not bogus SkyPatchID.
                            if (!is_corner) {spid1 = sp.shift(di, dj).id();}
                        }
                        // Write the new SkyPatchID to the spn array
                        spn[idx++] = spid1;
                    } // for over dj
                } // for over di
            } // for over j0
        } // for over i0
    } // for over f

    // Do *NOT* populate the neighbor distance table!
    // This is not always needed, so only build it on demand.
};

// *****************************************************************************
/// Need to manually delete two arrays that were allocated manually
SkyPatchNeighbor::~SkyPatchNeighbor()
{
    delete [] spn;
    delete [] spnd;
}

// *****************************************************************************
/// Need to manually delete two arrays that were allocated manually
void SkyPatchNeighbor::build_neighbor_distance()
{
    // Build SkyPatch table with precomputed SkyPatch positions.
    vector<SkyPatch> spt = make_SkyPatch_table();

    // Direction of the first and second skypatch; reused in inner loop.
    double u0[3];
    double u1[3];

    // Build a vector of all the sky patches indexed by SkyPatchID.
    // This will be reused in the neighbor distance calculation

    // Loop through the first SkyPatch
    for (int spid0=0; spid0<N_sp; spid0++)
    {
        // The first SkyPatch from the table
        SkyPatch sp0 = spt[spid0];
        // Direction of the first sky patch
        sp0.xyz(u0);       
        // The base index into the spn and spnd arrays for this row (neighbors of sp0)
        int idx = spid0*9;
        // Loop through the 9 neighbors of sp0
        for (int j=0; j<9; j++)
        {
            // The second SkyPatch ID on the table
            int spid1 = spn[idx];
            // Only build a SkyPatch instance for real entries; -1 signifies holes due to a corner
            if (spid1 >= 0)
            {
                // This is a real neighbor entry; get the SkyPatch instance from the table
                SkyPatch sp1 = spt[spid1];
                // Direction of the second sky patch
                sp1.xyz(u1);
                // Distance between the directions u0 and u1
                spnd[idx] = norm(u0, u1);
            }   // if real sky patch
            // Write the maximum distance of 2.0 for the holes
            else {spnd[idx] = 2.0;}
            // Advance the index into spn and spnd
            idx++;
        }   // for over j
    } // for over spid0
} // function
 
// *****************************************************************************
const int32_t* SkyPatchNeighbor::operator[](int32_t spid) const
    {return (spn + 9*spid);}

// *****************************************************************************
const double* SkyPatchNeighbor::neighbor_distance(int32_t spid) const
    {return (spnd + 9*spid);}
