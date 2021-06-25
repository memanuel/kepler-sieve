/*****************************************************************************
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/

// *****************************************************************************
// Included files
#include "sky_patch.h"

/******************************************************************************
Functions for working with SkyPatch objects
******************************************************************************/
namespace ks {

// *****************************************************************************
int fij2spid(int8_t f, int16_t i, int16_t j)
{
    return (M2*f) + (M*i) + j;
}

// *****************************************************************************
/** Initialize a SkyPatch from its integer ID.*/
SkyPatch SkyPatch_from_id(int32_t id)
{
    // First integer division; unpack id into cube face f and remainder x
    div_t qr = div(static_cast<int>(id), M2);
    int8_t f_ = static_cast<int8_t>(qr.quot);
    int32_t x = qr.rem;

    // Second integer division; unpack x into grid points i and j    
    qr = div(x, M);
    int16_t i_ = static_cast<int16_t>(qr.quot);
    int16_t j_ = static_cast<int16_t>(qr.rem);

    // Instantiate the sky patch
    return SkyPatch(f_, i_, j_);
}

// *****************************************************************************
int sky_patch_count()
{
    return static_cast<int>(N_sp);
}

// *****************************************************************************
spn_type make_sky_patch_neighbor_table()
{
    // Allocate an array of size N_spn (number of sky patch neighbors) to hold the 9 neighbors of each patch
    spn_type spn = new int32_t[N_spn];
    // Loop through the starting sky patch, with ID sky_patch_id_1
    for (int8_t f=0; f<6; f++)
    {
        // The offset for this grid face is M2*f
        int32_t idx_f = M2*f;
        for (int16_t i0=0; i0<M; i0++)
        {
            for (int16_t j0=0; j0<M; j0++)
            {
                // The starting SkyPatch and its ID
                SkyPatch sp = SkyPatch(f, i0, j0);
                int32_t spid0 = sp.id();

                // The starting index into the neighbors table
                size_t idx = spid0*9;

                // Get grid coordinates of 9 candidate neighbors.
                // First iterate over the shift in the i index, di               
                for (int16_t di=-1; di<=1; di++)
                {
                    // The grid entry i1 for the three candidate neighbors in this row
                    int16_t i1 = i0+di;
                    // The additive term to the index into the sky patch table for this row
                    int32_t idx_i1 = M*i1;
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
                            if (!is_corner)
                            {
                                spid1 = sp.shift(di, dj).id();
                            }
                        }
                        // Write the new SkyPatchID to the spn array
                        spn[idx++] = spid1;
                    } // for over dj
                } // for over di
            } // for over j0
        } // for over i0
    } // for over f
    // Return the populated array
    return spn;
}

// *****************************************************************************
spnd_type make_sky_patch_neighbor_dist_table(const spn_type spn)
{
    // Allocate an array of size N_spn (number of sky patch neighbors) to hold the distance
    // to the 9 neighbors of each patch.
    spnd_type spnd = new double[N_spn];
    // Direction of the first and second skypatch; reused in inner loop.
    double u0[3];
    double u1[3];

    // Loop through the first SkyPatch
    for (int spid0=0; spid0<N_sp; spid0++)
    {
        // The first SkyPatch
        SkyPatch sp0 = SkyPatch_from_id(spid0);
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
                // This is a real neighbor entry; build the SkyPatch instance
                SkyPatch sp1 = SkyPatch_from_id(spid1);
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
    // Return the populated array
    return spnd;
}
// *****************************************************************************
}; // namespace
