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
    div_t qr = div(id, M2);
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
    return static_cast<int>(N_spc);
}

// *****************************************************************************
void write_sky_patch_neighbor_table(int32_t* spn)
{
    // Loop through the starting sky patch, with ID sky_patch_id_1
    for (int8_t f=0; f<6; f++)
    {
        // The offset for this grid face is M2*f
        int32_t idx_f = M2*f;
        for (int16_t i0=0; i0<M; i0++)
        // for (int16_t i0=0; i0<10; i0++)
        {
            for (int16_t j0=0; j0<M; j0++)
            {
                // The starting SkyPatchID
                int32_t spid0 = idx_f + (M*i0) + j0;
                // The starting index
                size_t idx = spid0*9;
                
                // Get grid coordinates of 9 candidate neighbors.
                for (int16_t di=-1; di<=1; di++)
                {
                    // The grid entry i1 for the three candidate neighbors in this row
                    int16_t i1 = i0+di;
                    int32_t idx_i1 = M*i1;
                    for (int16_t dj=-1; dj<=1; dj++)
                    {
                        // The grid entry j1 for this candidate neighbor
                        int16_t j1 = j0+dj;
                        // Write the candidate neighbor to the array if it is on the same face.
                        // Don't worry about wrapping around edges and corners.
                        // Write the dummy value -1 if the neighbor is not a real grid point
                        bool is_on_grid = (0 <= i1) && (i1 < M) && (0 <= j1) && (j1 < M);
                        spn[idx++] = is_on_grid ? (idx_f+idx_i1+j1) : -1;
                    }
                }
            }
        }
    }
}

// *****************************************************************************
}; // namespace
