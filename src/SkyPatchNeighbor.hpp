/** @file   SkyPatchNeighbor.hpp
 *  @brief  Class to calculate and cache neighbors of each SkyPatch.
 *   
 *  @author Michael S. Emanuel
 *  @date   2021-06-28
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies

// Local dependencies
#include "SkyPatch.hpp"
    using ks::sky_patch::M;
    using ks::sky_patch::M2;
    using ks::sky_patch::N_sp;
    using ks::SkyPatch_from_id;

// *****************************************************************************
namespace ks {

// *****************************************************************************
namespace sky_patch{
//*The number of SkyPatch neighbors: each cell has 9 neighbors except for the 8 corners.
constexpr int N_spn = N_sp*9;
} // namespace ks::sky_patch

// *****************************************************************************
class SkyPatchNeighbor
{
public:
    /// Initialize a SkyPatchNeighbor table
    SkyPatchNeighbor();
    /// Destructor for SkyPatchNeighbor.
    ~SkyPatchNeighbor();

    /// Build the neighbor distance table
    void build_neighbor_distance();

    /// Get the neighbors of a skypatch given its integer ID
    const int32_t* operator[](int32_t spid) const;

    /// Get the neighbor distances of a skypatch given its integer ID
    const double* neighbor_distance(int32_t spid) const;

private:
    // Data
    /**Array of N_spn SkyPatchIDs; 9 entries indexed[spid*9, (spid+1)*9) 
       are the 9 neighbors of the sky patch with ID spid */
    int32_t* spn;
    /// Array of N_spn distances from spid to its neighbors
    double* spnd;
};

// *****************************************************************************
} // namespace ks
