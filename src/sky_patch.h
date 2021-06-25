/*****************************************************************************
 * Functions for working with SkyPatch and related classes.
 * Build a SkyPatch from an integer ID.
 * Assemble a tables of sky patch neighbors and their distances.
 * 
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/
#pragma once

// *****************************************************************************
// Included libraries
#include <string>
#include <algorithm>
#include <stdexcept>
#include <boost/format.hpp>

// Local dependencies
#include "CubeFace.h"
#include "SkyPatch.h"
#include "utils.h"
#include "astro_utils.h"

// *****************************************************************************
// Standard library and boost class names used
using std::string;
using std::range_error;

// Local names used
using ks::sky_patch::M;
using ks::sky_patch::M2;
using ks::sky_patch::N_sp;
using ks::norm;

// *****************************************************************************
// Put classes into the namespace ks (for Kepler Sieve)
namespace ks {

// *****************************************************************************
//* The number of SkyPatch neighbors: each cell has 9 neighbors except for the 8 corners.
namespace sky_patch {
constexpr int N_spn = N_sp*9;
}
using ks::sky_patch::N_spn;

// *****************************************************************************
//* Construct a SkyPatch from its integer ID, sky_patch_id.
SkyPatch SkyPatch_from_id(int32_t sky_patch_id);

//* Calculate a SkyPatchID from the CubeFaceID f and grid coordinates (i, j).
int fij2spid(int8_t f, int16_t i, int16_t j);

//* The number of SkyPatches for the selected grid size. Alias to access ks::sky_patch::N_sp.
int sky_patch_count();

using spn_type = int32_t*;
//* Construct a table of sky patch neighbors.
spn_type make_sky_patch_neighbor_table();

using spnd_type = double*;
//* Construct a table of sky patch neighbor distances.*/
spnd_type make_sky_patch_neighbor_dist_table(const spn_type spn);

// *****************************************************************************
}; // namespace
