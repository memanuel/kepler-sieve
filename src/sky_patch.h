#pragma once

// *****************************************************************************
// Included libraries
#include <string>
#include <stdexcept>
#include <boost/format.hpp>

// Local dependencies
#include "utils.h"
#include "CubeFace.h"
#include "SkyPatch.h"

// *****************************************************************************
// Standard library and boost class names used
using std::string;
using std::range_error;

// Local names used
using ks::sky_patch::M;
using ks::sky_patch::M2;
using ks::sky_patch::N_spc;

// *****************************************************************************
// Put classes into the namespace ks (for Kepler Sieve)
namespace ks {

// *****************************************************************************
//* Construct a SkyPatch from its integer ID, sky_patch_id.
SkyPatch SkyPatch_from_id(int32_t sky_patch_id);

//* Calculate a SkyPatchID from the CubeFaceID f and grid coordinates (i, j).
int fij2spid(int8_t f, int16_t i, int16_t j);

/** The number of SkyPatches for the selected grid size.*/
int sky_patch_count();

/** Construct a table of sky patch neighbors.*/
void write_sky_patch_neighbor_table(int32_t* spn);

// *****************************************************************************
}; // namespace
