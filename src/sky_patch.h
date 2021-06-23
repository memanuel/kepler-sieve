#pragma once

// *****************************************************************************
// Included libraries
#include <string>
#include <stdexcept>
#include <boost/format.hpp>

// Local dependencies
#include "utils.h"

// *****************************************************************************
// Standard library and boost class names used
using std::string;
using std::range_error;

// *****************************************************************************
// Put classes into the namespace ks (for Kepler Sieve)
namespace ks {

// *****************************************************************************
// Set the grid size for the sky patch at compile time
constexpr int N_sky_patch = 1024;

// *****************************************************************************
/** One face of a cube; used in SkyPatch. */
class CubeFace
{
    public:
        // Constructor and destructor
        CubeFace(int8_t id_);
        ~CubeFace();
        // Data
        const int8_t id;
        // Access attributes
        const string code() const;
        const int8_t i() const;
        const int8_t j1() const;
        const int8_t j2() const;
        const char alpha() const;
        const char beta() const;
        const char gamma() const;
        const double c() const;
        // Neighbors
        const CubeFace neighbor_i0() const;
        const CubeFace neighbor_i1() const;
        const CubeFace neighbor_j0() const;
        const CubeFace neighbor_j1() const;
};


// *****************************************************************************
/** One cell on the surface of the unit sphere, used for spatial indexing. */
class SkyPatch
{
    public:
        // Constructor and destructor
        SkyPatch(int8_t f_, int16_t i_, int16_t j_);
        ~SkyPatch();
        // Data
        const CubeFace f;
        const int16_t i;
        const int16_t j;
        // The integer ID of this sky patch
        const int32_t id() const;
        // Local coordinates of the center of this sky patch
        const double u() const;
        const double v() const;
        const double w() const;
        // Global coordinates of the center of this sky patch
        const double x() const;
        const double y() const;
        const double z() const;
        // Output string description
        const string str() const;
    private:
        // Distance to circumscribing cube
        const double r;
        // Local coordinates of the center of this sky patch, projected onto circumscribing cube
        const double a() const;
        const double b() const;
        const double c() const;
};

// *****************************************************************************
/** Construct a SkyPatch from its integer ID, sky_patch_id.*/
SkyPatch SkyPatch_from_id(int32_t sky_patch_id);

/** Calculate a SkyPatchID from the CubeFaceID f and grid coordinates (i, j).*/
int fij2spid(int8_t f, int16_t i, int16_t j);

/** The number of SkyPatches for the selected grid size.*/
int sky_patch_count();

/** Construct a table of sky patch neighbors.*/
void write_sky_patch_neighbor_table(int32_t* spn);

// *****************************************************************************
}; // namespace
