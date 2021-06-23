#pragma once

// *****************************************************************************
// Included files
#include <string>
#include <stdexcept>

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
        const string code();
        const int8_t i();
        const int8_t j1();
        const int8_t j2();
        const char alpha();
        const char beta();
        const char gamma();
        const double c() const;
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
        const int32_t id();
        // Local coordinates of the center of this sky patch
        const double u();
        const double v();
        const double w();
        // Global coordinates of the center of this sky patch
        const double x();
        const double y();
        const double z();
    private:
        // Distance to circumscribing cube
        const double r;
        // Local coordinates of the center of this sky patch, projected onto circumscribing cube
        const double a();
        const double b();
        const double c();
};

/** Construct a SkyPatch from its integer ID, sky_patch_id.*/
SkyPatch SkyPatch_from_id(int32_t sky_patch_id);

// *****************************************************************************
}; // namespace
