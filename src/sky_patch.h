#pragma once

// *****************************************************************************
// Included files
#include <string>

// *****************************************************************************
// Standard library and boost class names used
using std::string;

// *****************************************************************************
// Put classes into the namespace ks (for Kepler Sieve)
namespace ks {

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
};


// *****************************************************************************
}; // namespace
