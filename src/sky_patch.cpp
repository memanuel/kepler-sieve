// *****************************************************************************
// Included files
#include <stdexcept>
#include "sky_patch.h"

// *****************************************************************************
using std::range_error;
using ks::CubeFace;

/** Initialize a CubeFace from its ID; starts from 0.*/
CubeFace::CubeFace(int8_t id_) : 
id(id_) {}

/** Default destructor for CubeFace.*/
CubeFace::~CubeFace() {}

/** Two character description, e.g. "Z+"; see DB table KS.CubeFace. */
const string CubeFace::code()
{
    switch (id)
    {
        case 0:     return string("Z+");
        case 1:     return string("Y+");
        case 2:     return string("X+");
        case 3:     return string("X-");
        case 4:     return string("Y-");
        case 5:     return string("Z-");
        default:
            throw range_error("cube face id must be between 0 and 5, inclusive.\n");
    }
}
