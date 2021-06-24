#pragma once

// *****************************************************************************
// Included libraries
#include <cstdint>
#include <string>
#include <stdexcept>

// *****************************************************************************
// Standard library and boost class names used
using std::string;
using std::range_error;

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

        // String code and description
        //*String representation of this CubeFace, e.g. 'Z+. See DB table KS.CubeFace'
        const string str() const;   
        //*Long string description
        const string description() const;   
        // Integer index of the three axes on this cube face

        //*Index of axis indexed by i with local values u, e.g. 1 on 'Z+'
        const int k1() const;
        //*Index of axis indexed by j with local values v, e.g. 2 on 'Z+'
        const int k2() const;
        //*Index of largest axis, with local values w, e.g. 3 on 'Z+'
        const int k3() const;    

        //*The index corresponding to global 'X' axis
        const int index_x() const;
        //*The index corresponding to global 'Y' axis
        const int index_y() const;
        //*The index corresponding to global 'Z' axis
        const int index_z() const;

        // One letter labels for the three axes on this cube face
        //*Label of axis indexed by i with local values u, e.g. 'X' on 'Z+'
        const char alpha() const;
        //*Label of axis indexed by j with local values v, e.g. 'Y' on 'Z+'
        const char beta() const;
        //*Label of largest axis, with local values w, e.g. 'Z' on 'Z+'
        const char gamma() const;
        //*The shared value of c on circumscribed cube coordinates (a, b, c)
        const double c() const;

        // Neighbors
        //*Neighboring face when wrapping through i=0 (u=-1.0).
        const CubeFace neighbor_i0() const;
        //*Neighboring face when wrapping through i=M (u=+1.0).
        const CubeFace neighbor_i1() const;
        //*Neighboring face when wrapping through j=0 (v=-1.0).
        const CubeFace neighbor_j0() const;
        //*Neighboring face when wrapping through j=M (v=+1.0).
        const CubeFace neighbor_j1() const;     
        //*The opposite face
        const CubeFace opposite() const;
};

// *****************************************************************************
}; // namespace
