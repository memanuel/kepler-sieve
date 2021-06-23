// *****************************************************************************
// Included files
#include <string>

// *****************************************************************************
// Standard library and boost class names used
using std::string;

// *****************************************************************************
// Classes defined in sky_patch.cpp

/** One face of a cube; used in SkyPatch. */
class CubeFace
{
    private:
        int8_t id;
    public:
        CubeFace(int8_t id);
        const string code();
        const int8_t i();
        const int8_t j1();
        const int8_t j2();
        const char alpha();
        const char beta();
        const char gamma();
};
