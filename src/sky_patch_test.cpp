// *****************************************************************************
// Included libraries
#include <iostream>
#include <boost/format.hpp>

// Local dependencies
#include "sky_patch.h"

// *****************************************************************************
// Names used
using std::cout;
using boost::format;
using ks::CubeFace;

int main()
{
    // Initialize a CubeFace
    CubeFace cf = CubeFace(0);
    // Get some attributes of the cube face
    cout << ("CubeFace:\n");
    cout << format("id:   %d\n") % static_cast<int>(cf.id);
    cout << format("code: %s\n") % cf.code();

    // Normal program exit
    return 0;
}
