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
    // String attributes
    cout << format("code : %s\n") % cf.code();
    cout << format("alpha: %c\n") % cf.alpha();
    cout << format("beta : %c\n") % cf.beta();
    cout << format("gamma: %c\n") % cf.gamma();
    // Integer attributes
    cout << format("id   : %d\n") % static_cast<int>(cf.id);
    cout << format("i    : %d\n") % static_cast<int>(cf.i());
    cout << format("j1   : %d\n") % static_cast<int>(cf.j1());
    cout << format("j2   : %d\n") % static_cast<int>(cf.j2());
    cout << format("ci   :%+3.1f\n") % cf.ci();

    // Normal program exit
    return 0;
}
