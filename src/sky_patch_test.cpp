// *****************************************************************************
// Included libraries
#include <iostream>
#include <cmath>
#include <boost/format.hpp>

// Local dependencies
#include "sky_patch.h"

// *****************************************************************************
// Names used
using std::cout;
using boost::format;
using ks::CubeFace;
using ks::SkyPatch;
using ks::sqr;
constexpr int N = ks::N_sky_patch;

// *****************************************************************************
int main()
{
    // Initialize a CubeFace
    CubeFace cf = CubeFace(0);
    cout << "\nCubeFace:\n";
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
    cout << format("c    :%+3.1f\n") % cf.c();

    // Initialize a SkyPatch
    SkyPatch sp = SkyPatch(0, N/2-1, N/2+1);
    cout << "\nSkyPatch:\n";
    // Calculate sky patch integer attributes
    cout << format("f:   %d\n") % static_cast<int>(sp.f.id);
    cout << format("i:   %d\n") % static_cast<int>(sp.i);
    cout << format("j:   %d\n") % static_cast<int>(sp.j);
    // Calculate sky patch coordinates
    double u = sp.u();
    double v = sp.v();
    double w = sp.w();
    // Display local coordinates on unit sphere
    cout << format("u:   %+8.6f\n") % u;
    cout << format("v:   %+8.6f\n") % v;
    cout << format("w:   %+8.6f\n") % w;
    // Check that point really on the sphere
    double r_sph = sqrt(sqr(u) + sqr(v) + sqr(w));
    cout << "Recovered radius of point (u, v, w) on sphere:\n";
    cout << format("r    : %8.6f\n") % r_sph;
    double err = fabs(1.0 - r_sph);
    cout << format("Error: %8.6e\n") % err;

    // Normal program exit
    return 0;
}
