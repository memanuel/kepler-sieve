// *****************************************************************************
// Included libraries
#include <cmath>
#include <iostream>
#include <boost/format.hpp>

// Local dependencies
#include "utils.h"
#include "sky_patch.h"

// *****************************************************************************
// Names used
using std::cout;
using boost::format;
using ks::CubeFace;
using ks::SkyPatch;
using ks::SkyPatch_from_id;
using ks::sqr;
using ks::print_stars;
using ks::print_newline;
using ks::report_test;

// The grid size
constexpr int N = ks::N_sky_patch;
constexpr int M = 2*N;
constexpr int M2 = M*M;

// *****************************************************************************
int main()
{
    // Initialize a CubeFace
    CubeFace cf = CubeFace(0);
    // Report results
    print_stars();
    cout << "CubeFace:\n";
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

    // Set SkyPatch integer values for test
    int8_t f_ = 0;
    int16_t i = 512;
    int16_t j = 1024+512;
    int32_t sky_patch_id = (M2*f_ + M*i + j);
    // Initialize a SkyPatch
    SkyPatch sp = SkyPatch(f_, i, j);
    // Report results
    print_newline();
    print_stars();
    cout << "SkyPatch:\n";
    cout << sp.str();
    // Calculate sky patch integer attributes
    cout << "Integer attributes\n";
    cout << format("f:   %d\n") % static_cast<int>(sp.f.id);
    cout << format("i:   %d\n") % static_cast<int>(sp.i);
    cout << format("j:   %d\n") % static_cast<int>(sp.j);
    cout << format("id   %d\n") % sp.id();
    // Test that recovered ID matches
    bool is_ok_id = (sky_patch_id == sp.id());
    string test_name = (format("SkyPatch: recovered ID matches calculated ID (%s)") % sky_patch_id).str();
    report_test(test_name, is_ok_id);
    
    // Calculate local sky patch coordinates (u, v, w) on unit sphere
    double u = sp.u();
    double v = sp.v();
    double w = sp.w();
    // Display (u, v, w)
    cout << "\nCoordinates of midpoint on unit sphere\n";
    cout << format("u:   %+8.6f\n") % u;
    cout << format("v:   %+8.6f\n") % v;
    cout << format("w:   %+8.6f\n") % w;

    // Check that (u, v, w) is really on the sphere
    double r_sph = sqrt(sqr(u) + sqr(v) + sqr(w));
    cout << "\nRecovered radius of point (u, v, w) on sphere:\n";
    cout << format("r    : %8.6f\n") % r_sph;
    double err = fabs(1.0 - r_sph);
    cout << format("Error: %8.6e\n") % err;
    bool is_ok_sphere = (err < 1.0E-15);
    report_test("SkyPatch: point(u, v, w) on unit sphere", is_ok_sphere);

    // Create a second SkyPatch from the integer ID and check that it's the same as the first
    SkyPatch sp2 = SkyPatch_from_id(sky_patch_id);
    cout << format("\nSkyPatch built using this ID (%d):\n") % sky_patch_id;
    cout << format("(%+8.6f, %+8.6f, %+8.6f)\n") % sp2.u() % sp2.v() % sp2.w();
    cout << format("id:   %d\n") % sp2.id();
    bool is_ok_from_id = (sp2.id() == sky_patch_id);
    report_test("SkyPatch: instance built from ID matches input sky_patch_id", is_ok_from_id);

    // Combined test results
    bool is_ok = (is_ok_id && is_ok_sphere && is_ok_from_id);
    report_test("\nSkyPatch: overall test results", is_ok);

    // Normal program exit
    return 0;
}
