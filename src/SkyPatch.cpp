// *****************************************************************************
// Included files
#include "SkyPatch.h"

// *****************************************************************************
// Local names used
using ks::SkyPatch;
using ks::sky_patch::N;
using ks::sky_patch::M;
using ks::sky_patch::M2;

// *****************************************************************************
// Precompute 1.0/N for speed
constexpr double N_inv = 1.0 / static_cast<double>(N);

// *****************************************************************************
// Instantiate exceptions for a bad cube face f or grid entry (i or j)
namespace ks {
range_error err_sky_patch_f = range_error("SkyPatch cube face f must be between 0 and 6, inclusive.\n");
range_error err_sky_patch_grid_i = range_error("SkyPatch grid coordinate i must be between 0 and M=2048, exclusive.\n");
range_error err_sky_patch_grid_j = range_error("SkyPatch grid coordinate j must be between 0 and M=2048, exclusive.\n");
}
using ks::err_sky_patch_f, ks::err_sky_patch_grid_i, ks::err_sky_patch_grid_j;

/******************************************************************************
Implementation of the SkyPatch class
******************************************************************************/

// *****************************************************************************
SkyPatch::SkyPatch(int8_t f_, int16_t i_, int16_t j_) : 
f(CubeFace(f_)),
i(i_),
j(j_),
// The distance r to the point on the unit cube (where unit sphere is inscribed)
r(sqrt(1.0 + sqr(N_inv*(i+0.5)-1.0) + sqr(N_inv*(j+0.5)-1.0)))
{
    if ((f_ < 0) || (f_ > 5)) {throw err_sky_patch_f;}
    if ((i_ < 0) || (i_ > M)) {throw err_sky_patch_grid_i;}
    if ((j_ < 0) || (j_ > M)) {throw err_sky_patch_grid_j;}
}

// *****************************************************************************
SkyPatch::~SkyPatch() {}

// *****************************************************************************
/**Integer ID of this sky patch*/
const int32_t SkyPatch::id() const
{
    return (M2*f.id) + (M*i) + j;
}

// *****************************************************************************
/**Coordinate on unit cube of first varying axis (alpha).*/
const double SkyPatch::a() const
{
    return N_inv*(static_cast<double>(i)+0.5) - 1.0;
}

// *****************************************************************************
/**Coordinate on unit cube of second varying axis (beta).*/
const double SkyPatch::b() const
{
    return N_inv*(static_cast<double>(j)+0.5) - 1.0;
}

// *****************************************************************************
/**Coordinate on unit cube of third varying axis (gamma).*/
const double SkyPatch::c() const
{
    return f.c();
}

// *****************************************************************************
/**Coordinate on unit sphere of first varying axis (alpha).*/
const double SkyPatch::u() const
{
    return a() / r;
}

// *****************************************************************************
/**Coordinate on unit sphere of second varying axis (beta).*/
const double SkyPatch::v() const
{
    return b() / r;
}

// *****************************************************************************
/**Coordinate on unit sphere of third varying axis (gamma).*/
const double SkyPatch::w() const
{
    return c() / r;
}

// *****************************************************************************
/**Global x coordinate.*/
const double SkyPatch::x() const
{
    if (f.alpha()=='X') {return u();}
    if (f.beta()=='X') {return v();}
    if (f.gamma()=='X') {return w();}
    throw std::runtime_error("Bad cube face.  None of alpha, beta, gamma equal 'X'.");
}

// *****************************************************************************
/**Global y coordinate.*/
const double SkyPatch::y() const
{
    if (f.alpha()=='Y') {return u();}
    if (f.beta()=='Y') {return v();}
    if (f.gamma()=='Y') {return w();}
    throw std::runtime_error("Bad cube face.  None of alpha, beta, gamma equal 'Y'.");
}

// *****************************************************************************
/**Global z coordinate.*/
const double SkyPatch::z() const
{
    if (f.alpha()=='Z') {return u();}
    if (f.beta()=='Z') {return v();}
    if (f.gamma()=='Z') {return w();}
    throw std::runtime_error("Bad cube face.  None of alpha, beta, gamma equal 'Z'.");
}

// *****************************************************************************
const string SkyPatch::str() const
{
    // The Integer coordinates (f, i, j)
    int f_ = static_cast<int>(f.id);
    string fij = (format("(f, i, j) = (%d, %4d, %4d)") % f_ % i % j).str();
    // The midpoint local coordinates (u, v, w)
    // string uvw = (format("(u, v, w) = (%8.6f, %8.6f, %8.6f)") % u() % v() % w()).str();
    // The midpoint global coordinates (x, y, z)
    string xyz = (format("(x, y, z) = (%8.6f, %8.6f, %8.6f)") % x() % y() % z()).str();
    // Return one combined description
    return (format("spid = %8d. %s. %s.\n") % id() % fij % xyz).str();
}
