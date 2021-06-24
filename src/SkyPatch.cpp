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
//*Error condition for bad integer input for a CubeFace
range_error err_sky_patch_f = range_error("SkyPatch cube face f must be between 0 and 6, inclusive.\n");
//*Error condition for bad grid i coordinate
range_error err_sky_patch_grid_i = range_error("SkyPatch grid coordinate i must be between 0 and M=2048, exclusive.\n");
//*Error condition for bad grid j coordinate
range_error err_sky_patch_grid_j = range_error("SkyPatch grid coordinate j must be between 0 and M=2048, exclusive.\n");
}

// *****************************************************************************
// Various error conditions
/******************************************************************************
Implementation of the SkyPatch class
******************************************************************************/

// *****************************************************************************
//*Helper function used to initialize SkyPatch; calculate distance r to point on cube.
double r_initializer(int16_t i, int16_t j)
{
    // The u coordinate is (i+0.5)/N - 1.0; calculate its square
    double u2 = sqr(N_inv*(i+0.5) - 1.0);
    // The v coordinate is (j+0.5)/N - 1.0; calculate its square
    double v2 = sqr(N_inv*(j+0.5) - 1.0);
    // Calculate the distance r using u2 and v2 from above; w2 is always 1.0
    return sqrt(1.0 + u2 + v2);
}

// *****************************************************************************
SkyPatch::SkyPatch(CubeFace f_, int16_t i_, int16_t j_) : 
f(f_),
i(i_),
j(j_),
// The distance r to the point on the unit cube (where unit sphere is inscribed)
// r(sqrt(1.0 + sqr(N_inv*(i+0.5)-1.0) + sqr(N_inv*(j+0.5)-1.0)))
r(r_initializer(i, j))
{
    // Check that grid points were valid
    if ((i_ < 0) || (i_ > M)) 
    {
        int fi = static_cast<int>(f_.id);
        cout << format("Bad inputs to SkyPatch. (f, i, j) = (%d, %d, %d).\n") % fi % i_ % j_;
        throw err_sky_patch_grid_i;
    }
    if ((j_ < 0) || (j_ > M)) 
    {
        int fi = static_cast<int>(f_.id);
        cout << format("Bad inputs to SkyPatch. (f, i, j) = (%d, %d, %d).\n") % fi % i_ % j_;
        throw err_sky_patch_grid_j;
    }
}

// *****************************************************************************
// Delegate to constructor taking a CubeFace instance by building f from its integer ID
SkyPatch::SkyPatch(int8_t f_, int16_t i_, int16_t j_) : 
SkyPatch(CubeFace(f_), i_, j_) {};

// *****************************************************************************
SkyPatch::~SkyPatch() {}

// *****************************************************************************
const int32_t SkyPatch::id() const
{
    return (M2*f.id) + (M*i) + j;
}

// *****************************************************************************
const double SkyPatch::a() const
{
    return N_inv*(static_cast<double>(i)+0.5) - 1.0;
}

// *****************************************************************************
const double SkyPatch::b() const
{
    return N_inv*(static_cast<double>(j)+0.5) - 1.0;
}

// *****************************************************************************
const double SkyPatch::c() const
{
    return f.c();
}

// *****************************************************************************
const double SkyPatch::u() const
{
    return a() / r;
}

// *****************************************************************************
const double SkyPatch::v() const
{
    return b() / r;
}

// *****************************************************************************
const double SkyPatch::w() const
{
    return c() / r;
}

// *****************************************************************************
const double SkyPatch::x() const
{
   switch (f.index_x())
   {
        case 1:     return u();
        case 2:     return v();
        case 3:     return w();
        default:    throw err_sky_patch_f;
   }
}

// *****************************************************************************
/**Global y coordinate.*/
const double SkyPatch::y() const
{
   switch (f.index_y())
   {
        case 1:     return u();
        case 2:     return v();
        case 3:     return w();
        default:    throw err_sky_patch_f;
   }
}

// *****************************************************************************
/**Global z coordinate.*/
const double SkyPatch::z() const
{
   switch (f.index_z())
   {
        case 1:     return u();
        case 2:     return v();
        case 3:     return w();
        default:    throw err_sky_patch_f;
   }
}

// *****************************************************************************
const SkyPatch SkyPatch::shift_i(const int16_t di) const
{
    // Calculate candidate new i value
    int16_t i_ = i + di;

    // Check whether i has wrapped around either the low or high end of the cube face
    bool is_lo = (i_ < 0);
    bool is_hi = (i_ >= M);

    // Most common case by far is that i_ is on the same grid face
    if ((!is_lo) && (!is_hi))
    {
        return SkyPatch(f, i_, j);
    }

    // If we get here, we've wrapped around.
    // The sign of the cube face as an integer; e.g. Z+ has ci=1
    int ci = static_cast<int> (f.c());
    // Calculate the distance that i has gone past the end; this is a small positive integer
    int16_t wrap = is_lo ? -i_ : (i_-M);
    // Calculate the new grid index on what used to be the largest axis (gamma)
    int16_t k = N + (N-wrap)*ci;

    // Initialize j_ to an invalid value; both i_ and j_ will be populated correctly below
    int16_t j_ = -1;

    // Get the neighboring CubeFace
    CubeFace f_ = is_lo ? f.neighbor_i0() : f.neighbor_i1();

    // Assign the grid index of new i axis
    if (f_.k1() == f.k1()) {i_ = i;}    // New i is old i; shared axis
    if (f_.k1() == f.k2()) {i_ = j;}    // New i is old j; shared axis
    if (f_.k1() == f.k3()) {i_ = k;}    // New i is wrapped axis

    // Assign the grid index of new j axis
    if (f_.k2() == f.k1()) {j_ = i;}    // New j is old i; shared axis
    if (f_.k2() == f.k2()) {j_ = j;}    // New j is old j; shared axis
    if (f_.k2() == f.k3()) {j_ = k;}    // New j is wrapped axis

    // Now return new SkyPatch
    return SkyPatch(f_, i_, j_);
}

// *****************************************************************************
const SkyPatch SkyPatch::shift_j(const int16_t dj) const
{
    // Calculate candidate new j value
    int16_t j_ = j + dj;

    // Check whether j has wrapped around either the low or high end of the cube face
    bool is_lo = (j_ < 0);
    bool is_hi = (j_ >= M);

    // Most common case by far is that j_ is on the same grid face
    if ((!is_lo) && (!is_hi))
    {
        return SkyPatch(f, i, j_);
    }

    // If we get here, we've wrapped around.
    // The sign of the cube face as an integer; e.g. Z+ has ci=1
    int ci = static_cast<int> (f.c());
    // Calculate the distance that i has gone past the end; this is a small positive integer
    int16_t wrap = is_lo ? -j_ : (j_-M);
    // Calculate the new grid index on what used to be the largest axis (gamma)
    int16_t k = N + (N-wrap)*ci;

    // Initialize i_ to an invalid value; both i_ and j_ will be populated correctly below
    int16_t i_ = -1;

    // Get the neighboring CubeFace
    CubeFace f_ = is_lo ? f.neighbor_j0() : f.neighbor_j1();

    // Assign the grid index of new i axis
    if (f_.k1() == f.k1()) {i_ = i;}    // New i is old i; shared axis
    if (f_.k1() == f.k2()) {i_ = j;}    // New i is old j; shared axis
    if (f_.k1() == f.k3()) {i_ = k;}    // New i is wrapped axis

    // Assign the grid index of new j axis
    if (f_.k2() == f.k1()) {j_ = i;}    // New j is old i; shared axis
    if (f_.k2() == f.k2()) {j_ = j;}    // New j is old j; shared axis
    if (f_.k2() == f.k3()) {j_ = k;}    // New j is wrapped axis

    // Now return new SkyPatch
    return SkyPatch(f_, i_, j_);
}

// *****************************************************************************
const SkyPatch SkyPatch::shift(const int16_t di, const int16_t dj) const
{
    // Calculate candidate new i and j values
    int16_t i_ = i + di;
    int16_t j_ = j + dj;

    // Check whether i and j have wrapped around
    bool is_wrap_i = (i_ <0) || (i_ >M);
    bool is_wrap_j = (j_ <0) || (j_ >M);
    bool is_wrap = is_wrap_i || is_wrap_j;

    // Most common case: no wrap.  Then stay on the same CubeFace
    if (!is_wrap)
    {
        return SkyPatch(f, i_, j_);
    }

    // If we did not wrap in the j direction, delegate to shift_i
    if (!is_wrap_j)
    {
        return shift_i(di);
    }

    // If we did not wrap in the i direction, delegate to shift_j
    if (!is_wrap_i)
    {
        return shift_j(dj);
    }

    // If we get here, we are wrapping around twice.
    // Map these to the opposite face to make sure we have a valid SkyPatch, but one that is very far away.
    // This way it will be filtered out later when setting a maximum distance.
    CubeFace f_ = f.opposite();
    return SkyPatch(f_, i, j);
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
    string xyz = (format("(x, y, z) = (%+8.6f, %+8.6f, %+8.6f)") % x() % y() % z()).str();
    // Return one combined description
    return (format("spid = %8d. %s. %s.\n") % id() % fij % xyz).str();
}
