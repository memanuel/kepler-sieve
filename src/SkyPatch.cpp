/** @file SkyPatch.cpp
 *  @brief Implmentation of SkyPatch class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

// *****************************************************************************
// Local dependencies
#include "SkyPatch.hpp"
    using ks::SkyPatch;
    using ks::sky_patch::N;
    using ks::sky_patch::M;
    using ks::sky_patch::M2;
    using ks::sky_patch::N_sp;

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
} // namespace ks

/******************************************************************************
Implementation of the SkyPatch class
******************************************************************************/

// *****************************************************************************
int ks::fij2spid(int8_t f, int16_t i, int16_t j)
{
    return (M2*f) + (M*i) + j;
}

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
        print("Bad inputs to SkyPatch. (f, i, j) = ({:d}, {:d}, {:d}).\n", fi, i_, j_);
        throw err_sky_patch_grid_i;
    }
    if ((j_ < 0) || (j_ > M)) 
    {
        int fi = static_cast<int>(f_.id);
        print("Bad inputs to SkyPatch. (f, i, j) = ({:d}, {:d}, {:d}).\n", fi % i_ % j_);
        throw err_sky_patch_grid_j;
    }
}

// *****************************************************************************
// Delegate to constructor taking a CubeFace instance by building f from its integer ID
SkyPatch::SkyPatch(int8_t f_, int16_t i_, int16_t j_) : 
    SkyPatch(CubeFace(f_), i_, j_) {}

// *****************************************************************************
//*Fast constructor when r precomputed and f, i, j already validated. Used in copying.
SkyPatch::SkyPatch(CubeFace f_, int16_t i_, int16_t j_, double r_) : 
    f(f_),
    i(i_),
    j(j_),
    r(r_) {}

// *****************************************************************************
//*Use same approach as fast constructor for copy assignment
SkyPatch::SkyPatch(const ks::SkyPatch& sp) :
    f(sp.f),
    i(sp.i),
    j(sp.j),
    r(sp.r) {}

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
void SkyPatch::xyz(double *u)
{
    u[0] = x();
    u[1] = y();
    u[2] = z();
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
    int16_t wrap = is_lo ? -i_ : (i_-M+1);
    // Calculate the new grid index on what used to be the largest axis (gamma)
    int16_t k = N + (N-wrap)*ci;

    // Initialize j_ to an invalid value; both i_ and j_ will be populated correctly below
    int16_t j_ = -1;

    // Get the neighboring CubeFace
    CubeFace f_ = is_lo ? f.neighbor_i0() : f.neighbor_i1();

    // Assign the grid index of new i axis
    // All three cases shown, but new i = old i is impossible; adjacent cube faces always have different i axes.
    // if (f_.k1() == f.k1()) {i_ = i;}     // New i is old i; shared axis.
    if (f_.k1() == f.k2()) {i_ = j;}        // New i is old j; shared axis
    if (f_.k1() == f.k3()) {i_ = k;}        // New i is wrapped axis

    // Assign the grid index of new j axis
    // All three cases shown, but new j = old j is impossible; adjacent cube faces always have different i axes.
    if (f_.k2() == f.k1()) {j_ = i;}        // New j is old i; shared axis
    // if (f_.k2() == f.k2()) {j_ = j;}     // New j is old j; shared axis
    if (f_.k2() == f.k3()) {j_ = k;}        // New j is wrapped axis

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
    // Calculate the distance that j has gone past the end; this is a small positive integer
    int16_t wrap = is_lo ? -j_ : (j_-M+1);
    // Calculate the new grid index on what used to be the largest axis (gamma)
    int16_t k = N + (N-wrap)*ci;

    // Initialize i_ to an invalid value; both i_ and j_ will be populated correctly below
    int16_t i_ = -1;

    // Get the neighboring CubeFace
    CubeFace f_ = is_lo ? f.neighbor_j0() : f.neighbor_j1();

    // Assign the grid index of new i axis
    // All three cases shown, but new i = old i is impossible; adjacent cube faces always have different i axes.
    // if (f_.k1() == f.k1()) {i_ = i;}     // New i is old i; shared axis
    if (f_.k1() == f.k2()) {i_ = j;}        // New i is old j; shared axis
    if (f_.k1() == f.k3()) {i_ = k;}        // New i is wrapped axis

    // Assign the grid index of new j axis
    // All three cases shown, but new j = old j is impossible; adjacent cube faces always have different i axes.
    if (f_.k2() == f.k1()) {j_ = i;}        // New j is old i; shared axis
    // if (f_.k2() == f.k2()) {j_ = j;}     // New j is old j; shared axis
    if (f_.k2() == f.k3()) {j_ = k;}        // New j is wrapped axis

    // Now return new SkyPatch
    return SkyPatch(f_, i_, j_);
}

// *****************************************************************************
const SkyPatch SkyPatch::shift(const int16_t di, const int16_t dj) const
{
    // Calculate candidate new i and j values
    int16_t i_ = i + di;
    int16_t j_ = j + dj;

    // Determine if we're wrapping in the i and j directions
    bool is_on_grid_i = (0 <= i_) && (i_ < M);
    bool is_on_grid_j = (0 <= j_) && (j_ < M);
    // Is this the simple case where we are on the same cube face?
    bool is_on_grid = is_on_grid_i && is_on_grid_j;

    // Most common case: no wrap.  Then stay on the same CubeFace
    if (is_on_grid)
    {
        return SkyPatch(f, i_, j_);
    }

    // If we did not wrap in the i direction, delegate to shift_i first, then shift_j
    if (is_on_grid_i)
    {
        SkyPatch sp = shift_i(di);
        return sp.shift_j(dj);
    }

    // If we did not wrap in the j direction, delegate to shift_j first, then shift_i
    if (is_on_grid_j)
    {
        SkyPatch sp = shift_j(dj);
        return sp.shift_i(di);
    }

    // If we get here, we are wrapping around twice.
    // Map these to the opposite face to make sure we have a valid SkyPatch, but one that is very far away.
    // This way it will be filtered out later when setting a maximum distance.
    CubeFace f_ = f.opposite();
    // Wrap i_ and j_ mod M so they are legal
    i_ = (M+i_) % M;
    j_ = (M+j_) % M;
    return SkyPatch(f_, i_, j_);
}

// *****************************************************************************
const string SkyPatch::str() const
{
    // The Integer coordinates (f, i, j)
    int f_ = static_cast<int>(f.id);
    string fij = format("(f, i, j) = ({:d}, {:4d}, {:4d})", f_, i, j);
    // The midpoint global coordinates (x, y, z)
    string xyz = format("(x, y, z) = ({:+8.6f}, {:+8.6f}, {:+8.6f})", x(),  y(), z());
    // Return one combined description
    return format("spid = {:8d}. {:s}. {:s}.\n", id(), fij, xyz);
}

// *****************************************************************************
/** Initialize a SkyPatch from its integer ID.*/
SkyPatch ks::SkyPatch_from_id(int32_t id)
{
    // First integer division; unpack id into cube face f and remainder x
    div_t qr = div(static_cast<int>(id), M2);
    int8_t f_ = static_cast<int8_t>(qr.quot);
    int32_t x = qr.rem;

    // Second integer division; unpack x into grid points i and j    
    qr = div(x, M);
    int16_t i_ = static_cast<int16_t>(qr.quot);
    int16_t j_ = static_cast<int16_t>(qr.rem);

    // Instantiate the sky patch
    return SkyPatch(f_, i_, j_);
}

// *****************************************************************************
vector<SkyPatch> ks::make_SkyPatch_table()
{
    // Initialize an empty vector and reserve space for N_sp entries
    vector<SkyPatch> spt;
    spt.reserve(N_sp);

    // Efficiently loop through f, i, j to initialize the entries on the sky patch table
    for (int8_t f=0; f<6; f++)
    {
        for (int16_t i=0; i<M; i++)
        {
            for (int16_t j=0; j<M; j++)
            {
                // The SkyPatch and its ID
                SkyPatch sp = SkyPatch(f, i, j);
                // Save this sky patch to the table
                spt.push_back(sp);
            }
        } // for / j
    } // for / f
    return spt;
}
