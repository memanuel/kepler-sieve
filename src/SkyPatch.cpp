// *****************************************************************************
// Included files
#include <cmath>
#include "utils.h"
#include "sky_patch.h"

// *****************************************************************************
using std::div;
using std::div_t;
using std::range_error;
using boost::format;
using ks::CubeFace;
using ks::SkyPatch;
using ks::sqr;

// *****************************************************************************
// Set the grid size at compile time
constexpr int N = ks::N_sky_patch;
constexpr int M = 2*N;
constexpr int M2 = M*M;
constexpr double N_inv = 1.0 / static_cast<double>(N);
constexpr int N_spc = 6*M2;

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
/** Initialize a SkyPatch from its face ID and grid coordinates.*/
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
/** Default destructor for SkyPatch.*/
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
/**Coordinate on first varying axis (alpha).*/
const double SkyPatch::u() const
{
    return a() / r;
}

// *****************************************************************************
/**Coordinate on second varying axis (beta).*/
const double SkyPatch::v() const
{
    return b() / r;
}

// *****************************************************************************
/**Coordinate on third varying axis (gamma).*/
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

/******************************************************************************
Functions for working with SkyPatch objects
******************************************************************************/
namespace ks {

// *****************************************************************************
int fij2spid(int8_t f, int16_t i, int16_t j)
{
    return (M2*f) + (M*i) + j;
}

// *****************************************************************************
/** Initialize a SkyPatch from its integer ID.*/
SkyPatch SkyPatch_from_id(int32_t id)
{
    // First integer division; unpack id into cube face f and remainder x
    div_t qr = div(id, M2);
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
int sky_patch_count()
{
    return static_cast<int>(N_spc);
}

// *****************************************************************************
void write_sky_patch_neighbor_table(int32_t* spn)
{
    // Loop through the starting sky patch, with ID sky_patch_id_1
    for (int8_t f=0; f<6; f++)
    {
        // The offset for this grid face is M2*f
        int32_t idx_f = M2*f;
        for (int16_t i0=0; i0<M; i0++)
        // for (int16_t i0=0; i0<10; i0++)
        {
            for (int16_t j0=0; j0<M; j0++)
            {
                // The starting SkyPatchID
                int32_t spid0 = idx_f + (M*i0) + j0;
                // The starting index
                size_t idx = spid0*9;
                
                // Get grid coordinates of 9 candidate neighbors.
                for (int16_t di=-1; di<=1; di++)
                {
                    // The grid entry i1 for the three candidate neighbors in this row
                    int16_t i1 = i0+di;
                    int32_t idx_i1 = M*i1;
                    for (int16_t dj=-1; dj<=1; dj++)
                    {
                        // The grid entry j1 for this candidate neighbor
                        int16_t j1 = j0+dj;
                        // Write the candidate neighbor to the array if it is on the same face.
                        // Don't worry about wrapping around edges and corners.
                        // Write the dummy value -1 if the neighbor is not a real grid point
                        bool is_on_grid = (0 <= i1) && (i1 < M) && (0 <= j1) && (j1 < M);
                        spn[idx++] = is_on_grid ? (idx_f+idx_i1+j1) : -1;
                    }
                }
            }
        }
    }
}

// *****************************************************************************
}; // namespace
