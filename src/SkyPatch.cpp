// *****************************************************************************
// Included files
#include <cmath>
#include "sky_patch.h"

// *****************************************************************************
using std::div;
using std::div_t;
using std::range_error;
using ks::CubeFace;
using ks::SkyPatch;
using ks::sqr;

// *****************************************************************************
// Set the grid size at compile time
// constexpr int N=ks::N_sky_patch;
constexpr int N = ks::N_sky_patch;
constexpr int M=2*N;
constexpr int M2 = M*M;
// constexpr double Nf = static_cast<double>(N);
constexpr double N_inv = 1.0 / static_cast<double>(N);

// *****************************************************************************
// Instantiate exceptions for a bad cube face f or grid entry (i or j)
namespace ks {
range_error err_sky_patch_f = range_error("SkyPatch cube face f must be between 0 and 6, inclusive.\n");
range_error err_sky_patch_grid_i = range_error("SkyPatch grid coordinate i must be between 0 and M=2048, exclusive.\n");
range_error err_sky_patch_grid_j = range_error("SkyPatch grid coordinate j must be between 0 and M=2048, exclusive.\n");
}
using ks::err_sky_patch_f, ks::err_sky_patch_grid_i, ks::err_sky_patch_grid_j;

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
/** Default destructor for SkyPatch.*/
SkyPatch::~SkyPatch() {}

// *****************************************************************************
/**Integer ID of this sky patch*/
const int32_t SkyPatch::id()
{
    return (M2*f.id) + (M*i) + j;
}

// *****************************************************************************
/**Coordinate on unit cube of first varying axis (alpha).*/
const double SkyPatch::a()
{
    return N_inv*(static_cast<double>(i)+0.5) - 1.0;
}

// *****************************************************************************
/**Coordinate on unit cube of second varying axis (beta).*/
const double SkyPatch::b()
{
    return N_inv*(static_cast<double>(j)+0.5) - 1.0;
}

// *****************************************************************************
/**Coordinate on unit cube of third varying axis (gamma).*/
const double SkyPatch::c()
{
    return f.c();
}

// *****************************************************************************
/**Coordinate on first varying axis (alpha).*/
const double SkyPatch::u()
{
    return a() / r;
}

// *****************************************************************************
/**Coordinate on second varying axis (beta).*/
const double SkyPatch::v()
{
    return b() / r;
}

// *****************************************************************************
/**Coordinate on third varying axis (gamma).*/
const double SkyPatch::w()
{
    return c() / r;
}
