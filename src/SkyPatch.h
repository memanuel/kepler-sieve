/*****************************************************************************
 * Class encapsulating the concept of SkyPatch - one grid of the sky with
 * the unit sphere projected onto a circumscribing unit cube.
 * See DB tables SkyPatch, CubeFace and SkyPatchNeighbor.
 * 
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/
#pragma once

// *****************************************************************************
// Included libraries
#include <cmath>
#include <string>
#include <iostream>
#include <stdexcept>
#include <fmt/format.h>

// Local dependencies
#include "CubeFace.h"
#include "utils.h"

// *****************************************************************************
// Standard library class names used
using std::div;
using std::div_t;
using std::string;
using std::cout;
using std::range_error;
using fmt::print;
using fmt::format;

// *****************************************************************************
// Local names used
using ks::CubeFace;
using ks::sqr;

// *****************************************************************************
// Put classes into the namespace ks (for Kepler Sieve)
namespace ks {

// *****************************************************************************
// Set the grid size for the sky patch at compile time
namespace sky_patch{
    //*The multiplier for local coordinates; each cube face is on a 2Nx2N grid
    constexpr int N = 1024;
    //*The side length in dicrete grid units of each cube face
    constexpr int M = 2*N;
    //*The number of squares in each grid face
    constexpr int M2 = M*M;
    //*The number of SkyPatch entries; there are six faces on a cube
    constexpr int N_sp = 6*M2;
}

// *****************************************************************************
/** One cell on the surface of the unit sphere, used for spatial indexing. */
class SkyPatch
{
    public:
        // Constructor and destructor
        //*Initialize a SkyPatch from a CubeFace instance and grid coordinates.
        SkyPatch(CubeFace f_, int16_t i_, int16_t j_);
        //*Initialize a SkyPatch from its face ID and grid coordinates.
        SkyPatch(int8_t f_, int16_t i_, int16_t j_);
        //*Default destructor for SkyPatch.
        ~SkyPatch();

        // Data
        const CubeFace f;
        const int16_t i;
        const int16_t j;

        //* The integer ID of this sky patch
        const int32_t id() const;

        // Local coordinates of the center of this sky patch on the unit sphere
        //*The local coordinate indexed by i with label alpha
        const double u() const;
        //*The local coordinate indexed by j with label beta
        const double v() const;
        //*The largest local coordinate, with label gamma
        const double w() const;

        // Global coordinates of the center of this sky patch on the unit sphere
        //*Center of SkyPatch: global x
        const double x() const;
        //*Center of SkyPatch: global y
        const double y() const;
        //*Center of SkyPatch: global z
        const double z() const;

        //*Write the direction into an array of size 3
        void xyz(double *u);

        //*String description of this SkyPatch
        const string str() const;

        // Get another SkyPatch by shifting around the grid
        //*Nearby SkyPatch by shifting in the i direction
        const SkyPatch shift_i(const int16_t di) const;
        //*Nearby SkyPatch by shifting in the j direction
        const SkyPatch shift_j(const int16_t dj) const;
        //*Nearby SkyPatch by shifting both i and j
        const SkyPatch shift(const int16_t d, const int16_t dj) const;

    private:
        //*Distance to circumscribing cube, i.e. to point (a, b, c)
        const double r;

        // Local coordinates of the center of this sky patch, projected onto unit cube
        //*Local coordinate on cube face indexed by i (labeled alpha); corresponds to u.
        const double a() const;
        //*Local coordinate on cube face indexed by j (labeled beta); corresponds to v.
        const double b() const;
        //*Largest local coordinate on cube face; constant on the face; (labeled gamma); corresponds to w.
        const double c() const;
};

// *****************************************************************************
}; // namespace
