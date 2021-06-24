// *****************************************************************************
// Included files
#include "sky_patch.h"

// *****************************************************************************
using std::range_error;
using ks::CubeFace;

// *****************************************************************************
// Instantiate exception for a bad CubeFaceID
namespace ks{
range_error err_cube_face_id = range_error("CubeFace id must be between 0 and 5, inclusive.\n");
}
using ks::err_cube_face_id;

// *****************************************************************************
/** Initialize a CubeFace from its ID; starts from 0.*/
CubeFace::CubeFace(int8_t id_) : 
id(id_) 
{
    if ((id_ < 0) || (id_ > 5)) {throw err_cube_face_id;}
}

// *****************************************************************************
/** Default destructor for CubeFace.*/
CubeFace::~CubeFace() {}

// *****************************************************************************
const string CubeFace::str() const
{
    switch (id)
    {
        case 0:     return string("Z+");    // Face Z+; trio (u, v, w) = (X, Y, Z)
        case 1:     return string("Y+");    // Face Y+; trio (u, v, w) = (Z, X, Y)
        case 2:     return string("X+");    // Face X+; trio (u, v, w) = (Y, Z, X)
        case 3:     return string("X-");    // Face X-; trio (u, v, w) = (Y, Z, X)
        case 4:     return string("Y-");    // Face Y-; trio (u, v, w) = (Z, X, Y)
        case 5:     return string("Z-");    // Face Z-; trio (u, v, w) = (X, Y, Z)
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
const string CubeFace::description() const
{
    switch (id)
    {
        case 0:     return string("Face Z+; trio (u, v, w) = (X, Y, Z)");
        case 1:     return string("Face Y+; trio (u, v, w) = (Z, X, Y)");
        case 2:     return string("Face X+; trio (u, v, w) = (Y, Z, X)");
        case 3:     return string("Face X-; trio (u, v, w) = (Y, Z, X)");
        case 4:     return string("Face Y-; trio (u, v, w) = (Z, X, Y)");
        case 5:     return string("Face Z-; trio (u, v, w) = (X, Y, Z)");
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
const int CubeFace::k1() const
{
    switch (id)
    {
        case 0:     return 1;   // Face Z+; (X, Y, Z); alpha=X
        case 1:     return 3;   // Face Y+; (Z, X, Y); alpha=Z
        case 2:     return 2;   // Face X+; (Y, Z, X); alpha=Y
        case 3:     return 2;   // Face X-; (Y, Z, X); alpha=Y
        case 4:     return 3;   // Face Y-; (Z, X, Y); alpha=Z
        case 5:     return 1;   // Face Z-; (X, Y, Z); alpha=X
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
const int CubeFace::k2() const
{
    switch (id)
    {
        case 0:     return 2;   // Face Z+; (X, Y, Z); beta=Y
        case 1:     return 1;   // Face Y+; (Z, X, Y); beta=X
        case 2:     return 3;   // Face X+; (Y, Z, X); beta=Z
        case 3:     return 3;   // Face X-; (Y, Z, X); beta=Z
        case 4:     return 1;   // Face Y-; (Z, X, Y); beta=X
        case 5:     return 2;   // Face Z-; (X, Y, Z); beta=X
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Integer index of the largest axis, with local values w; e.g. for Z+, i=3 and gamma='Z'.*/
const int CubeFace::k3() const
{
    switch (id)
    {
        case 0:     return 3;   // Face Z+; (X, Y, Z); gamma=Z
        case 1:     return 2;   // Face Y+; (Z, X, Y); gamma=Y
        case 2:     return 1;   // Face X+; (Y, Z, X); gamma=X
        case 3:     return 1;   // Face X-; (Y, Z, X); gamma=X
        case 4:     return 2;   // Face Y-; (Z, X, Y); gamma=Y
        case 5:     return 3;   // Face Z-; (X, Y, Z); gamma=Z
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
// Get a string label from an integer index by adding the index to the letter before (X, Y, Z), i.e. 'W'
constexpr int label_base = static_cast<int>('W');

// *****************************************************************************
/**Name of the first varying axis, with local values 'u'; e.g. for Z+, j1=1 and alpha='X'.*/
const char CubeFace::alpha() const
{    
    return static_cast<const char>(label_base + k1());
}

// *****************************************************************************
/**Name of the second varying axis, with local values 'v'; e.g. for Z+, j2=2 and beta='Y'.*/
const char CubeFace::beta() const
{
    return static_cast<const char>(label_base + k2());
}

// *****************************************************************************
/**Name of the largest axis, with local values 'w'; e.g. for Z+, i=3 and gamma='Z'.*/
const char CubeFace::gamma() const
{
    return static_cast<const char>(label_base + k3());
}

// *****************************************************************************
/**Sign of the largest axis, with local values 'w'; e.g. for Z+, ci=1.0 .*/
const double CubeFace::c() const
{
    switch (id)
    {
        case 0:     return 1.0;     // Face Z+
        case 1:     return 1.0;     // Face Y+
        case 2:     return 1.0;     // Face X+
        case 3:     return -1.0;    // Face X-
        case 4:     return -1.0;    // Face Y_
        case 5:     return -1.0;    // Face Z_
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Neighbor when i=0 (u=-1.0).*/
const CubeFace CubeFace::neighbor_i0() const
{
    switch (id)
    {
        case 0: return 3;       // Z+; (X, Y, Z); neighbor X-
        case 1: return 5;       // Y+; (Z, X, Y); neighbor Z-
        case 2: return 4;       // X+; (Y, Z, X); neighbor Y-
        case 3: return 4;       // X-; (Y, Z, X); neighbor Y-
        case 4: return 5;       // Y-; (Z, X, Y); neighbor Z-
        case 5: return 3;       // Z-; (X, Y, Z); neighbor X-
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Neighbor when i=M (u=+1.0).*/
const CubeFace CubeFace::neighbor_i1() const
{
    switch (id)
    {
        case 0: return 2;       // Z+; (X, Y, Z); neighbor X+
        case 1: return 0;       // Y+; (Z, X, Y); neighbor Z+
        case 2: return 1;       // X+; (Y, Z, X); neighbor Y+
        case 3: return 1;       // X-; (Y, Z, X); neighbor Y+
        case 4: return 0;       // Y-; (Z, X, Y); neighbor Z+
        case 5: return 2;       // Z-; (X, Y, Z); neighbor X+
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Neighbor when j=0 (v=-1.0).*/
const CubeFace CubeFace::neighbor_j0() const
{
    switch (id)
    {
        case 0: return 4;       // Z+; (X, Y, Z); neighbor Y-
        case 1: return 3;       // Y+; (Z, X, Y); neighbor X-
        case 2: return 5;       // X+; (Y, Z, X); neighbor Z-
        case 3: return 5;       // X-; (Y, Z, X); neighbor Z-
        case 4: return 5;       // Y-; (Z, X, Y); neighbor X-
        case 5: return 4;       // Z-; (X, Y, Z); neighbor Y-
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Neighbor when j=1 (v=+1.0).*/
const CubeFace CubeFace::neighbor_j1() const
{
    switch (id)
    {
        case 0: return 1;       // Z+; (X, Y, Z); neighbor Y+
        case 1: return 2;       // Y+; (Z, X, Y); neighbor X+
        case 2: return 0;       // X+; (Y, Z, X); neighbor Z+
        case 3: return 0;       // X-; (Y, Z, X); neighbor Z+
        case 4: return 2;       // Y-; (Z, X, Y); neighbor X+
        case 5: return 1;       // Z-; (X, Y, Z); neighbor Y+
        default:    throw err_cube_face_id;
    }
}
