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
/** Two character description, e.g. "Z+"; see DB table KS.CubeFace. */
const string CubeFace::code()
{
    switch (id)
    {
        case 0:     return string("Z+");
        case 1:     return string("Y+");
        case 2:     return string("X+");
        case 3:     return string("X-");
        case 4:     return string("Y-");
        case 5:     return string("Z-");
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Integer index of the largest axis, with local values w; e.g. for Z+, i=3 and gamma='Z'.*/
const int8_t CubeFace::i()
{
    switch (id)
    {
        case 0:     return 3;
        case 1:     return 2;
        case 2:     return 1;
        case 3:     return 1;
        case 4:     return 2;
        case 5:     return 3;
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Integer index of the first varying axis, with local values 'u'; e.g. for Z+, j1=1 and alpha='X'.*/
const int8_t CubeFace::j1()
{
    switch (id)
    {
        case 0:     return 1;
        case 1:     return 3;
        case 2:     return 2;
        case 3:     return 2;
        case 4:     return 3;
        case 5:     return 1;
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Integer index of the second varying axis, with local values 'v'; e.g. for Z+, j2=2 and beta='Y'.*/
const int8_t CubeFace::j2()
{
    switch (id)
    {
        case 0:     return 2;
        case 1:     return 1;
        case 2:     return 3;
        case 3:     return 3;
        case 4:     return 1;
        case 5:     return 2;
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Name of the first varying axis, with local values 'u'; e.g. for Z+, j1=1 and alpha='X'.*/
const char CubeFace::alpha()
{
    switch (id)
    {
        case 0:     return 'X';
        case 1:     return 'Z';
        case 2:     return 'Y';
        case 3:     return 'Y';
        case 4:     return 'Z';
        case 5:     return 'X';
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Name of the second varying axis, with local values 'v'; e.g. for Z+, j2=2 and beta='Y'.*/
const char CubeFace::beta()
{
    switch (id)
    {
        case 0:     return 'Y';
        case 1:     return 'X';
        case 2:     return 'Z';
        case 3:     return 'Z';
        case 4:     return 'X';
        case 5:     return 'Y';
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Name of the largest axis, with local values 'w'; e.g. for Z+, i=3 and gamma='Z'.*/
const char CubeFace::gamma()
{
    switch (id)
    {
        case 0:     return 'Z';
        case 1:     return 'Y';
        case 2:     return 'X';
        case 3:     return 'X';
        case 4:     return 'Y';
        case 5:     return 'Z';
        default:    throw err_cube_face_id;
    }
}

// *****************************************************************************
/**Sign of the largest axis, with local values 'w'; e.g. for Z+, ci=1.0 .*/
const double CubeFace::c() const
{
    switch (id)
    {
        case 0:     return 1.0;
        case 1:     return 1.0;
        case 2:     return 1.0;
        case 3:     return -1.0;
        case 4:     return -1.0;
        case 5:     return -1.0;
        default:    throw err_cube_face_id;
    }
}
