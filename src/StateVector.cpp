/** @file   StateVector.cpp
 *  @brief  Implementation of class StateVector.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-02
 * 
 */

// *****************************************************************************
// Local dependencies
#include "StateVector.hpp"

// *****************************************************************************
namespace ks {

// *****************************************************************************
Position operator+ (const Position& p1, const Position& p2)
{
    // Add the two positions componentwise
    return Position
    {
        .qx = p1.qx + p2.qx,
        .qy = p1.qy + p2.qy,
        .qz = p1.qz + p2.qz
    };
}

// *****************************************************************************
Velocity operator+ (const Velocity& v1, const Velocity& v2)
{
    // Add the two positions componentwise
    return Velocity
    {
        .vx = v1.vx + v2.vx,
        .vy = v1.vy + v2.vy,
        .vz = v1.vz + v2.vz
    };
}

// *****************************************************************************
StateVector operator+ (const StateVector& v1, const StateVector& v2)
{
    // Add the two positions componentwise
    return StateVector
    {
        .qx = v1.qx + v2.qx,
        .qy = v1.qy + v2.qy,
        .qz = v1.qz + v2.qz,
        .vx = v1.vx + v2.vx,
        .vy = v1.vy + v2.vy,
        .vz = v1.vz + v2.vz
    };
}

// *****************************************************************************
Velocity operator* (const Velocity& v, const double alpha)
{
    return Velocity
    {
        .vx = v.vx*alpha,
        .vy = v.vy*alpha,
        .vz = v.vz*alpha
    };
}

// *****************************************************************************
Velocity operator* (const double alpha, const Velocity& v) 
{
    // Multiplication is commutative; delegate to to the previous definition with alpha on the right
    return operator* (v, alpha);
}

// *****************************************************************************
Position sv2pos(StateVector& sv)
{
    return Position {.qx=sv.qx, .qy=sv.qy, .qz=sv.qz};
}

// *****************************************************************************
Velocity sv2vel(StateVector& sv)
{
    return Velocity {.vx=sv.vx, .vy=sv.vy, .vz=sv.vz};
}

// *****************************************************************************
void print_state_vector(StateVector& sv, bool header)
{
    if (header)
    {print("{:9s} : {:9s} : {:9s} : {:9s} : {:9s} : {:9s}\n", 
        "qx", "qy", "qz", "vx", "vy", "vz");}
    else
    {print("{:+9.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f} : {:+9.6f}\n", 
        sv.qx, sv.qy, sv.qz, sv.vx, sv.vy, sv.vz);}
}

// *****************************************************************************
void print_state_vector_long(StateVector& sv)
{
    print("qx = {:+12.8f}\n", sv.qx);
    print("qy = {:+12.8f}\n", sv.qy);
    print("qz = {:+12.8f}\n", sv.qz);
    print("vx = {:+12.8f}\n", sv.vx);
    print("vy = {:+12.8f}\n", sv.vy);
    print("vz = {:+12.8f}\n", sv.vz);
}

// *****************************************************************************
double norm(Position p1, Position p2)
{
    double dq2 = sqr(p1.qx-p2.qx) + sqr(p1.qy-p2.qy) + sqr(p1.qz-p2.qz);
    return sqrt(dq2);
}

// *****************************************************************************
double norm(Velocity v1, Velocity v2)
{
    double dv2 = sqr(v1.vx-v2.vx) + sqr(v1.vy-v2.vy) + sqr(v1.vz-v2.vz);
    return sqrt(dv2);
}

// *****************************************************************************
double norm(StateVector v1, StateVector v2)
{
    double dq2 = sqr(v1.qx-v2.qx) + sqr(v1.qy-v2.qy) + sqr(v1.qz-v2.qz);
    double dv2 = sqr(v1.vx-v2.vx) + sqr(v1.vy-v2.vy) + sqr(v1.vz-v2.vz);
    return sqrt(dq2 + dv2);
}

// *****************************************************************************
} // namespace ks
