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
