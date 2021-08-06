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
// Add two vectors
// *****************************************************************************

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
    // Add the two velocities componentwise
    return Velocity
    {
        .vx = v1.vx + v2.vx,
        .vy = v1.vy + v2.vy,
        .vz = v1.vz + v2.vz
    };
}

// *****************************************************************************
StateVector operator+ (const StateVector& s1, const StateVector& s2)
{
    // Add the two state vectors componentwise
    return StateVector
    {
        .qx = s1.qx + s2.qx,
        .qy = s1.qy + s2.qy,
        .qz = s1.qz + s2.qz,
        .vx = s1.vx + s2.vx,
        .vy = s1.vy + s2.vy,
        .vz = s1.vz + s2.vz
    };
}

// *****************************************************************************
// Subtract two vectors
// *****************************************************************************

// *****************************************************************************
Position operator- (const Position& p1, const Position& p2)
{
    // Add the two positions componentwise
    return Position
    {
        .qx = p1.qx - p2.qx,
        .qy = p1.qy - p2.qy,
        .qz = p1.qz - p2.qz
    };
}

// *****************************************************************************
Velocity operator- (const Velocity& v1, const Velocity& v2)
{
    // Add the two velocities componentwise
    return Velocity
    {
        .vx = v1.vx - v2.vx,
        .vy = v1.vy - v2.vy,
        .vz = v1.vz - v2.vz
    };
}

// *****************************************************************************
StateVector operator- (const StateVector& s1, const StateVector& s2)
{
    // Add the two positions componentwise
    return StateVector
    {
        .qx = s1.qx - s2.qx,
        .qy = s1.qy - s2.qy,
        .qz = s1.qz - s2.qz,
        .vx = s1.vx - s2.vx,
        .vy = s1.vy - s2.vy,
        .vz = s1.vz - s2.vz
    };
}

// *****************************************************************************
// Multiply a velocity by a scalar
// *****************************************************************************

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
// Extract position and velocity from a state vector
// *****************************************************************************

// *****************************************************************************
Position sv2pos(const StateVector& s)
{
    return Position {.qx=s.qx, .qy=s.qy, .qz=s.qz};
}

// *****************************************************************************
Velocity sv2vel(const StateVector& s)
{
    return Velocity {.vx=s.vx, .vy=s.vy, .vz=s.vz};
}

// *****************************************************************************
// Norm and distance of positions, velocities and state vectors
// *****************************************************************************

// *****************************************************************************
double norm(const Position& p)
{
    return sqrt(sqr(p.qx) + sqr(p.qy) + sqr(p.qz));
}

// *****************************************************************************
double norm(const Velocity& v)
{
    return sqrt(sqr(v.vx) + sqr(v.vy) + sqr(v.vz));
}

// *****************************************************************************
double norm(const StateVector& s)
{
    return sqrt(sqr(s.qx) + sqr(s.qy) + sqr(s.qz) + sqr(s.vx) + sqr(s.vy) + sqr(s.vz));
}

// *****************************************************************************
double dist(const Position& p1, const Position& p2) {return norm(p2-p1);}

// *****************************************************************************
double dist(const Velocity& v1, const Velocity& v2) {return norm(v2-v1);}

// *****************************************************************************
double dist(const StateVector& s1, const StateVector& s2) {return norm(s2-s1);}

// *****************************************************************************
double dist_dq(const StateVector& s1, const StateVector& s2) {return dist(sv2pos(s1), sv2pos(s2));}

// *****************************************************************************
double dist_dv(const StateVector& s1, const StateVector& s2) {return dist(sv2vel(s1), sv2vel(s2));}

// *****************************************************************************
double dist(const StateVector& s1, const Position& q2)
{
    // Extract position component
    Position q1 = sv2pos(s1);
    return dist(q1, q2);
}

// *****************************************************************************
double dist(const Position& q1, const StateVector& s2)
{
    // Extract position component
    Position q2 = sv2pos(s2);
    return dist(q1, q2);
}

// *****************************************************************************
double dist(const StateVector& s1, const Velocity& v2)
{
    // Extract velocity component
    Velocity v1 = sv2vel(s1);
    return dist(v1, v2);
}

// *****************************************************************************
double dist(const Velocity& v1, const StateVector& s2)
{
    // Extract velocity component
    Velocity v2 = sv2vel(s2);
    return dist(v1, v2);
}

// *****************************************************************************
// Check if two vectors are close
// *****************************************************************************

// *****************************************************************************
bool is_close(const Position& p1, const Position& p2, double tol_dq)
    {return dist(p1, p2) < tol_dq;}

// *****************************************************************************
bool is_close(const Velocity& v1, const Velocity& v2, double tol_dv)
    {return dist(v1, v2) < tol_dv;}

// *****************************************************************************
bool is_close(const StateVector& s1, const StateVector& s2, double tol_dq, double tol_dv)
{
    // Extract position and velocity components
    Position q1 = sv2pos(s1);
    Position q2 = sv2pos(s2);
    Velocity v1 = sv2vel(s1);
    Velocity v2 = sv2vel(s2);
    return (dist(q1, q2) < tol_dq) && (dist(v1, v2) < tol_dv);
}

// *****************************************************************************
bool is_close(const StateVector& s1, const Position& q2, double tol_dq) 
    {return dist(s1, q2) < tol_dq;}

// *****************************************************************************
bool is_close(const Position& q1, const StateVector& s2, double tol_dq) 
    {return dist(q1, s2) < tol_dq;};

bool is_equal(const StateVector& s1, const StateVector& s2)
    {return (s1.qx==s2.qx) && (s1.qy==s2.qy) && (s1.qz==s2.qz) && 
            (s1.vx==s2.vx) && (s1.vy==s2.vy) && (s1.vz==s2.vz);}

// *****************************************************************************
// Print a position
// *****************************************************************************

// *****************************************************************************
void print_position_headers(const string prefix)
    {print("{:s}{:10s} : {:10s} : {:10s}\n", 
            prefix, " qx", " qy", " qz");}

// *****************************************************************************
void print_position(const Position& p, const string prefix)
    {print("{:s}{:+10.6f} : {:+10.6f} : {:+10.6f}\n", 
        prefix, p.qx, p.qy, p.qz);}

// *****************************************************************************
void print_position_sci(const Position& p, const string prefix)
    {print("{:s}{:+10.2e} : {:+10.2e} : {:+10.2e}\n", 
        prefix, p.qx, p.qy, p.qz);}

// *****************************************************************************
// Print a state vector
// *****************************************************************************

// *****************************************************************************
void print_state_vector_headers(const string prefix)
    {print("{:s}{:10s} : {:10s} : {:10s} : {:10s} : {:10s} : {:10s}\n", 
            prefix, " qx", " qy", " qz", " vx", " vy", " vz");}

// *****************************************************************************
void print_state_vector(const StateVector& s, const string prefix)
    {print("{:s}{:+10.6f} : {:+10.6f} : {:+10.6f} : {:+10.6f} : {:+10.6f} : {:+10.6f}\n", 
        prefix, s.qx, s.qy, s.qz, s.vx, s.vy, s.vz);}

// *****************************************************************************
void print_state_vector_sci(const StateVector& s, const string prefix)
    {print("{:s}{:+10.2e} : {:+10.2e} : {:+10.2e} : {:+10.2e} : {:+10.2e} : {:+10.2e}\n", 
        prefix, s.qx, s.qy, s.qz, s.vx, s.vy, s.vz);}

// *****************************************************************************
void print_state_vector_long(const StateVector& s)
{
    print("qx = {:+12.8f}\n", s.qx);
    print("qy = {:+12.8f}\n", s.qy);
    print("qz = {:+12.8f}\n", s.qz);
    print("vx = {:+12.8f}\n", s.vx);
    print("vy = {:+12.8f}\n", s.vy);
    print("vz = {:+12.8f}\n", s.vz);
}

// *****************************************************************************
} // namespace ks
