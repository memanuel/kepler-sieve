/** @file   Direction.cpp
 *  @brief  Implementation of class Direction and associated functions.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-08-06
 * 
 */

// *****************************************************************************
// Local dependencies
#include "Direction.hpp"

// *****************************************************************************
namespace ks {

// *****************************************************************************
// Add and subtract direction vectors
// *****************************************************************************

// *****************************************************************************
// Add two direction vectors; result no longer on unit sphere
Direction operator+ (const Direction& u1, const Direction& u2)
{
    // Add the two directions componentwise
    return Direction
    {
        .ux = u1.ux + u2.ux,
        .uy = u1.uy + u2.uy,
        .uz = u1.uz + u2.uz
    };
}

// *****************************************************************************
// Subtract two direction vectors; result no longer on unit sphere
Direction operator- (const Direction& u1, const Direction& u2)
{
    // Subtract the two directions componentwise
    return Direction
    {
        .ux = u1.ux - u2.ux,
        .uy = u1.uy - u2.uy,
        .uz = u1.uz - u2.uz
    };
}

// *****************************************************************************
// Norm and distance of direction vectors
// *****************************************************************************

// *****************************************************************************
// Return the norm of a direction vector
double norm(const Direction& u)
{
    return sqrt(sqr(u.ux) + sqr(u.uy) + sqr(u.uz));
}

// *****************************************************************************
// Return the distance between two direction vectors
double dist(const Direction& u1, const Direction& u2)
{
    return norm(u2-u1);
}

// *****************************************************************************
// Check if two directions are close
// *****************************************************************************

// *****************************************************************************
// Test if two position vectors are close within the given tolerance
bool is_close(const Direction& u1, const Direction& u2, double tol)
    {return dist(u1, u2) < tol;}

// *****************************************************************************
// Calculate a direction from position of observer and state vector of target
// *****************************************************************************

// *****************************************************************************
// Create a direction by normalizing a position (difference) vector
ObservationResult normalize_pos(const Position& pos)
{
    // First get the denominator
    double r = norm(pos);
    // Divide the components by r and wrap into a Direction
    Direction u
    {
        .ux = pos.qx / r,
        .uy = pos.qy / r,
        .uz = pos.qz / r
    };
    // Assemble into an ObservationResult
    return ObservationResult
    {
        .u = u,
        .r = r
    };
}

// *****************************************************************************
// Calculate observed direction from an object with Newtonian light model
ObservationResult observe(const Position& q_obs, const StateVector& s_tgt)
{
    // Calculate the instantaneous position difference
    Position dq 
    {
        .qx = s_tgt.qx - q_obs.qx,
        .qy = s_tgt.qy - q_obs.qy,
        .qz = s_tgt.qz - q_obs.qz
    };
    // Instantaneous distance and light time
    double r {norm(dq)};
    double lt {c_inv * r};
    // Adjust the displacement dq to account for the light time
    Position dq_adj 
    {
        .qx = dq.qx - s_tgt.vx * lt,
        .qy = dq.qy - s_tgt.vy * lt,
        .qz = dq.qz - s_tgt.vz * lt
    };
    // Wrap the adjusted displacement into an ObservationResult
    return normalize_pos(dq_adj);
}

// *****************************************************************************
// Calculate only the direction between two objects with Newtonian light model
Direction observe_dir(const Position& q_obs, const StateVector& s_tgt)
{
    // Calculate the instantaneous position difference
    Position dq 
    {
        .qx = s_tgt.qx - q_obs.qx,
        .qy = s_tgt.qy - q_obs.qy,
        .qz = s_tgt.qz - q_obs.qz
    };
    // Instantaneous distance and light time
    double r {norm(dq)};
    double lt {c_inv * r};
    // Adjust the displacement dq to account for the light time
    Position dq_adj 
    {
        .qx = dq.qx - s_tgt.vx * lt,
        .qy = dq.qy - s_tgt.vy * lt,
        .qz = dq.qz - s_tgt.vz * lt
    };
    // Adjusted distance
    double r_adj {norm(dq_adj)};
    // Adjusted direction
    return Direction
    {
        .ux = dq_adj.qx / r_adj,
        .uy = dq_adj.qy / r_adj,
        .uz = dq_adj.qz / r_adj
    };    
}

// *****************************************************************************
// Print a direction
// *****************************************************************************

// *****************************************************************************
/// Print a column headers for one line description of a direction
void print_direction_headers(const string prefix)
    {print("{:s}{:10s} : {:10s} : {:10s}\n", 
            prefix, " ux", " uy", " uz");}

// *****************************************************************************
// Print a one line description of a direction
void print_direction(const Direction& u, const string prefix)
    {print("{:s}{:+10.6f} : {:+10.6f} : {:+10.6f}\n", 
        prefix, u.ux, u.uy, u.uz);}

// *****************************************************************************
// Print a one line description of a direction in scientific notation
void print_direction_sci(const Direction& u, const string prefix)
    {print("{:s}{:+10.2e} : {:+10.2e} : {:+10.2e}\n", 
        prefix, u.ux, u.uy, u.uz);}

// *****************************************************************************
} // namespace ks
