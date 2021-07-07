/** @file CandidateElement.hpp
 *  @brief Class to manage one set of candidate orbital elements.
 *  Calculate its trajectory (exactly with rebound or approximately by Kepler approximation).
 *  Calculate direction at a DetectionTime.
 *  Assemble a list of detections within a threshold for gradient search to improve the elements.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-07
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <string>
    using std::string;
#include <vector>
    using std::vector;

// Local dependencies
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
// #include "BodyVector.hpp"
//     using ks::BodyVector;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

// *****************************************************************************
namespace ks {

// *****************************************************************************
class CandidateElement
{
public:
    // Build a CandidateElement from an OrbitalElement and a DetectionTime table
    CandidateElement(OrbitalElement elt, DetectionTimeTable dtt);
    /// Destructor for CandidateElement.
    ~CandidateElement();

    /// Calculate asteroid trajectory (positions and velocity)
    void calc_trajectory();
    /// Calculate direction from asteroid trajectory to observatory
    void calc_direction();

private:
    /// Number of detection times
    const int N;
    /// Array of mjds when detections taken; size N
    double *mjds;
    /// Array of positions of an asteroid with these candidate elements; size 3N
    double *q_ast;
    /// Array of velocities of an asteroid with these candidate elements; size 3N
    double *v_ast;
    /// Array of directions to an asteroid with these candidate elements; size 3N
    double *u_ast;
};

// *****************************************************************************
} // namespace ks
