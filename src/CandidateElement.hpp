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
#include <cstring>
    // memcpy

// Local dependencies
#include "constants.hpp"
    using ks::cs::mu_sun;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "StateVector.hpp"
    using ks::norm;
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::mean_motion;
    using ks::anomaly_M2f;
    using ks::elt2vec;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;

// *****************************************************************************
namespace ks {

// *****************************************************************************
class CandidateElement
{
public:
    // Build a CandidateElement from an OrbitalElement and a DetectionTime table
    CandidateElement(OrbitalElement& elt, DetectionTimeTable& dtt);
    /// Destructor for CandidateElement.
    ~CandidateElement();

    /// The underlying OrbitalElement
    OrbitalElement elt;
    /// Calculate asteroid trajectory (positions and velocity)
    void calc_trajectory();
    /// Calculate direction from asteroid trajectory to observatory
    void calc_direction();

    /// Read access to array of mjd when detections were observed (print time); size N
    double* get_mjd() const {return mjd;}
    /// Read access to array of positions of observatory; size 3N
    double* get_q_obs() const {return q_obs;}
    /// Read access to array of positions of an asteroid with these candidate elements; size 3N
    double* get_q_ast() const {return q_ast;}
    /// Read access to array of velocities of an asteroid with these candidate elements; size 3N
    double* get_v_ast() const {return v_ast;}
    /// Read access to array of directions of an asteroid with these candidate elements; size 3N
    double* get_u_ast() const {return u_ast;}

private:
    /// Number of detection times
    const int N;
    /// Array of mjd when detections were observed (print time); size N
    double* mjd;
    /// Array of positions of observatory; size 3N
    double* q_obs;
    /// Array of positions of an asteroid with these candidate elements; size 3N
    double* q_ast;
    /// Array of velocities of an asteroid with these candidate elements; size 3N
    double* v_ast;
    /// Array of directions to an asteroid with these candidate elements; size 3N
    double* u_ast;
};

// *****************************************************************************
} // namespace ks
