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
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
#include "BodyVector.hpp"
    using ks::BodyVector;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "Simulation.hpp"
    using ks::reb::Simulation;
    using ks::reb::make_sim_planets;

// *****************************************************************************
namespace ks {

// *****************************************************************************
class CandidateElement
{
public:
    // Build a CandidateElement from an OrbitalElement with a given size
    CandidateElement(OrbitalElement elt, int32_t candidate_id, int N_t);

    // Build a CandidateElement from an OrbitalElement at the desired output times
    CandidateElement(OrbitalElement elt, int32_t candidate_id, const double* mjd, int N_t);

    // Build a CandidateElement from an OrbitalElement using the shared DetectionTime table
    CandidateElement(OrbitalElement elt, int32_t candidate_id);

    /// Destructor for CandidateElement.
    ~CandidateElement();

    /// Initialize a newly constructed CandidateElement with desired output times
    void init(const double* mjd);

    /// The underlying OrbitalElement; mutable
    OrbitalElement elt;
    /// The integer ID of these candidate elements; plays the role of a body_id
    const int32_t candidate_id;

    /// Calibrate asteroid trajectory to a rebound integration
    void calibrate(const PlanetVector& pv);
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

    /// Extract a StateVector from the q_ast and v_ast arrays
    const StateVector state_vector(int i) const;

private:
    /// Initial value of element used to initialize this object
    const OrbitalElement elt0;

    /// One BodyVector object for Sun shared by all instances
    const inline static BodyVector bv_sun = BodyVector(SolarSystemBody_bv::sun);  
    /// One BodyVector object for Earth shared by all instances
    const inline static BodyVector bv_earth = BodyVector(SolarSystemBody_bv::earth);
    /// One DetectionTimeTable object shared by all instances
    const inline static DetectionTimeTable dtt = DetectionTimeTable();
    /// One DetectionTable object shared by all instances
    // const inline static DetectionTable dt = DetectionTable();

    /// Number of detection times
    const int N_t;
    /// Number of rows of data in spatial arrays for q and v
    const int N_row;
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
    /// Array of position shifts; includes (1) sun position (2) numerical calibration adjustment
    double* dq_ast;
    /// Array of velocity shifts; includes (1) sun velocity (2) numerical calibration adjustment
    double* dv_ast;
};

// *****************************************************************************
} // namespace ks
