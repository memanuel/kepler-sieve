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
#include "constants.hpp"
    using ks::cs::mu_sun;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "StateVector.hpp"
    using ks::norm;
#include "Direction.hpp"
    using ks::Direction;
    using ks::ObservationResult;
    using ks::observe;
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
// Constructor and destructor
public:
    // Build a CandidateElement from an OrbitalElement with a given size
    CandidateElement(OrbitalElement elt, int32_t candidate_id, int N_t);

    // Build a CandidateElement from an OrbitalElement at the desired output times
    CandidateElement(OrbitalElement elt, int32_t candidate_id, const double* const mjd_in, int N_t);

    // Build a CandidateElement from an OrbitalElement using the shared DetectionTime table
    CandidateElement(OrbitalElement elt, int32_t candidate_id);

    /// Destructor for CandidateElement.
    ~CandidateElement();

// Public interface
public:
    /// Calibrate asteroid trajectory to a rebound integration
    void calibrate();
    /// Calculate asteroid trajectory (positions and velocity)
    void calc_trajectory();
    /// Calculate direction from asteroid trajectory to observatory
    void calc_direction();

    /// Search for all rows where the direction within a tolerance of a detection in the DetectionTable.
    /// This is brute force version that iterates through every detection.
    vector<int> search_bf(double thresh_dist) const;
    /// Search for all rows where the direction within a tolerance of a detection in the DetectionTable
    /// This is "smart" version that uses SkyPatch indexing.
    vector<int> search(double thresh_dist) const;
    /// Construct a new CandidateElement by filtering down to a list of selected rows
    CandidateElement filter_rows(vector<int> ii);

    /// Extract a StateVector from the q_ast_ and v_ast_ arrays
    const StateVector state_vector(int i) const;
    /// Extract a Position of the observer from the q_obs_ array
    const Position observer_pos(int i) const;
    /// Extract a Direction to the asteroid from the u_ast_ array
    const Direction direction(int i) const;
    /// Extract a distance to the asteroid from the r_ast_ array
    const inline double distance(int i) const {return r_ast_[i];}
    /// Report time to build static objects
    double report_static_time() const;

// Static data members
private:
    /// One PlanetVector object shared by all instances
    const static PlanetVector pv;
    /// One BodyVector object for Sun shared by all instances
    const static BodyVector bv_sun;  
    /// One BodyVector object for Earth shared by all instances
    const static BodyVector bv_earth;
    /// Structure to save time required for initial construction
    const static std::map<string, double> time_profile;
public:    
    /// One DetectionTimeTable object shared by all instances
    const static DetectionTimeTable dtt;
    /// One DetectionTable object shared by all instances
    const static DetectionTable dt;

// Data members
public:
    /// The underlying OrbitalElement; mutable
    OrbitalElement elt;
    /// Initial value of element used to initialize this object
    const OrbitalElement elt0;
    /// The integer ID of these candidate elements; plays the role of a body_id
    const int32_t candidate_id;
private:
    /// Number of detection times
    const int N_t;
    /// Number of rows of data in spatial arrays for q and v
    const int N_row;
    /// Array of mjd when detections were observed (print time); size N
    double* mjd_;
    // TODO - add detection_time_id
    /// Array of detection_time_id, when applicable; otherwise 0
    /// Array of positions of observatory; size 3N
    double* q_obs_;
    /// Array of positions of an asteroid with these candidate elements; size 3N
    double* q_ast_;
    /// Array of velocities of an asteroid with these candidate elements; size 3N
    double* v_ast_;
    /// Array of directions to an asteroid with these candidate elements; size 3N
    double* u_ast_;
    /// Array of distances to an asteroid with these candidate elements; size N
    double* r_ast_;
    /// Array of position shifts; includes (1) sun position (2) numerical calibration adjustment
    double* dq_ast_;
    /// Array of velocity shifts; includes (1) sun velocity (2) numerical calibration adjustment
    double* dv_ast_;

// Read-only copies of data arrays    
public:
    /// Read-only copy of mjd
    const double* const mjd;    
    /// Read-only copy of u_ast
    const double* const u_ast;
    /// Read-only copy of r_ast
    const double* const r_ast;

// Implementation functions
private:    
    /// Initialize a newly constructed CandidateElement with desired output times
    void init(const double* mjd);
    // Calculate three array indices jx, jy, jz for spatial data from a time index i
    /// Calculate array index jx for spatial data from time index i (x component)
    inline const int i2jx(const int& i) const {return 3*i+0;}
    /// Calculate array index jy for spatial data from time index i (y component)
    inline const int i2jy(const int& i) const {return 3*i+1;}
    /// Calculate array index jz for spatial data from time index i (z component)
    inline const int i2jz(const int& i) const {return 3*i+2;}

};

// *****************************************************************************
} // namespace ks
