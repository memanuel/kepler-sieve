/** @file CandidateElement.cpp
 *  @brief Implmentation of CandidateElement class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-07
 */

// *****************************************************************************
// Local Dependencies
#include "CandidateElement.hpp"

// *****************************************************************************
// Local names used
using ks::CandidateElement;

// *****************************************************************************
// Constructor & destructor
// *****************************************************************************

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id, int N_t): 
    // Small data elements
    elt {elt},
    candidate_id {candidate_id},
    elt0 {elt},
    N_t {N_t},
    N_row {3*N_t},
    // Allocate arrays
    mjd_    {new double[N_t]},
    q_obs_  {new double[N_row]},
    q_ast_  {new double[N_row]},
    v_ast_  {new double[N_row]},
    u_ast_  {new double[N_row]},
    r_ast_  {new double[N_t]},
    dq_ast_ {new double[N_row]},
    dv_ast_ {new double[N_row]},
    /// Read-only copies of arrays in public interface
    mjd     {const_cast<double*> (mjd_) },
    u_ast   {const_cast<double*> (u_ast_) },
    r_ast   {const_cast<double*> (r_ast_) }
    {}

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id, const double* const mjd_in, int N_t):
    // Delegate to main constructor using mjd_ array from detection time table
    CandidateElement(elt, candidate_id, N_t)
{
    // Initialize mjd_ array with input times
    init(mjd_in);
}

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id):
    // Delegate to constructor taking an mjd array and size; use mjd array from detection time table
    CandidateElement(elt, candidate_id, dtt.mjd+1, dtt.N)
{}

// *****************************************************************************
// Delete manually allocated arrays
CandidateElement::~CandidateElement()
{
    delete [] mjd_;
    delete [] q_obs_;
    delete [] q_ast_;
    delete [] v_ast_;
    delete [] u_ast_;
    delete [] r_ast_;
    delete [] dq_ast_;
    delete [] dv_ast_;
}

// *****************************************************************************
void CandidateElement::init(const double* mjd_in)
{
    // DEBUG
    print("CandidateElement::init(elt, candidate_id).\n");
    print("N_t = {:d}.\n", N_t);

    // Copy from mjd_in to mjd_ on the candidate element
    for (int i=0; i<N_t; i++) {mjd_[i]=mjd_in[i];}

    // Populate q_cal_ and v_cal_ from BodyVector of the Sun
    // This will be a pretty accurate first pass before running an optional numerical integration
    for (int i=0; i<N_t; i++)
    {
        // The time of this entry
        double t = mjd_[i];
        // Interpolated position vector of the Earth at time i
        Position p_earth = bv_earth.interp_pos(t);
        // Interpolated state vector of the Sun at time i
        StateVector s_sun = bv_sun.interp_vec(t);
        // Index for arrays q_sun and v_sun
        const int jx {i2jx(i)};
        const int jy {i2jy(i)};
        const int jz {i2jz(i)};
        // Copy position components into q_obs
        q_obs_[jx] = p_earth.qx;
        q_obs_[jy] = p_earth.qy;
        q_obs_[jz] = p_earth.qz;
        // Copy position components into dq_ast
        dq_ast_[jx] = s_sun.qx;
        dq_ast_[jy] = s_sun.qy;
        dq_ast_[jz] = s_sun.qz;
        // Copy velocity components into dv_ast
        dv_ast_[jx] = s_sun.vx;
        dv_ast_[jy] = s_sun.vy;
        dv_ast_[jz] = s_sun.vz;
    }   // for / i
}   // CandidateElement::init

// *****************************************************************************
// Calculate and calibrate trajectory in Kepler model
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calc_trajectory()
{
    // Calculate the mean motion, n
    double n = mean_motion(elt.a, mu_sun);
    // Unpack the five fields of the orbital element that don't change
    double a = elt.a;
    double e = elt.e;
    double inc = elt.inc;
    double Omega = elt.Omega;
    double omega = elt.omega;
    // Get the reference date and mean anomaly for the elements
    double mjd0 = elt.mjd;
    double M0 = elt.M;

    // Iterate through the mjds in the array
    for (int i=0; i<N_t; i++)
    {
        // The time shift from the reference time
        double dt = mjd_[i] - mjd0;
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*dt;
        // Convert M to a true anomaly f
        double f = anomaly_M2f(M, e);
        // Convert orbital elements to a heliocentric position
        StateVector s = elt2vec(a, e, inc, Omega, omega, f, mu_sun);
        // Index for arrays q_ast_ and v_ast_
        const int jx {i2jx(i)};
        const int jy {i2jy(i)};
        const int jz {i2jz(i)};
        // Save into array q_ast_
        q_ast_[jx] = s.qx + dq_ast_[jx];
        q_ast_[jy] = s.qy + dq_ast_[jy];
        q_ast_[jz] = s.qz + dq_ast_[jz];
        // Save into array v_ast_
        v_ast_[jx] = s.vx + dv_ast_[jx];
        v_ast_[jy] = s.vy + dv_ast_[jy];
        v_ast_[jz] = s.vz + dv_ast_[jz];
    }
}

// *****************************************************************************
void CandidateElement::calibrate()
{
    // Get the reference time for the elements
    double mjd0 = elt.mjd;

    // Zero out old calibration
    for (int k=0; k<N_row; k++)
    {
        dq_ast_[k] = 0.0;
        dv_ast_[k] = 0.0;
    }

    // Calculate the trajectory without calibration
    calc_trajectory();

    // Build rebound simulation with planets at reference time
    Simulation sim = make_sim_planets(pv, mjd0);
    // Add asteroid to simulation with candidate elements
    sim.add_test_particle(elt, candidate_id);

    // Write the numerically integrated vectors from this simulation into q_cal and v_cal
    sim.write_vectors(dq_ast_, dv_ast_, mjd_, N_t, candidate_id);

    // Iterate through all the rows; subtract the Kepler component from the numerical component
    // This is equivalent to writing
    // dq_ast[k] = numerical_trajectory[k] - kepler_trajectory[k]
    // It just uses less memory by using dq_ast to store the numerical trajectory and subtracting in place
    for (int k=0; k<N_row; k++)
    {
        // double q_num = dq_ast[k];
        dq_ast_[k] -= q_ast_[k];
        dv_ast_[k] -= v_ast_[k];
        // print("k={:d}, q_num = {:+10.6f}, q_kep = {:+10.6f}, dq = {:+10.6f}.\n", k, q_num, q_ast_[k], dq_ast_[k]);
    }
}

// *****************************************************************************
// Get predicted state vectors and direction
// *****************************************************************************

// *****************************************************************************
const StateVector CandidateElement::state_vector(int i) const
{
    // The array index for the test date
    const int jx {i2jx(i)};
    const int jy {i2jy(i)};
    const int jz {i2jz(i)};

    // Wrap the predicted position and velocity of the asteroid on the test date into a StateVector
    return StateVector 
    {
        .qx = q_ast_[jx], 
        .qy = q_ast_[jy],
        .qz = q_ast_[jz],
        .vx = v_ast_[jx],
        .vy = v_ast_[jy],
        .vz = v_ast_[jz]
    };
}

// *****************************************************************************
const Position CandidateElement::observer_pos(int i) const
{
    // Wrap the observer position on the test date into a Position
    return Position 
    {
        .qx = q_obs_[i2jx(i)], 
        .qy = q_obs_[i2jy(i)],
        .qz = q_obs_[i2jz(i)]
    };
}

// Extract a Direction to the asteroid from the u_ast_ array
const Direction CandidateElement::direction(int i) const
{
    return Direction
    {
        .ux = u_ast_[i2jx(i)],
        .uy = u_ast_[i2jy(i)],
        .uz = u_ast_[i2jz(i)]
    };
}


// *****************************************************************************
// Calculate direction
// *****************************************************************************

// *****************************************************************************
// Calculate direction based on observer position q_obs and state vector q_ast, v_ast
// Caller is responsible to populate these before running calc_direction()
void CandidateElement::calc_direction()
{
    // Iterate over all N observation times
    for (int i=0; i<N_t; i++)
    {
        // Wrap the observer position
        Position q_obs {observer_pos(i)};
        // Wrap the asteroid state vector
        StateVector s_ast {state_vector(i)};
        // Calculate direction and distance in Newtonian light model
        ObservationResult res {observe(q_obs, s_ast)};      
        // Write to direction array
        u_ast_[i2jx(i)] = res.u.ux;
        u_ast_[i2jy(i)] = res.u.uy;
        u_ast_[i2jz(i)] = res.u.uz;
        // Write to distance array
        r_ast_[i]           = res.r;
    }
}

// *****************************************************************************
// Initialize the static member PlanetVector
// *****************************************************************************

// Shared body vector objects
BodyVector CandidateElement::bv_sun {BodyVector(SolarSystemBody_bv::sun)};
BodyVector CandidateElement::bv_earth {BodyVector(SolarSystemBody_bv::earth)};

// Build throwaway copy of DetectionTimeTable to get date range
DetectionTimeTable dtt {DetectionTimeTable()};

// Inputs to build a PlanetVector table suited to CandidateElement
constexpr int pad = 32;
const int mjd_first = static_cast<int>(floor(dtt.mjd_first()));
const int mjd_last = static_cast<int>(ceil(dtt.mjd_last()));
const int mjd0 = std::min(mjd_first, 58000) - pad;
const int mjd1 = std::max(mjd_last,  60000) + pad;
constexpr int dt_min = 60;
bool load = true;
PlanetVector CandidateElement::pv {PlanetVector(mjd0, mjd1, dt_min, load)};

// Shared copy of DetectionTimeTable
DetectionTimeTable CandidateElement::dtt {DetectionTimeTable()};

// DEBUG - build a small Detection table quickly for testing
// DetectionTable CandidateElement::dt = DetectionTable()
int d1 = 85000000;
int d2 = 86000000;
DetectionTable CandidateElement::dt = DetectionTable(d1, d2, load);
