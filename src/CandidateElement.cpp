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

// // One BodyVector for the Sun is shared
// BodyVector bv_sun {BodyVector("Sun")};

// // One BodyVector for the Earth is shared
// BodyVector bv_earth {BodyVector("Earth")};

// // One dectection time table is shared
// DetectionTimeTable dtt = DetectionTimeTable();


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
    mjd   {new double[N_t]},
    q_obs {new double[N_row]},
    q_ast {new double[N_row]},
    v_ast {new double[N_row]},
    u_ast {new double[N_row]},
    q_sun {new double[N_row]},
    v_sun {new double[N_row]},
    q_cal {new double[N_row]},
    v_cal {new double[N_row]}
    {}

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id, const double* mjd_, int N_t):
    // Delegate to main constructor using mjd array from detection time table
    CandidateElement(elt, candidate_id, N_t)
{
    // Initialize mjd array with input times
    init(mjd_);
}

// *****************************************************************************
CandidateElement::CandidateElement(OrbitalElement elt, int32_t candidate_id):
    // Delegate to main constructor using mjd array from detection time table
    CandidateElement(elt, candidate_id, dtt.get_mjd(), dtt.N())
{
    // Initialize mjd array from detection time table
    init(dtt.get_mjd());
}

// *****************************************************************************
// Delete manually allocated arrays
CandidateElement::~CandidateElement()
{
    delete [] mjd;
    delete [] q_obs;
    delete [] q_ast;
    delete [] v_ast;
    delete [] u_ast;
    delete [] q_sun;
    delete [] v_sun;
    delete [] q_cal;
    delete [] v_cal;
}

// *****************************************************************************
void CandidateElement::init(const double* mjd_)
{
    // Copy from mjd_ input to mjd on the candidate element
    for (int i=0; i<N_t; i++) {mjd[i]=mjd_[i];}

    // // Initialize q_obs with Earth center as a placeholder
    // // This will be overwritten with the exact observatory position later
    // for (int i=0; i<N_t; i++)
    // {
    //     // Interpolated state vector of the Sun at time i
    //     Position p = bv_earth.interp_pos(mjd[i]);
    //     // Array base into q and v
    //     int j = 3*i;
    //     // Copy position components
    //     q_obs[j+0] = p.qx;  q_obs[j+1] = p.qy;  q_obs[j+2] = p.qz;
    // }   // for / i
    // // DEBUG
    // print("CandidateElement constructor: populated q_obs array.\n");

    // // Populate q_sun and v_sun from BodyVector of the Sun
    // for (int i=0; i<N_t; i++)
    // {
    //     // Interpolated state vector of the Sun at time i
    //     StateVector s = bv_sun.interp_vec(mjd[i]);
    //     // Array base into q and v
    //     int j = 3*i;
    //     // Copy position components
    //     q_sun[j+0] = s.qx;  q_sun[j+1] = s.qy;  q_sun[j+2] = s.qz;
    //     // Copy velocity components
    //     v_sun[j+0] = s.vx;  v_sun[j+1] = s.vy;  v_sun[j+2] = s.vz;
    // }   // for / i
    // // DEBUG
    // print("CandidateElement constructor: populated q_sun, v_sun arrays.\n");

    // Intitialize calibration adjustments to zero
    for (int i=0; i<N_row; i++)
    {
        q_cal[i] = 0.0;
        v_cal[i] = 0.0;
    }   // for / i
    // DEBUG
    print("CandidateElement constructor: intialized q_cal, v_cal arrays.\n");

}   // CandidateElement::init

// *****************************************************************************
// Calculate and calibrate trajectory in Kepler model
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calc_trajectory(bool with_calibration)
{
    // Alias input flag
    bool cal = with_calibration;
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
        double t = mjd[i] - mjd0;
        // Mean anomaly at this time; changes at rate n and equal to elt.M at time elt.mjd
        double M = M0 + n*t;
        // Convert M to a true anomaly f
        double f = anomaly_M2f(M, e);
        // Convert orbital elements to a heliocentric position
        StateVector s = elt2vec(a, e, inc, Omega, omega, f, mu_sun);
        // Index for arrays q_ast and v_ast
        int j = 3*i;
        // Save into array q_ast with calibration adjustment if requested
        q_ast[j+0] = s.qx + (cal ? q_cal[j+0] : 0.0);
        q_ast[j+1] = s.qy + (cal ? q_cal[j+1] : 0.0);
        q_ast[j+2] = s.qz + (cal ? q_cal[j+2] : 0.0);
        // Save into array v_ast with calibration adjustment if requested
        v_ast[j+0] = s.vx + (cal ? v_cal[j+0] : 0.0);
        v_ast[j+1] = s.vy + (cal ? v_cal[j+1] : 0.0);
        v_ast[j+2] = s.vz + (cal ? v_cal[j+2] : 0.0);
        // DEBUG
        // print("calc_trajectory: i={:d}, mjd={:8.2f}, s= ", i, mjd[i]);
        // print_state_vector(s);
        // print("q_ast[j+0]={:8.2f}\n", q_ast[j+0]);
    }
}

// *****************************************************************************
void CandidateElement::calibrate(const PlanetVector& pv)
{
    // Get the reference time for the elements
    double mjd0 = elt.mjd;

    // Build rebound simulation with planets at reference time
    Simulation sim = make_sim_planets(pv, mjd0);
    // Add asteroid to simulation with candidate elements
    sim.add_test_particle(elt, candidate_id);
    // Write the numerically integrated vectors from this simulation into q_cal and v_cal
    // sim.write_vectors(mjd, candidate_id, N_t, q_cal, v_cal);

    // Calculate the trajectory without calibration
    calc_trajectory(false);

    // Iterate through all the rows; subtract the Kepler component from the numerical component
    // This is equivalent to writing
    // q_cal[k] = numerical_trajectory[k] - kepler_trajectory[k]
    // It just uses less memory by using q_cal to store the numerical trajectory and subtracting in place
    for (int k=0; k<N_row; k++)
    {
        q_cal[k] -= q_ast[k];
        v_cal[k] -= v_ast[k];
    }
}

// *****************************************************************************
// Get predicted state vectors
// *****************************************************************************

// *****************************************************************************
const StateVector CandidateElement::state_vector(int i) const
{
    // The array index for the test date
    int j = 3*i;

    // Wrap the predicted position and velocity of the asteroid on the test date into a StateVector
    return StateVector 
    {
        .qx = q_ast[j+0], 
        .qy = q_ast[j+1],
        .qz = q_ast[j+2],
        .vx = v_ast[j+0],
        .vy = v_ast[j+1],
        .vz = v_ast[j+2]
    };
}

// *****************************************************************************
// Calculate direction
// *****************************************************************************

// *****************************************************************************
void CandidateElement::calc_direction()
{
    ;
}

