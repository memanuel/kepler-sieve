/** @file   rebound_utils.cpp
 *  @brief  Implementation of rebound_utils.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-10
 * 
 */

// *****************************************************************************
// Local dependencies
#include "rebound_utils.hpp"

// *****************************************************************************
namespace ks {


// *****************************************************************************
Simulation* make_sim()
{
    // Create an empty simulation
    Simulation* s = reb_create_simulation();

    // Configure simulation
    // Unit system is (AU, Msun, day)

    // Set initial time step to 1.0 days
    s->dt = 1.0;            // in days

    // Set gravitational constant in this unit system, i.e. in AU^3 / Msun / day^2
    s->G  = G;

    // Use the IAS15 integrator
    s->integrator = reb_simulation::REB_INTEGRATOR_IAS15;

    // Don't set a heartbeat for now
    // s->heartbeat        = heartbeat;

    // Finish exactly at tmax in reb_integrate(). Default is already 1.
    s->exact_finish_time = 1;

    // Return the assembled simulation object
    return s;
}

// *****************************************************************************
} // Namespace ks
