/** @file   SimulationFactory.cpp
 *  @brief  Factory functions to build a rebound Simulation for the planets or DE435 collection.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-10
 * 
 */

// *****************************************************************************
// Local dependencies
#include "rebound_utils.hpp"
#include "Simulation.hpp"
    using ks::reb::Simulation;

// *****************************************************************************
namespace ks {
namespace reb {

// *****************************************************************************
// Factory functions to build a simulation with the planets
// *****************************************************************************

// *****************************************************************************
Simulation make_sim_planets(const PlanetVector& pv, double epoch)
{
    // Test input date vs. date range
    if ((epoch < mjd0_db) || (epoch > mjd1_db))
    {
        string msg = format("make_sim_planets(PlanetVector): bad input date {:8.2f}!. "
                            "Epoch must be bewteen {:d} and {:d}.\n", epoch, mjd0_db, mjd1_db);
        throw invalid_argument(msg);
    }

    // Start with an empty simulation
    Simulation sim;
    // Set the time field in the simulation to match the epoch
    sim.set_time(epoch);
    // Load table of masses
    MassiveBodyTable mbt = MassiveBodyTable();

    // Iterate through the bodies in this collection
    // Add a particle for each body with initial conditions based on interpolated state vector on PlanetVector
    for (int i=0; i<N_body_planets; i++)
    {
        // The body_id and primary_id of this body
        int32_t body_id = body_ids_planets[i];
        int32_t primary_body_id = get_primary_body_id(body_id);
        // Get the state vector of this particle at the epoch
        StateVector s = pv.interp_vec(body_id, epoch);
        // Look up the mass on the MassiveBodyTable
        double m = mbt.get_M(body_id);
        // Add this particle
        sim.add_particle(s, m, body_id, primary_body_id);
    }

    // Return the assembled simulation
    return sim;
}

// *****************************************************************************
Simulation make_sim_planets(const PlanetElement& pe, double epoch)
{
    // Test input date vs. date range
    if ((epoch < mjd0_db) || (epoch > mjd1_db))
    {
        string msg = format("make_sim_planets(PlanetElement): bad input date {:8.2f}!. "
                            "Epoch must be bewteen {:d} and {:d}.\n", epoch, mjd0_db, mjd1_db);
        throw invalid_argument(msg);
    }

    // Start with an empty simulation
    Simulation sim;
    // Set the time field in the simulation to match the epoch
    sim.set_time(epoch);
    // Load table of masses
    MassiveBodyTable mbt = MassiveBodyTable();
    // The body_ids; add them to the simulation in Jacobi order rather than sorted by body_id
    constexpr int32_t body_ids[] = {10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9};

    // Iterate through the bodies in this collection
    // Add a particle for each body with initial conditions based on interpolated state vector on PlanetElement
    for (int i=0; i<N_body_planets; i++)
    {
        // The body_id and primary_id of this body
        int32_t body_id = body_ids[i];
        int32_t primary_body_id = get_primary_body_id(body_id);
        // Get the state vector of this particle at the epoch
        StateVector s = pe.interp_vec(body_id, epoch);
        // Look up the mass on the MassiveBodyTable
        double m = mbt.get_M(body_id);
        // Add this particle
        sim.add_particle(s, m, body_id, primary_body_id);
    }

    // Return the assembled simulation
    return sim;
}

// *****************************************************************************
Simulation make_sim_horizons(db_conn_type& conn, const string collection, int epoch)
{
    // Test input date vs. date range
    if ((epoch < mjd0_hrzn) || (epoch > mjd1_hrzn))
    {
        string msg = format("make_sim_horizons: bad input date {:d}!. Epoch must be bewteen {:d} and {:d}.\n",
                                epoch, mjd0_hrzn, mjd1_hrzn);
        throw invalid_argument(msg);
    }
    // Test input date vs. stride
    if ((collection == "DE435") && (epoch%stride_hrzn_de435 != 0))
    {
        string msg = format(
            "make_sim_horizons: bad input date {:d}!. For collection DE435, epoch must be a multiple of {:d}.\n",
            epoch, stride_hrzn_de435);
        throw invalid_argument(msg);
    }

    // Start with an empty simulation
    Simulation sim;
    // Convert integer epoch to a floating point time
    double t = static_cast<double>(epoch);
    // Set the time field in the simulation to match the epoch
    sim.set_time(t);

    // Run the stored procedure to get state vectors of the planets from Horizons
    string sp_name = "JPL.GetHorizonsStateCollection";
    // Run the stored procedure to get JPL state vectors for the named collection
    vector<string> params = {wrap_string(collection), to_string(epoch)};
    // DEBUG
    print("{:s}({:s}, {:s}).\n", sp_name, params[0], params[1]);
    ResultSet* rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset
        // The body ID and name
        int32_t body_id = rs->getInt("BodyID");
        string body_name = rs->getString("BodyName");
        // The mass
        double m = rs->getDouble("m");
        // Six state vector components
        double qx = rs->getDouble("qx");
        double qy = rs->getDouble("qy");
        double qz = rs->getDouble("qz");
        double vx = rs->getDouble("vx");
        double vy = rs->getDouble("vy");
        double vz = rs->getDouble("vz");
        // Wrap state vector
        StateVector s {.qx=qx, .qy=qy, .qz=qz, .vx=vx, .vy=vy, .vz=vz};
        // Primary body of this particle
        int32_t primary_body_id = get_primary_body_id(body_id);
        // Add this particle
        sim.add_particle(s, m, body_id, primary_body_id);
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;

    // Return the assembled simulation
    return sim;
}

// *****************************************************************************
Simulation make_sim_planets_horizons(db_conn_type& conn, int epoch)
    {return make_sim_horizons(conn, "Planets", epoch);}
    
// *****************************************************************************
Simulation make_sim_de435_horizons(db_conn_type& conn, int epoch)
    {return make_sim_horizons(conn, "DE435", epoch);}

// *****************************************************************************
} // Namespace reb
} // Namespace ks
