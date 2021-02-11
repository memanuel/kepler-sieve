"""
Tools for working with asteroid orbital elements.

Michael S. Emanuel
2021-02-10
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Local imports
from db_utils import sp2df

# ********************************************************************************************************************* 
def get_ast_ref_elts_jpl(epoch: int) -> pd.DataFrame:
    """
    Get reference orbital elements, as quoted by JPL, on the given epoch.
    INPUTS:
        epoch:  The epoch as of which to get the reference elements
    OUTPUTS:
        elts:   DataFrame with orbital elements exactly as provided by JPL.
    """
    # Get the elements on this epoch from the database
    elts = sp2df(sp_name='JPL.GetAsteroidRefElements', params={'epoch':epoch})
    # Return assembled DataFrame
    return elts

# ********************************************************************************************************************* 
def get_ast_ref_elts(epoch: int) -> pd.DataFrame:
    """
    Get calculated reference orbital elements as of the given epoch.
    INPUTS:
        epoch:  The epoch as of which to get the reference elements
    OUTPUTS:
        elts:   DataFrame with orbital elements calculated directly from the JPL provided elements.
    """
    # Get the elements on this epoch from the database
    elts = sp2df(sp_name='KS.GetAsteroidRefElements', params={'epoch':epoch})
    # Return assembled DataFrame
    return elts

# ********************************************************************************************************************* 
def add_asteroid_elts(sim: rebound.Simulation, elts: pd.DataFrame) -> np.array:
    """
    Add asteroids with the provided orbital elements to a simulation.
    INPUTS:
        sim:    The simulation the asteroid will be added to.
        elts:   DataFrame of orbital elements
    OUTPUTS:
        asteroid_ids: Numpy array of asteroid IDs from KS.Asteroid table
    Modifies sim in place by adding the desired asteroids.
    """

    # Add the asteroids one at a time
    ast_count: int = elts.shape[0]
    for i in range(ast_count):
        # Unpack the orbital elements
        a = elts.a[i]
        e = elts.e[i]
        inc = elts.inc[i]
        Omega = elts.Omega[i]
        omega = elts.omega[i]
        M = elts.M[i]
        # name = elts.AsteroidName[i]
        # Set the primary to the sun (NOT the solar system barycenter!)
        primary = sim.particles['Sun']
        # Add the new asteroid
        sim.add(m=0.0, a=a, e=e, inc=inc, Omega=Omega, omega=omega, M=M, primary=primary)
        # Set the hash to the asteroid's AsteroidID
        sim.particles[-1].hash = int(elts.AsteroidID[i])
        # Save the name of this asteroid to the particle entry (hack)
        sim.particles[-1].name = elts.AsteroidName[i]

    # The corresponding list of asteroid IDs and names
    asteroid_ids = elts.AsteroidID.values

    # Return the new simulation including the asteroid IDs
    return asteroid_ids

# ********************************************************************************************************************* 
def update_asteroid_elements(sim, elts, epoch) -> None:
    """
    Get updated orbital elements from sim and apply them to elts DataFrame.
    INPUTS:
        sim:   Rebound simulation object
        elts:  DataFrame of orbital elements
        epoch: New epoch (after integration) as of which elements are quoted
    """
    # Calculate all orbits; primary is the Sun
    sim.calculate_orbits(primary=sim.particles[0])

    # Number of asteroids and active particles
    N: int = elts.shape[0]
    N_active: int = sim.N_active

    # The TimeID at this epoch
    TimeID = np.int32(np.rint(epoch * 24 * 60))

    # Loop through asteroids in the simulation
    for i in range(N):
        # The index j in the simulation is offset by number of active particles
        j: int = N_active + i
        # Updated time variables
        elts.loc[i, 'TimeID'] = TimeID
        elts.loc[i, 'epoch'] = epoch
        # Updated orbital elements
        elts.loc[i, 'a'] = sim.particles[j].a
        elts.loc[i, 'e'] = sim.particles[j].e
        elts.loc[i, 'inc'] = sim.particles[j].inc
        elts.loc[i, 'Omega'] = sim.particles[j].Omega
        elts.loc[i, 'omega'] = sim.particles[j].omega
        elts.loc[i, 'f'] = sim.particles[j].f
        elts.loc[i, 'M'] = sim.particles[j].M

