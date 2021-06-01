"""
Tools for working with asteroid orbital elements.

Functions in this module:
get_asteroids(key_to_body_id)
get_ast_ref_elts_jpl(epoch)
get_ast_ref_elts(epoch, n0, n1, missing)
add_asteroid_elts(sim, elts)
update_asteroid_elements(sim, elts, epoch)
make_sim_asteroids_ref(epoch, n0, n1, missing)
get_ast_elts(n0, n1, epoch)
get_ast_elts_ts(n0, n1, mjd0, mjd1)
get_ast_data(n0, n1, epoch)
make_sim_asteroids(epoch, n0, n1)

Michael S. Emanuel
2021-02-10
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Local imports
from rebound_sim import make_sim_planets
from db_utils import sp2df

# ********************************************************************************************************************* 
def get_asteroids(key_to_body_id: bool=False) -> pd.DataFrame:
    """
    Return list of known asteroid names and IDs
    INPUTS:
        key_to_body_id: When true, reindex this DataFrame to BodyID.  Default (false) keys to AsteroidID.
    """
    ast: pd.DataFrame = sp2df(sp_name='KS.GetAsteroids')
    if key_to_body_id:
        ast.set_index(keys='BodyID', drop=False, inplace=True)
    else:
        ast.set_index(keys='AsteroidID', drop=False, inplace=True)
    return ast

# ********************************************************************************************************************* 
# Work with reference orbital elements for asteroids as of one date provided by JPL
# Enrich them to calculated orbital elements on many dates that are saved in the database
# ********************************************************************************************************************* 

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
def get_ast_ref_elts(epoch: int, n0: int, n1: int, missing: bool) -> pd.DataFrame:
    """
    Get calculated reference orbital elements as of the given epoch.
    INPUTS:
        epoch:      The epoch as of which to get the reference elements
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        missing:    Flag - true (default) only returns elements for asteroids that haven't already been integrated.  
                    False returns all available elements.
    OUTPUTS:
        elts:   DataFrame with orbital elements calculated directly from the JPL provided elements.
    """
    # Select SP name based on the missing flag
    sp_name = 'KS.GetAsteroidRefElementsMissing' if missing else 'KS.GetAsteroidRefElements'
    # Assemble input parameters (same for both stored procedures)
    params = {
        'epoch': epoch,
        'n0': n0,
        'n1': n1
    }
    # Get the elements on this epoch from the database
    elts = sp2df(sp_name=sp_name, params=params)    
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
    # Does the elements DataFrame include the optional AsteroidName field?
    # If not, put in dummy asteroid names
    if 'AsteroidName' not in elts.columns:
        elts['AsteroidName'] = 'Asteroid.'
        elts['AsteroidName'] = elts.AsteroidName + elts.AsteroidID.astype(str)

    # Does the elements DataFrame include the optional BodyID field?
    # If not, populate it
    if 'BodyID' not in elts.columns:
        elts['BodyID'] = elts.AsteroidID + 1000000

    # Does the elements DataFrame include the optional BodyName field?
    # If not, put in dummy body names
    if 'BodyName' not in elts.columns:
        elts['BodyName'] = 'SB.'
        elts['BodyName'] = elts.BodyName + elts.AsteroidName

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
        # Save the name of this asteroid to the particle entry (hack)
        sim.particles[-1].name = elts.AsteroidName[i]

    # Get the ID columns: AsteroidID, BodyID and BodyName
    asteroid_ids = elts.AsteroidID.values
    body_ids = elts.BodyID.values
    body_names = elts.BodyName.values

    # Save the asteroid IDs
    sim.asteroid_ids = asteroid_ids
    # Append the body IDs and names
    sim.body_ids = np.concatenate([sim.body_ids, body_ids])
    sim.body_names = np.concatenate([sim.body_names, body_names])

    # Return the new simulation including the asteroid IDs
    return asteroid_ids

# ********************************************************************************************************************* 
def update_asteroid_elements(sim: rebound.Simulation, elts: pd.DataFrame, epoch: int) -> None:
    """
    Get updated orbital elements from sim and apply them to elts DataFrame.
    INPUTS:
        sim:   Rebound simulation object
        elts:  DataFrame of orbital elements
        epoch: New epoch (after integration) as of which elements are quoted
    """
    # Calculate all orbits; primary is the Sun
    orbits = sim.calculate_orbits(primary=sim.particles[0])

    # Number of asteroids and active particles
    N: int = elts.shape[0]
    N_active: int = sim.N_active

    # The TimeID at this epoch
    TimeID = np.int32(np.rint(epoch * 24 * 60))

    # Loop through asteroids in the simulation
    for i in range(N):
        # The index j in the orbits is offset by number of active particles-1; 
        # the minus 1 is because the primary (sun) has no elements
        j: int = N_active + i - 1
        # The orbit of this asteroid
        orb = orbits[j]
        # Updated time variables
        elts.loc[i, 'TimeID'] = TimeID
        elts.loc[i, 'epoch'] = epoch
        # Updated orbital elements
        elts.loc[i, 'a'] = orb.a
        elts.loc[i, 'e'] = orb.e
        elts.loc[i, 'inc'] = orb.inc
        elts.loc[i, 'Omega'] = orb.Omega
        elts.loc[i, 'omega'] = orb.omega
        elts.loc[i, 'f'] = orb.f
        elts.loc[i, 'M'] = orb.M

# ********************************************************************************************************************* 
def make_sim_asteroids_ref(epoch: int, n0: int, n1: int, missing: bool):
    """
    Create a simulation with the planets and a block of asteroids as of an epoch.
    INPUTS:
        epoch:  Epoch as of which the planets and asteroids are initialized
        n0:     n0 - AsteroidID of the first asteroid to add (inclusive)
        n1:     n1 - AsteroidID of the last asteroid to add (exclusive)
        missing: When true, only add elements for asteroids missing from AsteroidVectors or AsteroidElements
    """
    # Initialize a planets simulation as of epoch
    sim = make_sim_planets(epoch=epoch, load_file=True)
    sim.t = np.float64(epoch)
    sim.N_active = sim.N

    # Get the orbital elements of this slice of asteroids; only request m
    elts = get_ast_ref_elts(epoch=epoch, n0=n0, n1=n1, missing=missing)
    # Add these asteroids with their orbital elements
    asteroid_ids = add_asteroid_elts(sim=sim, elts=elts)

    return sim
    
# ********************************************************************************************************************* 
# Get calculated orbital elements and wrap them into a simulation
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def get_ast_elts(n0: int, n1: int, epoch: int) -> pd.DataFrame:
    """
    Get calculated reference orbital elements as of the given epoch.
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        epoch:      Date on which elements are requested
    OUTPUTS:
        elts:       DataFrame with orbital elements from the saved MSE integration
    """
    # The back end stored procedure
    sp_name = 'KS.GetAsteroidElements'
    # Convert the epoch to a date range that will return just the one date if it's available
    mjd0: int = epoch
    mjd1: int = epoch+1

    # Assemble input parameters (same for both stored procedures)
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    # Get the elements on this epoch from the database
    elts = sp2df(sp_name=sp_name, params=params)    
    return elts

# ********************************************************************************************************************* 
def get_ast_elts_ts(n0: int, n1: int, mjd0: int, mjd1: int) -> pd.DataFrame:
    """
    Get calculated orbital elements as of the given epoch.
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        mjd0:       First date on which to return data (inclusive)
        mjd1:       Last date on which to return data (exclusive)
    OUTPUTS:
        elts:       DataFrame with orbital elements from the saved MSE integration
    """
    # The back end stored procedure
    sp_name = 'KS.GetAsteroidElements'
    # Assemble input parameters (same for both stored procedures)
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    # Get the elements in this date range from the database
    elts = sp2df(sp_name=sp_name, params=params)    
    return elts

# ********************************************************************************************************************* 
def get_ast_data(n0: int, n1: int, epoch: int) -> pd.DataFrame:
    """
    Get both calculated state vectors and orbital elements as of the given epoch.
    INPUTS:
        n0:         First asteroid number to return (inclusive)
        n1:         Last asteroid number to return (exclusive)
        epoch:      Date on which elements are requested
    OUTPUTS:
        elts:       DataFrame with both state vectors and orbital elements from the saved MSE integration
    """
    # The back end stored procedure
    sp_name = 'KS.GetAsteroidData'
    # Convert the epoch to a date range that will return just the one date if it's available
    mjd0: int = epoch
    mjd1: int = epoch+1

    # Assemble input parameters (same for both stored procedures)
    params = {
        'n0': n0,
        'n1': n1,
        'mjd0': mjd0,
        'mjd1': mjd1,
    }
    # Get the elements on this epoch from the database
    elts = sp2df(sp_name=sp_name, params=params)    
    return elts

# ********************************************************************************************************************* 
def make_sim_asteroids(epoch: int, n0: int, n1: int):
    """
    Create a simulation with the planets and a block of asteroids as of an epoch.
    Available dates are in the range mjd= 48000-63000 at multiples of 4.
    INPUTS:
        epoch:  Epoch as of which the planets and asteroids are initialized
        n0:     n0 - AsteroidID of the first asteroid to add (inclusive)
        n1:     n1 - AsteroidID of the last asteroid to add (exclusive)
    """
    # Initialize a planets simulation as of epoch
    sim = make_sim_planets(epoch=epoch, load_file=True)
    sim.t = np.float64(epoch)
    sim.N_active = sim.N

    # Get the orbital elements of this slice of asteroids on the epoch
    elts = get_ast_elts(n0=n0, n1=n1, epoch=epoch)
    # Add these asteroids with their orbital elements
    asteroid_ids = add_asteroid_elts(sim=sim, elts=elts) 

    return sim
