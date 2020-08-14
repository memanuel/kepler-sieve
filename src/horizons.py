"""
Harvard IACS Masters Thesis
Utilities for building Rebound simulations using NASA JPL Horizons data

Michael S. Emanuel
06-Aug-2020
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Utility
from datetime import datetime
import collections

# MSE
from astro_utils import mjd_to_jd
from db_config import db_engine

# Typing
from typing import List, Dict, Set

# ********************************************************************************************************************* 
# Gravitational constant from JPL
G_jpl = 2.95912208285591095E-04

horizon_entry = collections.namedtuple('horizon_entry', 
    ['BodyID', 'BodyName', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz',])    

# ********************************************************************************************************************* 
def quote(s):
    """Wrap a string s in single quotes"""
    return f"\'{s}\'"

# ********************************************************************************************************************* 
def get_hrzn_state_coll(body_collection: str, epoch: int) -> pd.DataFrame:
    """
    Get state vectors from Horizons for all bodies in the named collection as of the given epoch.
    INPUTS:
        body_collection: String name of a collection of bodies on DB table KS.BodyCollection, e.g.
                         'Planets', 'DE-435', 'DE-435-Top-16'
        epoch:           MJD as of which state vectors are taken.
                         Must be an INTEGER in the supported date range 40400 (1969-06-28) to 77600 (2071-05-04)
    RETURNS:
        states:          Pandas DataFrame including columns BodyID, BodyName, m, qx, qy, qz, vx, vy, vz
    """
    # Assemble SQL to call the stored procedure JPL.GetHorizonsStateCollection
    body_collection = quote(body_collection)
    sql = f"CALL JPL.GetHorizonsStateCollection({body_collection}, {epoch});"

    # Run the SP and wrap results as a DataFrame
    with db_engine.connect() as conn:
        states = pd.read_sql(sql, con=conn)

    return states

# ********************************************************************************************************************* 
def make_sim_horizons(body_collection: str, epoch: int) -> rebound.Simulation:
    """
    Create a new rebound simulation with initial data from the NASA Horizons system
    INPUTS:
        body_collection: String name of a collection of bodies on DB table KS.BodyCollection, e.g.
                         'Planets', 'DE435', 'DE435-Top-16'
        epoch:           MJD as of which state vectors are taken.
                         Must be an INTEGER in the supported date range 40400 (1969-06-28) to 77600 (2071-05-04)
    RETURNS:
        sim:             Rebound simulation initialized with these state vectors
    """
    # Create an empty simulation
    sim = rebound.Simulation()
    
    # Set units
    sim.units = ('day', 'AU', 'Msun')
    
    # Set gravitational constant
    sim.G = G_jpl

    # Get state vectors
    states = get_hrzn_state_coll(body_collection=body_collection, epoch=epoch)

    # When initializing a horizons simulation from a body collection, we ALWAYS treat bodies as massive
    add_as_test: bool = False

    # Set additional attributes on sim to store BodyID and BodyName field for each particle
    sim.body_ids = np.array([], dtype=np.int64)
    sim.body_names = np.array([], dtype='object')

    # Add bodies in this collection to the empty simulation
    add_hrzn_bodies(sim=sim, states=states, add_as_test=add_as_test)

    # Set extra attribute with epoch
    sim.epoch = epoch

    return sim

# ********************************************************************************************************************* 
def add_hrzn_bodies(sim: rebound.Simulation, states: pd.DataFrame, add_as_test: bool) -> None:
    """
    Add all the bodies in a collection of state vectors to a Rebound simulation.
    INPUTS:
        sim:         Rebound simulation these bodies will be added to
        states:      DataFrame including columns BodyName, m, qx, qy, qz, vx, vy, vz
                     Units are day, AU, Msun
        add_as_test: Flag indicating whether bodies added as massless test particles or not
    RETURNS:
        None.  Modifies the input simulation, sim, in place by adding the bodies.
    """

    # Number of bodies in the collection
    n = states.shape[0]

    # Set mass multiplier based on add_as_test
    mass_multiplier: float = 0.0 if add_as_test else 1.0

    # Generate set of object hashes present before we start; used to avoid adding duplicates
    hashes_present: Set[int] = set(p.hash.value for p in sim.particles) if sim.particles else set()
    
    # Mask indicating whether each body is added (true for non-duplicates, false for duplicates)
    mask = np.ones(shape=n, dtype=bool)

    # Iterate over each body in the collection
    for i in range(n):
        # Get the ith body's state vector
        s = states.iloc[i]

        # Particle name and rebound hash of that name
        body_name: str = s.BodyName
        hash_id: int = rebound.hash(body_name).value

        # If this particle is already present, don't add a duplicate
        if hash_id in hashes_present:
            mask[i] = False
            continue

        # Mass to use depends on add_as_test
        m: float = s.m * mass_multiplier

        # Add this particle with the given state vector
        sim.add(m=m, x=s.qx, y=s.qy, z=s.qz, vx=s.vx, vy=s.vy, vz=s.vz, hash=hash_id)

        # Add the names and IDs of this bodies
        np.append(sim.body_ids,  s.BodyID)
        np.append(sim.body_names,  s.BodyName)

    # # Add the names and IDs of only those bodies that were not duplicates
    sim.body_ids = np.concatenate([sim.body_ids, states.BodyID.values[mask]])
    sim.body_names = np.concatenate([sim.body_names, states.BodyName.values[mask]])

# ********************************************************************************************************************* 
def extend_sim_horizons(sim: rebound.Simulation, body_names: List[str], add_as_test: bool) -> None:
    """
    Extend a rebound simulation with initial data from the NASA Horizons system
    INPUTS:
        sim:         A rebound simulation object
        body_names:  List of string body names; references DB table KS.Body, field name BodyName
        add_as_test: Flag indicating whether bodies added as massless test particles or not
    RETURNS:
        None:       Modifies sim in place
    """
    # Get the epoch from the simulation
    epoch: int = sim.epoch

    # Generate list of missing object names
    hashes_present: Set[int] = set(p.hash.value for p in sim.particles)
    bodies_missing: List[str] = [nm for nm in body_names if rebound.hash(nm).value not in hashes_present]
    is_body_missing: bool = (len(bodies_missing) > 0)
    
    # Quit early if there are no missing bodies
    if not is_body_missing:
        return

    # Current number of total particles
    N_active = sim.N

    # Add missing objects one at a time if not already present
    with db_engine.connect() as conn:
        for body_name in bodies_missing:
            # Look up the state vector using stored procedure JPL.GetHorizonsState()
            body_name_quoted = quote(body_name)
            sql = f"CALL JPL.GetHorizonsState({body_name_quoted}, {epoch});"
            state = pd.read_sql(sql, con=conn)
            # Add this body to the simulation
            add_hrzn_bodies(sim, states=state, add_as_test=add_as_test)

    # Update number of active particles if necessary
    if add_as_test:
        sim.N_active = N_active

# ********************************************************************************************************************* 
def extend_sim_horizons_ast(sim: rebound.Simulation, n_ast: int, add_as_test: bool) -> None:
    """
    Extend a rebound simulation with initial data from the NASA Horizons system
    INPUTS:
        sim:         A rebound simulation object
        n_ast:       The number of asteroids to be added; the first n_ast asteroids are added to the simulation
        add_as_test: Flag indicating whether bodies added as massless test particles or not
    RETURNS:
        None:       Modifies sim in place
    """
    # Get the epoch from the simulation
    epoch: int = sim.epoch

    # Current number of total particles
    N_active = sim.N

    # Add missing objects one at a time if not already present
    with db_engine.connect() as conn:
        # Look up the state vector using stored procedure JPL.GetHorizonsStateAsteroids()
        sql = f"CALL JPL.GetHorizonsStateAsteroids({n_ast}, {epoch});"
        states = pd.read_sql(sql, con=conn)
        # Add this body to the simulation
        add_hrzn_bodies(sim, states=states, add_as_test=add_as_test)

    # Update number of active particles if necessary
    if add_as_test:
        sim.N_active = N_active
