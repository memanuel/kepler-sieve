"""
Harvard IACS Masters Thesis
Utilities for building Rebound simulations using NASA JPL Horizons data

Michael S. Emanuel
06-Aug-2020
"""

# Library imports
import pandas as pd
import rebound
from datetime import datetime
import collections
from typing import List, Dict, Set

# Local imports
from astro_utils import mjd_to_jd
from db_config import db_engine

# ********************************************************************************************************************* 
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
def add_hrzn_bodies(sim: rebound.Simulation, states: pd.DataFrame):
    """
    Add all the bodies in a collection of state vectors to a Rebound simulation.
    INPUTS:
        sim: Rebound simulation these particles will be added to
        states: DataFrame including columns BodyName, m, qx, qy, qz, vx, vy, vz
                Units are day, AU, Msun
    """

    # Number of bodies in the collection
    n = states.shape[0]
    
    # Iterate over each body in the collection
    for i in range(n):
        # Get the ith body's state vector
        s = states.iloc[i]

        # Particle name and rebound hash of that name
        body_name = s.BodyName
        hash_id = rebound.hash(body_name)

        # Add this particle with the given state vector
        sim.add(m=s.m, x=s.qx, y=s.qy, z=s.qz, vx=s.vx, vy=s.vy, vz=s.vz, hash=hash_id)

# ********************************************************************************************************************* 
def make_sim_horizons(body_collection: str, epoch: int) -> rebound.Simulation:
    """
    Create a new rebound simulation with initial data from the NASA Horizons system
    INPUTS:
        body_collection: String name of a collection of bodies on DB table KS.BodyCollection, e.g.
                         'Planets', 'DE-435', 'DE-435-Top-16'
        epoch:           MJD as of which state vectors are taken.
                         Must be an INTEGER in the supported date range 40400 (1969-06-28) to 77600 (2071-05-04)
    RETURNS:
        sim:             Rebound simulation initialized with these state vectors
    """
    # Create an empty simulation
    sim = rebound.Simulation()
    
    # Set units
    sim.units = ('day', 'AU', 'Msun')
    
    # Get state vectors
    states = get_hrzn_state_coll(body_collection=body_collection, epoch=epoch)

    # Add bodies in this collection to the empty simulation
    add_hrzn_bodies(sim, states)

    # Move to center of mass
    sim.move_to_com()
    
    return sim

# ********************************************************************************************************************* 
# def extend_sim_horizons(sim: rebound.Simulation, body_names: List[str], epoch: int) -> None:
#     """Extend a rebound simulation with initial data from the NASA Horizons system"""
#     # Generate list of missing object names
#     hashes_present: Set[int] = set(p.hash.value for p in sim.particles)
#     bodies_missing: List[str] = [nm for nm in body_names if rebound.hash(nm).value not in hashes_present]
    
#     # Add missing objects one at a time if not already present
#     for body_name in bodies_missing:
#         add_one_object_hrzn(sim=sim, object_name=object_name, epoch_dt=epoch_dt, hrzn=hrzn)
    
#     # Move to center of mass
#     sim.move_to_com()
