"""
Assemble orbital elements for known asteroids Pandas DataFrame.
Main interface used by consumers is load_ast_elt()

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Library imports
import numpy as np
import pandas as pd
import rebound

# Utility
from datetime import datetime
from tqdm import tqdm as tqdm_console

# MSE
from astro_utils import anomaly_M2E, anomaly_E2f

# Types
from typing import List, Tuple, Dict, Optional

# Radians in a circle
tau = 2.0 * np.pi

# Gravitational constant in units of day, AU, Msun
# sim = rebound.Simulation()
# sim.units = ('day', 'AU', 'Msun')
# G_ = sim.G
G_ = 0.00029591220828559104
mu = G_

# ********************************************************************************************************************* 
def load_data_numbered() -> pd.DataFrame:
    """Load the asteroid data for numbered asteroids into a Pandas DataFrame"""
    # The source for this file is at https://ssd.jpl.nasa.gov/?sb_elem
    fname: str = '../data/jpl/orbital_elements/asteroid_numbered.txt'

    # The field names in the JPL file and their column positions
    names: List[str] = ['Num', 'Name', 'Epoch', 'a', 'e', 'i', 'w', 'Node', 'M', 'H', 'G', 'Ref']
    colspec_tbl: Dict[str, Tuple[int, int]] = {
        'Num': (0,6), 
        'Name': (7, 25), 
        'Epoch': (25, 30), 
        'a': (31, 41), 
        'e': (42, 52), 
        'i': (54, 62), 
        'w': (63, 72),
        'Node': (73, 82),
        'M': (83, 94),
        'H': (95, 100),
        'G': (101, 105),
        'Ref': (106, 113),
    }
    
    # Other arguments for Pandas file import
    colspecs: List[Tuple[int, int]] = [colspec_tbl[nm] for nm in names]
    header: int = 0
    skiprows: List[int] = [1]
    dtype: Dict[str, int] = {
        'Num': int,
        'Name': str,
        'Epoch': float,
        'a': float,
        'e': float,
        'i': float,
        'w': float,
        'Node': float,
        'M': float,
        'H': float,
        'G': float,
        'Ref': str,
    }

    # Read the DataFrame
    df: pd.DataFrame = pd.read_fwf(fname, colspecs=colspecs, header=header, names=names, skiprows=skiprows, dtype=dtype)
    # Set the asteroid number field to be the index
    df.set_index(keys=['Num'], drop=False, inplace=True)
    return df

# ********************************************************************************************************************* 
def load_data_unnumbered() -> pd.DataFrame:
    """Load the asteroid data for unnumbered asteroids into a Pandas DataFrame"""

    fname: str = '../data/jpl/orbital_elements/asteroid_unnumbered.txt'

    # The field names in the JPL file and their column positions
    names: List[str] = ['Name', 'Epoch', 'a', 'e', 'i', 'w', 'Node', 'M', 'H', 'G', 'Ref']
    colspec_tbl: Dict[str, Tuple[int, int]] = {
        'Name': (0, 12), 
        'Epoch': (12, 17), 
        'a': (18, 28), 
        'e': (29, 39), 
        'i': (40, 50), 
        'w': (51, 60),
        'Node': (61, 70),
        'M': (71, 82),
        'H': (83, 88),
        'G': (89, 93),
        'Ref': (94, 104),
    }

    # Other arguments for Pandas file import
    colspecs: List[Tuple[int, int]] = [colspec_tbl[nm] for nm in names]
    header: int = 0
    skiprows: List[int] = [1]
    dtype: Dict[str, int] = {
        'Name': str,
        'Epoch': float,
        'a': float,
        'e': float,
        'i': float,
        'w': float,
        'Node': float,
        'M': float,
        'H': float,
        'G': float,
        'Ref': str,
    }

    # Read the DataFrame
    df: pd.DataFrame = pd.read_fwf(fname, colspecs=colspecs, header=header, names=names, skiprows=skiprows, dtype=dtype)
    # Populate the asteroid_num field and add it to DataFrame
    ast_num_offset = 1000001
    ast_num = df.index.values + ast_num_offset
    df.insert(loc=0, column='Num', value=ast_num)
    # Set the asteroid number field to be the index
    df.set_index(keys=['Num'], drop=False, inplace=True)
    return df

# ********************************************************************************************************************* 
def load_data_impl() -> pd.DataFrame:
    """Load the combined asteroid data into a Pandas DataFrame"""
    # Load main data file
    df1 = load_data_numbered()
    # Load auxiliary data file
    df2 = load_data_unnumbered()
    # Return combined DataFrame
    df = pd.concat([df1, df2])
    # Field indicating whether asteroid is IAU numbered or not
    df['IsNumberedAsteroid'] = (df.Num <= 1000000)
    # Rename columns for asteroid number and name to MSE spec
    mapper = {
        'Num': 'AsteroidNumber',
        'Name': 'AsteroidName',
    }
    df.rename(columns=mapper, inplace=True)
    # Filter to remove invalid entries with eccentricity outside of [0, 1)
    mask = (0 <= df.e) & (df.e < 1.0)
    return df[mask]

# ********************************************************************************************************************* 
def convert_data(df_in: pd.DataFrame, epoch: Optional[float]=None) -> pd.DataFrame:
    """
    Convert data from the JPL format to be friendly to rebound integrator and matching selected epoch
    INPUTS:
        df_in: DataFrame with orbital elements in JPL format
        epoch: Optional epoch as an MJD; used to filter for only matching epochs.
               If not specified, take the whole DataFrame, which will have different epochs
    """
    # Create a mask with only the matching rows if epoch was specified
    mask = (df_in.Epoch == epoch) if epoch is not None else np.ones_like(df_in.Epoch, dtype=bool)
    
    # Initialize Dataframe with asteroid numbers, names, and IsNumberedAsteroid flag
    columns = ['AsteroidNumber', 'AsteroidName', 'IsNumberedAsteroid']
    df = pd.DataFrame(data=df_in[mask][columns])

    # Add fields one at a time
    df['epoch'] = df_in.Epoch[mask]
    df['a'] = df_in.a[mask]
    df['e'] = df_in.e[mask]
    df['inc'] = np.radians(df_in.i[mask]) % tau
    df['Omega'] = np.radians(df_in.Node[mask]) % tau
    df['omega'] = np.radians(df_in.w[mask]) % tau
    df['M'] = np.radians(df_in.M[mask]) % tau
    df['H'] = df_in.H[mask]
    df['G'] = df_in.G[mask]
    df['Ref'] = df_in.Ref[mask]

    # Set the asteroid number field to be the index
    df.set_index(keys=['AsteroidNumber'], drop=False, inplace=True)

    # Return the newly assembled DataFrame
    return df

# ********************************************************************************************************************* 
def ast_data_add_calc_elements_reb(ast_elt) -> pd.DataFrame:
    """
    Add the true anomaly and other calculated orbital elements to the asteroid DataFrame
    Older version - delegates all computations to Rebound.
    """   
    # Number of asteroids
    N: int = len(ast_elt)

    # Initialize empty arrays for computed orbital elements
    f = np.zeros(N)
    P = np.zeros(N)
    mean_motion = np.zeros(N)
    long = np.zeros(N)
    theta = np.zeros(N)
    pomega = np.zeros(N)
    T_peri = np.zeros(N)

    # Base Rebound simulation with just the Sun
    # We are NOT integrating this simulation, only using it to convert orbital elements
    # Therefore we don't need the planets, or the initial configuration of the sun; just its mass.
    sim = rebound.Simulation()
    sim.units = ('day', 'AU', 'Msun')
    sim.add(m=1.0)
    sim.N_active = 1
        
    # All the available asteroid numbers; wrap as a tqdm iterator for progress bar
    # nums = ast_elt.index
    nums = ast_elt.Num.values
    iters = tqdm_console(nums)

    # Make a gigantic simulation with all these asteroids
    print(f'Making big simulation with all {N} asteroids...')
    for num in iters:
        # Unpack the orbital elements
        a = ast_elt.a[num]
        e = ast_elt.e[num]
        inc = ast_elt.inc[num]
        Omega = ast_elt.Omega[num]
        omega = ast_elt.omega[num]
        M = ast_elt.M[num]
        # Set the primary to the sun (NOT the solar system barycenter!)
        primary = sim.particles[0]
        # Add the asteroid to the simulation as a massless test particle.  Just want the elements!
        sim.add(m=0.0, a=a, e=e, inc=inc, Omega=Omega, omega=omega, M=M, primary=primary)
       
    # Calculate orbital elements for all particles; must specify primary = Sun!!!
    print(f'Computing orbital elements...')
    orbits = sim.calculate_orbits(primary=sim.particles[0])
    
    # Iterate over all the asteroids in the simulation
    print(f'Copying additional orbital elements to DataFrame...')
    iters = list(enumerate(nums))
    for i, num in tqdm_console(iters):
        # Look up the orbit of asteroid i
        orb = orbits[i]
        # Unpack the additional (calculated) orbital elements
        f[i] = orb.f
        P[i] = orb.P
        mean_motion[i] = orb.n
        long[i] = orb.l
        theta[i] = orb.theta
        pomega[i] = orb.pomega
        T_peri[i] = orb.T

    # Save computed orbital elements to the DataFrame
    ast_elt['f'] = f % tau
    ast_elt['P'] = P
    ast_elt['n'] = mean_motion
    ast_elt['long'] = long % tau
    ast_elt['theta'] = theta % tau
    ast_elt['pomega'] = pomega % tau
    ast_elt['T_peri'] = T_peri

    # Return the updated DataFrame 
    return ast_elt

# ********************************************************************************************************************* 
def ast_data_add_calc_elements(ast_elt) -> pd.DataFrame:
    """
    Add the true anomaly and other calculated orbital elements to the asteroid DataFrame
    New version.  Uses fast, mathematical calculations only; no dependence on Rebound library.
    """   
    # Extract input orbital elements: mean anomaly eccentricity, semi-major axis, 
    M = ast_elt.M.values
    e = ast_elt.e.values
    a = ast_elt.a.values
    Omega = ast_elt.Omega.values
    omega = ast_elt.omega.values    

    # Compute eccentric anomaly E from mean anomaly M
    E = anomaly_M2E(M=M, e=e)

    # Compute true anomaly f from eccentric anomaly E
    f = anomaly_E2f(E=E, e=e)

    # Compute the orbital period P in days
    P = tau * np.sqrt(a**3 / mu)

    # Mean motion in radians per day
    mean_motion = tau / P

    # Mean longitude
    long = Omega + omega + M

    # True longitude 
    theta = Omega + omega + f

    # Longitude of pericenter
    # pomega = 0.0

    # Save computed orbital elements to the DataFrame
    ast_elt['eccentric_anomaly'] = E
    ast_elt['f'] = f
    ast_elt['period'] = P
    ast_elt['mean_motion'] = mean_motion
    ast_elt['long'] = long % tau
    ast_elt['theta'] = theta % tau

    # Return the updated DataFrame 
    return ast_elt

# ********************************************************************************************************************* 
def load_ast_elt() -> pd.DataFrame:
    """Load the asteroid orbital elements data into a Pandas Dataframe"""
    # The name for the saved DataFrame
    fname: str = '../data/jpl/orbital_elements/orb_elements_asteroid.h5'
    
    # Try to load from disk if available
    ast_elt: pd.DataFrame
    try:
        ast_elt = pd.read_hdf(fname, key='ast_elt')
    except:
        # Load data from JPL asteroids file
        df_in = load_data_impl()
        # Convert data to rebound format
        ast_elt = convert_data(df_in=df_in)
        # Add the calculated orbital elements
        ast_elt = ast_data_add_calc_elements(ast_elt)
        # Add the row number field
        ast_elt['row_num'] = np.arange(ast_elt.shape[0], dtype=np.int32)
        # Intuitive column order
        cols = ['AsteroidNumber', 'AsteroidName', 'IsNumberedAsteroid', 
        'epoch', 'a', 'e', 'inc', 'Omega', 'omega', 
        'f', 'M', 'eccentric_anomaly', 'period', 'mean_motion',
        'H', 'G', 'Ref', 'row_num']
        ast_elt = ast_elt[cols]
        # Save orbital elements dataframe to h5
        ast_elt.to_hdf(fname, key='ast_elt', mode='w')
    
    return ast_elt
