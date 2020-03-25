"""
Harvard IACS Masters Thesis
candidate_elements.py: 
Generate candidate orbital elements for asteroid search

Michael S. Emanuel
Wed Mar 25 09:55 2020
"""

# Core
import numpy as np
import pandas as pd

# Astronomy
import rebound

# Local
from asteroid_integrate import load_ast_elt
from planets import make_sim_planets
from astro_utils import mjd_to_datetime

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
def orbital_element_batch(ast_nums: np.ndarray) -> pd.DataFrame:
    """
    Return a batch of orbital elements as a DataFrame
    INPUTS:
        ast_nums: Numpy array of asteroid numbers to include in batch
    OUTPUTS:
        elts: DataFrame with columns for a, e, inc, Omega, omega, f, epoch
    """
    # The orbital elements and epoch
    dtype = np.float32
    a = ast_elt.a[ast_nums].values.astype(dtype)
    e = ast_elt.e[ast_nums].values.astype(dtype)
    inc = ast_elt.inc[ast_nums].values.astype(dtype)
    Omega = ast_elt.Omega[ast_nums].values.astype(dtype)
    omega = ast_elt.omega[ast_nums].values.astype(dtype)
    f = ast_elt.f[ast_nums].values.astype(dtype)
    epoch = ast_elt.epoch_mjd[ast_nums].to_numpy().astype(dtype)
    
    # Wrap into dictionary
    elts_dict = {
        # 'asteroid_num': ast_nums,
        'element_id': ast_nums,
        'a': a,
        'e': e,
        'inc': inc,
        'Omega': Omega,
        'omega': omega,
        'f': f,
        'epoch': epoch
    }

    # Convert dict to DataFrame
    elts = pd.DataFrame(elts_dict)
    return elts

# ********************************************************************************************************************* 
def orbital_element_batch_by_ast_num(n0: int, batch_size: int=64):
    """
    Return a batch of orbital elements for asteroids in a batch of consecutive asteroid numbers.
    DEPRECATED.
    INPUTS:
        n0: first asteroid number, e.g. 1
        batch_size: number of asteroids in batch
    OUTPUTS:
        elts: Dictionary with seven keys for a, e, inc, Omega, omega, f, epoch
    """
    # Get start and end index location of this asteroid number
    i0: int = ast_elt.index.get_loc(n0)
    i1: int = i0 + batch_size
    ast_nums = np.arange(i0, i1+1, dytpe=np.int32)

    # Delegate to orbital_element_batch
    return orbital_element_batch(ast_nums)

# ********************************************************************************************************************* 
def perturb_elts(elts: pd.DataFrame, sigma_a=0.05, sigma_e=0.10, 
                 sigma_f_deg=5.0, sigma_Omega_deg=0.0, sigma_omega_deg=0.0,
                 mask_pert=None, random_seed: int = 42):
    """Apply perturbations to orbital elements"""
    # Copy the elements
    elts_new = elts.copy()

    # Default for mask_pert is all elements
    if mask_pert is None:
        mask_pert = np.ones_like(elts['a'], dtype=bool)

    # Number of elements to perturb
    num_shift = np.sum(mask_pert)

    # Set random seed
    np.random.seed(seed=random_seed)

    # Apply shift log(a)
    log_a = np.log(elts['a'])
    log_a[mask_pert] += np.random.normal(scale=sigma_a, size=num_shift)
    elts_new['a'] = np.exp(log_a)
    
    # Apply shift to log(e)
    log_e = np.log(elts['e'])
    log_e[mask_pert] += np.random.normal(scale=sigma_e, size=num_shift)
    elts_new['e'] = np.exp(log_e)
    
    # Apply shift directly to true anomaly f
    f = elts['f']
    sigma_f = np.deg2rad(sigma_f_deg)
    f[mask_pert] += np.random.normal(scale=sigma_f, size=num_shift)
    elts_new['f'] = f
    
    # Apply shift directly to Omega
    Omega = elts['Omega']
    sigma_Omega = np.deg2rad(sigma_Omega_deg)
    Omega[mask_pert] += np.random.normal(scale=sigma_Omega, size=num_shift)
    elts_new['Omega'] = Omega

    # Apply shift directly to omega
    omega = elts['omega']
    sigma_omega = np.deg2rad(sigma_omega_deg)
    omega[mask_pert] += np.random.normal(scale=sigma_omega, size=num_shift)
    elts_new['omega'] = omega

    return elts_new

# ********************************************************************************************************************* 
def random_elts(element_id_start: np.int32 = 0, 
                size: np.int32 = 64, 
                random_seed: np.int32 = 42):
    """
    Generate a DataFrame of random orbital elements
    INPUTS:
        element_id_start: First element_id used to label these elements
        size:             Number of elements to draw
        random_seed:      Random seed for the elements
    OUTPUTS:
        elts:        DataFrame of random orbital elemnents.
                     Columns include element_id, a, e, inc, Omega, omega, f, epoch
    """
    # Set random state
    np.random.seed(random_seed)

    # Randomly sample a, e, inc, Omega from empirical observations
    a = np.random.choice(ast_elt.a, size=size)
    e = np.random.choice(ast_elt.e, size=size)
    inc = np.random.choice(ast_elt.inc, size=size)
    Omega = np.random.choice(ast_elt.Omega, size=size)

    # Sample mean anomaly M and omega randomly
    two_pi = 2.0*np.pi
    M = np.random.uniform(low=0.0, high=two_pi, size=size)
    omega = np.random.uniform(low=0.0, high=two_pi, size=size)

    # Allocate array for the true anomaly f
    f = np.zeros(size)

    # The epoch
    epoch_mjd = ast_elt.epoch_mjd.values[0]

    # Epoch as a datetime
    epoch_dt = mjd_to_datetime(epoch_mjd)

    # Base Rebound simulation of the planets and moons on this date
    sim = make_sim_planets(epoch=epoch_dt)
    # Set the number of active particles to the base simulation
    sim.N_active = sim.N

    # Iterate over candidate elements
    for i in range(size):
            # Set the primary to the sun (NOT the solar system barycenter!)
            # Put this inside the loop b/c not guaranteed to remain constant as particles are added
            primary = sim.particles['Sun']
            # Add the new asteroid
            sim.add(m=0.0, a=a[i], e=e[i], inc=inc[i], Omega=Omega[i], omega=omega[i], M=M[i], primary=primary)
            # The asteroid that was just added
            ast = sim.particles[-1]
            # Extract the true anomaly f
            f[i] = ast.f
            # Set the hash to the asteroid's number in this batch
            ast.hash = rebound.hash(f'{i}')

    # The element_id and epoch arrays
    element_id = np.arange(element_id_start, size, dtype=np.int32)
    epoch = np.full(shape=size, fill_value=epoch_mjd)

    # Elements as a Python dict with required columns in order
    elts_dict = {
        'element_id': element_id,
        'a': a,
        'e': e,
        'f': f,
        'inc': inc,
        'Omega': Omega,
        'omega': omega,
        'epoch': epoch
    }

    # Convert these arrays to a DataFrame
    elts = pd.DataFrame(elts_dict)

    return elts
