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
from scipy.interpolate import PchipInterpolator
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF

# Astronomy
import rebound

# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local
from asteroid_element import load_ast_elt
from planets import make_sim_planets
from astro_utils import mjd_to_datetime, deg2dist, dist2deg

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
# Default data type
dtype = np.float64

# ********************************************************************************************************************* 
# Convert between different representations for orbital elements: Numpy table, dictionary and DataFrame
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def elts_np2dict(elts_np):
    """Convert an Nx7 array of orbital elements into a dict"""
    # Dictionary
    elts_dict = {
        'a': elts_np[:,0],
        'e': elts_np[:,1],
        'inc': elts_np[:,2],
        'Omega': elts_np[:,3],
        'omega': elts_np[:,4],
        'f': elts_np[:,5],
        'epoch': elts_np[:,6],
    }
    return elts_dict

# ********************************************************************************************************************* 
def elts_np2df(elts_np):
    """Convert an Nx7 array of orbital elements into a DataFrame"""
    # Dictionary
    elts_dict = elts_np2dict(elts_np=elts_np)
    # Return a DataFrame
    return pd.DataFrame(elts_dict)

# ********************************************************************************************************************* 
def elts_df2dict(elts_df):
    """Convert a DataFrame of orbital elements into a dict (built-in Dataframe.to_dict() method fails)"""
    # Columns in the elements DataFrame
    cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
    # Return a dict
    return {col: elts_df[col] for col in cols_elt}

# ********************************************************************************************************************* 
# Generate candidate elements from known asteroids (including perturbations)
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def asteroid_elts(ast_nums: np.ndarray, dtype = dtype) -> pd.DataFrame:
    """
    Return a batch of orbital elements as a DataFrame
    INPUTS:
        ast_nums: Numpy array of asteroid numbers to include in batch
        score_by_elt: DataFrame generated by ztf_score_by_element in ztf_element.py
        R_deg:    Initial value of resolution parameter (in degrees) in mixture model
        dtype:    Data type for the DataFrame.
    OUTPUTS:
        elts:     DataFrame with columns for a, e, inc, Omega, omega, f, epoch.
                  Also includes h and R for mixture model
    """
    # The orbital elements and epoch
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
def perturb_elts(elts: pd.DataFrame, 
                 sigma_a=0.00, sigma_e=0.00, sigma_inc_deg=0.0,
                 sigma_f_deg=1.0, sigma_Omega_deg=0.0, sigma_omega_deg=0.0,
                 mask_pert=None, random_seed: int = 42):
    """Apply perturbations to orbital elements"""
    # Copy the elements; overwrite the original reference to prevent accidentally changing them
    # This can happen due to separate handles to the same numpy array! Caused an obscure bug.
    elts = elts.copy()

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
    elts['a'] = np.exp(log_a)
    
    # Apply shift to log(e)
    log_e = np.log(elts['e'])
    log_e[mask_pert] += np.random.normal(scale=sigma_e, size=num_shift)
    elts['e'] = np.exp(log_e)
    
    # Apply shift directly to inclination inc
    inc = elts['inc']
    sigma_inc = np.deg2rad(sigma_inc_deg)
    inc[mask_pert] += np.random.normal(scale=sigma_inc, size=num_shift)
    elts['inc'] = inc
    
    # Apply shift directly to true anomaly f
    f = elts['f']
    sigma_f = np.deg2rad(sigma_f_deg)
    f[mask_pert] += np.random.normal(scale=sigma_f, size=num_shift)
    elts['f'] = f
    
    # Apply shift directly to Omega
    Omega = elts['Omega']
    sigma_Omega = np.deg2rad(sigma_Omega_deg)
    Omega[mask_pert] += np.random.normal(scale=sigma_Omega, size=num_shift)
    elts['Omega'] = Omega

    # Apply shift directly to omega
    omega = elts['omega']
    sigma_omega = np.deg2rad(sigma_omega_deg)
    omega[mask_pert] += np.random.normal(scale=sigma_omega, size=num_shift)
    elts['omega'] = omega

    return elts

# ********************************************************************************************************************* 
# Generate candidate elements by randomly sampling elements of known asteroids
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def random_elts(element_id_start: np.int32 = 0, 
                size: np.int32 = 64, 
                random_seed: np.int32 = 42,
                dtype = dtype):
    """
    Generate a DataFrame of random orbital elements
    INPUTS:
        element_id_start: First element_id used to label these elements
        size:             Number of elements to draw
        random_seed:      Random seed for the elements
        dtype:            Data type for the DataFrame.
    OUTPUTS:
        elts:        DataFrame of random orbital elemnents.
                     Columns include element_id, a, e, inc, Omega, omega, f, epoch.
                     Also includes h and R for mixture model.
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
        'a': a.astype(dtype),
        'e': e.astype(dtype),
        'f': f.astype(dtype),
        'inc': inc.astype(dtype),
        'Omega': Omega.astype(dtype),
        'omega': omega.astype(dtype),
        'epoch': epoch.astype(dtype)
    }

    # Convert these arrays to a DataFrame
    elts = pd.DataFrame(elts_dict)

    return elts

# ********************************************************************************************************************* 
# Augment candidate orbital elements with mixture model parameters
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def elts_add_hit_rate(elts: pd.DataFrame, score_by_elt: pd.DataFrame, num_hits: int = 10):
    """
    Populate the hit rate h by guessing each element has num_hits hits
    INPUTS:
        elts: DataFrame with columns a, e, inc, Omega, omega, f, epoch
        score_by_elt: DataFrame generated by ztf_score_by_element in ztf_element.py
        num_hits:     Number of hits to assume for computing the initial guess on hit rate h
    OUTPUTS:
        Modifies elts in place
    """
    # Filter down to the intersection
    is_match = elts.element_id.isin(score_by_elt.index)
    # elts = elts[is_match]
    idx_missing = elts.index[~is_match]
    elts.drop(idx_missing, inplace=True)
    
    # Update column h
    elts['h'] = num_hits / score_by_elt.num_obs.values

# ********************************************************************************************************************* 
def elts_add_R_deg(elts: pd.DataFrame, R_deg: float, dtype=dtype):
    """
    Add two columns to a DataFrame of orbital elements for the mixture parameters h and R.
    INPUTS:
        elts: DataFrame with columns a, e, inc, Omega, omega, f, epoch
        R:    Resolution parameter (Cartesian distance)
    OUTPUTS:
        Modifies elts in place
    """
    # Number of asteroids in this batch
    N_ast = elts.shape[0]

    # Convert R from degrees to Cartesian
    R = deg2dist(R_deg)

    # Add column for R
    elts['R'] = np.full(fill_value=R, shape=N_ast, dtype=dtype)

# ********************************************************************************************************************* 
def elts_add_mixture_params(elts: pd.DataFrame, score_by_elt: pd.DataFrame, num_hits: int, R_deg: float, dtype=dtype):
    """
    Add two columns to a DataFrame of orbital elements for the mixture parameters h and R.
    INPUTS:
        elts: DataFrame with columns a, e, inc, Omega, omega, f, epoch
        score_by_elt: DataFrame generated by ztf_score_by_element in ztf_element.py
        num_hits:     Number of hits to assume for computing the initial guess on hit rate h
        R_deg:        Resolution parameter in degrees
    OUTPUTS:
        Modifies elts in place
    """
    # Number of asteroids in this batch
    N_ast = elts.shape[0]

    # Convert R from degrees to Cartesian
    R = deg2dist(R_deg)

    # Add columns for h and R
    elts_add_hit_rate(elts=elts, score_by_elt=score_by_elt, num_hits=num_hits)
    elts_add_R_deg(elts=elts, R_deg=R_deg)

# ********************************************************************************************************************* 
# Transform orbital elements to U = Unif[0, 1] or Z = Norm(0, 1)
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_interp_x(x: np.ndarray, cdf_samp: np.ndarray, output_type: str='z'):
    """
    Build a PCHIP interpolator that takes as inputs x and returns a variable z that is standard normal.
    INPUTS:
        x: Values whose CDF is to be taken, e.g. log(a)
        cdf_samp: Desired CDF levels on the spline
        output_type: Either 'z' for standard normal or 'u' or standard uniform 
    OUTPUT:
        x_to_z: a PCHIP interpolator mapping x to standard normal z
    """
    # Build empirical CDF for x
    ecdf_x = ECDF(x)
    
    # Unique observations of x
    x_unq = np.unique(x)
    
    # CDF of unique observations of x
    x_cdf = ecdf_x(x_unq)
    
    # Interpolator from CDF to x
    cdf_to_x = PchipInterpolator(x=x_cdf, y=x_unq)
    
    # x at the sampled CDF levels
    x_samp = cdf_to_x(cdf_samp)
    
    # Take unique items again for edge case with e.g. sin(Omega)
    x_samp_unq, idx = np.unique(x_samp, return_index=True)

    # z at the sampled CDF levels
    u_samp_unq = cdf_samp[idx]
    z_samp_unq = norm.ppf(u_samp_unq)

    # The selected output
    output_tbl = {
        'u': u_samp_unq,
        'z': z_samp_unq, 
    }
    y_samp_unq = output_tbl[output_type]
    
    # Interpolator from x to z
    x_to_z = PchipInterpolator(x=x_samp_unq, y=y_samp_unq)
    
    return x_to_z

# ********************************************************************************************************************* 
def ast_elt_transform(ast_elt: pd.DataFrame):
    """
    Build interpolators from orbital elements to uniform z and populate DataFrame of transformed elements
    """

    # Set number of sample points
    N_samp_u: int = 2**16
    z_range: float = 6.0
    N_samp_z: int = int(200*z_range + 1)

    # Sample CDF levels: N_samp evenly spaced by U
    cdf_samp_u = (np.arange(N_samp_u) + 0.5) / N_samp_u

    # Sample CDF levels: N_samp_z points evenly spaced by Z
    z_samp = np.linspace(-z_range, z_range, N_samp_z)
    cdf_samp_z = norm.cdf(z_samp)
    pdf_samp_z = norm.pdf(z_samp)

    # Combine the two sets of sample points
    cdf_samp = np.unique(np.hstack([cdf_samp_u, cdf_samp_z]))    
    
    # Build interpolator from a to z
    log_a = np.log(ast_elt.a.values)
    log_a_to_z = make_interp_x(x=log_a, cdf_samp=cdf_samp, output_type='z')
    # Compute transformed values
    log_a_z = log_a_to_z(log_a)
    
    # Build interpolator from e to z
    e = ast_elt.e.values
    e_to_z = make_interp_x(x=e, cdf_samp=cdf_samp, output_type='z')
    # Compute transformed values
    e_z = e_to_z(e)
    
    # Build interpolator from sin(inc) to z
    sin_inc = np.sin(ast_elt.inc.values)
    sin_inc_to_z = make_interp_x(x=sin_inc, cdf_samp=cdf_samp, output_type='z')
    # Compute transformed values
    sin_inc_z = sin_inc_to_z(sin_inc)
    
    # Build interpolator from sin(Omega) and cos(Omega) to z
    sin_Omega = np.sin(ast_elt.Omega.values)
    cos_Omega = np.cos(ast_elt.Omega.values)
    sin_Omega_to_z = make_interp_x(x=sin_Omega, cdf_samp=cdf_samp, output_type='z')
    cos_Omega_to_z = make_interp_x(x=cos_Omega, cdf_samp=cdf_samp, output_type='z')
    # Compute transformed values
    sin_Omega_z = sin_Omega_to_z(sin_Omega)
    cos_Omega_z = cos_Omega_to_z(cos_Omega)
    
    # Build interpolator from sin(omega) and cos(omega) to z
    sin_omega = np.sin(ast_elt.omega.values)
    cos_omega = np.cos(ast_elt.omega.values)
    sin_omega_to_z = make_interp_x(x=sin_omega, cdf_samp=cdf_samp, output_type='z')
    cos_omega_to_z = make_interp_x(x=cos_omega, cdf_samp=cdf_samp, output_type='z')    
    # Compute transformed values
    sin_omega_z = sin_omega_to_z(sin_omega)
    cos_omega_z = cos_omega_to_z(cos_omega)
    
    # Build interpolator from sin(f) and cos(f) to z
    sin_f = np.sin(ast_elt.Omega.values)
    cos_f = np.cos(ast_elt.Omega.values)
    sin_f_to_z = make_interp_x(x=sin_f, cdf_samp=cdf_samp, output_type='z')
    cos_f_to_z = make_interp_x(x=cos_f, cdf_samp=cdf_samp, output_type='z')
    # Compute transformed values
    sin_f_z = sin_f_to_z(sin_f)
    cos_f_z = cos_f_to_z(cos_f)
    
    # Dictionary of interpolators
    interp_tbl = {
        'log_a': log_a_to_z,
        'e': e_to_z,
        'sin_inc': sin_inc_to_z,
        'sin_Omega': sin_Omega_to_z,
        'cos_Omega': cos_Omega_to_z,
        'sin_omega': sin_omega_to_z,
        'cos_omega': cos_omega_to_z,
        'sin_f': sin_f_to_z,
        'cos_f': cos_f_to_z,
    }
    
    # Create table of transformed elements
    cols_key = ['Num', 'Name']
    ast_elt_xf = ast_elt[cols_key].copy()
    
    # Original orbital elements
    ast_elt_xf['a'] = ast_elt.a
    ast_elt_xf['e'] = ast_elt.e
    ast_elt_xf['inc'] = ast_elt.inc
    ast_elt_xf['Omega'] = ast_elt.Omega
    ast_elt_xf['omega'] = ast_elt.omega
    ast_elt_xf['f'] = ast_elt.f
    ast_elt_xf['epoch_mjd'] = ast_elt.epoch_mjd

    # Mathematical transforms of the original elements
    ast_elt_xf['log_a'] = log_a
    ast_elt_xf['sin_inc'] = sin_inc
    ast_elt_xf['sin_Omega'] = sin_Omega
    ast_elt_xf['cos_Omega'] = cos_Omega
    ast_elt_xf['sin_omega'] = sin_omega
    ast_elt_xf['cos_omega'] = cos_omega
    ast_elt_xf['sin_f'] = sin_f
    ast_elt_xf['cos_f'] = cos_f
    
    # Transformations to z
    ast_elt_xf['log_a_z'] = log_a_z
    ast_elt_xf['e_z'] = e_z
    ast_elt_xf['sin_inc_z'] = sin_inc_z
    ast_elt_xf['sin_Omega_z'] = sin_Omega_z
    ast_elt_xf['cos_Omega_z'] = cos_Omega_z
    ast_elt_xf['sin_omega_z'] = sin_omega_z
    ast_elt_xf['cos_omega_z'] = cos_omega_z
    ast_elt_xf['sin_f_z'] = sin_f_z
    ast_elt_xf['cos_f_z'] = cos_f_z  

    return ast_elt_xf, interp_tbl

# ********************************************************************************************************************* 
def plot_elt_transform_pdf(ast_elt_xf: pd.DataFrame, elt_name: str):
    """Plot transformed orbital elements"""
    # Set number of sample points
    z_range: float = 6.0
    N_samp_z: int = int(200*z_range + 1)

    # Sample CDF levels: z_samp points evenly spaced by Z
    z_samp = np.linspace(-z_range, z_range, N_samp_z)
    cdf_samp_z = norm.cdf(z_samp)
    pdf_samp_z = norm.pdf(z_samp)
    
    # Transformed element
    elt_name_z = f'{elt_name}_z'
    x = ast_elt_xf[elt_name_z].values

    # Plot transformed element
    fig, ax = plt.subplots()
    ax.set_title(f'Transformed Orbital Element Histogram: {elt_name}')
    ax.set_xlabel(f'{elt_name}_z: Transform of {elt_name}')
    ax.set_ylabel('Density')
    n, bins, patches = ax.hist(x=x, bins=1000, density=True, cumulative=False, label='data', color='blue')
    ax.plot(z_samp, pdf_samp_z, label='N(0,1)', color='red')
    ax.set_xlim([-3.0, 3.0])
    ax.legend()
    ax.grid()

# ********************************************************************************************************************* 
def plot_elt_transform_map(ast_elt_xf: pd.DataFrame, elt_name: str):
    """Plot transformed orbital elements"""
    # Set number of sample points
    z_range: float = 6.0
    N_samp_z: int = int(200*z_range + 1)

    # Sample CDF levels: z_samp points evenly spaced by Z
    z_samp = np.linspace(-z_range, z_range, N_samp_z)
    cdf_samp_z = norm.cdf(z_samp)
    pdf_samp_z = norm.pdf(z_samp)
    
    # Name of transformed element
    elt_name_xf = f'{elt_name}_z'
    
    # Values for plot
    idx = np.argsort(ast_elt_xf[elt_name_xf].values)
    x_plot = ast_elt_xf[elt_name_xf].values[idx]
    y_plot = ast_elt_xf[elt_name].values[idx]

    # Plot transformed element
    fig, ax = plt.subplots()
    ax.set_title(f'Transformed Orbital Element Histogram: {elt_name}')
    ax.set_xlabel(f'{elt_name_xf}: Transform of {elt_name}')    
    ax.set_ylabel(f'{elt_name}')
    ax.set_xlim([-3.0, 3.0])
    ax.plot(x_plot, y_plot, label=elt_name, color='blue')
    ax.legend(loc='lower right')
    ax.grid()    