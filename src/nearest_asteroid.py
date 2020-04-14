"""
Harvard IACS Masters Thesis
nearest_asteroid.py: 
Search known asteroids for the one nearest to input orbital elements

Michael S. Emanuel
Wed Mar 25 09:55 2020
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF

# Machine learning
import tensorflow as tf

# Utility
import os
from datetime import datetime

# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# MSE imports
from asteroid_element import load_ast_elt
from asteroid_model import make_model_ast_pos
from astro_utils import datetime_to_mjd

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
# Set time array for known asteroids
t0 = datetime_to_mjd(datetime(2010,1,1))
t1 = datetime_to_mjd(datetime(2030,1,1))
dt = 31
ts = np.arange(t0, t1, dt)

# ********************************************************************************************************************* 
# Module level constants defined below after manufacturing functions defined:
# q_known_ast = calc_elt_pos(elts=elts_ast, ts=ts)
# X_known_ast = keras.backend.constant(value=q_known_ast, dtype=dtype)
# ast_elt_xf, interp_tbl = ast_elt_transform(ast_elt)
# beta, X_beta = calc_beta(ast_elts_xf)

# ********************************************************************************************************************* 
# Default data type
dtype = np.float64

# ********************************************************************************************************************* 
# Search for known asteroid with closest elements by comparing orbital trajectories in Cartesian space
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def calc_elt_pos(elts, ts):
    """
    Calculate position of orbital elements on a fixed array of times
    """    
    # Number of asteroids in elts DataFrame
    N_ast = elts.shape[0]
    
    # Number of time points
    N_t = ts.size
    
    # Number of spatial dimensions
    space_dims = 3

    # Tile the time array N_ast times
    ts_np = np.tile(ts, N_ast).astype(np.float32)    

    # Build row_lengths; all the same
    row_lengths_np = np.repeat(N_t, N_ast).astype(np.int32)

    # Build asteroid position model
    model_pos = make_model_ast_pos(ts_np=ts_np, row_lengths_np=row_lengths_np)

    # Extract orbital elements from elts DataFrame
    a = elts.a.values
    e = elts.e.values
    inc = elts.inc.values
    Omega = elts.Omega.values
    omega = elts.omega.values
    f = elts.f.values
    epoch = elts.epoch.values
    
    # Position according to position model
    q_ast, v_ast = model_pos([a, e, inc, Omega, omega, f, epoch])
    
    # Reshape q_ast to fixed size
    q_ast = tf.reshape(q_ast.values, (N_ast, N_t, space_dims)).numpy()
    
    return q_ast

# ********************************************************************************************************************* 
def load_known_ast_pos():
    """
    Calculate position of asteroids on a fixed array of times
    """
    # Filename
    file_path = f'../data/candidate_elt/known_ast_pos.npz'
    
    # DataFrame of known asteroid positions to be returned
    ast_pos: pd.DataFrame

    try:
        q_ast = np.load(file_path)['q_ast']
        # print(f'Loaded {file_path} from disk.')
    except FileNotFoundError:
        print(f'Unable to load {file_path}, generating it...')
        # Calculate the asteroid position and save to disk in chunks of 100,000 to avoid exhausting GPU memory
        N_ast = ast_elt.shape[0]
        N_t = ts.size
        space_dims = 3
        # Set up chunks
        chunk_size = 100000
        num_chunks = N_ast // chunk_size
        # Preallocate q_ast
        q_ast = np.zeros((N_ast, N_t, space_dims))        
        # Calculate and save chunks one at a time
        for i in range(num_chunks):
            r0 = i * chunk_size
            r1 = r0 + chunk_size
            q_ast_chunk = calc_elt_pos(elts=ast_elt.iloc[r0:r1], ts=ts)
            q_ast[r0:r1, :, :] = q_ast_chunk
        # Save the big file
        np.savez(file=f'{file_path}', q_ast=q_ast)
    
    return q_ast

# ********************************************************************************************************************* 
# Load the known asteroid positions
q_known_ast = load_known_ast_pos()
# Convert to a float32 tensor
X_known_ast = tf.constant(q_known_ast, dtype=tf.float32)

# ********************************************************************************************************************* 
def nearest_ast_elt_cart(elts, X_known_ast=X_known_ast):
    """
    Search for nearest asteroid element based on Cartesian orbits
    """
    # Extract element_id
    element_id = elts.element_id.values

    # Calculate position of element
    q_elt = calc_elt_pos(elts, ts)
    # Convert to a tensor; use float32 to save memory
    Y = tf.constant(q_elt, dtype=tf.float32)

    # Number of asteroids and test elements
    N_ast = X_known_ast.shape[0]
    N_elt = Y.shape[0]
    # Sqrt of N_ast used to get RMS distance
    sqrt_N_ast = np.sqrt(N_ast).astype(np.float32)

    # Need to do this one candidate at a time b/c run out of memory when doing all at once
    idx = np.zeros(N_elt, dtype=np.int32)
    ast_num = np.zeros(N_elt, dtype=np.int32)
    dist = np.zeros(N_elt, dtype=np.float32)

    # Iterate over candidates    
    for i in range(N_elt):
        # dist_i is the distance between element i and all N_ast asteroids; shape [N_ast,]
        dist_i = tf.linalg.norm(X_known_ast - Y[i], axis=(-2, -1))
        # The index of the closest asteroid
        j = tf.argmin(dist_i)
        # Save the asteroid number and distance of the closest asteroid
        idx[i] = j
        dist[i] = dist_i[j] / sqrt_N_ast
        # Save the asteroid number of the closest asteroid
        ast_num[i] = ast_elt.Num.values[j]

    # Columns from ast_elt we want in ast_elt_nearest
    cols = ['Num', 'Name', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
    # Select the closest asteroids and these columns
    ast_elt_nearest = ast_elt.loc[ast_num, cols]
    # Use naming conventions for candidate elements
    col_rename = {
        'Num': 'nearest_ast_num',
        'Name': 'nearest_ast_name',
    }
    ast_elt_nearest.rename(columns=col_rename, inplace=True)        

    # Add column with the element_id this asteroid is closest to
    ast_elt_nearest.insert(loc=0, column='element_id', value=element_id)

    # Add columns to ast_elt_nearest showing the element_id and distance
    ast_elt_nearest.insert(loc=3, column='nearest_ast_dist', value=dist)

    # Index the nearest asteroid frame by row number so it aligns with 
    ast_elt_nearest.reset_index(drop=True, inplace=True)

    # Add two extra columns to elts showing the asteroid number and distance
    if 'nearest_ast_num' not in elts.columns:
        elts.insert(loc=elts.columns.size, column='nearest_ast_num', value=ast_num)
    else:
        elts['nearest_ast_num'] = ast_num
    if 'nearest_ast_dist' not in elts.columns:
        elts.insert(loc=elts.columns.size, column='nearest_ast_dist', value=dist)
    else:
        elts['nearest_ast_dist'] = dist

    return ast_elt_nearest

# ********************************************************************************************************************* 
# Transform orbital elements to approximate Z = Norm(0, 1) using mean, standard deviation
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_interp_x(x: np.ndarray):
    """
    Build a linear function that takes as inputs x and returns a variable z that is approximately standard normal.
    INPUTS:
        x: Values whose CDF is to be taken, e.g. log(a)
    OUTPUT:
        x_to_z: a lambda function that standardizes x based on its first 2 moments
    """
    # Mean of x
    mu = np.mean(x)

    # Standard deviation of x
    sigma = np.std(x)    
    
    # Interpolator from x to z
    x_to_z = lambda x: (x - mu) / sigma
    
    return x_to_z

# ********************************************************************************************************************* 
def trig_to_z(x):
    """Convert a trigonometric variable to a standard uniform"""
    tiny = 2.0**-20
    u = (0.5+np.arcsin(x)/np.pi)*(1-tiny)
    return norm.ppf(u)

# ********************************************************************************************************************* 
def ast_elt_transform(ast_elt: pd.DataFrame):
    """
    Build interpolators from orbital elements to uniform z and populate DataFrame of transformed elements
    """

    # Build interpolator from a to z
    log_a = np.log(ast_elt.a.values)
    log_a_to_z = make_interp_x(x=log_a)
    # Compute transformed values
    log_a_z = log_a_to_z(log_a)
    
    # Build interpolator from e to z
    e = ast_elt.e.values
    e_to_z = make_interp_x(x=e)
    # Compute transformed values
    e_z = e_to_z(e)
    
    # Build interpolator from sin(inc) to z
    sin_inc = np.sin(ast_elt.inc.values)
    sin_inc_to_z = make_interp_x(x=sin_inc)
    # Compute transformed values
    sin_inc_z = sin_inc_to_z(sin_inc)
    
    # Build interpolator from sin(Omega) and cos(Omega) to z
    sin_Omega = np.sin(ast_elt.Omega.values)
    cos_Omega = np.cos(ast_elt.Omega.values)
    # sin_Omega_to_z = make_interp_x(x=sin_Omega)
    # cos_Omega_to_z = make_interp_x(x=cos_Omega)
    # Compute transformed values
    # sin_Omega_z = sin_Omega_to_z(sin_Omega)
    # cos_Omega_z = cos_Omega_to_z(cos_Omega)
    sin_Omega_z = trig_to_z(sin_Omega)
    cos_Omega_z = trig_to_z(cos_Omega)
    
    # Build interpolator from sin(omega) and cos(omega) to z
    sin_omega = np.sin(ast_elt.omega.values)
    cos_omega = np.cos(ast_elt.omega.values)
    sin_omega_to_z = make_interp_x(x=sin_omega)
    cos_omega_to_z = make_interp_x(x=cos_omega)
    # Compute transformed values
    # sin_omega_z = sin_omega_to_z(sin_omega)
    # cos_omega_z = cos_omega_to_z(cos_omega)
    sin_omega_z = trig_to_z(sin_omega)
    cos_omega_z = trig_to_z(cos_omega)
    
    # Build interpolator from sin(f) and cos(f) to z
    sin_f = np.sin(ast_elt.f.values)
    cos_f = np.cos(ast_elt.f.values)
    # sin_f_to_z = make_interp_x(x=sin_f)
    # cos_f_to_z = make_interp_x(x=cos_f)
    # Compute transformed values
    # sin_f_z = sin_f_to_z(sin_f)
    # cos_f_z = cos_f_to_z(cos_f)
    sin_f_z = trig_to_z(sin_f)
    cos_f_z = trig_to_z(cos_f)

    
    # Dictionary of interpolators
    interp_tbl = {
        'log_a': log_a_to_z,
        'e': e_to_z,
        'sin_inc': sin_inc_to_z,
        # 'sin_Omega': sin_Omega_to_z,
        # 'cos_Omega': cos_Omega_to_z,
        # 'sin_omega': sin_omega_to_z,
        # 'cos_omega': cos_omega_to_z,
        # 'sin_f': sin_f_to_z,
        # 'cos_f': cos_f_to_z,
        'sin_Omega': trig_to_z,
        'cos_Omega': trig_to_z,
        'sin_omega': trig_to_z,
        'cos_omega': trig_to_z,
        'sin_f': trig_to_z,
        'cos_f': trig_to_z,
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
    ast_elt_xf['epoch'] = ast_elt.epoch

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
# Build transformed elements and interpolators
ast_elt_xf, interp_tbl = ast_elt_transform(ast_elt)

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

# ********************************************************************************************************************* 
def elts_to_X_cov(elts):
    """Transform elements to X where each dimension is standard normal"""
    # Copy index only of elts
    elts_xf = elts.iloc[:, 0:0].copy()

    # Add transformed values of input elements
    elts_xf['log_a'] = np.log(elts.a)
    elts_xf['e'] = elts.e
    elts_xf['sin_inc'] = np.sin(elts.inc)
    elts_xf['sin_Omega'] = np.sin(elts.Omega)
    elts_xf['cos_Omega'] = np.cos(elts.Omega)
    elts_xf['sin_omega'] = np.sin(elts.omega)
    elts_xf['cos_omega'] = np.cos(elts.omega)
    elts_xf['sin_f'] = np.sin(elts.f)
    elts_xf['cos_f'] = np.cos(elts.f)

    # Transform to normally distributed Z
    elts_xf['log_a_z'] = interp_tbl['log_a'](elts_xf.log_a)
    elts_xf['e_z'] = interp_tbl['e'](elts_xf.e)
    elts_xf['sin_inc_z'] = interp_tbl['sin_inc'](elts_xf.sin_inc)
    elts_xf['sin_Omega_z'] = interp_tbl['sin_Omega'](elts_xf.sin_Omega)
    elts_xf['cos_Omega_z'] = interp_tbl['cos_Omega'](elts_xf.cos_Omega)
    elts_xf['sin_omega_z'] = interp_tbl['sin_omega'](elts_xf.sin_omega)
    elts_xf['cos_omega_z'] = interp_tbl['cos_omega'](elts_xf.cos_omega)
    elts_xf['sin_f_z'] = interp_tbl['sin_f'](elts_xf.sin_f)
    elts_xf['cos_f_z'] = interp_tbl['cos_f'](elts_xf.cos_f)

    # Relevant columns
    cols_xf = ['log_a_z', 'e_z', 'sin_inc_z', 'sin_Omega_z', 'cos_Omega_z', 'sin_omega_z', 'cos_omega_z', 'sin_f_z', 'cos_f_z',]
    # Importance weights for the columns
    importance = np.array([1.0, 1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,])
    # Normalize importance so it sums to 1.0
    importance = importance / np.sum(importance)
    # We multiply each column by the square root of its importance
    sqrt_importance = np.sqrt(importance)

    # Nx9 matrix of transformed elements
    X = elts_xf[cols_xf].values
    # Scale columns by importance weights
    # (need to represent them as sin, cos pair but don't want to overweight them)
    X = np.dot(X, np.diag(sqrt_importance))
    return X

# ********************************************************************************************************************* 
def calc_beta(elts):
    """
    Calculate the matrix beta such that (X * beta) has covariance matrix I_n, i.e.
    (X*beta)^T (X*beta) = I_n
    """
    # Build X
    X = elts_to_X_cov(elts)

    # Covariance matrix of these X's
    Q = np.cov(m=X, rowvar=False)
    # Eigenvalue decomposition
    lam, P = np.linalg.eig(Q)
    # Filter out tiny imaginary components due to roundoff
    lam = np.real(lam)
    # Square root of diagonal matrix
    d = np.diag(np.sqrt(lam))
    d_inv = np.diag(1.0 / np.sqrt(lam))
    # The beta matrix
    beta = np.dot(P, d_inv)
    # The product X_beta
    X_beta = np.dot(X, beta)
    
    return beta, X_beta

# ********************************************************************************************************************* 
# Build tranformation matrix beta and X_beta for computing distance to orbital elements
beta, X_beta = calc_beta(ast_elt)

# ********************************************************************************************************************* 
def elt_q_norm(elts: pd.DataFrame, ast_num: np.ndarray):
    """
    Compute Q-norm style distance between orbital elements and specific asteroids
    INPUTS:
        elts: DataFrame of orbital elements
        ast_num: Numpy array of asteroid numbers
    OUTPUTS:
        q_norm: The covariance style distance between elts and these asteroids
    """
    # Transform of inputs Y; shape [batch_size, 9,]
    Y = elts_to_X_cov(elts)
    # Multiply by beta
    Y_beta = np.dot(Y, beta)
    
    # Row numbers of the input asteroid numbers; numpy array of shape [batch_size,]
    row_num = ast_elt.loc[ast_num].row_num.values
    # X_beta for these asteroid numbers
    X_beta_ast = X_beta[row_num]
    
    # Difference in the transformed space
    dU = Y_beta - X_beta_ast
    
    return np.linalg.norm(dU, axis=1)
    
# ********************************************************************************************************************* 
def nearest_ast_elt_cov(elts):
    """
    Search the known asteroid elements for the nearest one to the input elements.
    Calculation based on covariance of trasnformed orbital elements.
    INPUTS:
        elts: DataFrame of elements to search
    """

    # Extract element_id
    element_id = elts.element_id.values

    # Copy index only of elts
    elts_xf = elts.iloc[:, 0:0].copy()

    # Add transformed values of input elements
    Y = elts_to_X_cov(elts)
    # Multiply by beta
    Y_beta = np.dot(Y, beta)

    # Distance from Y_beta to X_beta; shape [N_ast, N_elt] e.g. [733489, 64]
    # Use numpy broadcasting trick to avoid an expensive for loop
    Q_norm_tbl = np.linalg.norm(X_beta.reshape(-1, 1, 9) - Y_beta.reshape(1, -1, 9), axis=-1)

    # Row number of nearest asteroid elements
    row_idx = np.argmin(Q_norm_tbl, axis=0)
    # Asteroid numbers of nearest asteroid elements
    ast_num = ast_elt.Num.values[row_idx]

    # Q_norm to nearest asteroid element
    col_idx = np.arange(row_idx.size, dtype=np.int32)
    Q_norm = Q_norm_tbl[row_idx, col_idx]

    # The closest asteroid elements
    # Columns from ast_elt we want in ast_elt_nearest
    cols = ['Num', 'Name', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
    # Select the closest asteroids and these columns
    ast_elt_nearest = ast_elt.loc[ast_num, cols]
    # Use naming conventions for candidate elements
    col_rename = {
        'Num': 'nearest_ast_num',
        'Name': 'nearest_ast_name',
    }
    ast_elt_nearest.rename(columns=col_rename, inplace=True)        
    
    # Add column with the element_id this asteroid is closest to
    ast_elt_nearest.insert(loc=0, column='element_id', value=element_id)

    # Add columns to ast_elt_nearest showing the element_id and distance
    ast_elt_nearest.insert(loc=3, column='nearest_ast_Q_norm', value=Q_norm)

    # Index the nearest asteroid frame by row number so it aligns with 
    ast_elt_nearest.reset_index(drop=True, inplace=True)

    # Add two extra columns to elts showing the asteroid number and Q_norm
    if 'nearest_ast_num_cov' not in elts.columns:
        elts.insert(loc=elts.columns.size, column='nearest_ast_num_cov', value=ast_num)
    else:
        elts['nearest_ast_num_cov'] = ast_num
    if 'nearest_ast_Q_norm' not in elts.columns:
        elts.insert(loc=elts.columns.size, column='nearest_ast_Q_norm', value=Q_norm)
    else:
        elts['nearest_ast_Q_norm'] = Q_norm

    return ast_elt_nearest
