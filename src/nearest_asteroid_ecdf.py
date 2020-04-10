"""
Harvard IACS Masters Thesis
nearest_asteroid_ecdf.py: 
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
# q_ast = calc_elt_pos(elts=elts_ast, ts=ts)
# ast_elt_xf, interp_tbl = ast_elt_transform(ast_elt)
# beta, X_beta = calc_beta(ast_elts_xf)

# ********************************************************************************************************************* 
# Default data type
dtype = np.float64

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
    sin_f = np.sin(ast_elt.f.values)
    cos_f = np.cos(ast_elt.f.values)
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
def calc_beta(ast_elt_xf):
    """
    Calculate the matrix beta such that (X * beta) has covariance matrix I_n, i.e.
    (X*beta)^T (X*beta) = I_n
    """
    # Relevant columns
    cols_xf = ['log_a_z', 'e_z', 'sin_inc_z', 'sin_Omega_z', 'cos_Omega_z', 'sin_omega_z', 'cos_omega_z', 'sin_f_z', 'cos_f_z',]
    # Nx9 matrix of transformed elements
    X = ast_elt_xf[cols_xf].values
    # Scale columns 3:9 by sqrt(1/2) so they are not weighted 2x (need to represent them as sin, cos pair but don't want to overweight them)
    X[:, 3:9] *= np.sqrt(0.5)

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
beta, X_beta = calc_beta(ast_elt_xf)

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
    # Nx9 matrix of transformed elements
    X = elts_xf[cols_xf].values
    # Scale columns 3:9 by sqrt(1/2) so they are not weighted 2x 
    # (need to represent them as sin, cos pair but don't want to overweight them)
    X[:, 3:9] *= np.sqrt(0.5)
    return X

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
    row_num = ast_elt.row_num.loc[ast_num].values
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
    dist = np.linalg.norm(X_beta.reshape(-1, 1, 9) - Y_beta.reshape(1, -1, 9), axis=-1)

    # Row number of nearest asteroid elements
    row_idx = np.argmin(dist, axis=0)
    # Asteroid numbers of nearest asteroid elements
    ast_num = ast_elt.Num.values[row_idx]

    # Q_norm to nearest asteroid element
    col_idx = np.arange(row_idx.size, dtype=np.int32)
    Q_norm = dist[row_idx, col_idx]

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
        elts['nearest_ast_dist'] = Q_norm

    return ast_elt_nearest
