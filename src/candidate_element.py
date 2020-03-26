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
from scipy.special import expit as sigmoid, logit

# Astronomy
import rebound

# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local
from asteroid_integrate import load_ast_elt
from planets import make_sim_planets
from astro_utils import mjd_to_datetime

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
# Set plot style variables
mpl.rcParams['figure.figsize'] = [16.0, 10.0]
mpl.rcParams['font.size'] = 16

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

# ********************************************************************************************************************* 
def score_by_elt(ztf_elt):
    """Score the ztf observations by element"""
    # Add log(v) column
    ztf_elt['log_v'] = np.log(ztf_elt.v)

    # Group by element_id
    elt_score = ztf_elt['log_v'].groupby(ztf_elt.element_id).agg(['count', 'mean', 'std'])

    # Rename columns
    col_name_tbl = {
        'count': 'n_obs',
        'mean': 'log_v_mean',
        'std': 'log_v_std',
    }
    elt_score.rename(columns=col_name_tbl, inplace=True)
    # Add t-score column
    elt_score['t'] = elt_score.log_v_mean * np.sqrt(elt_score.n_obs)

    # Return score grouped by element
    return elt_score

# ********************************************************************************************************************* 
def score_batches(elt_score_tbl: dict):
    """
    Summarize the score over a batch
    INPUTS:
        elt_score_tbl: Dict. key = batch_name, value= elt_score DataFrame
    """
    # Initialize batch_score frame
    batch_score = pd.DataFrame(columns=['t_mean', 't_std']) 
    # Iterate through scores in the table
    for batch_name, elt_score in elt_score_tbl.items():
        t_mean, t_std = elt_score['t'].agg(['mean', 'std'])
        batch_score.loc[batch_name] = [t_mean, t_std]
        
    return batch_score

# ********************************************************************************************************************* 
def report_v(ztf_elt: pd.DataFrame, elt_name: str):
    """Summary statistics for v and log(v)"""

    v: np.ndarray = ztf_elt.v
    
    # Summary statistics for v
    # E[V] = 0.5, Var[V] = 1/12
    v_mean = np.mean(v)
    v_std = np.std(v)
    v_min = np.min(v)
    v_max = np.max(v)
    
    # Number of observations
    n = v.size

    # Summary statistics for log(v)
    # E[log(V)] = -1.0, Var[log(V)] = 1.0
    log_v_mean = np.mean(np.log(v))
    log_v_std = np.std(np.log(v))
    # t-score; theoretical standard deviation for 
    # log_v_t = (-log_v_mean - 1.0) / (log_v_std * np.sqrt(n))
    log_v_t = (-log_v_mean - 1.0) * np.sqrt(n)
    
    # Report
    print(f'v = (1-z) / (1-z_thresh) for {elt_name} orbital elements:')
    print(f'Mean: {v_mean:10.6f}')
    print(f'Std : {v_std:10.6f}')
    # print(f'Min : {v_min:10.6f}')
    # print(f'Max : {v_max:10.6f}')
    print(f'\nlog(v):')
    print(f'Mean: {log_v_mean:10.6f}')
    print(f'Std : {log_v_std:10.6f}')
    print(f't   : {log_v_t:10.6f}')

# ********************************************************************************************************************* 
def plot_v(ztf_elt: pd.DataFrame, elt_name: str, v_max = 1.0, is_cum: bool = False):
    """Plot histogram of v"""

    # Extract v
    v = ztf_elt.v
    
    # Set range
    v_range = [0.0, v_max]
    
    # Set y_label
    y_label = 'Cumulative Probability' if is_cum else 'Density'

    # Plot density of v
    fig, ax = plt.subplots()
    ax.set_title(f'Score Histogram for {elt_name}: v')
    ax.set_xlabel('v')
    ax.set_ylabel(y_label)
    n, bins, patches = ax.hist(x=v, bins=101, range=v_range, density=True, cumulative=is_cum, color='blue')
    if is_cum:
        ax.plot(v, v, color='red')
    else:
        ax.axhline(1.0, color='red')
    ax.set_xlim(v_range)
    # ax.legend()
    ax.grid()
    # fig.savefig('../figs/elts/elt_hist_a.png', bbox_inches='tight')
    plt.show()

# ********************************************************************************************************************* 
def mixture_log_like(h, lam, x):
    """
    Log-Likelihood function for mixture model.
    h:   Fraction of hits
    lam: Exponential decay rate for hit probability as v increases
    x:   Array of observed v for observations
    Probability density for hits is exponential: 
        x|hit ~ Expo(lambda)
        g(x) = lambda exp(-lambda x) / N
        N normalization factor 1-exp(-lam) for truncated normal on [0, 1]
    Probability density for misses is uniform:
        x|miss ~ Unif(0, 1)
    The mixture has probability h of being a hit, and (1-h) of being a miss, so
        f(x) = h g(x) + (1-h)
    """
    # Exponential of -lambda x
    emlx = np.exp(-lam * x)
    # Exponential distribution before normalization
    u = lam * emlx
    # Normalization for the denominator of the exponential
    N = (1.0 - np.exp(-lam))
    
    # Conditional probability on the exponential distribution
    g = u / N
    # Mixture probability
    f = h * g + (1.0 - h)
    
    # Log likelihood
    log_f = np.log(f)
       
    # Total of log likelihood
    L_tot = np.sum(log_f)
    
    # First derivative g'(lam)
    du = emlx * (1.0 - lam * x)
    dN = 1.0 - N
    N2 = N**2
    dg = (du * N - u * dN) / N2

    # First derivatives of f
    df_dh = g - 1.0
    df_dl = h * dg
    # First derivatives of log likelihood L = log(f)
    dL_dh = df_dh / f
    dL_dl = df_dl / f
    
    # Second derivative of g
    d2u = emlx * (lam * x - 2.0) * x
    d2g = (d2u - 2*dN * dg) / N + (u*dN) / N2
    
    # Second derivatives of f
    d2f_dh2 = np.zeros_like(x)
    d2f_dl2 = h * d2g
    d2f_dhdl = dg

    # Second derivatives of L   
    d2L_dh2 = d2f_dh2 / f - (df_dh/f)**2
    d2L_dl2 = d2f_dl2 / f - (df_dl/f)**2
    d2L_dhdl = d2f_dhdl / f - (df_dh * df_dl)/f**2
    
    # Gradient of log likelihood w.r.t. lambda
    dL_dh_tot = np.sum(dL_dh)
    dL_dl_tot = np.sum(dL_dl)
    
    # Assemble gradient
    grad = np.array([dL_dh_tot, dL_dl_tot])
    
    # Assemble Hessian
    d2L_dh2_tot = np.sum(d2L_dh2)
    d2L_dl2_tot = np.sum(d2L_dl2)
    d2L_dhdl_tot = np.sum(d2L_dhdl)
    hess = np.array([[d2L_dh2_tot, d2L_dhdl_tot], [d2L_dhdl_tot, d2L_dl2_tot]])
    
    return L_tot, grad, hess

# ********************************************************************************************************************* 
def mixture_objective_all(x, v):
    """
    Optimization function for unconstrained minimization
    Apply coordinate transformations
    h = sigmoid(x[0])
    lam = exp(x[1])
    for an unconstrained minimization problem.
    Return the negative of the log likelihood so it is a minimization
    """
    # Unpack h, lam from x and apply transformations
    h = sigmoid(x[0])
    lam = np.exp(x[1])

    # Delegate to log_like
    L, grad_L, hess_L = mixture_log_like(h=h, lam=lam, x=v)
    # Unpack gradient
    dL_dh, dL_dl = grad_L
    # Unpack hessian
    d2L_dh2 = hess_L[0, 0]
    d2L_dl2 = hess_L[1, 1]
    d2L_dhdl = hess_L[0 ,1]
       
    # Apply chain rule and flip signs
    F = -L
    dh_dx = h * (1.0 - h)
    dl_dy = lam
    d2h_dx2 = dh_dx * (1.0 - 2*h)
    dF_dx = -dL_dh * dh_dx
    dF_dy = -dL_dl * dl_dy
    grad_F = np.array([dF_dx, dF_dy])
    
    # Second derivative of F    
    d2F_dx2 = -(d2L_dh2 * dh_dx**2 + dL_dh * d2h_dx2)
    d2F_dy2 = -(d2L_dl2 * lam**2 + dL_dl * lam)
    d2F_dxdy = -(d2L_dhdl * lam * dh_dx)
    hess_F = np.array([[d2F_dx2, d2F_dxdy], [d2F_dxdy, d2F_dy2]])

    # Return F and its gradient
    return F, grad_F, hess_F

# ********************************************************************************************************************* 
def mixture_objective(x, v):
    """Objective function; value only, no derivatives"""
    F, grad_F, hess_F = mixture_objective_all(x=x, v=v)
    return F

# ********************************************************************************************************************* 
def mixture_jacobian(x, v):
    """Derivative of objective function"""
    F, grad_F, hess_F = mixture_objective_all(x=x, v=v)
    return grad_F

# ********************************************************************************************************************* 
def mixture_hessian(x, v):
    """Hessian of objective function"""
    F, grad_F, hess_F = mixture_objective_all(x=x, v=v)
    return hess_F