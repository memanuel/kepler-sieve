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
from scipy.optimize import minimize

# Utility
from tqdm.auto import tqdm

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
# Set plot style variables
mpl.rcParams['figure.figsize'] = [16.0, 10.0]
mpl.rcParams['font.size'] = 16

# Default data type
dtype = np.float64

# ********************************************************************************************************************* 
# Report summary attributes of ZTF elements generated from different candidates (EDA of elements and ZTF interaction)
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def score_by_elt(ztf_elt, thresh_deg=None, fit_mixture: bool = False):
    """Score the ztf observations by element"""   
    # Add log(v) column
    ztf_elt['log_v'] = np.log(ztf_elt.v)

    # Group by element_id
    elt_score = ztf_elt['log_v'].groupby(ztf_elt.element_id).agg(['count', 'mean', 'std'])
    hit_count = ztf_elt['is_hit'].groupby(ztf_elt.element_id.astype(np.int32)).agg(['sum']).rename(columns={'sum':'hits'})

    # Rename columns
    col_name_tbl = {
        'count': 'n_obs',
        'mean': 'log_v_mean',
        'std': 'log_v_std',
    }
    elt_score.rename(columns=col_name_tbl, inplace=True)

    # Add t-score column
    elt_score['t'] = -elt_score.log_v_mean * np.sqrt(elt_score.n_obs)
    # Add column with hit count
    elt_score['hits'] = hit_count

    # Add columns with output of mixture model fit: log_like, h, lambda
    mixture_cols = ['log_like', 'h']
    if fit_mixture:
        # Distinct element_ids
        element_ids = np.unique(ztf_elt.element_id)

        # Array for lambdas
        lams = np.zeros(element_ids.size)
        # Create placeholder columns
        for col in mixture_cols:
            elt_score[col] = 0.0

        # Set low iteration count
        maxiter = 10
        # Threshold distance
        thresh_dist = deg2dist(thresh_deg)

        # Iterate through unique elements
        iterates = tqdm(list(enumerate(element_ids)))
        for i, element_id in iterates:
            # Extract v for this element
            mask = ztf_elt.element_id == element_id
            v = ztf_elt[mask].v
            # Fit mixture model
            log_like, h, lam = mixture_fit(v, maxiter=maxiter)
            # Save mixture fit to elt_score frame
            elt_score.loc[element_id, mixture_cols] = np.array([log_like, h])
            lams[i] = lam

        # Fit width
        mu_v = thresh_dist / lams
        # Convert from z to s
        mu_s = np.sqrt(2.0*mu_v)
        # Convert from s to arc seconds
        mu_deg = dist2deg(mu_s)
        mu_sec = 3600.0 * mu_deg
        # Save to elt_score frame
        elt_score['lambda'] = lams
        elt_score['mu_sec'] = mu_sec

    return elt_score

# ********************************************************************************************************************* 
def score_batches(elt_score_tbl: dict):
    """
    Summarize the score over a batch
    INPUTS:
        elt_score_tbl: Dict. key = batch_name, value= elt_score DataFrame
    """
    
    # Column collections
    cols_log_like = ['log_like_sum', 'log_like_mean', 'log_like_std']
    cols_t = ['t_mean', 't_std']
    # Initialize batch_score frame
    cols_all = cols_log_like + cols_t
    batch_score = pd.DataFrame(columns=cols_all) 
    
    # Iterate through scores in the table
    for batch_name, elt_score in elt_score_tbl.items():
        # Was the mixture model fit?
        has_mixture: bool = 'log_like' in elt_score.columns
        # The mean and standard of the t-score
        t_mean, t_std = elt_score['t'].agg(['mean', 'std'])
        batch_score.loc[batch_name, cols_t] = [t_mean, t_std]
        # The mixture columns if available
        if has_mixture:
            log_like_sum, log_like_mean, log_like_std = elt_score['log_like'].agg(['sum', 'mean', 'std'])
            batch_score.loc[batch_name, cols_log_like] = [log_like_sum, log_like_mean, log_like_std]
        
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
def hl2x(h, lam):
    """Convert h and lam to array logit(h), log(lam) for unconstrained optimization"""
    x = np.array([logit(h), np.log(lam)])
    return x

# ********************************************************************************************************************* 
def x2hl(x):
    """Convert x array to h and lam after unconstrained optimization"""
    tiny = 2.0**-50
    h = sigmoid(x[0]) + tiny
    lam = np.exp(x[1])
    return h, lam
    
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
    # h = sigmoid(x[0])
    # lam = np.exp(x[1])
    h, lam = x2hl(x)

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

# ********************************************************************************************************************* 
def mixture_fit(v, maxiter: int = 20):
    """
    Fit a mixture distribution to observed v
    """
    # Initial guess
    h0 = 0.5
    lam0 = 1.0
    # Convert to array
    x0 = hl2x(h0, lam0)

    # Alias inputs to minimize
    fun = mixture_objective
    jac = mixture_jacobian
    hess = mixture_hessian

    # Additional options for minimize
    options = {
        'maxiter': maxiter,
        'disp': False,
    }

    # Minimize objective with scipy.minimize, Newton Conjugate Gradient method
    res = minimize(fun=fun, x0=x0, args=v, method='Newton-CG', jac=jac, hess=hess, options=options)

    # Extract the log likelihood
    log_like = -res.fun
    # Extract result and convert to h, lam
    h, lam = x2hl(res.x)
    # Return log likelihood and parameters
    return (log_like, h, lam,)
