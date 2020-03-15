"""
Harvard IACS Masters Thesis
ZTF Data Visualization
Visualize exploratory data analysis for ZTF2 data and nearest asteroids.

Michael S. Emanuel
27-Feb-2020
"""

# Standard libraries
import numpy as np
import pandas as pd

# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython.display import Image

# ********************************************************************************************************************* 
def ztf_obs_by_month(ztf):
    """Generate a chart summarizing observations by month in ZTF data"""
    # Extract the year-month tuple for each observation for summarizing
    tt = Time(ztf.mjd, format='mjd')
    isotimes = tt.iso
    ym = np.array([isotime[0:7] for isotime in isotimes])
    ym_ser = pd.Series(data=ym, index=ztf.index)

    # Group data by month for monthly summary
    obs_monthly = ztf.groupby(ym_ser)
    obs_monthly_count = obs_monthly.size()

    # Calculations for plot
    month_strs = obs_monthly_count.index.values
    x_values = np.arange(obs_monthly_count.size)
    x_dates = [date(int(x[0:4]), int(x[5:7]), 1) for x in month_strs]
    y_values = obs_monthly_count.values

    # Plot the number of observations by month
    fig, ax = plt.subplots()
    ax.set_title('Alerce Asteroid Observations by Month')
    ax.set_xlabel('Month')
    ax.set_ylabel('Asteroid Observations')
    # ax.bar(x=x_values, height=y_values, tick_label=month_strs, color='blue')
    ax.bar(x=x_values, height=y_values, color='blue')
    ax.set_xticks(x_values[::3])
    ax.set_xticks(x_values, minor=True)
    ax.set_xticklabels(month_strs[::3], minor=False)
    # ax.legend()
    ax.grid()
    fig.savefig('../figs/ztf/alerce_ast_per_month.png', bbox_inches='tight')
    plt.show()

# ********************************************************************************************************************* 
def dist2rad(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to radians in [0, pi]"""
    x_rad = np.arcsin(0.5 * dist) * 2.0
    return x_rad

# ********************************************************************************************************************* 
def rad2dist(x_rad):
    """Convert a distance on unit sphere from radians in [0, pi] to cartesian distance in [0, 2]"""
    return np.sin(0.5 * x_rad) * 2.0

# ********************************************************************************************************************* 
def dist2deg(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to degrees in [0, 180]"""
    x_rad = dist2rad(dist)
    return np.rad2deg(x_rad)

# ********************************************************************************************************************* 
def deg2dist(x_deg):
    """Convert a distance on unit sphere from degrees in [0, 180] to cartesian distance in [0, 2]"""
    x_rad = np.deg2rad(x_deg)
    return rad2dist(x_rad)

# ********************************************************************************************************************* 
def cdf_nearest_dist(dist: np.ndarray, n: int, thresh_deg: float = 180.0):
    """
    Compute the theoretical CDF of the nearest distance to n points that are distributed uniformly at random on the sphere.
    INPUTS:
        dist: Distances in cartesian space (so dist ranges from 0 to 2.0)
        n:    Number of asteroids compared to each observation, e.g. 16 or 100
        thresh_deg: Threshold in degrees; only distances smaller than threshold are included
    OUTPUTS:
        p:    CDF; these will be distributed as Unif[0, 1] for random data
    """
    # Convert dist to radians; result will be in the interval [0, 2pi]
    x = dist2rad(dist)
    # The theoretical CDF of the min of n distances
    alpha = 0.5*(1.0 + np.cos(x))
    p = 1.0 - alpha**n
    # Compute the CDF at the threshold
    thresh = np.deg2rad(thresh_deg)
    alpha_thresh = 0.5*(1.0 + np.cos(thresh))
    p_thresh = 1.0 - alpha_thresh**n
    # return the conditional probability distribution
    p_cond = p / p_thresh
    return p_cond

# ********************************************************************************************************************* 
def angle_to_str(angle_deg) -> str:
    """
    Convert an angle in degrees to a user friendly description in degrees, arc minutes or arc seconds
    INPUTS:
        angle_deg: an angle in degrees
    OUTPUTS:
        angle_caption: a string suitable for a chart caption, e.g. '1.0 Arc Seconds'
        angle_str:     a terse string without spaces, e.g. '1_arc_sec'
        angle_unit:    units for quoting this angle; one of 'deg', 'min', 'sec'
    """
    # String description of threshold
    angle_min = 60 * angle_deg
    angle_sec = 60 * angle_min
    if 1.0 <= angle_deg:
        angle_caption = f'{angle_deg:.0f} Degrees'
        angle_str = f'{angle_deg:.0f}_deg'
        angle_unit = 'deg'
    elif 1.0 <= angle_min:
        angle_caption = f'{angle_min:.0f} Arc Minutes'
        angle_str = f'{angle_min:.0f}_arc_min'
        angle_unit = 'min'
    elif 1.0 <= angle_sec:
        angle_caption = f'{angle_sec:.0f} Arc Seconds'
        angle_str = f'{angle_sec:.0f}_arc_sec'
        angle_unit = 'sec'
    return angle_caption, angle_str, angle_unit
    
# ********************************************************************************************************************* 
def plot_cdf_uncond(ztf: pd.DataFrame, n: int, bins=100):
    """
    Generate plot of two CDFs of asteroid distance.
    INPUTS:
        ztf: DataFrame with distance to nearest asteroid
        n:   Number of asteroids, e.g. 16
    """
    # Number of rows in data
    N_obs = ztf.shape[0]
    # Histogram: unconditional probability
    p = cdf_nearest_dist(dist=ztf.nearest_ast_dist.values, n=n, thresh_deg=180.0)

    # Plot unconditional histogram
    fig, ax = plt.subplots()
    ax.set_title(f'Distance to Nearest Asteroid for n={n} (No Threshold)')
    ax.set_xlabel('Percentile of Random Distribution')
    ax.set_ylabel('Observed Frequency')
    freq, bins_np, patches = ax.hist(x=p, bins=bins, color='blue', label='observed')
    bin_count = bins_np.size - 1
    random_freq = N_obs / bin_count
    ax.axhline(y=random_freq, color='red', label='random')
    ax.legend()
    ax.grid()
    fig.savefig(f'../figs/ztf/nearest_ast_hist_n={n}.png', bbox_inches='tight')
    plt.show()
    return fig, ax

# ********************************************************************************************************************* 
def plot_cdf_cond(ztf, n: int, thresh_deg: float = 1.0, bins=20):
    """
    Generate plot of two CDFs of asteroid distance.
    INPUTS:
        ztf: DataFrame with distance to nearest asteroid
        n:   Number of asteroids, e.g. 16
        thresh_deg: Threshold in degrees for close observations
    """
    # Number of rows in data
    N_obs = ztf.shape[0]

    # Threshold distance and flag indicating whether observations are within threshold
    thresh_dist = deg2dist(thresh_deg)
    is_close = ztf.nearest_ast_dist < thresh_dist

    # Probability of being close at this threshold on a random observation 
    prob_close = cdf_nearest_dist(dist=thresh_dist, n=n)

    # CDF 
    p = cdf_nearest_dist(dist=ztf[is_close].nearest_ast_dist.values, n=n, thresh_deg=1.0)

    # String description of threshold
    thresh_caption, thresh_str, angle_unit = angle_to_str(angle_deg=thresh_deg)

    # Plot unconditional histogram
    fig, ax = plt.subplots()
    ax.set_title(f'Distance to Nearest Asteroid for n={n} (Threshold {thresh_caption} )')
    ax.set_xlabel('Percentile of Random Distribution')
    ax.set_ylabel('Observed Frequency')
    freq, bins_np, patches = ax.hist(x=p, bins=bins, color='blue', label='observed')
    bin_count = bins_np.size - 1
    random_freq = N_obs * prob_close/ bin_count
    ax.axhline(y=random_freq, color='red', label='random')
    ax.legend()
    ax.grid()
    fig.savefig(f'../figs/ztf/nearest_ast_hist_cdf_n={n}_thresh={thresh_str}.png', bbox_inches='tight')
    plt.show()
    return fig, ax

# ********************************************************************************************************************* 
def plot_dist(ztf, n: int, thresh_deg: float = 1.0, bins=20):
    """
    Generate plot with histogram of distance to the nearest asteroid vs. random baseline.
    INPUTS:
        ztf: DataFrame with distance to nearest asteroid
        n:   Number of asteroids, e.g. 16
        thresh_deg: Threshold in degrees for close observations
    """
    # Number of rows in data
    N_obs = ztf.shape[0]

    # Threshold distance and flag indicating whether observations are within threshold
    thresh_dist = deg2dist(thresh_deg)
    is_close = ztf.nearest_ast_dist < thresh_dist

    # Cartesian distance of close observations
    dist = ztf[is_close].nearest_ast_dist.values
    # Distance in degrees and arc seconds
    dist_deg = dist2deg(dist)
    # dist_sec = 3600.0 * dist_deg

    # String description of threshold
    thresh_caption, thresh_str, angle_unit = angle_to_str(angle_deg=thresh_deg)

    # data for x axis; scale according to selected units
    # key = angle_unit, value = (plot data, threshold)
    plot_x_tbl = {
        'deg': (dist_deg, thresh_deg,),
        'min': (dist_deg * 60.0, thresh_deg * 60.0,),
        'sec': (dist_deg * 3600.0, thresh_deg * 3600.0),
    }
    plot_x, thresh_x = plot_x_tbl[angle_unit]

    # x-axis label depends on unit choice
    # key = angle_unit, value = axis label
    axis_label_tbl = {
        'deg': 'Distance in Degrees',
        'min': 'Distance in Arc Minutes',
        'sec': 'Distance in Arc Seconds',
    }

    # Plot frequency histogram vs. distance to nearest asteroid
    fig, ax = plt.subplots()
    ax.set_title(f'Distance to Nearest Asteroid for n={n} (Threshold {thresh_caption} )')
    ax.set_xlabel(axis_label_tbl[angle_unit])
    ax.set_ylabel('Observed Frequency')

    # bins in plot units and cartesian distance
    bins_x = np.linspace(0.0, thresh_x, bins+1)
    bins_dist = np.linspace(0.0, thresh_dist, bins+1)

    # Predicted number of points in each bin
    N_pred_cum = cdf_nearest_dist(dist=bins_dist, n=n, thresh_deg=180.0) * N_obs
    N_pred = np.diff(N_pred_cum)    
    # Prepare (x_ran, y_ran) to make a high resolution plot-style chart of this
    dx_ran = thresh_x / 1000
    x_ran = np.arange(dx_ran, thresh_x, dx_ran)
    idx = np.searchsorted(bins_x, x_ran, side='right') - 1
    y_ran = N_pred[idx]
    
    # Histogram of observed frequency in each bin
    ax.hist(x=plot_x, bins=bins_x, color='blue', label='observed')
    # Line chart with what the histogram would look like in the baseline with uniformly random distribution
    ax.plot(x_ran, y_ran, color='red', label='random')

    ax.legend()
    ax.grid()
    fig.savefig(f'../figs/ztf/nearest_ast_hist_dist_n={n}_thresh={thresh_str}.png', bbox_inches='tight')
    plt.show()

    return fig, ax

# ********************************************************************************************************************* 
def plot_rel_density(ztf, n: int, thresh_deg: float = 1.0, bins=20, chart_type: str = 'rel'):
    """
    Generate two CDFs of asteroid distance.
    INPUTS:
        ztf: DataFrame with distance to nearest asteroid
        n:   Number of asteroids, e.g. 16
        thresh_deg: Threshold in degrees for close observations
        chart_type: Type of density to plot.  Either 'rel' or 'abs'.
                    'rel' : relative density to theoretical
                    'abs' : absolute density in hits per square degree
    """
    # Number of rows in data
    N_obs = ztf.shape[0]

    # Threshold distance and flag indicating whether observations are within threshold
    thresh_dist = deg2dist(thresh_deg)
    is_close = ztf.nearest_ast_dist < thresh_dist

    # Cartesian distance of close observations
    dist = ztf[is_close].nearest_ast_dist.values
    # Distance in degrees and arc seconds
    dist_deg = dist2deg(dist)
    # dist_sec = 3600.0 * dist_deg

    # String description of threshold
    thresh_caption, thresh_str, angle_unit = angle_to_str(angle_deg=thresh_deg)

    # data for x axis; scale according to selected units
    # key = angle_unit, value = (plot data, threshold)
    plot_x_tbl = {
        'deg': (dist_deg, thresh_deg,),
        'min': (dist_deg * 60.0, thresh_deg * 60.0,),
        'sec': (dist_deg * 3600.0, thresh_deg * 3600.0),
    }
    plot_x, thresh_x = plot_x_tbl[angle_unit]

    # x-axis label depends on unit choice
    # key = angle_unit, value = axis label
    axis_label_tbl = {
        'deg': 'Distance in Degrees',
        'min': 'Distance in Arc Minutes',
        'sec': 'Distance in Arc Seconds',
    }

    # Set caption and y_label based on chart type
    caption = {
        'rel': 'Relative Density',
        'abs': 'Absolute Density',
    }[chart_type]

    y_label = {
        'rel': 'Relative Density (Observed / Random)',
        'abs': 'Absolute Density (Hit Prob. per Square Degree)',
    }[chart_type]

    # Plot frequency histogram vs. distance to nearest asteroid
    fig, ax = plt.subplots()
    ax.set_title(f'{caption}: Distance to Nearest Asteroid for n={n}')
    ax.set_xlabel(axis_label_tbl[angle_unit])
    ax.set_ylabel(y_label)

    # bins in plot units and cartesian distance
    bins_x = np.linspace(0.0, thresh_x, bins+1)
    bins_dist = np.linspace(0.0, thresh_dist, bins+1)

    # Predicted number of points in each bin
    N_pred_cum = cdf_nearest_dist(dist=bins_dist, n=n, thresh_deg=180.0) * N_obs
    N_pred = np.diff(N_pred_cum) 

    # Relative density in each bin
    hist, edges = np.histogram(a=plot_x, bins=bins_x)
    dens_rel = hist / N_pred

    # Absolute density; hits per deg^2
    sq_deg_in_sky = 4.0 * np.pi * np.rad2deg(1.0)**2
    area_cum = cdf_nearest_dist(dist=bins_dist, n=1)*sq_deg_in_sky
    area_bin = np.diff(area_cum)
    dens_abs = hist / area_bin / N_obs
    # log_dens_abs = np.log(dens_abs)

    # Histogram of observed frequency in each bin
    x_bar = bins_x[0:bins]
    width = thresh_x / bins

    # Plot relative or absolute
    if chart_type == 'rel':
        ax.bar(x=x_bar, height=dens_rel, width=width, align='edge', color='blue', label='relative')
    elif chart_type == 'abs':
        ax.bar(x=x_bar, height=dens_abs, width=width, align='edge', color='blue', label='absolute')
        ax.set_yscale('log')

    # ax.legend()
    ax.grid()
    fig.savefig(f'../figs/ztf/nearest_ast_hist_dens_{chart_type}_n={n}_thresh={thresh_str}.png', bbox_inches='tight')
    plt.show()

    return fig, ax

# ********************************************************************************************************************* 
def plot_close_freq(ztf, n: int, thresh_deg: float, is_cum=True, bins=100):
    """
    Generate plot with histogram of distance to the nearest asteroid vs. random baseline.
    INPUTS:
        ztf: DataFrame with distance to nearest asteroid
        n:   Number of asteroids, e.g. 16000
        thresh_deg: Threshold in degrees for close observations
    """
    # Threshold distance and flag indicating whether observations are within threshold
    thresh_dist = deg2dist(thresh_deg)
    is_close = ztf.nearest_ast_dist < thresh_dist
    # View of ztf limited to close observations
    ztfc = ztf[is_close]

    # Group close observations by asteroid number
    close_by_ast = ztfc.groupby(ztfc.nearest_ast_num)
    close_by_ast_count = close_by_ast.size()
    close_ast_num = close_by_ast_count.index.values
    close_ast_count = close_by_ast_count.values

    # String description of threshold
    thresh_caption, thresh_str, angle_unit = angle_to_str(angle_deg=thresh_deg)

    # Title and axis labels
    cum_prefix = f'Cumulative ' if is_cum else ''
    chart_title = f'{cum_prefix}Frequency of Close Observations by Asteroid (Threshold {thresh_caption})'
    chart_x_label = f'Number of Close Observations < {thresh_caption} for one Asteroid'
    chart_y_label = f'{cum_prefix}Frequency (Number of Asteroids)'

    # Plot unconditional histogram
    fig, ax = plt.subplots()
    ax.set_title(chart_title)
    ax.set_xlabel(chart_x_label)
    ax.set_ylabel(chart_y_label)
    ax.set_xticks(np.arange(0, bins+1, 10))
    bins_np = np.arange(bins+1)
    ax.hist(close_ast_count, color='blue', bins=bins_np, cumulative=is_cum)
    # ax.legend()
    ax.grid()
    cum_str = f'_cum' if is_cum else ''
    fig.savefig(f'../figs/ztf/nearest_ast_count_by_dist{cum_str}_n={n}_thresh={thresh_str}.png', bbox_inches='tight')
    plt.show()
