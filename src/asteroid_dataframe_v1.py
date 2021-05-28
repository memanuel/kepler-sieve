"""
Harvard IACS Masters Thesis
Load integrated asteroid trajectories as Pandas DataFrames

Michael S. Emanuel
Sat Sep 21 10:38:38 2019
"""

# Core
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

# Astronomy
import astropy
from astropy.units import au, day

# Utility
from datetime import datetime
from tqdm.auto import tqdm

# Local imports
from utils import range_inc
from astro_utils import datetime_to_mjd
from asteroid_element import load_ast_elt
from asteroid_integrate import calc_ast_pos_all
from rebound_utils import load_sim_np
from ra_dec import qv2dir, dir2radec, calc_topos, astrometric_dir, direction_diff

# Type names
from typing import Optional, Tuple, Dict

# ********************************************************************************************************************* 
# DataFrame of asteroid snapshots
ast_elt = load_ast_elt()
# Size of the blocks in asteroid integrator files
ast_block_size: int = 1000
# space_dims = 3



# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# ********************************************************************************************************************* 
def load_ast_data_block(block: int, 
                        mjd0: Optional[float]=None, 
                        mjd1: Optional[float]=None) -> \
                        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the MSE asteroid integrations for this range of asteroids.
    INPUTS:
        block: Block of asteroid data, e.g. 0 for asteroids [0, 1000)
               Integrator saves them in blocks of 1000 in data\asteroids
        mjd0:  Start modfified julian date used to filter output.
        mjd1:  Last modified julian date used to filter output.
               Default for mjd0 and mjd1 is None; then return all available time steps
    OUTPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
        df_earth: Position & velocity of earth in barycentric frame; heliocentric orbital elements
        df_sun:   Position & velocity of sun in barycentric frame
    """

    # First asteroid in the block
    n0: int = ast_block_size * (block + 0)
    n1: int = ast_block_size * (block + 1)

    # Name of the numpy archive
    fname_np: str = f'../data/asteroids/sim_asteroids_n_{n0:06}_{n1:06}.npz'

    # The full array of positions and velocities
    q, v, elts, catalog = load_sim_np(fname_np=fname_np)

    # The object names
    object_names = catalog['object_names']

    # The snapshot times; offset to start time t0=0
    mjd_all = catalog['ts']
    # Convert ts from relative time vs. t0 to MJD
    dt0 = datetime(2000, 1, 1)
    t_offset = datetime_to_mjd(dt0)
    mjd_all += t_offset

    # mask for selected time interval
    mjd0 = 0.0 if mjd0 is None else mjd0
    mjd1 = 1E9 if mjd1 is None else mjd1
    mask_t = (mjd0 <= mjd_all) & (mjd_all < mjd1)
    # Apply filter to find mjd that will be sent back
    mjd = mjd_all[mask_t]

    # Times as astropy Time class
    ts_ap = astropy.time.Time(mjd, format='mjd')
    # Times as julian dates 
    jd = ts_ap.jd
    # Times as integer key counting the number of hours in MJD era
    # time_key = np.int32(np.round(mjd*24*3600))

    # mask for selected asteroids
    mask_ast = (n0 <= ast_elt.Num) & (ast_elt.Num < n1)

    # count of selected asteroids
    asteroid_name = ast_elt.Name[mask_ast].values
    asteroid_num = ast_elt.Num[mask_ast].values
    N_ast: int = np.sum(mask_ast)
    # offset for indexing into asteroids; the first [10] objects are sun and planets
    ast_offset: int = len(object_names) - N_ast

    # filter for only selected times
    q = q[mask_t, :, :]
    v = v[mask_t, :, :]

    # Position of the sun in barycentric coordinates
    sun_idx = 0
    q_sun = q[:, sun_idx, :]
    v_sun = v[:, sun_idx, :]
    # Position of the earth in barycentric coordinates
    earth_idx = 3
    q_earth = q[:, earth_idx, :]
    v_earth = v[:, earth_idx, :]

    # Position of the asteroids in barycentric coordinates
    q_ast = q[:, ast_offset:, :]
    v_ast = v[:, ast_offset:, :]

    # Swap axes so they will be indexed (asteroid_num, time_step, space_dim)
    # initial order is (time_step, asteroid_number, space_dim)
    q_ast = np.swapaxes(q_ast, 0, 1)
    v_ast = np.swapaxes(v_ast, 0, 1)    

    # Assemble concatenated 3d arrays for splining
    # qv_ast = np.concatenate([q_ast, v_ast], axis=2)

    # number of times
    N_t = mjd.size
    # number of rows in df_ast
    N_row = N_ast * N_t

    # Arrays for orbital elements; restrict to selected times, index into asteroids
    orb_a = elts.a[mask_t, ast_offset:]
    orb_e = elts.e[mask_t, ast_offset:]
    orb_inc = elts.inc[mask_t, ast_offset:]
    orb_Omega = elts.Omega[mask_t, ast_offset:]
    orb_omega = elts.omega[mask_t, ast_offset:]
    orb_f = elts.f[mask_t, ast_offset:]

    # asteroid DataFrame
    ast_dict = {
        # indexing : asteroid number and time
        'asteroid_num': np.repeat(asteroid_num, N_t),
        'mjd': np.tile(mjd, N_ast),
        # 'time_key' : np.tile(time_key, N_ast),
        # cartesian coordinates
        'qx': q_ast[:, :, 0].reshape(N_row),
        'qy': q_ast[:, :, 1].reshape(N_row),
        'qz': q_ast[:, :, 2].reshape(N_row),
        'vx': v_ast[:, :, 0].reshape(N_row),
        'vy': v_ast[:, :, 1].reshape(N_row),
        'vz': v_ast[:, :, 2].reshape(N_row),
        # orbital elements
        'a':     orb_a.swapaxes(0,1).reshape(N_row),
        'e':     orb_e.swapaxes(0,1).reshape(N_row),
        'inc':   orb_inc.swapaxes(0,1).reshape(N_row),
        'Omega': orb_Omega.swapaxes(0,1).reshape(N_row),
        'omega': orb_omega.swapaxes(0,1).reshape(N_row),
        'f':     orb_f.swapaxes(0,1).reshape(N_row),    
    }
    df_ast = pd.DataFrame(ast_dict)

    # earth DataFrame
    earth_dict = {
        # indexing
        'mjd': mjd,
        # 'time_key' : time_key,
        # cartesian coordinate
        'qx': q_earth[:, 0],
        'qy': q_earth[:, 1],
        'qz': q_earth[:, 2],
        'vx': v_earth[:, 0],
        'vy': v_earth[:, 1],
        'vz': v_earth[:, 2],
        # orbital elements
        'a':     elts.a[mask_t, earth_idx],
        'e':     elts.e[mask_t, earth_idx],
        'inc':   elts.inc[mask_t, earth_idx],
        'Omega': elts.Omega[mask_t, earth_idx],
        'omega': elts.omega[mask_t, earth_idx],
        'f':     elts.f[mask_t, earth_idx],
    }
    df_earth = pd.DataFrame(earth_dict)

    # sun DataFrame
    sun_dict = {
        # indexing
        'mjd': mjd,
        # 'time_key' : time_key,
        # cartesian coordinate
        'qx': q_sun[:, 0],
        'qy': q_sun[:, 1],
        'qz': q_sun[:, 2],
        'vx': v_sun[:, 0],
        'vy': v_sun[:, 1],
        'vz': v_sun[:, 2],
    }
    df_sun = pd.DataFrame(sun_dict)

    return df_ast, df_earth, df_sun

# ********************************************************************************************************************* 
def load_ast_data(n0: int, n1: int, 
                  mjd0: Optional[float]=None, 
                  mjd1: Optional[float]=None,
                  progbar: bool = True) -> \
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the MSE asteroid integrations for this range of asteroids.
    INPUTS:
        n0:  First asteroid to load, e.g. 0
        n1:  Last asteroid to load, (exclusive) e.g. 1000
        mjd0: Start modfified julian date used to filter output.
        mjd1: Last modified julian date used to filter output.
              Default for mjd0 and mjd1 is None; then return all available time steps
        progbar: Whether to print a progress bar to the console
    OUTPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
        df_earth: Position & velocity of earth in barycentric frame; heliocentric orbital elements
        df_sun:   Position & velocity of sun in barycentric frame
    """
    # Range of blocks required
    block_min = n0 // ast_block_size
    block_max = (n1-1) // ast_block_size
    num_blocks = block_max - block_min + 1

    # List of asteroid frames for each data block
    dfs_ast = np.empty(num_blocks, dtype=object)

    # Set up iterates and prepare progress bar if requested
    iterates = list(enumerate(range_inc(block_min, block_max)))
    if progbar:
        print(f'Loading asteroid data from n0={n0} to n1={n1} in {num_blocks} blocks...')
        iterates = tqdm(iterates)

    # Iterate through the blocks
    for i, block in iterates:
        # Load the data for this block
        df_ast_i, df_earth_i, df_sun_i = load_ast_data_block(block=block, mjd0=mjd0, mjd1=mjd1)
        # If it's the first or last block, filter it necessary
        if i in (0, num_blocks-1):
            mask = (n0 <= df_ast_i.asteroid_num) & (df_ast_i.asteroid_num < n1)
            df_ast_i = df_ast_i[mask]
        # Save asteroid frame to list of frames
        dfs_ast[i] = df_ast_i

    # Concatenate frames for asteroids only! Earth and sun just one copy on the distinct dates
    df_ast = pd.concat(dfs_ast)
    # Just copy the last frame for earth and sun
    df_earth = df_earth_i
    df_sun = df_sun_i
    
    return df_ast, df_earth, df_sun
    
# ********************************************************************************************************************* 
def calc_ast_data(elts: Dict[str, np.ndarray], 
                  mjd0: np.float,
                  mjd1: np.float,
                  element_id: Optional[np.ndarray]=None,
                  progbar: bool = True) -> \
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build asteroid integration for this block orbital elements on the fly.
    INPUTS:
        elts:       Dictionary of orbital elements as numpy arrays
        mjd0:       Start date as an mjd
        mjd1:       End date an mjd
        element_id: Optional array of element_ids for describing / saving elements
        progbar:    Whether to print a progress bar to the console
    OUTPUTS:
        df_ast:     Position & velocity of asteroids in barycentric frame
        df_earth:   Position & velocity of earth in barycentric frame
        df_sun:     Position & velocity of sun in barycentric frame
    """
    # Extract the epocha
    epoch = elts['epoch'].iloc[0]

    # Daily snapshots in range of data
    ts = np.arange(mjd0, mjd1, dtype=np.float64)

    # Number of asteroids in this batch of elements
    N_elt = elts['a'].size

    # Number of times and number of rows in flattened DataFrame
    N_t = ts.size
    N_row = N_elt * N_t

    # Element ID in this batch; counts from e.g. 0 to 63
    if element_id is None:
        element_id = np.arange(N_elt, dtype=np.int32)
    # Whether or not element_id was provided, it needs to be repeated to N_t times to line up with flattened frame.
    element_id = np.repeat(element_id, ts.size)

    # The mjd column of the DataFrame; times must be tiled N_elt times for asteroids
    mjd_ast = np.tile(ts, N_elt)    
    mjd_earth = ts
    mjd_sun = ts

    # Calculate asteroid positions on batch
    pos_tbl = calc_ast_pos_all(elts=elts, epoch=epoch, ts=ts)

    # Extract asteroid position and velocity
    q_ast = pos_tbl['q_ast']
    v_ast = pos_tbl['v_ast']
    # Reshape q and v; now flat data frames, indexing element_id first, then timestamp_id
    q_ast = q_ast.reshape(N_row, 3)
    v_ast = v_ast.reshape(N_row, 3)

    # Extract asteroid position and velocity of earth and sun
    q_earth = pos_tbl['q_earth']
    v_earth = pos_tbl['v_earth']
    q_sun = pos_tbl['q_sun']
    v_sun = pos_tbl['v_sun']

    # Build asteroid DataFrame
    df_ast = pd.DataFrame({
        'element_id':element_id, 
        'mjd': mjd_ast,
    })

    # Build earth and sun DataFrame
    df_earth = pd.DataFrame({
        'mjd': mjd_earth,
    })
    df_sun = pd.DataFrame({
        'mjd': mjd_earth,
    })
    
    # Columns for position and velocity
    cols_q = ['qx', 'qy', 'qz']
    cols_v = ['vx', 'vy', 'vz']

    # Save position columns
    for k, col in enumerate(cols_q):
        df_ast[col] = q_ast[:, k]
        df_earth[col] = q_earth[:, k] 
        df_sun[col] = q_sun[:, k]

    # Save velocity columns
    for k, col in enumerate(cols_v):
        df_ast[col] = v_ast[:, k] 
        df_earth[col] = v_earth[:, k]
        df_sun[col] = v_sun[:, k] 

    return df_ast, df_earth, df_sun

# ********************************************************************************************************************* 
def spline_ast_vec_df(df_ast: pd.DataFrame, df_earth: pd.DataFrame, df_sun: pd.DataFrame, mjd: np.ndarray,
                      include_elts: bool = True, progbar: bool = True) \
                      -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the MSE asteroid integrations for this range of asteroids.
    INPUTS:
        df_ast:    DataFrame with asteroid position and velocity
        df_earth:  DataFrame with earth position and velocity
        df_sun:    DataFrame with earth position and velocity
        mjd:       Array of times at which splined output is desired
        include_elts: Whether to include orbital elements for the asteroids and earth
        progbar:   Whether to print a progress bar to the console
    OUTPUTS:
        df_ast_out:   Position & velocity of asteroids in barycentric frame
        df_earth_out: Position & velocity of earth in barycentric frame
        df_sun_out:   Position & velocity of sun in barycentric frame
    """
    # Time key from mjd at spline points
    # time_key = np.int32(np.round(mjd*24*3600))

    # Distinct asteroids or elements; get the name of the relevant column and array of distinct values
    id_col = 'asteroid_num' if 'asteroid_num' in df_ast.columns else 'element_id'
    id_val_unq = np.unique(df_ast[id_col].values)

    # Number of asteroids or elements
    N_ast = id_val_unq.size
    # Number of times in input and splined output
    N_t_in = df_earth.mjd.size
    N_t_out = mjd.size
    # Number of rows in input and splined output
    N_row_in = N_ast * N_t_in
    N_row_out = N_ast * N_t_out

    # Data to be splined: x axis is time
    x_spline = df_earth.mjd

    # Desired output columns to be splined for asteroids, earth and sun
    cols_cart = ['qx', 'qy', 'qz', 'vx', 'vy', 'vz']
    cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f']
    cols_ast = cols_cart + cols_elt if include_elts else cols_cart
    cols_earth = cols_cart + cols_elt if include_elts else cols_cart
    cols_sun = cols_cart

    # the indexing of y_spline_ast is (ast_num, time_step, data_dim) with shape e.g. (16, 3653, 6)
    y_spline_ast = df_ast[cols_ast].values.reshape(N_ast, N_t_in, -1)
    # Earth and sun data to be splined don't need to be reshaped
    y_spline_earth = df_earth[cols_earth].values
    y_spline_sun = df_sun[cols_sun].values

    # splining functions for asteroids, earth, and sun
    spline_func_ast = CubicSpline(x=x_spline, y=y_spline_ast, axis=1)
    spline_func_earth = CubicSpline(x=x_spline, y=y_spline_earth, axis=0)
    spline_func_sun = CubicSpline(x=x_spline, y=y_spline_sun, axis=0)

    # splined output for asteroids, earth and sun
    spline_data_ast = spline_func_ast(mjd)
    spline_data_earth = spline_func_earth(mjd)
    spline_data_sun = spline_func_sun(mjd)

    # ID column (asteroid_num or element_id) corresponding to splined output
    id_val = np.repeat(id_val_unq, N_t_out)

    # asteroid DataFrame
    ast_dict_keys = {
        id_col: id_val,
        'mjd': np.tile(mjd, N_ast),
        # 'time_key' : np.tile(time_key, N_ast)
    }
    ast_dict_data = {col:spline_data_ast[:,:,k].reshape(N_row_out) for k, col in enumerate(cols_ast)}
    ast_dict = dict(**ast_dict_keys, **ast_dict_data)
    df_ast_out = pd.DataFrame(ast_dict)

    # earth DataFrame
    earth_dict_keys = {
        'mjd': mjd,
        # 'time_key': time_key
    }
    earth_dict_data = {col:spline_data_earth[:,k] for k, col in enumerate(cols_earth)}
    earth_dict = dict(**earth_dict_keys, **earth_dict_data)
    df_earth_out = pd.DataFrame(earth_dict)

    # sun DataFrame
    sun_dict_keys = earth_dict_keys
    sun_dict_data = {col:spline_data_sun[:,k] for k, col in enumerate(cols_sun)}
    sun_dict = dict(**sun_dict_keys, **sun_dict_data)
    df_sun_out = pd.DataFrame(sun_dict)
    
    return df_ast_out, df_earth_out, df_sun_out
    
# ********************************************************************************************************************* 
def spline_ast_vec(n0: int, n1: int, mjd: np.ndarray, 
                   include_elts: bool = True, progbar: bool = True) \
                   -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the MSE asteroid integrations for this range of asteroids.
    INPUTS:
        n0:  First asteroid to load, inclusive, e.g. 0
        n1:  Last asteroid to load, exclusive, e.g. 1000
        mjd: Array of modified julian dates on which splined output is desired
        progbar: Whether to print a progress bar to the console
    OUTPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
        df_earth: Position & velocity of earth in barycentric frame; heliocentric orbital elements
        df_sun:   Position & velocity of sun in barycentric frame
    """
    # Compute mjd0 and mjd1 from mjd
    mjd0 = np.min(mjd) - 1
    mjd1 = np.max(mjd) + 1

    # Load data in this date range
    df_ast, df_earth, df_sun = load_ast_data(n0=n0, n1=n1, mjd0=mjd0, mjd1=mjd1, progbar=progbar)

    # Delegate to spline_ast_vec_df
    return spline_ast_vec_df(df_ast=df_ast, df_earth=df_earth, df_sun=df_sun, mjd=mjd, 
                             include_elts=include_elts, progbar=progbar)

# ********************************************************************************************************************* 
def calc_ast_dir(df_ast: pd.DataFrame, df_earth: pd.DataFrame, site_name: str) -> pd.DataFrame:
    """
    Generate a DataFrame of predicted directions given positions of asteroids and earth, and a site name.
    INPUTS:
        df_ast:    Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
        df_earth:  Position & velocity of earth in barycentric frame; heliocentric orbital elements
        site_name: Name of observatory site, e.g. 'geocenter' or 'palomar'
    """
    # distinct asteroids or elements; get the name of the relevant column and array of distinct values
    id_col = 'asteroid_num' if 'asteroid_num' in df_ast.columns else 'element_id'
    id_val = df_ast[id_col].values

    # observation times
    mjd_ast = df_ast.mjd.values
    mjd_earth = df_earth.mjd.values

    # extract position and velocity of asteroids
    cols_q = ['qx', 'qy', 'qz']
    cols_v = ['vx', 'vy', 'vz']
    q_ast = df_ast[cols_q].values * au
    v_ast = df_ast[cols_v].values * au / day

    # number of asteroids
    # asteroid_num = df_ast.asteroid_num.values
    # asteroid_num_unq = np.unique(asteroid_num)
    # N_ast = asteroid_num_unq.size
    id_val_unq = np.unique(id_val)
    N_ast = id_val_unq.size

    # position of earth and topos adjustment
    q_earth_once = df_earth[cols_q].values * au
    # Topos adjustment
    dq_topos = calc_topos(obstime_mjd=mjd_earth, site_name=site_name)
    # position of the observatory in ecliptic frame
    q_obs_once = q_earth_once + dq_topos

    # tile q_obs so it can be subtracted from q_ast
    q_obs = np.tile(q_obs_once, (N_ast,1))    

    # calculate the astrometric direction with the observer position
    u_ast, delta = astrometric_dir(q_body=q_ast, v_body=v_ast, q_obs=q_obs)

    # extract components from u
    ux = u_ast[:, 0]
    uy = u_ast[:, 1]
    uz = u_ast[:, 2]

    # compute RA/DEC from direction
    ra, dec = dir2radec(u=u_ast, obstime_mjd=mjd_ast)

    # build the observation DataFrame
    dir_dict = {
        id_col: id_val,
        'mjd' : mjd_ast,
        # 'time_key': np.int32(np.round(mjd_ast*24)),
        'ra': ra.value,
        'dec': dec.value,
        'ux': ux,
        'uy': uy,
        'uz': uz,
        'delta': delta
    }
    df_dir = pd.DataFrame(dir_dict)
    
    return df_dir

# ********************************************************************************************************************* 
def spline_ast_vec_dir(n0: int, n1: int, mjd: np.ndarray, 
                       site_name: str = 'geocenter', progbar: bool = True) -> \
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Convenience API to (1) obtain splined vectors (2) observed angles
    INPUTS:
        n0:  First asteroid to load, e.g. 1
        n1:  Last asteroid to load, e.g. 64
        mjd: Array of modified julian dates on which splined output is desired
        site_name: Name of observatory; defaults to 'geocenter'
        progbar: Whether to print a progress bar to the console
    OUTPUTS:
        df_ast:   Position & velocity of asteroids in barycentric frame; heliocentric orbital elements
        df_earth: Position & velocity of earth in barycentric frame; heliocentric orbital elements
        df_obs:   Projected observation of asteroid from earth at given observatory
    """

    # Calculate splined vectors from spline_ast_data
    print(f'Splining asteroid and earth vectors for asteroid numbers {n0} to {n1}...')
    df_ast, df_earth, df_sun = spline_ast_vec(n0=n0, n1=n1, mjd=mjd, progbar=progbar)

    # Calculate predicted RA / DEC and direction with spline_ast_obs
    print(f'Computing astrometric asteroid directions for asteroid numbers {n0} to {n1}...')
    df_dir = calc_ast_dir(df_ast=df_ast, df_earth=df_earth, site_name=site_name)

    # Return earth and asteroid vectors plus aligned observations
    return df_ast, df_earth, df_dir

# ********************************************************************************************************************* 
def compare_df_vec(df_mse, df_jpl, name: str):
    """Compare DataFrames with MSE vs. JPL vectors"""
    # Columnwise mean absolute error
    err_mjd = np.mean(np.abs(df_mse.mjd - df_jpl.mjd))
    # err_jd = np.mean(np.abs(df_mse.jd - df_jpl.JulianDate))

    # norm of position error
    q_mse = df_mse[['qx', 'qy', 'qz']].values
    q_jpl = df_jpl[['X', 'Y', 'Z']].values
    err_q = np.linalg.norm(q_mse-q_jpl, axis=1)
    err_q_mean = np.mean(err_q)
    err_q_max = np.max(err_q)

    # norm of velocity error
    v_mse = df_mse[['vx', 'vy', 'vz']].values
    v_jpl = df_jpl[['VX', 'VY', 'VZ']].values
    err_v = np.linalg.norm(v_mse-v_jpl, axis=1)
    err_v_mean = np.mean(err_v)
    v_jpl_mean = np.mean(np.linalg.norm(v_jpl, axis=1))
    err_v_rel = err_v_mean / np.mean(np.linalg.norm(v_jpl, axis=1))

    print(f'Mean absolute error for df_{name}_mse vs. df_{name}_jpl:')
    print(f'mjd: {err_mjd:6.2e} days')
    # print(f' jd: {err_jd:6.2e} days')
    print(f'  q: {err_q_mean:6.2e} AU     (max {err_q_max:6.2e})')
    print(f'  v: {err_v_mean:6.2e} AU/day (rel {err_v_rel:6.2e})')

# ********************************************************************************************************************* 
def compare_df_dir(df_mse, df_jpl, name: str):
    """Compare DataFrames with MSE vs. JPL observations"""
    # Columnwise mean absolute error
    err_mjd = np.mean(np.abs(df_mse.mjd - df_jpl.mjd))

    # direction from both frames
    u_mse = df_mse[['ux', 'uy', 'uz']].values
    u_jpl = df_jpl[['ux_jpl', 'uy_jpl', 'uz_jpl']].values

    print(f'Mean absolute error for {name} observations: MSE vs. JPL')
    print(f'mjd: {err_mjd:6.2e} days')

    direction_diff(name1='mse', name2='jpl', u1=u_mse, u2=u_jpl, verbose=True)