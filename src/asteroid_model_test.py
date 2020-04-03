"""
Harvard IACS Masters Thesis
Solar Asteroid Model: Predict the movement of a test particle (e.g. asteroid) in the solar system
using the Kepler approximation with the sun as a fixed central attractor.

Michael S. Emanuel
Sun Oct 13 11:56:50 2019
"""

# Core
import tensorflow as tf
import numpy as np

# Utility
import time

# Local imports
from asteroid_model import make_model_ast_pos, make_model_ast_dir
# from asteroid_model import AsteroidPosition
from asteroid_data import make_dataset_ast_pos, make_dataset_ast_dir, make_dataset_ast_dir_spline
from asteroid_integrate import calc_ast_pos
from astro_utils import dist2deg, dist2sec

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
def test_ast_pos() -> bool:
    """Test asteroid position model"""
    # Load data for the first 1000 asteroids
    ds: tf.data.Dataset = make_dataset_ast_pos(n0=0, num_files=1, include_vel=True)
    # Get reference times
    batch_in, batch_out = list(ds.take(1))[0]
    ts = batch_in['ts'][0]
    batch_size = 64
    # Create the model to predict asteroid trajectories
    model: keras.Model = make_model_ast_pos(ts=ts, batch_size=batch_size)
    # Compile with MSE (mean squared error) loss
    model.compile(loss='MSE')
    # Evaluate this model
    mse_qv, mse_q, mse_v = model.evaluate(ds)
    # RMS errors for q and v
    rmse_q: float = np.sqrt(mse_q)
    rmse_v: float = np.sqrt(mse_v)
    # Threshold for passing
    thresh_q: float = 0.125
    thresh_v: float = 0.125
    isOK_mse: bool = (rmse_q < thresh_q) and (rmse_v < thresh_v)
    # Report results
    msg: str = 'PASS' if isOK_mse else 'FAIL'
    print(f'\nRoot MSE for asteroid model on first 1000 asteroids:')
    print(f'q in AU     = {rmse_q:8.6f}')
    print(f'v in AU/day = {rmse_v:8.6f}')
    print(f'***** {msg} *****')

    # Evaluate q, v on first batch before adjustments
    q_true = batch_out['q']
    v_true = batch_out['v']
    q_pred, v_pred = model.predict(batch_in)
    err_q_pre = np.mean(np.linalg.norm(q_pred - q_true, axis=2))
    err_v_pre = np.mean(np.linalg.norm(v_pred - v_true, axis=2))

    # Assemble orbital elements
    elts = {
        'a': batch_in['a'].numpy(),
        'e': batch_in['e'].numpy(),
        'inc': batch_in['inc'].numpy(),
        'Omega': batch_in['Omega'].numpy(),
        'omega': batch_in['omega'].numpy(),
        'f': batch_in['f'].numpy(),
        'epoch': batch_in['epoch'].numpy(),
        }
    # Epoch must be constant in a batch
    epoch = batch_in['epoch'][0].numpy()
    
    # Compute numerical orbit for calibration
    t0 = time.time()
    q_ast, q_earth, v_ast = calc_ast_pos(elts=elts, epoch=epoch, ts=ts)
    t1 = time.time()
    calc_time = t1 - t0
    print(f'\nNumerical integration of orbits took {calc_time:5.3f} seconds.')
    
    # Calibrate the position model
    model.ast_pos_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)
    # Evaluate error after calibration
    q_pred, v_pred = model.predict(batch_in)
    err_q_post = np.mean(np.linalg.norm(q_pred - q_true, axis=2))
    err_v_post = np.mean(np.linalg.norm(v_pred - v_true, axis=2))

    # Relative errors
    err_q_pre_rel = err_q_pre / np.mean(np.linalg.norm(q_true, axis=2))
    err_v_pre_rel = err_v_pre / np.mean(np.linalg.norm(v_true, axis=2))
    err_q_post_rel = err_q_post / np.mean(np.linalg.norm(q_true, axis=2))
    err_v_post_rel = err_v_post / np.mean(np.linalg.norm(v_true, axis=2))

    # Report results for position after calibration
    thresh_q = 2.0E-6
    isOK_q: bool = (err_q_post < thresh_q)
    msg = 'PASS' if isOK_q else 'FAIL'
    print(f'\nMean position error on first batch of 64 asteroids in AU:')
    print(f'Before calibration: {err_q_pre:5.3e} ({err_q_pre_rel:5.3e} relative)')
    print(f'After calibration:  {err_q_post:5.3e} ({err_q_post_rel:5.3e} relative)')
    print(f'***** {msg} *****')

    # Report results for velocity after calibration
    thresh_v = 1.0E-8
    isOK_v: bool = (err_v_post < thresh_v)
    msg = 'PASS' if isOK_v else 'FAIL'
    print(f'\nMean velocity error on first batch of 64 asteroids in AU/day:')
    print(f'Before calibration: {err_v_pre:5.3e} ({err_v_pre_rel:5.3e} relative)')
    print(f'After calibration:  {err_v_post:5.3e} ({err_v_post_rel:5.3e} relative)')
    print(f'***** {msg} *****')

    return (isOK_mse and isOK_q and isOK_v)

# ********************************************************************************************************************* 
def test_ast_dir(use_spline: bool = False, site_name: str = 'geocenter') -> bool:
    """Test the asteroid direction model"""
    # Set number of time steps and batch_size
    N_t: int = 1000
    batch_size: int = 64

    # If the source is not splined, we can only test at the geocenter
    if not use_spline:
        site_name = 'geocenter'

    # Load data for the first 1000 asteroids
    # Source depends on whether to use the spline or direct load from the integration output files
    ds: tf.data.Dataset
    if use_spline:
        ds = make_dataset_ast_dir_spline(n0=1, n1=65, site_name=site_name, N_t=N_t, batch_size=batch_size)
    else:
        ds = make_dataset_ast_dir(n0=0, num_files=1)

    # Get reference times
    batch_in, batch_out = list(ds.take(1))[0]
    ts = batch_in['ts'][0]
    # Create the model to predict asteroid trajectories
    model: keras.Model = make_model_ast_dir(ts=ts, site_name=site_name, batch_size=batch_size)
    # num_ast: int = 1000
    # Compile with MSE (mean squared error) loss
    model.compile(loss='MSE')
    # Evaluate this model
    mse_ur, mse_u, mse_r = model.evaluate(ds)
    rmse_u: float = np.sqrt(mse_u)
    rmse_r: float = np.sqrt(mse_r)
    # Convert error from unit vector to angle
    rmse_deg = dist2deg(rmse_u)
    rmse_sec = rmse_deg * 3600.0
    # Threshold for passing
    thresh_deg_pre: float = 1.0
    isOK_mse: bool = (rmse_deg < thresh_deg_pre)

    # Report results
    test_source = 'spline' if use_spline else 'file'
    print(f'\nTesting asteroid direction against {test_source} at observer site {site_name}.')
    msg: str = 'PASS' if isOK_mse else 'FAIL'
    print(f'\nMSE for asteroid model on first 1000 asteroids = {mse_u:8.6f}')
    print(f'Angle error = {rmse_u:5.3e} cart / {rmse_deg:8.6f} degrees / {rmse_sec:6.2f} arc seconds')
    print(f'***** {msg} *****')

    # Evaluate on first batch before adjustments
    u_true = batch_out['u']
    u_pred, r_pred = model.predict(batch_in)
    # Error in cartesian distance and arc seconds, before calibration
    err_dist_pre = np.mean(np.linalg.norm(u_pred - u_true, axis=2))
    err_deg_pre = dist2deg(err_dist_pre)
    err_sec_pre = err_deg_pre * 3600.0
    
    # Assemble orbital elements for calibration
    elts = {
        'a': batch_in['a'].numpy(),
        'e': batch_in['e'].numpy(),
        'inc': batch_in['inc'].numpy(),
        'Omega': batch_in['Omega'].numpy(),
        'omega': batch_in['omega'].numpy(),
        'f': batch_in['f'].numpy(),
        'epoch': batch_in['epoch'].numpy(),
        }
    # Epoch must be common in a batch to allow for numerical integration
    epoch = elts['epoch'][0]
    
    # Compute correction factor dq; also time this operation
    t0 = time.time()
    q_ast, q_earth, v_ast = calc_ast_pos(elts=elts, epoch=epoch, ts=ts)
    t1 = time.time()
    calc_time = t1 - t0
    print(f'\nNumerical integration of orbits took {calc_time:5.3f} seconds.')
    
    # Evaluate error after correction
    model.ast_dir_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)
    u_pred, r_pred = model.predict(batch_in)
    err_dist_post = np.mean(np.linalg.norm(u_pred - u_true, axis=2))
    err_deg_post = dist2deg(err_dist_post)
    err_sec_post = err_deg_post * 3600.0
    
    # Report results
    thresh_sec_post = 1.0
    isOK_post: bool = (err_sec_post < thresh_sec_post)
    msg = 'PASS' if isOK_post else 'FAIL'
    print(f'Mean angle error on first batch of 64 asteroids in degrees:')
    print(f'Before calibration: {err_deg_pre:5.3e} degrees ({err_sec_pre:8.3f} arc-seconds)' )
    print(f'After calibration:  {err_deg_post:5.3e} degrees ({err_sec_post:8.3f} arc-seconds)')
    print(f'***** {msg} *****')
    
    return (isOK_mse and isOK_post)

# ********************************************************************************************************************* 
def main():
    test_ast_pos()
    test_ast_dir(use_spline=False, site_name='geocenter')
    test_ast_dir(use_spline=True, site_name='palomar')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
