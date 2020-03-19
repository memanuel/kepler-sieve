"""
Harvard IACS Masters Thesis
asteroid_search.py: Search for orbital elements of asteroids given observational data.

Michael S. Emanuel
Thu Oct 17 15:24:10 2019
"""

# Core
import numpy as np
import pandas as pd

# Tensorflow / ML
import tensorflow as tf
from tensorflow.python.keras import backend as K

# Utility
import time
from datetime import timedelta

# Local imports
from asteroid_search_model import make_model_asteroid_search
from ztf_data import load_ztf_easy_batch
from asteroid_data import make_ztf_dataset, orbital_element_batch
from asteroid_integrate import calc_ast_pos
from asteroid_search_report import report_model, report_training_progress
from utils import print_header

# Typing
from typing import Dict

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
def perturb_elts(elts, sigma_a=0.05, sigma_e=0.10, sigma_f_deg=5.0, mask=None):
    """Apply perturbations to orbital elements"""
    # Copy the elements
    elts_new = elts.copy()

    # Default for mask is all elements
    if mask is None:
        mask = np.ones_like(elts['a'], dtype=bool)

    # Number of elements to perturb
    num_shift = np.sum(mask)
    
    # Apply shift log(a)
    log_a = np.log(elts['a'])
    log_a[mask] += np.random.normal(scale=sigma_a, size=num_shift)
    elts_new['a'] = np.exp(log_a)
    
    # Apply shift to log(e)
    log_e = np.log(elts['e'])
    log_e[mask] += np.random.normal(scale=sigma_e, size=num_shift)
    elts_new['e'] = np.exp(log_e)
    
    # Apply shift directly to true anomaly f
    f = elts['f']
    sigma_f = np.deg2rad(sigma_f_deg)
    f[mask] += np.random.normal(scale=sigma_f, size=num_shift)
    elts_new['f'] = f
    
    return elts_new

# ********************************************************************************************************************* 
def test_easy_batch(time_batch_size: int = 128, elt_batch_size: int = 64, epochs=20):
    """
    Test asteroid search on ZTF easy batch
    """
    # Load all ZTF data with nearest asteroid calculations
    ztf, elts = load_ztf_easy_batch(batch_size=elt_batch_size)

    # Build a TensorFlow DataSet from ZTF DataFrame
    ds, ts, row_len = make_ztf_dataset(ztf=ztf, batch_size=time_batch_size)

    # Trajectory size and steps per batch
    traj_size = ts.shape[0]
    # If time_batch_size was None, use all the time points in each batch
    if time_batch_size is None:
        time_batch_size = traj_size
    steps = int(np.ceil(traj_size / time_batch_size))

    # Resolution and threshold in degrees
    R_deg: float = 0.5
    thresh_deg: float = 1.0

    # Pop asteroid number from and epoch from elts DataFrame
    ast_nums = elts.pop('ast_num')    
    epoch = elts.pop('epoch')[0]

    # The correct orbital elements as an array of shape Nx6
    elts_true = elts.values

    # Batch of perturbed orbital elements for asteroid model
    # ast_nums = np.unique(ztf.nearest_ast_num)
    # elts_np = orbital_element_batch(ast_nums)
    # epoch = elts_np['epoch'][0]

    # Get example batch
    batch_in, batch_out = list(ds.take(1))[0]
    # Contents of this batch
    t = batch_in['t']
    idx = batch_in['idx']
    row_len = batch_in['row_len']
    u_obs = batch_in['u_obs']

    # Get max_obs and number of observations
    max_obs: int = u_obs.shape[1]
    # The number of observations is the TOTAL FOR THE ZTF DATA SET!
    # It's not just the size of this easy batch, b/c the easy batch has been harvested to be close!
    # num_obs: float = np.sum(row_len, dtype=np.float32)
    num_obs: float = 5.7E6

    # The correct orbital elements as an array
    # elts_true = np.array([elts_np['a'], elts_np['e'], elts_np['inc'], elts_np['Omega'], 
    #                      elts_np['omega'], elts_np['f'], elts_np['epoch']]).transpose()

    # Mask where data expected vs not
    mask_good = np.arange(elt_batch_size) < (elt_batch_size//2)
    mask_bad = ~mask_good
    # Perturb second half of orbital elements
    # elts_np2 = perturb_elts(elts_np, sigma_a=0.00, sigma_e=0.00, sigma_f_deg=0.0, mask=mask_bad)
    elts_np2 = perturb_elts(elts_np, mask=mask_bad)

    # Orbits for calibration
    if 'q_cal' not in globals():
        print(f'Numerically integrating calibration trajectories q_cal...')
        q_cal = calc_ast_pos(elts=elts_np2, epoch=epoch, ts=ts)
    # q_cal = None

    # Set calibration flag
    use_calibration: bool = True

    # Alpha and beta parameters for the objective function
    alpha = 1.0
    beta = 0.0

    # Build functional model for asteroid score
    model = make_model_asteroid_search(\
        ts=ts, elts_np=elts_np2, max_obs=max_obs, num_obs=num_obs,
        elt_batch_size=elt_batch_size, time_batch_size=time_batch_size,
        R_deg=R_deg, thresh_deg=thresh_deg, alpha=alpha, beta=beta, 
        q_cal=q_cal, use_calibration=use_calibration)

    # Use Adam optimizer with gradient clipping
    # learning_rate = 2.0e-5
    learning_rate = 2.0E-5  # default 1.0E-3
    beta_1 = 0.900          # default 0.900
    beta_2 = 0.999          # default 0.999
    epsilon = 1.0E-7        # default 1.0E-7
    amsgrad = False         # default False
    clipvalue = 5.0         # default not used
    opt = keras.optimizers.Adam(learning_rate=learning_rate, 
                                beta_1=beta_1, 
                                beta_2=beta_2,
                                epsilon=epsilon, 
                                amsgrad=amsgrad,
                                # clipvalue=clipvalue, 
                                )
    model.compile(optimizer=opt)
    # Whether to display results
    display = False

    # Report losses before training
    print(f'Processed easy batch with {elt_batch_size} asteroids. Perturbed second half.')
    if display:
        print_header('Model Before Training:')
    pred0 = model.predict_on_batch(ds)
    elts0, R0, u_pred0, z0, _ = pred0
    scores0, traj_err0, elt_err0 = \
        report_model(model=model, ds=ds, R_deg=R_deg, mask_good=mask_good, 
                     batch_size=elt_batch_size, steps=steps, elts_true=elts_true, display=display)

    # Get intital gradients on entire data set
    with tf.GradientTape(persistent=True) as gt:
        gt.watch([model.elements.e_, model.elements.inc_, model.elements.R_])
        pred = model.predict_on_batch(ds.take(traj_size))
        elts, R, u_pred, z, scores = pred
        # unpack elements
        a = elts[:,0]
        e = elts[:,1]
        inc = elts[:,2]
        Omega = elts[:,3]
        omega = elts[:,4]
        f = elts[:,5]
        epoch = elts[:,6]
        # unpack scores
        raw_score = scores[:, 0]
        mu = scores[:, 1]
        sigma2 = scores[:, 2]
        objective = scores[:, 3]
        # total loss function
        loss = tf.reduce_sum(-objective)

    #    # Derivatives of elements w.r.t. control variables
    #    da_da_ = gt.gradient(a, model.elements.a_)
    #    de_de_ = gt.gradient(e, model.elements.e_)
    #    dinc_dinc_ = gt.gradient(inc, model.elements.inc_)
    #    dOmega_dOmega_ = gt.gradient(Omega, model.elements.Omega_)
    #    domega_domega_ = gt.gradient(omega, model.elements.omega_)
    #    df_df_ = gt.gradient(f, model.elements.f_)

    # Derivatives of loss w.r.t. control variables for elements and R
    dL_da = gt.gradient(loss, model.elements.a_) / steps
    dL_de = gt.gradient(loss, model.elements.e_) / steps
    dL_dinc = gt.gradient(loss, model.elements.inc_) / steps
    dL_dOmega = gt.gradient(loss, model.elements.Omega_) / steps
    dL_domega = gt.gradient(loss, model.elements.omega_) / steps
    dL_df = gt.gradient(loss, model.elements.f_) / steps
    dL_dR = gt.gradient(loss, model.elements.R_) / steps
    del gt

    # Train model
    step_multiplier = 5
    steps_per_epoch = steps*step_multiplier
    rows_per_epoch = steps_per_epoch * time_batch_size
    print_header(f'Training for {epochs} Epochs of Size {rows_per_epoch} Observation Days...')
    print(f'alpha:         {alpha:5.1f}')
    print(f'beta:          {beta:5.1f}')
    print(f'R (degrees):   {R_deg:5.1f}')
    print(f'learning_rate:   {learning_rate:5.2e}')
    print(f'clipvalue:      {clipvalue:5.2f}')

    train_time_0 = time.time()
    hist = model.fit(ds, epochs=epochs, steps_per_epoch=steps_per_epoch)
    train_time_1 = time.time()
    train_time = train_time_1 - train_time_0
    print(f'Elapsed Time: {str(timedelta(seconds=train_time))}')

    # Report results
    if display:
        print_header('Model After Training:')
    pred1 = model.predict_on_batch(ds)
    elts1, R1, u_pred1, z1, _ = pred1
    scores1, traj_err1, elt_err1 = \
        report_model(model=model, ds=ds, R_deg=R_deg, mask_good=mask_good, 
                     batch_size=elt_batch_size, steps=steps, elts_true=elts_true, display=display)

    # Unpack scores after training
    raw_score = scores1[:,0]
    mu = scores1[:,1]
    sigma2 = scores1[:,2]
    objective = scores1[:,3]
    sigma = scores1[:, 4]
    t_score = scores1[:, 5]
    eff_obs = raw_score - mu

    ## Change in scores
    d_scores = scores1 - scores0
    d_traj_err = traj_err1 - traj_err0
    d_elt_err = elt_err1 - elt_err0
    d_R = R1 - R0

    # Report training progress: scores, orbital element errors, and resolution
    print_header('Progress During Training:')
    # report_training_progress(d_scores, d_traj_err, d_elt_err, d_R)
    scores_01 = (scores0, scores1)
    traj_err_01 = (traj_err0, traj_err1)
    elt_err_01 = (elt_err0, elt_err1)
    R_01 = (R0, R1)
    # report_training_progress(d_scores, d_traj_err, d_elt_err, d_R)
    report_training_progress(scores_01, traj_err_01, elt_err_01, R_01, mask_good)

# ********************************************************************************************************************* 
if __name__ == '__main__':
    test_easy_batch(time_batch_size=None, elt_batch_size=64, epochs=5)

