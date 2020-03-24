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
import argparse

# Local imports
from asteroid_search_model import make_model_asteroid_search
from ztf_data import load_ztf_easy_batch, make_ztf_batch, report_ztf_score
from asteroid_data import make_ztf_dataset, orbital_element_batch
from asteroid_integrate import calc_ast_pos, load_ast_elt
from asteroid_search_report import report_model, report_training_progress
from utils import print_header
from tf_utils import tf_quiet, gpu_grow_memory, get_gpu_device

# Typing
from typing import Dict

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Tensorflow config
# Run TF quietly
tf_quiet()

# Configure TensorFlow to use GPU memory variably; this is done in asteroid_model
# gpu_grow_memory(verbose=True)

# ********************************************************************************************************************* 
# Constants
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
def random_elts(size: int = 64, random_seed: int = 42):
    """Generate a DataFrame of random orbital elements"""
    pass

# ********************************************************************************************************************* 
def perturb_elts(elts: pd.DataFrame, sigma_a=0.05, sigma_e=0.10, sigma_f_deg=5.0, mask_pert=None, random_seed: int = 42):
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
    
    return elts_new

# ********************************************************************************************************************* 
def report_initial_scores(elts_pert: pd.DataFrame):
    """Report the initial scores of the perturbed elements"""
    # ZTF batch of perturbed elements
    file_path = '../data/ztf/ztf-easy-batch-pert.h5'
    try:
        ztfp = pd.read_hdf(file_path)
    except:
        ztfp = make_ztf_batch(elts=elts_pert, thresh_deg=thresh_deg, near_ast=True)
        ztfp.to_hdf(file_path, key='ztfp', mode='w')
    
    # The orbital element_id's
    element_id = np.unique(ztfp.element_id)
    elt_batch_size = element_id.size

    # True and perturbed half of ztfp
    elt_id_pert = element_id[elt_batch_size//2]  # start of perturbed element_ids
    mask_true = (ztfp.element_id < elt_id_pert)
    mask_pert = ~mask_true
    # ZTF including the scores on the true and perturbed part
    ztf_t = ztfp[mask_true]
    ztf_p = ztfp[mask_pert]

    # Report the scores on each part
    print(f'Unperturbed elements:')
    report_ztf_score(ztf_t)
    print(f'\nPerturbed elements:')
    report_ztf_score(ztf_p)

# ********************************************************************************************************************* 
def test_easy_batch(R_deg: float = 1.0, 
                    thresh_deg: float = 1.0,
                    R_is_trainable: bool = True,
                    elt_batch_size: int = 64, 
                    time_batch_size: int = None, 
                    epochs: int = 10,
                    cycles_per_epoch: int = 10,
                    use_calibration: bool = True):
    """
    Test asteroid search on ZTF easy batch
    INPUTS:
        R_deg:            resolution in degrees
        thresh_deg:       threshold in degrees
        R_is_trainable:   whether resolution R is trainable
        time_batch_size:  size of time_batch; default None uses all time snaps
        epochs:           number of epochs to train
        cycles_per_epoch: number of times to cycle through all data points before declaring an epoch
        use_calibration:  whether to calibrate trajectories with numerically integrated one
    OUTPUTS:
       scores_01, traj_err_01, elt_orr_01, R_01, mask_true
    """
    # Load all ZTF data with nearest asteroid calculations
    ztf, elts = load_ztf_easy_batch(batch_size=elt_batch_size, thresh_deg=thresh_deg)

    # Build a TensorFlow DataSet from ZTF DataFrame
    ds, ts, row_len = make_ztf_dataset(ztf=ztf, batch_size=time_batch_size)

    # Trajectory size and steps per batch
    traj_size = ts.shape[0]
    # If time_batch_size was None, use all the time points in each batch
    if time_batch_size is None:
        time_batch_size = traj_size
    steps = int(np.ceil(traj_size / time_batch_size))

    # The element_id is the same as the asteroid number on the easy batch
    element_id = elts.element_id.values
    # The epoch (scalar)
    epoch = elts['epoch'][0]

    # The correct orbital elements as an array of shape Nx6
    elts_true = elts.copy()

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
    num_obs: float = 5.69E6

    # Mask where data perturbed vs not
    mask_true = np.arange(elt_batch_size) < (elt_batch_size//2)
    mask_pert = ~mask_true

    # Perturb second half of orbital elements
    sigma_a = 0.05 
    sigma_e = 0.01 
    sigma_f_deg = 5.0
    elts_pert = perturb_elts(elts, sigma_a=sigma_a, sigma_e=sigma_e, sigma_f_deg=sigma_f_deg, mask_pert=mask_pert)

    # Report element shifts
    print(f'Perturbing elements:')
    print(f'sigma_a: {sigma_a:5.06f} (log)')
    print(f'sigma_e: {sigma_e:5.06f} (log)')
    print(f'sigma_f: {sigma_f_deg:5.06f} (degrees)\n')

    # Report initial scores
    report_initial_scores(elts_pert=elts_pert)

    # Orbits for calibration
    if use_calibration:
        print(f'Numerically integrating calibration trajectories q_cal...')
        q_cal = calc_ast_pos(elts=elts_pert, epoch=epoch, ts=ts)
    else:
        q_cal = None

    # The observation site for ZTF data is Palomar Mountain.
    site_name = 'palomar'

    # Alpha and beta parameters for the objective function
    # objective = score - alpha * mu - beta * sigma
    alpha = 1.0
    beta = 1.0

    # Build functional model for asteroid score
    model = make_model_asteroid_search(\
        ts=ts, elts=elts_pert, max_obs=max_obs, num_obs=num_obs, site_name=site_name,
        elt_batch_size=elt_batch_size, time_batch_size=time_batch_size,
        R_deg=R_deg, thresh_deg=thresh_deg, R_is_trainable=R_is_trainable, alpha=alpha, beta=beta, 
        q_cal=q_cal, use_calibration=use_calibration)

    # Use Adam optimizer with gradient clipping
    # learning_rate = 2.0e-5
    learning_rate = 2.0E-5  # default 1.0E-3
    beta_1 = 0.900          # default 0.900
    beta_2 = 0.999          # default 0.999
    epsilon = 1.0E-7        # default 1.0E-7
    amsgrad = False         # default False
    clipvalue = None        # default not used
    # Optimizer arguments
    opt_args = {
        'learning_rate': learning_rate,
        'beta_1': beta_1,
        'beta_2': beta_2,
        'epsilon': epsilon,
        'amsgrad': amsgrad,
    }
    # Add the clip value if it was set
    if clipvalue is not None:
        opt_args['clipvalue'] = clipvalue
    # Build the optimizer
    opt = keras.optimizers.Adam(**opt_args)
    # Compile the model
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
        report_model(model=model, ds=ds, R_deg=R_deg, mask_true=mask_true, 
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
    steps_per_epoch = steps*cycles_per_epoch
    rows_per_epoch = steps_per_epoch * time_batch_size
    print_header(f'Training for {epochs} Epochs of Size {rows_per_epoch} ' 
                 f'Observation Days (Cycles/Epoch = {cycles_per_epoch})...')
    print(f'R (degrees):      {R_deg:8.6f}')
    print(f'thresh (deg):     {thresh_deg:8.6f}')
    print(f'R_is_trainable:   {R_is_trainable}')
    print(f'alpha:          {alpha:5.1f}')
    print(f'beta:           {beta:5.1f}')
    print(f'learning_rate:    {learning_rate:7.2e}')
    if clipvalue is not None:
        print(f'clipvalue:      {clipvalue:5.2f}')
    print()

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
        report_model(model=model, ds=ds, R_deg=R_deg, mask_true=mask_true, 
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
    scores_01 = (scores0, scores1)
    traj_err_01 = (traj_err0, traj_err1)
    elt_err_01 = (elt_err0, elt_err1)
    R_01 = (R0, R1)
    report_training_progress(scores_01, traj_err_01, elt_err_01, R_01, mask_true)

    # Return training results
    return scores_01, traj_err_01, elt_err_01, R_01, mask_true

# ********************************************************************************************************************* 
if __name__ == '__main__':    
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Search for asteroids.')
    parser.add_argument('-R_deg', nargs='?', metavar='R_deg', type=float, default=0.2,
                        help='initial value of resolution factor, R, in degrees')
    parser.add_argument('-thresh_deg', nargs='?', metavar='thresh_deg', type=float, default=1.0,
                        help='threshold in degrees; observations further from candidate elements than this are ignored')
    parser.add_argument('-R_is_trainable', default=True, action='store_true',
                        help='whether resolution R is trainable')
    parser.add_argument('-elt_batch_size', nargs='?', metavar='elt_batch_size', type=int, default=64,
                        help='the number of candidate orbital elements per batch')
    parser.add_argument('-epochs', nargs='?', metavar='elt_batch_size', type=int, default=10,
                        help='the number of epochs to train')
    parser.add_argument('-cycles_per_epoch', nargs='?', metavar='elt_batch_size', type=int, default=50,
                        help='number of times to cycle through all data points before declaring an epoch complete')
    parser.add_argument('-gpu_num', nargs='?', metavar='gpu_num', type=int, default=1,
                        help='the GPU to use; defaults to 1 to avoid clash with Jupyter sessions')
    parser.add_argument('--test', default=True, action='store_true',
                        help='test on easy batch')
    args = parser.parse_args()

    # Alias arguments
    R_deg = args.R_deg
    thresh_deg = args.thresh_deg
    R_is_trainable = args.R_is_trainable
    elt_batch_size = args.elt_batch_size
    epochs = args.epochs
    cycles_per_epoch = args.cycles_per_epoch
    gpu_num = args.gpu_num

    # The selected gpu device
    gpu_device = get_gpu_device(gpu_num)

    # Set time_batch_size to None always to use all available snapshots (with threshold this is way better)
    time_batch_size = None

    # Set use_calibration; should be True unless doing quick code testing
    use_calibration = True

    # If run in test mode, test on the easy batch
    if args.test:
        with gpu_device:
            print(f'Running test_easy_batch with gpu {gpu_num}.')
            test_easy_batch(R_deg=R_deg, 
                            thresh_deg=thresh_deg, 
                            R_is_trainable=R_is_trainable,
                            elt_batch_size=elt_batch_size,
                            time_batch_size=time_batch_size,
                            epochs=epochs,
                            cycles_per_epoch=cycles_per_epoch,
                            use_calibration=use_calibration)
