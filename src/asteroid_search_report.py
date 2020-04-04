"""
Harvard IACS Masters Thesis
asteroid_search_report.py: Functions for reporting results during asteroid search.

Michael S. Emanuel
Thu Oct 17 15:24:10 2019
"""

# Core
import numpy as np
import pandas as pd
import tensorflow as tf

# Local imports
from candidate_element import elts_np2df, elts_df2dict
from asteroid_model import make_model_ast_pos


# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
def traj_diff(elts0: pd.DataFrame, elts1: pd.DataFrame, model_pos: keras.Model):
    """
    Compute difference between two orbital elements trajectories
    INPUTS:
        elts0: First set of orbital elements as DataFrame
        elts1: Second set of orbital elements as DataFrame
        Element DataFrames must have columns: a, e, inc, Omega, omega, f, epoch
    """
    # Convert DataFrame to dict 
    # TensorFlow interprets a DataFrame as a rectangular array, not like a dict
    elts0_dict = elts_df2dict(elts0)
    elts1_dict = elts_df2dict(elts1)

    # Predict position and velocity
    q0, v0 = model_pos.predict(elts0_dict)
    q1, v1 = model_pos.predict(elts1_dict)
    
    # Displacement between predicted trajsectories; shape (batch_size, traj_size, 3)
    dq = q1 - q0
    # Distance between predicted trajectories; shape (batch_size, traj_size)
    distance = tf.reduce_sum(dq**2, axis=-1)
    # Return norm of the difference by row
    return tf.reduce_mean(distance, axis=-1)
    
# ********************************************************************************************************************* 
def report_model_attribute(att: np.array, mask_true: np.ndarray, att_name: str):
    """Report mean and stdev of a model attribute on true and perturbed masks"""
    # Complementary mask
    mask_pert = ~mask_true

    # Attribute on masks
    att_g = att[mask_true]
    att_b = att[mask_pert]
    mean_g = np.mean(att_g, axis=0)
    mean_b = np.mean(att_b, axis=0)
    mean_a = np.mean(att, axis=0)
    std_g = np.std(att_g, axis=0)
    std_b = np.std(att_b, axis=0)
    std_a = np.std(att, axis=0)

    print(f'\nMean & Std {att_name} by Category:')
    print(f'True: {mean_g:8.2f} +/- {std_g:8.2f}')
    print(f'Pert  {mean_b:8.2f} +/- {std_b:8.2f}')
    print(f'All:  {mean_a:8.2f} +/- {std_a:8.2f}')
    
# ********************************************************************************************************************* 
def report_model_attribute_change(att0: np.array, att1: np.array, mask_true: np.array, 
                                  att_name: str, dp: int = 6, dp0: int = None):
    """Report mean and stdev of a model attribute on true and perturbed masks"""
    # Default arguments
    if dp0 is None:
        dp0 = dp

    # Complementary mask
    mask_pert = ~mask_true
    
    # Change
    chng = att1 - att0
    
    # Attribute on masks
    att_g0 = att0[mask_true]
    att_b0 = att0[mask_pert]
    att_g1 = att1[mask_true]
    att_b1 = att1[mask_pert]
    chng_g = chng[mask_true]
    chng_b = chng[mask_pert]

    # Mean of attribute
    mean_g0 = np.mean(att_g0, axis=0)
    mean_b0 = np.mean(att_b0, axis=0)
    mean_a0 = np.mean(att0, axis=0)
    mean_g1 = np.mean(att_g1, axis=0)
    mean_b1 = np.mean(att_b1, axis=0)
    mean_a1 = np.mean(att1, axis=0)

    # Mean change
    mean_chng_g = np.mean(chng_g, axis=0)
    mean_chng_b = np.mean(chng_b, axis=0)
    mean_chng_a = np.mean(chng, axis=0)

    # Std change
    std_chng_g = np.std(chng_g, axis=0)
    std_chng_b = np.std(chng_b, axis=0)
    std_chng_a = np.std(chng, axis=0)
    
    # Percentage change
    pct_chng_g = mean_chng_g / mean_g0
    pct_chng_b = mean_chng_b / mean_b0
    pct_chng_a = mean_chng_a / mean_a0

    print(f'\nMean & Std Change in {att_name} by Category:')
    print(f'True: {mean_chng_g:+8.{dp}f} +/- {std_chng_g:8.{dp}f} -- from {mean_g0:8.{dp0}f} to {mean_g1:8.{dp0}f}'
          f'     ({pct_chng_g:+3.1%})')
    print(f'Pert: {mean_chng_b:+8.{dp}f} +/- {std_chng_b:8.{dp}f} -- from {mean_b0:8.{dp0}f} to {mean_b1:8.{dp0}f}'
          f'     ({pct_chng_b:+4.1%})')
    print(f'All:  {mean_chng_a:+8.{dp}f} +/- {std_chng_a:8.{dp}f} -- from {mean_a0:8.{dp0}f} to {mean_a1:8.{dp0}f}'
          f'     ({pct_chng_a:+4.1%})')
   
# ********************************************************************************************************************* 
def report_model(model: keras.Model, ds: tf.data.Dataset, R_deg: float, 
                 mask_true: np.ndarray, batch_size: int, steps: int, 
                 elts_true: pd.DataFrame, display: bool =True):
    """
    Report summary of model on true and perturbed
    INPUTS:
        model:          Keras model to be evaluated
        ds:             Tensorflow DataSet
        R_deg:          Resolution in degrees
        mask_true:      Mask for true (unperturbed) elements
        batch_size:     Number of orbital elements in a batch, e.g. 64
        steps:          Number of batches to work through tf.Dataset once
        elts_true:      True elements as DataFrame
        display:        Whether to display additional details
    """

    # Mask for perturbed ("bad") elements
    mask_pert = ~mask_true
    
    # Get scores on the whole data set
    pred = model.predict(ds, steps=steps)
    elts, R, u_pred, z, scores_all = pred

    # Consolidate results to batch_size rows
    num_rows = elts.shape[0]
    score_cols = scores_all.shape[1]
    row_idx = np.arange(num_rows, dtype=np.int32) % batch_size
    elts = elts[0:batch_size]
    R = R[0:batch_size]
    u_pred = u_pred[0:batch_size]

    # Output elements from TF model is a numpy array; convert to a DataFrame
    elts_out = elts_np2df(elts)

    # Model to predict position from elements
    ts = model.direction.ts
    model_pos = make_model_ast_pos(ts=ts, batch_size=batch_size)

    # Difference between predicted and true trajectory in Kepler 2 body model
    traj_err = traj_diff(elts0=elts_out, elts1=elts_true, model_pos=model_pos)

    # Consolidate the scores; create 2 extra columns for sigma and t_score
    scores = np.zeros((batch_size, score_cols+2))
    for batch_idx in range(batch_size):
        mask = (row_idx == batch_idx)
        scores[batch_idx, 0:score_cols] = scores_all[mask].sum(axis=0)

    # Unpock scores
    raw_score = scores[:,0]
    mu = scores[:,1]
    sigma2 = scores[:,2]
    objective = scores[:,3]

    # Compute derived scores after aggregation
    sigma = np.sqrt(sigma2)
    eff_obs = raw_score - mu
    t_score = eff_obs / sigma
    
    # Pack sigma and t_score at the end of scores
    scores[:, 4] = sigma
    scores[:, 5] = t_score
    
    # Error in orbital elements; express as an Nx6 numpy array, extracted from DataFrames
    cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f']
    elt_err = np.abs(elts_out[cols_elt].values - elts_true[cols_elt].values)
    # Convert angles from radians to degrees
    elt_err[:, 2:6] = np.rad2deg(elt_err[:, 2:6])
    # Mean element error on true and perturbed masks
    elt_err_g = elt_err[mask_true]
    elt_err_b = elt_err[mask_pert]
    mean_err_g = np.mean(elt_err_g[0:6], axis=0)
    mean_err_b = np.mean(elt_err_b[0:6], axis=0)

    if display:
        # Report trajectory error
        report_model_attribute(traj_err, mask_true, 'Trajectory Error (AU) vs. True Elements')

        # Report errors in orbital elements
        print('\nError in orbital elements:')
        print(f'(Angles shown in degrees)')
        print('      a          e         inc       Omega      omega     f')
        print(f'True: {mean_err_g[0]:8.6f},  {mean_err_g[1]:8.6f}, {mean_err_g[2]:8.6f}, '
              f'{mean_err_g[3]:8.6f},  {mean_err_g[4]:8.6f}, {mean_err_g[5]:8.6f}, ')
        print(f'Pert: {mean_err_b[0]:8.6f},  {mean_err_b[1]:8.6f}, {mean_err_b[2]:8.6f}, '
              f'{mean_err_b[3]:8.6f},  {mean_err_b[4]:8.6f}, {mean_err_b[5]:8.6f}, ')
    
        # Report effective observations, mu, sigma, and t_score    
        report_model_attribute(raw_score, mask_true, 'Raw Score')
        report_model_attribute(mu, mask_true, 'Mu')
        report_model_attribute(eff_obs, mask_true, 'Effective Observations')
        report_model_attribute(sigma, mask_true, 'Sigma')
        report_model_attribute(t_score, mask_true, 't_score')
        report_model_attribute(objective, mask_true, 'Objective Function')
        report_model_attribute(R, mask_true, 'Resolution R')

    return scores, traj_err, elt_err

# ********************************************************************************************************************* 
def report_training_progress(scores_01, traj_err_01, elt_err_01, R_01, mask_true, report_elt_change: bool = False):
    """Report progress while model trained"""
    
    # Unpack inputs pairs
    scores0, scores1 = scores_01
    traj_err0, traj_err1 = traj_err_01
    elt_err0, elt_err1 = elt_err_01
    R0, R1, = R_01

    # Unpack score components
    raw_score0 = scores0[:,0]
    raw_score1 = scores1[:,0]
    mu0 = scores0[:,1]
    mu1 = scores1[:,1]
    sigma2_0 = scores0[:,2]
    sigma2_1 = scores1[:,2]
    objective0 = scores0[:,3]
    objective1 = scores1[:,3]
    sigma0 = scores0[:,4]
    sigma1 = scores1[:,4]
    t_score0 = scores0[:,5]
    t_score1 = scores1[:,5]
    
    # Calculations
    eff_obs0 = raw_score0 - mu0
    eff_obs1 = raw_score1 - mu1

    # Changes in trajectory errors
    report_model_attribute_change(traj_err0, traj_err1, mask_true, 'Trajectory Error (AU)', dp=6)
    report_model_attribute_change(raw_score0, raw_score1, mask_true, 'Raw Score', dp=0)
    report_model_attribute_change(mu0, mu1, mask_true, 'Mu', dp=0)
    report_model_attribute_change(eff_obs0, eff_obs1, mask_true, 'Effective Observations', dp=0)
    report_model_attribute_change(sigma0, sigma1, mask_true, 'Sigma', dp=0)
    report_model_attribute_change(t_score0, t_score1, mask_true, 't-score', dp=1)
    report_model_attribute_change(objective0, objective1, mask_true, 'Objective Function', dp=0)
    report_model_attribute_change(np.rad2deg(R0), np.rad2deg(R1), mask_true, 'Resolution R (degrees)', dp=4)
    
    # Unpack element errors
    err_a0 = elt_err0[:,0]
    err_a1 = elt_err1[:,0]
    err_e0 = elt_err0[:,1]
    err_e1 = elt_err1[:,1]
    err_inc0 = elt_err0[:,2]
    err_inc1 = elt_err1[:,2]
    # err_Omega0 = elt_err0[:,3]
    # err_Omega1 = elt_err1[:,3]
    # err_omega0 = elt_err0[:,4]
    # err_omega1 = elt_err1[:,4]
    err_f0 = elt_err0[:,5]
    err_f1 = elt_err1[:,5]
    
    # Report element errors if requested
    if report_elt_change:
        report_model_attribute_change(err_a0, err_a1, mask_true, 'Error in semi-major axis, a', dp=6)
        report_model_attribute_change(err_e0, err_e1, mask_true, 'Error in eccentricity, e', dp=6)
        report_model_attribute_change(err_inc0, err_inc1, mask_true, 'Error in inclination, inc', dp=6)
        report_model_attribute_change(err_f0, err_f1, mask_true, 'Error in true anomaly, f', dp=6)
    
    # Changes in element errors and R
    d_elt_err = elt_err1 - elt_err0
    d_R = R1 - R0
    mask_pert = ~mask_true
    d_elt_err_g = np.mean(d_elt_err[mask_true], axis=0)
    d_elt_err_b = np.mean(d_elt_err[mask_pert], axis=0)
    d_R_g = np.mean(d_R[mask_true])
    d_R_b = np.mean(d_R[mask_pert])

    print(f'\nChange in Orbital Element error by Category:')
    print(f'(Angles shown in degrees)')
    print('          a           e         inc         Omega       omega      f          R')
    print(f'd_err_g: {d_elt_err_g[0]:+8.6f},  {d_elt_err_g[1]:+8.6f}, {d_elt_err_g[2]:+8.6f}, '
          f'{d_elt_err_g[3]:+8.6f},  {d_elt_err_g[4]:+8.6f}, {d_elt_err_g[5]:+8.6f}, {d_R_g:+8.6f}')
    print(f'd_err_b: {d_elt_err_b[0]:+8.6f},  {d_elt_err_b[1]:+8.6f}, {d_elt_err_b[2]:+8.6f}, '
          f'{d_elt_err_b[3]:+8.6f},  {d_elt_err_b[4]:+8.6f}, {d_elt_err_b[5]:+8.6f}, {d_R_b:+8.6f}')
