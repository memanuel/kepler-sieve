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
from ztf_data import make_ztf_easy_batch
# from ztf_data import load_ztf_det_all, load_ztf_nearest_ast, ztf_nearest_ast, calc_hit_freq
from asteroid_data import make_ztf_dataset, orbital_element_batch
from asteroid_integrate import load_data as load_data_asteroids, calc_ast_pos
from asteroid_model import AsteroidDirection
from search_score_functions import score_mean, score_var, score_mean_2d, score_var_2d
from asteroid_search_report import report_model, report_training_progress
from astro_utils import deg2dist
from utils import print_header
from tf_utils import tf_quiet, gpu_grow_memory, Identity

# Typing
from typing import Dict

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Run TF quietly
tf_quiet()
# Configure TensorFlow to use GPU memory variably; done automatically by importing asteroid_model
# gpu_grow_memory(verbose=True)

# Constants
space_dims = 3

# Load asteroid names and orbital elements
ast_elt = load_data_asteroids()

# ********************************************************************************************************************* 
# Transformations of orbital elements for search
# Range for a
a_min_: float = 0.5
a_max_: float = 32.0

# Range for e
e_min_: float = 0.0
e_max_: float = 1.0 - 2.0**-10

# Range for resolution parameter R
R_min_ = deg2dist(1.0/3600)
R_max_ = deg2dist(10.0)
log_R_min_ = np.log(R_min_)
log_R_max_ = np.log(R_max_)

# ********************************************************************************************************************* 
# Custom Layers
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
class OrbitalElements(keras.layers.Layer):
    """Custom layer to maintain state of candidate orbital elements and resolutions."""

    def __init__(self, elts_np: dict, batch_size: int, R_deg: float, **kwargs):
        super(OrbitalElements, self).__init__(**kwargs)
        
        # Configuration for serialization
        self.cfg = {
            'elts_np': elts_np,
            'batch_size': batch_size,
            'R_deg': R_deg,
        }
        # Save batch size, orbital elements as numpy array, and resolution in degrees
        self.batch_size = batch_size
        self.elts_np = elts_np
        self.R_deg = R_deg
        
        # Control over a_, in range 0.0 to 1.0
        self.a_min = tf.constant(a_min_, dtype=tf.float32)
        self.log_a_range = tf.constant(tf.math.log(a_max_) - tf.math.log(a_min_), dtype=tf.float32)
        self.a_ = tf.Variable(initial_value=self.inverse_a(elts_np['a']), trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='a_')
        
        # Control over e_, in range e_min to e_max
        self.e_min = tf.constant(e_min_, dtype=tf.float32)
        self.e_max = tf.constant(e_max_, dtype=tf.float32)
        # self.e_range = tf.constant(e_max_ - e_min_, dtype=tf.float32)
        # self.e_ = tf.Variable(initial_value=self.inverse_e(elts_np['e']), trainable=True, 
        #                      constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='e_')
        self.e_ = tf.Variable(initial_value=elts_np['e'], trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, self.e_min, self.e_max), name='e_')
        
        # Control over inc_, in range -pi/2 to pi/2
        self.inc_max = tf.constant(np.pi/2*(1-2**-20), dtype=tf.float32)
        self.inc_min = -self.inc_max
        self.inc_range = tf.constant(self.inc_max - self.inc_min, dtype=tf.float32)
        self.inc_ = tf.Variable(initial_value=self.inverse_inc(elts_np['inc']), trainable=True, 
                                constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='inc_')
        
        # Scale factor for unconstrained angles is 2*pi
        self.two_pi = tf.constant(2*np.pi, dtype=tf.float32)
        
        self.Omega_ = tf.Variable(initial_value=self.inverse_angle(elts_np['Omega']), trainable=True, name='Omega_')
        self.omega_ = tf.Variable(initial_value=self.inverse_angle(elts_np['omega']), trainable=True, name='omega_')
        self.f_ = tf.Variable(initial_value=self.inverse_angle(elts_np['f']), trainable=True, name='f_')

        # The epoch is not trainable
        self.epoch = tf.Variable(initial_value=elts_np['epoch'], trainable=False, name='epoch')
        
        # Control of the resolution factor R_, in range 0.0 to 1.0
        R_init = np.deg2rad(R_deg) * np.ones_like(elts_np['a'])
        # log_R_init  = np.log(R_init)
        self.R_min = tf.constant(R_min_, dtype=tf.float32)
        self.log_R_range = tf.constant(log_R_max_ - log_R_min_, dtype=tf.float32)
        self.R_ = tf.Variable(initial_value=self.inverse_R(R_init), trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='R_')
        
    def get_a(self):
        """Transformed value of a"""
        return self.a_min * tf.exp(self.a_ * self.log_a_range)

    def inverse_a(self, a):
        """Inverse transform value of a"""
        return tf.math.log(a / self.a_min) / self.log_a_range

    def get_e(self):
        """Transformed value of e"""
        # return self.e_min + self.e_ * self.e_range
        return self.e_

    # def inverse_e(self, e):
    #    """Inverse transform value of e"""
    #    return (e - self.e_min) / self.e_range

    def get_inc(self):
        """Transformed value of inc"""
        return self.inc_min + self.inc_ * self.inc_range

    def inverse_inc(self, inc):
        """Inverse transform value of inc"""
        return (inc - self.inc_min) / self.inc_range

    def get_angle(self, angle_):
        """Forward transform of an unconstrained angle variable (Omega, omega, f)"""
        return self.two_pi * angle_

    def inverse_angle(self, angle):
        """Forward transform of an unconstrained angle variable (Omega, omega, f)"""
        return angle / self.two_pi

    def get_Omega(self):
        return self.get_angle(self.Omega_)

    def get_omega(self):
        return self.get_angle(self.omega_)

    def get_f(self):
        return self.get_angle(self.f_)

    def get_R(self):
        """Transformed value of R"""
        return self.R_min * tf.exp(self.R_ * self.log_R_range)

    def inverse_R(self, R):
        """Inverse transform value of R"""
        return tf.math.log(R / self.R_min) / self.log_R_range

    def call(self, inputs):
        """Return the current orbital elements and resolution"""
        # print(f'type(inputs)={type(inputs)}.')
        # Transform a, e, and R from log to linear
        a = self.get_a()
        e = self.get_e()
        inc = self.get_inc()
        Omega = self.get_Omega()
        omega = self.get_omega()
        f = self.get_f()
        R = self.get_R()
        return a, e, inc, Omega, omega, f, self.epoch, R

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class DirectionDifference(keras.layers.Layer):
    """Compute the difference in direction between observed and predicted directions"""
    def __init__(self, batch_size: int, traj_size: int, max_obs: int, **kwargs):
        super(DirectionDifference, self).__init__(**kwargs)

        # Configuration for serialization
        self.cfg = {
            'batch_size': batch_size,
            'traj_size': traj_size,
            'max_obs': max_obs,
        }

        # Save sizes
        self.batch_size = batch_size
        self.traj_size = traj_size
        self.max_obs = max_obs
    
    def call(self, u_obs, u_pred, idx):
        """
        INPUTS:
            u_obs: observed directions, PADDED to a regular tensor; shape (traj_size, max_obs, 3,)
            u_pred: predicted directions; shape (batch_size, traj_size, 3,)
        """
        # Get sizes
        batch_size = self.batch_size
        traj_size = self.traj_size
        max_obs = self.max_obs

        # Slice of the full trajectory
        i0 = idx[0]
        i1 = idx[-1] + 1
        u_pred_slice = u_pred[:,i0:i1]
        # Manually set the shapes to work around documented bug on slices losing shape info
        u_slice_shape = (batch_size, traj_size, 3)
        u_pred_slice.set_shape(u_slice_shape)

        # Debug
        # print(f'u_obs.shape = {u_obs.shape}')
        # print(f'u_pred.shape = {u_pred.shape}')
        # print(f'u_pred_slice.shape = {u_pred_slice.shape}')
        # print(f'batch_size={batch_size}, traj_size={traj_size}, max_obs={max_obs}.')
        # print(f'i0={i0}, i1={i1}.')

        # The observations; broadcast to shape (1, traj_size, max_obs, 3)
        y = tf.broadcast_to(u_obs, (1, traj_size, max_obs, space_dims))
        # print(f'y.shape = {y.shape}')
        # The predicted directions; reshape to (batch_size, traj_size, 1, 3)
        x = tf.reshape(u_pred_slice, (batch_size, traj_size, 1, space_dims))
        # print(f'x.shape = {x.shape}')
        
        # The difference in directions; size (batch_size, traj_size, max_obs, 3)
        z = tf.subtract(y, x, name='z')
        # print(f'z.shape = {z.shape}')

        return z
    
    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class TrajectoryScore(keras.layers.Layer):
    """Score candidate trajectories"""
    def __init__(self, batch_size: int, alpha: float, beta: float, **kwargs):
        """
        INPUTS:
            batch_size: this is element_batch_size, the number of orbital elements per batch
            alpha: multiplicative factor on mu in objective function
            beta: multiplicative factor on sigma2 in objective function
        """
        super(TrajectoryScore, self).__init__(**kwargs)

        # Configuration for seralization
        self.cfg = {
            'batch_size': batch_size,
            'alpha': alpha,
            'beta': beta,
        }

        # Save sizes and parameters
        self.batch_size = batch_size
        self.alpha = alpha
        self.beta = beta
        
        # Max values of mu and sigma2
        A_min = 1.0 / R_max_**2
        self.mu_max = score_mean(A_min)
        self.sigma2_max = score_var(A_min)

    def call(self, z: tf.Tensor, R: tf.Tensor, num_obs: float):
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            z: difference in direction between u_pred and u_obs
            R: resolution factor in radians for score function
            num_obs: total number of observations (real, not padded!)
        """
        # The scaling coefficient for scores; score = exp(-1/2 A epsilon^2)
        A = 1.0 / R**2
        
        # The coefficient that multiplies epsilon^2
        B = tf.reshape(-0.5 * A, (self.batch_size, 1, 1,))
        # print(f'B.shape = {B.shape}')
        
        # Argument to the exponential
        # arg = tf.multiply(B, tf.linalg.norm(z, axis=-1))
        # z2 = tf.square(tf.linalg.norm(z, axis=-1))
        z2 = K.sum(tf.square(z), axis=(-1))
        arg = tf.multiply(B, z2)
        # print(f'arg.shape = {arg.shape}')
        
        # Filter to only exponentiate where the argument is close enough to matter
        is_close = arg > -10
        close_idx = tf.where(is_close)
        raw_scores = tf.scatter_nd(indices=close_idx, updates=tf.exp(arg[is_close]), shape=arg.shape)
        
        # The score function
        # raw_score = K.sum(tf.exp(arg), axis=(1,2))
        raw_score = K.sum(raw_scores, axis=(1,2))
        # print(f'raw_score.shape = {raw_score.shape}')
        
        # The expected score
        mu_per_obs = score_mean(A)
        mu = tf.multiply(num_obs, mu_per_obs, name='mu')
        
        # The expected variance
        var_per_obs = score_var(A)
        sigma2 = tf.multiply(num_obs, var_per_obs, name='sigma2')
        # sigma = tf.sqrt(sigma2, name='sigma')

        # Assemble the objective function to be maximized
        objective = raw_score - self.alpha * mu - self.beta * sigma2
        
        # Minimum of objective
        # objective_min = num_obs * tf.dtypes.cast(-self.alpha * self.mu_max - self.beta * self.sigma2_max, tf.float32)
        
        # Return both the raw and t scores
        return raw_score, mu, sigma2, objective

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
# Functional API model
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_model_asteroid_search(ts: tf.Tensor,
                               elts_np: Dict,
                               max_obs: int,
                               num_obs: float,
                               site_name: str='geocenter',
                               elt_batch_size: int=64, 
                               time_batch_size: int=None,
                               R_deg: float = 5.0,
                               alpha: float = 2.0,
                               beta: float = 1.0,
                               q_cal = None,
                               use_calibration: bool = True):
    """Make functional API model for scoring elements"""

    # The full trajectory size
    traj_size: int = ts.shape[0]
    # Default for time_batch_size is full trajectory size
    if time_batch_size is None:
        time_batch_size = traj_size

    # Inputs
    t = keras.Input(shape=(), batch_size=time_batch_size, dtype=tf.float32, name='t' )
    idx = keras.Input(shape=(), batch_size=time_batch_size, dtype=tf.int32, name='idx')
    row_len = keras.Input(shape=(), batch_size=time_batch_size, dtype=tf.int32, name='row_len')
    u_obs = keras.Input(shape=(max_obs, space_dims), batch_size=time_batch_size, dtype=tf.float32, name='u_obs')
    
    # Output times are a constant
    ts = keras.backend.constant(ts, name='ts')

    # Set of trainable weights with candidate
    elements_layer = OrbitalElements(elts_np=elts_np, batch_size=elt_batch_size, R_deg=R_deg, name='candidates')
    a, e, inc, Omega, omega, f, epoch, R = elements_layer(idx)
    
    # Alias the orbital elements; a, e, inc, Omega, omega, and f are trainable; epoch is fixed
    a = Identity(name='a')(a)
    e = Identity(name='e')(e)
    inc = Identity(name='inc')(inc)
    Omega = Identity(name='Omega')(Omega)
    omega = Identity(name='omega')(omega)
    f = Identity(name='f')(f)
    epoch = Identity(name='epoch')(epoch)

    # Alias the resolution output
    R = Identity(name='R')(R)

    # The orbital elements; stack to shape (elt_batch_size, 7)
    elts = tf.stack(values=[a, e, inc, Omega, omega, f, epoch], axis=1, name='elts')

    # The predicted direction
    direction_layer = AsteroidDirection(ts=ts, site_name=site_name, batch_size=elt_batch_size, name='u_pred')

    # Compute numerical orbits for calibration if necessary
    if use_calibration:
        if q_cal is not None:
            q_ast, q_earth, v_ast = q_cal
        else:
            print(f'Numerically integrating orbits for calibration...')
            epoch0 = elts_np['epoch'][0]
            q_ast, q_earth, v_ast = calc_ast_pos(elts=elts_np, epoch=epoch0, ts=ts)

    # Calibrate the direction prediction layer
    if use_calibration:
        direction_layer.calibrate(elts=elts_np, q_ast=q_ast, v_ast=v_ast)
    # Tensor of predicted directions
    u_pred, r_pred = direction_layer(a, e, inc, Omega, omega, f, epoch)

    # Difference in direction between u_obs and u_pred
    dir_diff_layer = DirectionDifference(batch_size=elt_batch_size, 
                                         traj_size=time_batch_size, 
                                         max_obs=max_obs, 
                                         name='z')
    z = dir_diff_layer(u_obs, u_pred, idx)

    # Calculate score compoments
    score_layer = TrajectoryScore(batch_size=elt_batch_size, alpha=alpha, beta=beta)
    raw_score,  mu, sigma2, objective  = score_layer(z, R, num_obs)
    # Stack the scores
    scores = tf.stack(values=[raw_score, mu, sigma2, objective], axis=1, name='scores')

    # Wrap inputs and outputs
    inputs = (t, idx, row_len, u_obs)
    outputs = (elts, R, u_pred, z, scores)

    # Create model with functional API
    model = keras.Model(inputs=inputs, outputs=outputs)

    # Bind the custom layers to model
    model.elements = elements_layer
    model.direction = direction_layer
    model.dir_diff = dir_diff_layer
    model.score = score_layer
    
    # Log transform
    # objective_min = 0.0

    # Add the loss function
    model.add_loss(-tf.reduce_sum(objective))
    # model.add_loss(-tf.reduce_sum(tf.math.log(objective-objective_min)))

    return model

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
def main(time_batch_size: int = 128, elt_batch_size: int = 64, epochs=20):
    # Load all ZTF data with nearest asteroid calculations
    ztf = make_ztf_easy_batch(batch_size=elt_batch_size)

    # Build a TensorFlow DataSet from ZTF DataFrame
    ds, ts, row_len = make_ztf_dataset(ztf=ztf, batch_size=time_batch_size)

    # Trajectory size and steps per batch
    traj_size = ts.shape[0]
    steps = int(np.ceil(traj_size / time_batch_size))

    # Batch of perturbed orbital elements for asteroid model
    R_deg: float = 2.0
    ast_nums = np.unique(ztf.nearest_ast_num)
    elts_np = orbital_element_batch(ast_nums)
    epoch = elts_np['epoch'][0]

    # Get example batch
    batch_in, batch_out = list(ds.take(1))[0]
    # Contents of this batch
    t = batch_in['t']
    idx = batch_in['idx']
    row_len = batch_in['row_len']
    u_obs = batch_in['u_obs']

    # Get max_obs and number of observations
    max_obs: int = u_obs.shape[1]
    num_obs: float = np.sum(row_len, dtype=np.float32)

    # The correct orbital elements as an array
    elts_true = np.array([elts_np['a'], elts_np['e'], elts_np['inc'], elts_np['Omega'], 
                          elts_np['omega'], elts_np['f'], elts_np['epoch']]).transpose()

    # Mask where data expected vs not
    mask_good = np.arange(elt_batch_size) < (elt_batch_size//2)
    mask_bad = ~mask_good
    # Perturb second half of orbital elements
    elts_np2 = perturb_elts(elts_np, sigma_a=0.00, sigma_e=0.00, sigma_f_deg=0.0, mask=mask_bad)

    # Orbits for calibration
    if 'q_cal' not in globals():
        print(f'Numerically integrating calibration trajectories q_cal...')
        q_cal = calc_ast_pos(elts=elts_np2, epoch=epoch, ts=ts)
    # q_cal = None

    # Set calibration flag
    use_calibration: bool = True

    # Alpha and beta parameters for the objective function
    alpha = 8.0
    beta = 20.0

    # Build functional model for asteroid score
    model = make_model_asteroid_search(\
        ts=ts, elts_np=elts_np2, max_obs=max_obs, num_obs=num_obs,
        elt_batch_size=elt_batch_size, time_batch_size=time_batch_size,
        R_deg=R_deg, alpha=alpha, beta=beta, q_cal=q_cal, use_calibration=use_calibration)

    # Use Adam optimizer with gradient clipping
    learning_rate = 2.0e-5
    clipvalue = 5.0
    opt = keras.optimizers.Adam(learning_rate=learning_rate, 
                                beta_1=0.900, 
                                beta_2=0.999,
                                epsilon=1.0E-7, 
                                clipvalue=clipvalue, 
                                amsgrad=False)
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
    main(time_batch_size=128, elt_batch_size=64, epochs=5)

