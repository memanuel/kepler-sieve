"""
Harvard IACS Masters Thesis
asteroid_search_model.py: Tensorflow layers and models used to search for asteroid orbital elements.

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
from asteroid_model import AsteroidDirection, elts_np2df
from asteroid_integrate import calc_ast_pos
from search_score_functions import score_mean, score_var, score_mean_var
from candidate_element import perturb_elts
from astro_utils import deg2dist
from tf_utils import tf_quiet, Identity

# Typing
from typing import Optional

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Run TF quietly
tf_quiet()

# Constants
space_dims = 3

# ********************************************************************************************************************* 
# Transformations of orbital elements for search
# Range for a
a_min_: float = 0.5
a_max_: float = 32.0

# Range for e
e_min_: float = 0.0
e_max_: float = 1.0 - 2.0**-10

# Range for hit rate h
h_min_ = 2.0**-12
h_max_ = 1.0 - h_min_

# Range for exponential decay parameter lambda
lam_min_ = 2.0**-6
lam_max_ = 2.0**24
log_lam_min_ = np.log(lam_min_)
log_lam_max_ = np.log(lam_max_)

# Range for resolution parameter R
R_min_ = deg2dist(1.0/3600)
R_max_ = deg2dist(1.0)
log_R_min_ = np.log(R_min_)
log_R_max_ = np.log(R_max_)

# ********************************************************************************************************************* 
# Custom Layers
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
class OrbitalElements(keras.layers.Layer):
    """Custom layer to maintain state of candidate orbital elements and resolutions."""

    def __init__(self, 
                elts: pd.DataFrame, 
                batch_size: int, 
                h: float = 0.5,
                lam: float = 1.0,
                thresh_deg: float = 1.0,
                R_deg: float = 1.0, 
                R_is_trainable: bool = True, 
                **kwargs):
        super(OrbitalElements, self).__init__(**kwargs)
        
        # Configuration for serialization
        self.cfg = {
            'elts': elts,
            'batch_size': batch_size,
            'h': h,
            'lam': lam,
            'thresh_deg': thresh_deg,
            'R_deg': R_deg,
            'R_is_trainable': R_is_trainable,
        }
        # Save batch size, orbital elements as numpy array
        self.batch_size = batch_size
        self.elts = elts
        
        # Control over a_, in range 0.0 to 1.0
        self.a_min = tf.constant(a_min_, dtype=tf.float32)
        self.log_a_range = tf.constant(tf.math.log(a_max_) - tf.math.log(a_min_), dtype=tf.float32)
        self.a_ = tf.Variable(initial_value=self.inverse_a(elts['a']), trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='a_')
        
        # Control over e_, in range e_min to e_max
        self.e_min = tf.constant(e_min_, dtype=tf.float32)
        self.e_max = tf.constant(e_max_, dtype=tf.float32)
        # self.e_range = tf.constant(e_max_ - e_min_, dtype=tf.float32)
        # self.e_ = tf.Variable(initial_value=self.inverse_e(elts['e']), trainable=True, 
        #                      constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='e_')
        self.e_ = tf.Variable(initial_value=elts['e'], trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, self.e_min, self.e_max), name='e_')
        
        # Control over inc_, in range -pi/2 to pi/2
        self.inc_max = tf.constant(np.pi/2*(1-2**-20), dtype=tf.float32)
        self.inc_min = -self.inc_max
        self.inc_range = tf.constant(self.inc_max - self.inc_min, dtype=tf.float32)
        self.inc_ = tf.Variable(initial_value=self.inverse_inc(elts['inc']), trainable=True, 
                                constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='inc_')
        
        # Scale factor for unconstrained angles is 2*pi
        self.two_pi = tf.constant(2*np.pi, dtype=tf.float32)
        
        # Angle variables Omega, omega and f
        self.Omega_ = tf.Variable(initial_value=self.inverse_angle(elts['Omega']), trainable=True, name='Omega_')
        self.omega_ = tf.Variable(initial_value=self.inverse_angle(elts['omega']), trainable=True, name='omega_')
        self.f_ = tf.Variable(initial_value=self.inverse_angle(elts['f']), trainable=True, name='f_')

        # The epoch is not trainable
        self.epoch = tf.Variable(initial_value=elts['epoch'], trainable=False, name='epoch')
        
        # The threshold in degrees is not trainable
        self.trhesh_deg = tf.Variable(initial_value=thresh_deg, trainable=False, name='thresh_deg')

        # Control of the hit rate h_, in range h_min to h_max
        self.h_min = tf.constant(h_min_, dtype=tf.float32)
        self.h_max = tf.constant(h_max_, dtype=tf.float32)
        self.h_ = tf.Variable(initial_value=h, trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, self.h_min, self.h_max), name='h_')

        # Control of the exponential decay paramater lam_, in range lam_min to lam_max
        lam_init = lam * np.ones_like(elts['a'])
        self.lam_min = tf.constant(h_min_, dtype=tf.float32)
        self.log_lam_range = tf.constant(log_lam_max_ - log_lam_min_, dtype=tf.float32)
        self.lam_ = tf.Variable(initial_value=self.inverse_lam(lam_init), trainable=True, 
                                constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='h_')

        # Control of the resolution factor R_, in range 0.0 to 1.0
        R_init = np.deg2rad(R_deg) * np.ones_like(elts['a'])
        self.R_min = tf.constant(R_min_, dtype=tf.float32)
        self.log_R_range = tf.constant(log_R_max_ - log_R_min_, dtype=tf.float32)
        self.R_ = tf.Variable(initial_value=self.inverse_R(R_init), trainable=R_is_trainable, 
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

    def get_h(self):
        """Transformed value of h"""
        return self.h_

    def get_lam(self):
        """Transformed value of lambda"""
        return self.lam_min * tf.exp(self.lam_ * self.log_lam_range)

    def inverse_lam(self, lam):
        """Inverse transform value of lam"""
        return tf.math.log(lam / self.lam_min) / self.log_lam_range

    def get_R(self):
        """Transformed value of R"""
        return self.R_min * tf.exp(self.R_ * self.log_R_range)

    def inverse_R(self, R):
        """Inverse transform value of R"""
        return tf.math.log(R / self.R_min) / self.log_R_range

    def call(self, inputs):
        """Return the current orbital elements and resolution"""
        # print(f'type(inputs)={type(inputs)}.')
        # Transform a, e from log to linear
        a = self.get_a()
        e = self.get_e()
        # Angles transformed by factor of 2*pi
        inc = self.get_inc()
        Omega = self.get_Omega()
        omega = self.get_omega()
        f = self.get_f()
        # Transform search parameters h and lambda
        h = self.get_h()
        lam = self.get_lam()
        R = self.get_R()
        return a, e, inc, Omega, omega, f, self.epoch, h, lam, R

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
    def __init__(self, batch_size: int, thresh_deg: float, 
                 alpha: float, beta: float, **kwargs):
        """
        INPUTS:
            batch_size: this is element_batch_size, the number of orbital elements per batch
            thresh_deg: threshold in degrees for observations to be included in the batch.  impacts mu, sigma2
            alpha: multiplicative factor on mu in objective function
            beta: multiplicative factor on sigma2 in objective function
        """
        super(TrajectoryScore, self).__init__(**kwargs)

        # Configuration for seralization
        self.cfg = {
            'batch_size': batch_size,
            'alpha': alpha,
            'beta': beta,
            'thresh_deg': thresh_deg
        }

        # Save sizes and parameters
        self.batch_size = batch_size
        # Convert threshold from degrees to Cartesian distance
        self.thresh = deg2dist(thresh_deg)
        self.thresh2 = self.thresh**2
        # Tuning factors (get rid of these)
        self.alpha = alpha
        self.beta = beta
        
    def call(self, z: tf.Tensor, R: tf.Tensor, num_obs: float):
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            z: difference in direction between u_pred and u_obs
            R: resolution factor in Cartesian distance for score function
            num_obs: total number of observations (real, not padded!)
        """
        # The scaling coefficient for scores; score = exp(-1/2 A epsilon^2)
        A = 1.0 / R**2
        
        # The coefficient that multiplies epsilon^2
        B = tf.reshape(-0.5 * A, (self.batch_size, 1, 1,))
        # print(f'B.shape = {B.shape}')
        
        # Argument to the exponential
        z2 = K.sum(tf.square(z), axis=(-1))
        arg = tf.multiply(B, z2)
        # print(f'arg.shape = {arg.shape}')
        
        # Filter to only include terms where z2 is within the threshold distance^2
        is_close = z2 < self.thresh2
        close_idx = tf.where(is_close)
        raw_scores = tf.scatter_nd(indices=close_idx, updates=tf.exp(arg[is_close]), shape=arg.shape)
        
        # The score function
        raw_score = K.sum(raw_scores, axis=(1,2))
        # print(f'raw_score.shape = {raw_score.shape}')
        
        # The expected score and variance per observation
        mu_per_obs, sigma2_per_obs = score_mean_var(A, thresh=self.thresh)

        # The expected score
        mu = tf.multiply(num_obs, mu_per_obs, name='mu')
        
        # The expected variance
        # Note; variance adds up over batches, not std dev.  
        # However, if the batch size is all of the observations, then std dev can be meaninfully averaged
        sigma2 = tf.multiply(num_obs, sigma2_per_obs, name='sigma2')
        sigma = tf.sqrt(sigma2)

        # The t-statistic
        # t_stat = (raw_score - mu) / sigma

        # Assemble the objective function to be maximized
        # objective = raw_score - self.alpha * mu - self.beta * sigma2
        objective = (raw_score - self.alpha * mu - self.beta + sigma)
       
        # Return all the interesting outputs
        return raw_score, mu, sigma2, objective

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
# Functional API model
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_model_asteroid_search(ts: tf.Tensor,
                               elts: pd.DataFrame,
                               max_obs: int,
                               num_obs: float,
                               site_name: str='geocenter',
                               elt_batch_size: int=64, 
                               time_batch_size: Optional[int]=None,
                               h: float = 0.5,
                               lam: float = 1.0,
                               R_deg: float = 5.0,
                               thresh_deg: float = 1.0,
                               R_is_trainable: bool = True,
                               alpha: float = 2.0,
                               beta: float = 0.0,
                               q_cal = None,
                               use_calibration: bool = True):
    """Make functional API model for scoring elements"""

    # The full trajectory size
    traj_size: int = ts.shape[0]
    # If time_batch_size was not specifiued, use full trajectory size
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
    elements_layer = OrbitalElements(elts=elts, batch_size=elt_batch_size, 
                                     h=h, lam=lam,
                                     R_deg=R_deg, R_is_trainable=R_is_trainable, 
                                     name='candidates')
    a, e, inc, Omega, omega, f, epoch, h, lam, R = elements_layer(idx)
    
    # Alias the orbital elements; a, e, inc, Omega, omega, and f are trainable; epoch is fixed
    a = Identity(name='a')(a)
    e = Identity(name='e')(e)
    inc = Identity(name='inc')(inc)
    Omega = Identity(name='Omega')(Omega)
    omega = Identity(name='omega')(omega)
    f = Identity(name='f')(f)
    epoch = Identity(name='epoch')(epoch)

    # Alias the resolution output
    h = Identity(name='h')(h)
    lam = Identity(name='lam')(lam)
    R = Identity(name='R')(R)

    # The orbital elements; stack to shape (elt_batch_size, 7)
    elts_tf = tf.stack(values=[a, e, inc, Omega, omega, f, epoch], axis=1, name='elts')

    # The predicted direction
    direction_layer = AsteroidDirection(ts=ts, site_name=site_name, batch_size=elt_batch_size, name='u_pred')

    # Compute numerical orbits for calibration if necessary
    if use_calibration:
        if q_cal is not None:
            q_ast, q_earth, v_ast = q_cal
        else:
            print(f'Numerically integrating orbits for calibration...')
            epoch0 = elts['epoch'][0]
            q_ast, q_earth, v_ast = calc_ast_pos(elts=elts, epoch=epoch0, ts=ts)

    # Calibrate the direction prediction layer
    if use_calibration:
        direction_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)
    # Tensor of predicted directions
    u_pred, r_pred = direction_layer(a, e, inc, Omega, omega, f, epoch)

    # Difference in direction between u_obs and u_pred
    dir_diff_layer = DirectionDifference(batch_size=elt_batch_size, 
                                         traj_size=time_batch_size, 
                                         max_obs=max_obs, 
                                         name='z')
    z = dir_diff_layer(u_obs, u_pred, idx)

    # Calculate score compoments
    score_layer = TrajectoryScore(batch_size=elt_batch_size, alpha=alpha, beta=beta, thresh_deg=thresh_deg)
    raw_score,  mu, sigma2, objective  = score_layer(z, R, num_obs)
    # Stack the scores
    scores = tf.stack(values=[raw_score, mu, sigma2, objective], axis=1, name='scores')

    # Wrap inputs and outputs
    inputs = (t, idx, row_len, u_obs)
    outputs = (elts_tf, R, u_pred, z, scores)

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

    return model

# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
