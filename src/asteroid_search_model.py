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
from astro_utils import deg2dist, dist2deg
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

# Data type
dtype = tf.float32
dtype_np = np.float32

# ********************************************************************************************************************* 
# Transformations of orbital elements for search
# Range for a
a_min_: float = 0.5
a_max_: float = 32.0

# Range for e
e_min_: float = 0.0
e_max_: float = 1.0 - 2.0**-10

# Range for hit rate h
h_min_ = 2.0**-8
h_max_ = 1.0 - h_min_

# Range for exponential decay parameter lambda
lam_min_ = 2.0**-2
lam_max_ = 2.0**24
log_lam_min_ = np.log(lam_min_)
log_lam_max_ = np.log(lam_max_)

# # Range for resolution parameter R
# R_min_ = deg2dist(1.0/3600)
# R_max_ = deg2dist(1.0)
# log_R_min_ = np.log(R_min_)
# log_R_max_ = np.log(R_max_)

# ********************************************************************************************************************* 
# Custom Layers
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
class OrbitalElements(keras.layers.Layer):
    """Custom layer to maintain state of candidate orbital elements and resolutions."""

    def __init__(self, 
                elts: pd.DataFrame, 
                elt_batch_size: int, 
                h: float = 0.125,
                lam: float = 1.0,
                **kwargs):
        super(OrbitalElements, self).__init__(**kwargs)
        
        # Configuration for serialization
        self.cfg = {
            'elts': elts,
            'elt_batch_size': elt_batch_size,
            'h': h,
            'lam': lam,
        }

        # Alias elt_batch_size to batch_size for legibility
        batch_size = elt_batch_size

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
        
        # Control of the hit rate h_, in range h_min to h_max
        h_init = h * np.ones_like(elts['a'])
        self.h_min = tf.constant(h_min_, dtype=tf.float32)
        self.h_max = tf.constant(h_max_, dtype=tf.float32)
        self.h_ = tf.Variable(initial_value=h_init, trainable=True, dtype=dtype,
                              constraint=lambda t: tf.clip_by_value(t, self.h_min, self.h_max), name='h_')

        # Control of the exponential decay paramater lam_, in range lam_min to lam_max
        lam_init = lam * np.ones_like(elts['a'])
        self.lam_min = tf.constant(h_min_, dtype=tf.float32)
        self.log_lam_range = tf.constant(log_lam_max_ - log_lam_min_, dtype=tf.float32)
        self.lam_ = tf.Variable(initial_value=self.inverse_lam(lam_init), trainable=True, dtype=dtype,
                                constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='lam_')

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

    def call(self, inputs=None):
        """Return the current orbital elements and mixture model parameters"""
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
        return a, e, inc, Omega, omega, f, self.epoch, h, lam

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class TrajectoryScore(keras.layers.Layer):
    """Score candidate trajectories"""
    def __init__(self, row_lengths: tf.Tensor, thresh_deg: float, **kwargs):
        """
        INPUTS:
            row_lengths: Number of observations for each element; shape [elt_batch_size]
            thresh_deg: threshold in degrees for observations to be included
        """
        super(TrajectoryScore, self).__init__(**kwargs)

        # Configuration for seralization
        self.cfg = {
            'row_lengths': row_lengths,
            'thresh_deg': thresh_deg
        }

        # Sizes of the data as keras constants
        self.elt_batch_size = keras.backend.constant(value=row_lengths.shape[0], dtype=tf.int32)
        self.data_size = keras.backend.constant(value=tf.reduce_sum(row_lengths), dtype=tf.int32)
        self.row_lengths = keras.backend.constant(value=row_lengths, shape=row_lengths.shape, dtype=tf.int32)
        u_shape = (self.data_size, space_dims,)        

        # Save the observed directions as a keras constant
        # self.u_obs = keras.backend.constant(value=u_obs, shape=u_shape, dtype=dtype)

        # The threshold distance and its square
        self.thresh_s = keras.backend.constant(value=deg2dist(thresh_deg), dtype=dtype, name='thresh_s')
        self.thresh_s2 = keras.backend.constant(value=self.thresh_s**2, dtype=dtype, name='thresh_s2')
        
    def call(self, u_pred: tf.Tensor, u_obs: tf.Tensor, h: tf.Tensor, lam: tf.Tensor):
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            u_pred: predicted directions with current batch of elements
            u_obs:  observed directions
            h:      fraction of hits in mixture model
            lam:    exponential decay parameter in mixture model
        """
        # Difference between actual and predicted directions
        du = keras.layers.subtract(inputs=[u_pred, u_obs], name='du')
        # Squared distance bewteen predicted and observed directions
        s2 = tf.reduce_sum(tf.square(du), axis=(-1), name='s2')
        
        # Filter to only include terms where z2 is within the threshold distance^2
        # is_close = s2 < self.thresh_s2
        is_close = tf.math.less(s2, self.thresh_s2, name='is_close')
        
        # Relative distance v on data inside threshold
        # v = tf.divide(s2[is_close], self.thresh_s2, name='v')
        v = tf.divide(tf.boolean_mask(tensor=s2, mask=is_close), self.thresh_s2, name='v')

        # Row_lengths, for close observations only
        # is_close_r = tf.RaggedTensor.from_row_lengths(values=is_close, row_lengths=self.row_lengths, name='is_close_r')
        ragged_map_func = lambda x : tf.RaggedTensor.from_row_lengths(values=x, row_lengths=self.row_lengths)
        is_close_r = tf.keras.layers.Lambda(function=ragged_map_func, name='is_close_r')(is_close)
        row_lengths_close = tf.reduce_sum(tf.cast(is_close_r, tf.int32), axis=1, name='row_lengths_close')

        # Shape of parameters
        close_size = tf.reduce_sum(row_lengths_close)
        param_shape = (close_size,)

        # Upsample h and lambda
        h_rep = tf.repeat(input=h, repeats=row_lengths_close, name='h_rep')
        h_vec = tf.reshape(tensor=h_rep, shape=param_shape, name='h_vec')
        lam_rep = tf.repeat(input=lam, repeats=row_lengths_close, name='lam_rep')
        lam_vec = tf.reshape(tensor=lam_rep, shape=param_shape, name='lam_vec')

        # Probability according to mixture model
        # p = h_vec * tf.exp(-lam_vec * v) + (1.0 - h_vec)
        p_hit_cond = tf.exp(-lam_vec * v, name='p_hit_cond') 
        p_hit = tf.multiply(h_vec, p_hit_cond, name='p_hit')
        p_miss = tf.subtract(1.0, h_vec, name='p_miss')
        p = tf.add(p_hit, p_miss, name='p')
        log_p_flat = keras.layers.Activation(tf.math.log, name='log_p_flat')(p)

        # Rearrange to ragged tensors
        log_p = tf.RaggedTensor.from_row_lengths(values=log_p_flat, row_lengths=row_lengths_close, name='log_p')
        # ragged_map_func_close = lambda x : tf.RaggedTensor.from_row_lengths(values=x, row_lengths=row_lengths_close)
        # log_p = tf.keras.layers.Lambda(function=ragged_map_func_close, name='log_p')(log_p_flat)

        # Log likelihood by element
        log_like = tf.reduce_sum(log_p, axis=1, name='log_like')

        # Return the log likelihood by element
        return log_like

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
# Functional API model
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_model_asteroid_search(elts: pd.DataFrame,
                               ztf_elt: pd.DataFrame,
                               site_name: str='geocenter',
                               h: float = 0.5,
                               lam: float = 1.0,
                               thresh_deg: float = 1.0):
    """
    Make functional API model for scoring elements
    INPUTS:
        elts:       DataFrame with initial guess for orbital elements.
                    Columns: element_id, a, e, inc, Omega, omega, f, epoch
                    Output of orbital_element_batch, perturb_elts or random_elts
        ztf_elt:    DataFrame with ZTF observations within thresh_deg degrees of
                    of the orbits predicted by these elements.
                    Output of make_ztf_batch or load_ztf_batch
        site_name:  Used for topos adjustment, e.g. 'geocenter' or 'palomar'
        h:          Initial value of hit probability in mixture model
        lam:        Initial value of exponential decay parameter in mixture model
    """

    # Element batch size comes from elts
    elt_batch_size = elts.shape[0]

    # Numpy array and tensor of observation times; flat, shape (data_size,)
    ts_np = ztf_elt.mjd.values.astype(dtype_np)
    ts = keras.backend.constant(value=ts_np, shape=ts_np.shape, dtype=dtype, name='ts')

    # Get observation count per element
    row_lengths_np = ztf_elt.element_id.groupby(ztf_elt.element_id).count()
    row_lengths = keras.backend.constant(value=row_lengths_np, shape=row_lengths_np.shape, 
                                         dtype=tf.int32, name='row_lengths')

    # Shape of the observed trajectories
    data_size = ztf_elt.shape[0]
    traj_shape = (data_size, space_dims)

    # Input is observed directions as a ragged tensor of shape (elt_batch_size, (num_obs), 3,)
    u_obs = keras.Input(shape=(3,), batch_size=data_size, name='u_obs')

    # # Flatten the ragged tensor of observations
    # ragged_map_func = lambda x : tf.RaggedTensor.from_row_lengths(values=x, row_lengths=row_lengths)
    # u_obs_r = tf.keras.layers.Lambda(function=ragged_map_func, name='u_obs_r')(u_obs)
    # u_obs_flat = u_obs_r.values

    # Dummy inputs: weights on the elements
    elt_wts = keras.Input(shape=(), batch_size=elt_batch_size, name='elt_wt')

    # # Observed directions; extract from ztf_elt DataFrame
    # cols_u_obs = ['ux', 'uy', 'uz']
    # u_obs_np = ztf_elt[cols_u_obs].values.astype(dtype_np)
    # u_obs = keras.backend.constant(value=u_obs_np, dtype=dtype, name='u_obs')

    # Set of trainable weights with candidate orbital elements; initialize according to elts
    elements_layer = OrbitalElements(elts=elts, elt_batch_size=elt_batch_size, h=h, lam=lam, name='candidates')
    
    # Extract the candidate elements and mixture parameters; pass dummy inputs to satisfy keras Layer API
    a, e, inc, Omega, omega, f, epoch, h, lam = elements_layer(inputs=None)
    
    # Alias the orbital elements; a, e, inc, Omega, omega, and f are trainable; epoch is fixed
    a = Identity(name='a')(a)
    e = Identity(name='e')(e)
    inc = Identity(name='inc')(inc)
    Omega = Identity(name='Omega')(Omega)
    omega = Identity(name='omega')(omega)
    f = Identity(name='f')(f)
    epoch = Identity(name='epoch')(epoch)

    # Alias the mixture model parameters
    h = Identity(name='h')(h)
    lam = Identity(name='lam')(lam)

    # The orbital elements; stack to shape (elt_batch_size, 7)
    elts_tf = tf.stack(values=[a, e, inc, Omega, omega, f, epoch], axis=1, name='elts')

    # The predicted direction
    direction_layer = AsteroidDirection(ts=ts, row_lengths=row_lengths, site_name=site_name, name='direction_layer')

    # Calibration arrays (flat)
    cols_q_ast = ['qx', 'qy', 'qz']
    cols_v_ast = ['vx', 'vy', 'vz']
    q_ast = ztf_elt[cols_q_ast].values.astype(dtype_np)
    v_ast = ztf_elt[cols_v_ast].values.astype(dtype_np)

    # Run calibration
    direction_layer.q_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)

    # Tensor of predicted directions
    u_pred, r_pred = direction_layer(a, e, inc, Omega, omega, f, epoch)

    # Score layer for these observations
    score_layer = TrajectoryScore(row_lengths=row_lengths, thresh_deg=thresh_deg, name='score_layer')

    # Compute the log likelihood by element from the predicted direction and mixture model parameters
    log_like = score_layer(u_pred, u_obs=u_obs, h=h, lam=lam)

    # Compute the weighted sum log likelihood so inputs and outputs are connected
    # log_like_wtd = tf.reduce_sum(tf.multiply(elt_wts, log_like), name='log_like_wtd')

    # Wrap inputs and outputs
    inputs = (u_obs)
    outputs = (log_like, elts_tf, u_pred)

    # Create model with functional API
    model = keras.Model(inputs=inputs, outputs=outputs)

    # Bind the custom layers to model
    model.elements = elements_layer
    model.direction = direction_layer
    model.score = score_layer
    
    # Add the loss function - the NEGATIVE of the log likelihood
    # (Take negative b/c TensorFlow minimizes the loss function)
    model.add_loss(-tf.reduce_sum(log_like))

    return model

# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
