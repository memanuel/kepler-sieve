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
from asteroid_model import AsteroidDirection
from asteroid_integrate import calc_ast_pos
from candidate_element import elts_np2df, perturb_elts
from astro_utils import deg2dist, dist2deg
from tf_utils import tf_quiet, Identity

# Typing
from typing import Optional, Union

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

# # Range for exponential decay parameter lambda
# lam_min_ = 1.0
# lam_max_ = 2.0**24
# log_lam_min_ = np.log(lam_min_)
# log_lam_max_ = np.log(lam_max_)

# Range for resolution R: 1.0 arc second to 2.0 degrees
R_min_sec_ = 1.0
R_max_deg_ = 2.0
R_min_ = deg2dist(R_min_sec_ / 3600.0)
R_max_ = deg2dist(R_max_deg_)
log_R_min_ = np.log(R_min_)
log_R_max_ = np.log(R_max_)

# ********************************************************************************************************************* 
def R2lam(R, thresh_s):
    """Convert a resolution and threshold in degrees to a decay parameter lambda"""
    # Squares
    R2 = R**2
    thresh_s2 = thresh_s**2
    lam = thresh_s2 / (2.0 * R2)
    return lam

def Rdeg2lam(R_deg, thresh_deg):
    """Convert a resolution and threshold in degrees to a decay parameter lambda"""
    # Convert from degrees to Cartesian distance
    R = deg2dist(R_deg)
    thresh_s = deg2dist(thresh_deg)
    return R2lam(R=R, thresh_s=thresh_s)

# ********************************************************************************************************************* 
# Custom Layers
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
class CandidateElements(keras.layers.Layer):
    """Custom layer to maintain state of candidate orbital elements and resolutions."""

    def __init__(self, 
                 elts: pd.DataFrame, 
                 h: Union[float, np.ndarray] = 0.125,
                 R: Union[float, np.ndarray] = deg2dist(1.0),
                 thresh_deg: float = 1.0,
                 **kwargs):
        super(CandidateElements, self).__init__(**kwargs)
        
        # Configuration for serialization
        self.cfg = {
            'elts': elts,
            'h': h,
            'R': R,
        }

        # Infer batch size from elts DataFrame
        batch_size = elts.shape[0]
        elt_shape = (batch_size,)

        # Threshold distance; 0.5 thresh_s^2 used to convert R to lambda
        self.half_thresh_s2 = keras.backend.constant(0.5 * deg2dist(thresh_deg)**2)

        # If h was passed as a scalar, promote it to an array
        h_init = h if isinstance(h, np.ndarray) else np.full(shape=elt_shape, fill_value=h)
        # If lam was passed as a scalar, promote it to an array
        # lam_init = lam if isinstance(lam, np.ndarray) else np.full(shape=elt_shape, fill_value=lam)
        # If R was passed as a scalar, promote it to an array
        R_init = R if isinstance(R, np.ndarray) else np.full(shape=elt_shape, fill_value=R)

        # Save batch size, orbital elements as numpy array
        self.batch_size = batch_size
        self.elts = elts
        
        # Control over a_, in range 0.0 to 1.0
        self.a_min = keras.backend.constant(a_min_, dtype=dtype)
        self.log_a_range = keras.backend.constant(tf.math.log(a_max_) - tf.math.log(a_min_), dtype=dtype)
        self.a_ = tf.Variable(initial_value=self.inverse_a(elts['a']), trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='a_')
        
        # Control over e_, in range e_min to e_max
        self.e_min = keras.backend.constant(e_min_, dtype=dtype)
        self.e_max = keras.backend.constant(e_max_, dtype=dtype)
        # self.e_range = tf.constant(e_max_ - e_min_, dtype=dtype)
        # self.e_ = tf.Variable(initial_value=self.inverse_e(elts['e']), trainable=True, 
        #                      constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='e_')
        self.e_ = tf.Variable(initial_value=elts['e'], trainable=True, 
                              constraint=lambda t: tf.clip_by_value(t, self.e_min, self.e_max), name='e_')
        
        # Control over inc_, in range -pi/2 to pi/2
        self.inc_max = keras.backend.constant(np.pi/2*(1-2**-20), dtype=dtype)
        self.inc_min = -self.inc_max
        self.inc_range = keras.backend.constant(self.inc_max - self.inc_min, dtype=dtype)
        self.inc_ = tf.Variable(initial_value=self.inverse_inc(elts['inc']), trainable=True, 
                                constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='inc_')
        
        # Scale factor for unconstrained angles is 2*pi
        self.two_pi = keras.backend.constant(2*np.pi, dtype=dtype)
        
        # Angle variables Omega, omega and f
        self.Omega_ = tf.Variable(initial_value=self.inverse_angle(elts['Omega']), trainable=True, name='Omega_')
        self.omega_ = tf.Variable(initial_value=self.inverse_angle(elts['omega']), trainable=True, name='omega_')
        self.f_ = tf.Variable(initial_value=self.inverse_angle(elts['f']), trainable=True, name='f_')

        # The epoch is not trainable
        self.epoch = tf.Variable(initial_value=elts['epoch'], trainable=False, name='epoch')
        
        # Control of the hit rate h_, in range h_min to h_max
        self.h_min = keras.backend.constant(h_min_, dtype=dtype)
        self.h_max = keras.backend.constant(h_max_, dtype=dtype)
        self.h_ = tf.Variable(initial_value=h_init, trainable=True, dtype=dtype,
                              constraint=lambda t: tf.clip_by_value(t, self.h_min, self.h_max), name='h_')

        # # Control of the exponential decay paramater lam_, in range lam_min to lam_max
        # This is controlled indirectly, by controlling the resolution R, then converting it to a decay rate lambda

        # Control of the resolution paramater R_, in range R_min to R_max
        self.R_min = keras.backend.constant(R_min_, dtype=dtype)
        self.log_R_range = keras.backend.constant(log_R_max_ - log_R_min_, dtype=dtype)
        self.R_ = tf.Variable(initial_value=self.inverse_R(R_init), trainable=True, dtype=dtype,
                                constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='R_')
    @tf.function
    def get_a(self):
        """Transformed value of a"""
        return self.a_min * tf.exp(self.a_ * self.log_a_range)

    def inverse_a(self, a):
        """Inverse transform value of a"""
        return tf.math.log(a / self.a_min) / self.log_a_range

    @tf.function
    def get_e(self):
        """Transformed value of e"""
        # return self.e_min + self.e_ * self.e_range
        return self.e_

    # def inverse_e(self, e):
    #    """Inverse transform value of e"""
    #    return (e - self.e_min) / self.e_range

    @tf.function
    def get_inc(self):
        """Transformed value of inc"""
        return self.inc_min + self.inc_ * self.inc_range

    def inverse_inc(self, inc):
        """Inverse transform value of inc"""
        return (inc - self.inc_min) / self.inc_range

    @tf.function
    def get_angle(self, angle_):
        """Forward transform of an unconstrained angle variable (Omega, omega, f)"""
        return tf.multiply(self.two_pi, angle_)

    def inverse_angle(self, angle):
        """Forward transform of an unconstrained angle variable (Omega, omega, f)"""
        return tf.divide(angle, self.two_pi)

    @tf.function
    def get_Omega(self):
        return self.get_angle(self.Omega_)

    @tf.function
    def get_omega(self):
        return self.get_angle(self.omega_)

    @tf.function
    def get_f(self):
        return self.get_angle(self.f_)

    @tf.function
    def get_h(self):
        """Transformed value of h"""
        return self.h_

    @tf.function
    def get_R(self):
        """Transformed value of R"""
        return self.R_min * tf.exp(self.R_ * self.log_R_range)

    def inverse_R(self, R):
        """Inverse transform value of R"""
        return tf.math.log(R / self.R_min) / self.log_R_range

    @tf.function
    def get_R_deg(self):
        """Transformed value of R in degrees"""
        R = self.get_R()        
        return dist2deg(R)

    @tf.function
    def R_to_lam(self, R):
        """Convert a resolution R to an exponential decay term lambda"""
        return tf.divide(self.half_thresh_s2, tf.square(R))

    @tf.function
    def get_lam(self):
        """Transformed value of lambda"""
        # return self.lam_min * tf.exp(self.lam_ * self.log_lam_range)
        R = self.get_R()
        return self.R_to_lam(R)

    def inverse_lam(self, lam):
        """Inverse transform value of lam"""
        return tf.math.log(lam / self.lam_min) / self.log_lam_range

    @tf.function
    def call(self, inputs=None):
        """Return the current candidate orbital elements and mixture model parameters"""
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
        return (a, e, inc, Omega, omega, f, self.epoch, h, lam, R,)

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class TrajectoryScore(keras.layers.Layer):
    """Score candidate trajectories"""
    def __init__(self, u_obs_np, row_lengths_np: np.ndarray, thresh_deg: float, **kwargs):
        """
        INPUTS:
            u_obs_np:       Observed positions; shape [data_size, 3]
            row_lengths_np: Number of observations for each element; shape [elt_batch_size]
            thresh_deg: threshold in degrees for observations to be included
        """
        super(TrajectoryScore, self).__init__(**kwargs)

        # Configuration for seralization
        self.cfg = {
            'u_obs_np': u_obs_np,
            'row_lengths_np': row_lengths_np,
            'thresh_deg': thresh_deg
        }

        # Sizes of the data as keras constants
        self.elt_batch_size = keras.backend.constant(value=row_lengths_np.shape[0], dtype=tf.int32)
        self.data_size = keras.backend.constant(value=tf.reduce_sum(row_lengths_np), dtype=tf.int32)
        self.row_lengths = keras.backend.constant(value=row_lengths_np, shape=row_lengths_np.shape, dtype=tf.int32)
        u_shape = (self.data_size, space_dims,)        

        # Save the observed directions as a keras constant
        self.u_obs = keras.backend.constant(value=u_obs_np, shape=u_shape, dtype=dtype)

        # The threshold distance and its square
        self.thresh_s = keras.backend.constant(value=deg2dist(thresh_deg), dtype=dtype, name='thresh_s')
        self.thresh_s2 = keras.backend.constant(value=self.thresh_s**2, dtype=dtype, name='thresh_s2')

    @tf.function        
    def call(self, u_pred: tf.Tensor, h: tf.Tensor, lam: tf.Tensor):
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            u_pred: predicted directions with current batch of elements
            u_obs:  observed directions
            h:      fraction of hits in mixture model
            lam:    exponential decay parameter in mixture model
        """
        # Difference between actual and predicted directions
        du = keras.layers.subtract(inputs=[u_pred, self.u_obs], name='du')
        # Squared distance bewteen predicted and observed directions
        s2 = tf.reduce_sum(tf.square(du), axis=(-1), name='s2')
        
        # Filter to only include terms where z2 is within the threshold distance^2
        is_close = tf.math.less(s2, self.thresh_s2, name='is_close')
        
        # Relative distance v on data inside threshold
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
        emlx = tf.exp(-lam_vec * v, name='emlx')
        p_hit_cond_num = tf.multiply(emlx, lam_vec, name='p_hit_cond_num')
        p_hit_cond_den = tf.subtract(1.0, tf.exp(-lam_vec), name='p_hit_cond_den')
        p_hit_cond = tf.divide(p_hit_cond_num, p_hit_cond_den, name='p_hit_cond')
        p_hit = tf.multiply(h_vec, p_hit_cond, name='p_hit')
        p_miss = tf.subtract(1.0, h_vec, name='p_miss')
        p = tf.add(p_hit, p_miss, name='p')
        log_p_flat = keras.layers.Activation(tf.math.log, name='log_p_flat')(p)

        # Rearrange to ragged tensors
        # log_p = tf.RaggedTensor.from_row_lengths(values=log_p_flat, row_lengths=row_lengths_close, name='log_p')
        ragged_map_func_close = lambda x : tf.RaggedTensor.from_row_lengths(values=x, row_lengths=row_lengths_close)
        log_p = tf.keras.layers.Lambda(function=ragged_map_func_close, name='log_p')(log_p_flat)

        # Log likelihood by element
        log_like = tf.reduce_sum(log_p, axis=1, name='log_like')

        # Return the log likelihood by element
        return log_like

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
