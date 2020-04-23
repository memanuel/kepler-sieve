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

# Utility
import time
from datetime import timedelta

# Local imports
from asteroid_model import AsteroidDirection, sigma_mag_normalizer_
from asteroid_integrate import calc_ast_pos
from candidate_element import elts_np2df, perturb_elts
from astro_utils import deg2dist, dist2deg
from tf_utils import tf_quiet, gpu_grow_memory, Identity

# Typing
from typing import List, Tuple, Dict, Optional, Union

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Run TF quietly
tf_quiet()

# Configure TensorFlow to use GPU memory variably
# gpu_grow_memory(verbose=True)

# ********************************************************************************************************************* 
# Constants
space_dims: int = 3

# Data type
dtype = tf.float32
dtype_np = np.float32

# Tensorflow constants for true and false
tf_true = tf.constant(True, dtype=tf.bool)
tf_false = tf.constant(False, dtype=tf.bool)

# ********************************************************************************************************************* 
# Transformations of orbital elements for search
# Range for a
a_min_: float = 0.5
a_max_: float = 32.0

# Range for e
e_min_: float = 0.0
# e_max_: float = 1.0 - 2.0**-10
e_max_: float = 1.0 - 2.0**-6

# Range for inc
sin_inc_max_: float = 1.0 - 2.0**-8
inc_max_: float = np.arcsin(sin_inc_max_)

# Range for num_hits
num_hits_min_: float = 6.0
num_hits_max_: float = 1024.0

# Minimum resolution R: 1.0 arc second to 1.0 degrees
R_min_sec_: float = 1.0
R_min_: float = deg2dist(R_min_sec_ / 3600.0)

# ********************************************************************************************************************* 
def R2lam(R: Union[tf.Tensor, np.ndarray], 
          thresh_s: Union[tf.Tensor, np.ndarray]):
    """
    Convert a resolution and threshold in degrees to a decay parameter lambda
    INPUTS:
        R: Resolution as Cartesian distance; tf.Tensor or np.ndarray
        thresh_s: Threshold as Cartesian distance; tf.Tensor or np.ndarray
        R and thresh_s must have the same size
    """
    # Squares
    R2 = R**2
    thresh_s2 = thresh_s**2
    lam = thresh_s2 / (2.0 * R2)
    return lam

def Rdeg2lam(R_deg: Union[tf.Tensor, np.ndarray], 
             thresh_deg: Union[tf.Tensor, np.ndarray]):
    """
    Convert a resolution and threshold in degrees to a decay parameter lambda
    INPUTS:
        R_deg:      Resolution in degrees; tf.Tensor or np.ndarray
        thresh_deg: Threshold in degrees; tf.Tensor or np.ndarray
        R_deg and thresh_deg must have the same size    
    """
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

    def __init__(self, elts: pd.DataFrame, **kwargs):
        super(CandidateElements, self).__init__(**kwargs)
        
        # Configuration for serialization
        self.cfg: Dict = {
            'elts': elts,
        }

        # Infer batch size from elts DataFrame
        batch_size: int = elts.shape[0]
        # elt_shape: Tuple = (batch_size,)

        # Save batch size, orbital elements as numpy array
        self.batch_size: int = batch_size
        self.elts: pd.DataFrame = elts
        
        # Control over a_, in range 0.0 to 1.0
        self.a_min: keras.backend.constant = keras.backend.constant(a_min_, dtype=dtype)
        self.log_a_range: keras.backend.constant = \
            keras.backend.constant(tf.math.log(a_max_) - tf.math.log(a_min_), dtype=dtype)
        self.a_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_a(elts['a']), trainable=True, dtype=dtype,
                        constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='a_')
        
        # Control over e_, in range e_min to e_max
        self.e_min: keras.backend.constant = keras.backend.constant(e_min_, dtype=dtype)
        self.e_max: keras.backend.constant = keras.backend.constant(e_max_, dtype=dtype)
        self.e_: tf.Variable = \
            tf.Variable(initial_value=elts['e'], trainable=True, dtype=dtype,
                        constraint=lambda t: tf.clip_by_value(t, self.e_min, self.e_max), name='e_')
        
        # Control over inc_, in range -pi/2 to pi/2
        self.inc_max = keras.backend.constant(inc_max_, dtype=dtype)
        self.inc_min = -self.inc_max
        self.inc_range = keras.backend.constant(self.inc_max - self.inc_min, dtype=dtype)
        self.inc_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_inc(elts['inc']), trainable=True, dtype=dtype,
                        constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='inc_')
        
        # Scale factor for unconstrained angles is 2*pi
        self.two_pi: keras.backend.constant = keras.backend.constant(2.0*np.pi, dtype=dtype)
        
        # Angle variables Omega, omega and f
        self.Omega_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_angle(elts['Omega']), 
                        trainable=True, dtype=dtype, name='Omega_')
        self.omega_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_angle(elts['omega']), 
                        trainable=True, dtype=dtype, name='omega_')
        self.f_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_angle(elts['f']), 
                        trainable=True, dtype=dtype, name='f_')

        # The epoch is not trainable
        self.epoch: tf.Variable = \
            tf.Variable(initial_value=elts['epoch'], 
                        trainable=False, dtype=dtype, name='epoch')

    @tf.function
    def get_a(self) -> tf.Tensor:
        """Transformed value of a"""
        return self.a_min * tf.exp(self.a_ * self.log_a_range)

    def inverse_a(self, a) -> tf.Tensor:
        """Inverse transform value of a"""
        return tf.math.log(a / self.a_min) / self.log_a_range

    @tf.function
    def get_e(self) -> tf.Tensor:
        """Transformed value of e"""
        return self.e_

    def inverse_e(self, e) -> tf.Tensor:
        """Inverse transform value of e"""
        return e

    @tf.function
    def get_inc(self) -> tf.Tensor:
        """Transformed value of inc"""
        return self.inc_min + self.inc_ * self.inc_range

    def inverse_inc(self, inc) -> tf.Tensor:
        """Inverse transform value of inc"""
        return (inc - self.inc_min) / self.inc_range

    @tf.function
    def get_angle(self, angle_) -> tf.Tensor:
        """Forward transform of an unconstrained angle variable (Omega, omega, f)"""
        return tf.multiply(self.two_pi, angle_)

    def inverse_angle(self, angle) -> tf.Tensor:
        """Forward transform of an unconstrained angle variable (Omega, omega, f)"""
        return tf.divide(angle, self.two_pi)

    @tf.function
    def get_Omega(self) -> tf.Tensor:
        return self.get_angle(self.Omega_)

    @tf.function
    def get_omega(self) -> tf.Tensor:
        return self.get_angle(self.omega_)

    @tf.function
    def get_f(self) -> tf.Tensor:
        return self.get_angle(self.f_)

    @tf.function
    def call(self, inputs=None) -> Tuple[tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor,]:
        """
        Return the current candidate orbital elements and mixture model parameters
        INPUTS:
            inputs - a dummy input to satisfy the keras API.  Any tensor of shape [batch_size,].
        OUTPUTS:
            a, e, inc, Omega, omega, f, epoch
            Tuple of the seven orbital elements including epoch
        """
        # print(f'type(inputs)={type(inputs)}.')
        # Transform a, e from log to linear
        a: tf.Tensor = self.get_a()
        e: tf.Tensor = self.get_e()
        # Angles transformed by factor of 2*pi
        inc: tf.Tensor = self.get_inc()
        Omega: tf.Tensor = self.get_Omega()
        omega: tf.Tensor = self.get_omega()
        f: tf.Tensor = self.get_f()
        # Return the elements
        return (a, e, inc, Omega, omega, f, self.epoch,)

    def load(self, elts):
        """Load values from elts DataFrame"""
        # Compute inverse transforms of the orbital elements
        a_ = self.inverse_a(elts['a'].values)
        e_ = self.inverse_e(elts['e'].values)
        inc_ = self.inverse_inc(elts['inc'].values)
        Omega_ = self.inverse_angle(elts['Omega'].values)
        omega_ = self.inverse_angle(elts['omega'].values)
        f_ = self.inverse_angle(elts['f'].values)
        epoch = elts['epoch'].values
        # Assign values to the weight variables
        self.a_.assign(a_)
        self.e_.assign(e_)
        self.inc_.assign(inc_)
        self.Omega_.assign(Omega_)
        self.omega_.assign(omega_)
        self.f_.assign(f_)
        self.epoch.assign(epoch)

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class MixtureParameters(keras.layers.Layer):
    """Custom layer to maintain state of mixture parameters."""

    def __init__(self, elts: pd.DataFrame, **kwargs):
        super(MixtureParameters, self).__init__(**kwargs)
        
        # Configuration for serialization
        self.cfg: Dict = {
            'elts': elts,
        }

        # Infer batch size from elts DataFrame
        batch_size: int = elts.shape[0]
        # elt_shape = (batch_size,)

        # Save batch size, orbital elements as numpy array
        self.batch_size: int = batch_size
        self.elts: pd.DataFrame = elts
        
        # Control of the number of hits num_hits with num_hits_
        # Min and max of num_hits are static; controlled gloablly
        self.num_hits_min: keras.backend.constant = keras.backend.constant(value=num_hits_min_, dtype=dtype)
        self.num_hits_max: keras.backend.constant = keras.backend.constant(value=num_hits_max_, dtype=dtype)
        # Dynamic range of num_hits
        self.log_num_hits_range: keras.backend.constant = \
            keras.backend.constant(value=np.log(num_hits_max_ / num_hits_min_), dtype=dtype)
        # Tensor with control variable for num_hits
        self.num_hits_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_num_hits(elts['num_hits']), trainable=True, dtype=dtype,
                        constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='num_hits_')

        # Control of the resolution parameter R_, in range R_min to R_max
        # R_min is set globally; distance corresponding to 1.0 arc second
        R_min: float = R_min_
        thresh_s: np.ndarray = elts['thresh_s'].values
        # Set R_max to half the initial threshold distance 
        R_max: np.ndarray = 0.5 * thresh_s
        # Dynamic range of resolution R
        log_R_range: np.ndarray = np.log(R_max / R_min)
        # Save these as keras constants or variables as appropriate
        self.R_min: keras.backend.constant = keras.backend.constant(R_min, dtype=dtype)
        self.log_R_range: tf.Variable = \
                tf.Variable(initial_value=log_R_range, trainable=False, dtype=dtype, name='log_R_range')
        # Tensor with control variable for R; shape [batch_size,]
        self.R_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_R(elts['R']), trainable=True, dtype=dtype,
                        constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0), name='R_')

    @tf.function
    def get_num_hits(self) -> tf.Tensor:
        """Transformed value of h"""
        return self.num_hits_min * tf.exp(self.num_hits_ * self.log_num_hits_range)

    def inverse_num_hits(self, num_hits) -> tf.Tensor:
        """Inverse transform value of num_hits"""
        return tf.math.log(num_hits / self.num_hits_min) / self.log_num_hits_range

    def set_num_hits(self, num_hits: np.ndarray) -> tf.Tensor:
        """Set num_hits parameter"""
        self.num_hits_.assign(self.inverse_num_hits(num_hits))

    @tf.function
    def get_R(self) -> tf.Tensor:
        """Transformed value of R"""
        return self.R_min * tf.exp(self.R_ * self.log_R_range)

    def inverse_R(self, R) -> tf.Tensor:
        """Inverse transform value of R"""
        return tf.math.log(R / self.R_min) / self.log_R_range

    def set_R(self, R: np.ndarray) -> None:
        """Set the resolution parameter"""
        self.R_.assign(self.inverse_R(R))

    @tf.function
    def get_R_deg(self) -> tf.Tensor:
        """Transformed value of R in degrees"""
        R = self.get_R()        
        return dist2deg(R)
    
    @tf.function
    def get_R_max(self) -> tf.Tensor:
        """Maximum value of R"""
        return self.R_min * tf.exp(self.log_R_range)

    def set_R_max(self, R_max: np.ndarray) -> None:
        """Set the maximum of the resolution parameter"""
        # Get old values of R
        R_old = self.get_R()
        # Apply the constraint to the current values; this is the new value
        R = np.minimum(R_old, R_max)
        # Update the variable with the dynamic range of log_R
        log_R_range = np.log(R_max / self.R_min.numpy())
        self.log_R_range.assign(log_R_range)
        # Assign the updated R back to the layer
        self.set_R(R)

    def set_R_deg_max(self, R_deg_max: np.ndarray) -> None:
        """Set the maximum of the resolution parameter in degrees"""
        # Convert R_deg_max to distance
        R_max = deg2dist(R_deg_max)
        # Delegate to set_R_max
        self.set_R_max(R_max)

    @tf.function
    def call(self, inputs=None) -> None:
        """Return the current mixture model parameters"""
        # Transform search parameters num_hits, R and R_max
        num_hits = self.get_num_hits()
        R = self.get_R()
        R_max = self.get_R_max()
        return (num_hits, R, R_max)

    def load(self, elts):
        """Load values from elts DataFrame"""
        self.set_num_hits(elts.num_hits.values)
        self.set_R_max(elts.R_max.values)
        self.set_R(elts.R.values)

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class TrajectoryScore(keras.layers.Layer):
    """Score candidate trajectories"""
    def __init__(self, 
                 u_obs_np: np.ndarray, 
                 mag_app_np: np.ndarray,
                 row_lengths_np: np.ndarray, 
                 thresh_deg: np.ndarray, 
                 **kwargs):
        """
        INPUTS:
            u_obs_np:         Observed positions; shape [data_size, 3,]
            mag_app_np:       Observed apparent magnitude; shape [data_size, ]
            row_lengths_np:   Number of observations for each element; shape [batch_size,]
            thresh_deg:       Initial setting for the threshold in degrees for observations to be included;
                              Not the same as the threshold in degrees for the original observation data!
                              The threshold for the data never changes, but the live value of the scoring
                              threshold parameter is updated as the model trains.
        """
        super(TrajectoryScore, self).__init__(**kwargs)

        # Configuration for seralization
        self.cfg : Dict = {
            'u_obs_np': u_obs_np,
            'mag_app_np': mag_app_np,
            'row_lengths_np': row_lengths_np,
            'thresh_deg': thresh_deg
        }

        # Sizes of the data as keras constants
        self.batch_size: keras.backend.constant = \
            keras.backend.constant(value=row_lengths_np.shape[0], dtype=tf.int32)
        self.data_size: keras.backend.constant = \
            keras.backend.constant(value=tf.reduce_sum(row_lengths_np), dtype=tf.int32)
        self.row_lengths: keras.backend.constant = \
            keras.backend.constant(value=row_lengths_np, shape=row_lengths_np.shape, dtype=tf.int32)
        u_shape: Tuple[int, int] = (int(self.data_size), space_dims,)
        mag_shape: Tuple[int] = (int(self.data_size), )
        
        # Save the observed directions as a keras constant
        self.u_obs: keras.backend.constant = \
            keras.backend.constant(value=u_obs_np, shape=u_shape, dtype=dtype)

        # Save the observed apparent magnitude as a keras constant
        self.mag_obs: keras.backend.constant = \
            keras.backend.constant(value=mag_app_np, shape=mag_shape, dtype=dtype)

        # Normalizer for sigma_mag; multiplty magnitude PDF by this so it has mean ~1
        self.sigma_mag_normalizer: keras.backend.constant = \
               keras.backend.constant(value=np.exp(1.0) * sigma_mag_normalizer_, dtype=dtype)

        # Mathematical constant for the normal PDF
        self.inv_root_2pi: keras.backend.constant = \
                keras.backend.constant(value=1.0 / np.sqrt(2.0 * np.pi), name='inv_root_2pi')

        # The threshold distance and its square; numpy arrays of shape [batch_size,]
        thresh_s: np.ndarray = deg2dist(thresh_deg)
        thresh_s2: np.ndarray = thresh_s**2

        # Support for making thresh_deg trainable
        # Allowed range for thresh_s2: from 10 arc seconds to the initial threshold
        thresh_s2_min: dtype_np = deg2dist(10.0/3600)**2
        thresh_s2_max: np.ndarray = thresh_s2
        log_thresh_s2_range: np.ndarray = np.log(thresh_s2_max / thresh_s2_min)
        # Save these as keras constants or variables as appropriate
        self.thresh_s2_min: keras.backend.constant = \
            keras.backend.constant(value=thresh_s2_min, dtype=dtype)
        self.log_thresh_s2_range: tf.Variable = \
            tf.Variable(initial_value=log_thresh_s2_range, trainable=False, 
                        dtype=dtype, name='log_thresh_s2_range')

        # Tensor with control variable for thresh_s2; shape [batch_size]
        self.thresh_s2_: tf.Variable = \
            tf.Variable(initial_value=self.inverse_thresh_s2(thresh_s2), trainable=True, 
                        constraint=lambda t: tf.clip_by_value(t, 0.0, 1.0),
                        dtype=dtype, name='thresh_s2_')

        # Threshold distance for counting a hit; 10.0 arc seconds
        thresh_hit_s: float = deg2dist(10.0/3600)
        thresh_hit_s2: float = thresh_hit_s**2
        self.thresh_hit_s2: keras.backend.constant = \
                keras.backend.constant(value=thresh_hit_s2, dtype=dtype)

        # Threshold posterior hit probability for counting a hit
        self.thresh_hit_prob_post: keras.backend.constant = \
            keras.backend.constant(value=0.95, dtype=dtype, name='thresh_hit_prob_post')

        # Maintain state of row_lengths_close to avoid needlessly passing it around between methods
        shape: Tuple[int] = (int(self.batch_size), )
        self.row_lengths_close: tf.Variable = \
            tf.Variable(initial_value=self.row_lengths, shape=shape, 
                        dtype=tf.int32, trainable=False)

        # Variables for calibrating magnitude
        self.mag_is_good: tf.Variable = \
            tf.Variable(initial_value=np.full(shape, False), shape=shape, dtype=dtype, trainable=False)
        self.delta_H: tf.Variable = \
            tf.Variable(initial_value=np.zeros(shape), shape=shape, dtype=dtype, trainable=False)
        self.sigma_hat: tf.Variable = \
            tf.Variable(initial_value=np.ones(shape), shape=shape, dtype=dtype, trainable=False)

    # *******************************************************************************************************
    # Get and set parameters on score layer
    # *******************************************************************************************************

    def get_thresh_s2(self) -> tf.Tensor:
        """Transformed value of thresh_s2"""
        return self.thresh_s2_min * tf.exp(self.thresh_s2_ * self.log_thresh_s2_range)
    
    def inverse_thresh_s2(self, thresh_s2) -> tf.Tensor:
        """Inverse transform value of thresh_s2"""
        return tf.math.log(thresh_s2 / self.thresh_s2_min) / self.log_thresh_s2_range   

    def set_thresh_s2(self, thresh_s2: np.ndarray) -> None:
        """
        Set the threshold parameter for which observations are included 
        in the conditional probability score calculations
        """        
        # The threshold distance and its square; numpy arrays of shape [batch_size,]
        self.thresh_s2_.assign(self.inverse_thresh_s2(thresh_s2))

    def get_thresh_deg(self) -> np.ndarray:
        """Return the current threshold in degrees"""
        thresh_s2: tf.Tensor = self.get_thresh_s2()
        thresh_s: tf.Tensor = tf.math.sqrt(thresh_s2)
        thresh_s_deg: np.ndarray = dist2deg(thresh_s.numpy())
        return thresh_s_deg
    
    def set_thresh_deg(self, thresh_deg: np.ndarray) -> None:
        """
        Set the threshold parameter for which observations are included 
        in the conditional probability score calculations
        INPUTS:
            thresh_deg: Threshold in degrees for including an observation; one for each candidate element
                        (In practice these may be the same.)
        """        
        # The threshold distance and its square; numpy arrays of shape [batch_size,]
        thresh_s: np.ndarray = deg2dist(thresh_deg)
        thresh_s2: np.ndarray = thresh_s**2
        self.set_thresh_s2(thresh_s2)

    def set_thresh_sec(self, thresh_sec: np.ndarray) -> None:
        thresh_deg = thresh_sec / 3600.0
        self.set_thresh_deg(thresh_deg)

    def set_thresh_deg_max(self, thresh_deg_max: np.ndarray) -> None:
        """Set the maximum of the thresh_deg parameter"""
        # Get old values of the threshold
        thresh_s2_old: tf.Tensor = self.get_thresh_s2()
        # Convert the constraint from degrees into s2
        thresh_s_max: np.ndarray = deg2dist(thresh_deg_max)
        thresh_s2_max: np.ndarray = thresh_s_max**2
        # Apply the constraint to the current values
        thresh_s2: np.ndarray = np.minimum(thresh_s2_old, thresh_s2_max)

        # Update the variable with the dynamic range of log_thresh_s2
        log_thresh_s2_range: np.ndarray = np.log(thresh_s2_max / self.thresh_s2_min.numpy())
        self.log_thresh_s2_range.assign(log_thresh_s2_range)
        # Assign the updated thresh_s2 back to the layer
        self.set_thresh_s2(thresh_s2)

    def set_thresh_sec_max(self, thresh_sec_max: np.ndarray) -> None:
        """Set the maximum of the thresh_deg parameter, specified in arc seconds"""
        thresh_deg_max = thresh_sec_max / 3600.0
        self.set_thresh_deg_max(thresh_deg_max)

    # *******************************************************************************************************
    # Calculate conditional probabilities; helper for mixture probability model
    # *******************************************************************************************************

    # *******************************************************************************************************
    @tf.function        
    def calc_dist(self, u_pred: tf.Tensor, thresh_s2_elt: tf.Tensor) \
             -> Tuple[tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor]:
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            u_pred:     predicted directions with current batch of elements
        OUTPUTS:
            s2:         squared distance; only on close elements
            v:          s^2 / thresh^2; only on close elements; in [0, 1], uniformly distributed
            is_close:   mask indicating which observations are within the threshold
        """
        # Difference between actual and predicted directions
        # du: tf.Tensor = keras.layers.subtract(inputs=[u_pred, self.u_obs], name='du')
        du: tf.Tensor = tf.subtract(u_pred, self.u_obs, name='du')
        # Squared distance bewteen predicted and observed directions
        s2_all: tf.Tensor = tf.reduce_sum(tf.square(du), axis=(-1), name='s2')
        
        # Upsample thresh_s2 so it matches input shape
        data_shape: Tuple[int] = (self.data_size,)
        thresh_s2_rep: tf.Tensor = \
            tf.repeat(input=thresh_s2_elt, repeats=self.row_lengths, name='thresh_s2_rep')
        thresh_s2_all: tf.Tensor = \
            tf.reshape(tensor=thresh_s2_rep, shape=data_shape, name='thresh_s2_all')

        # Filter to only include terms where z2 is within the threshold distance^2
        is_close: tf.Tensor = tf.math.less(s2_all, thresh_s2_all, name='is_close')

        # Relative distance v on data inside threshold; shape is [num_close,]
        s2: tf.Tensor = tf.boolean_mask(tensor=s2_all, mask=is_close, name='s2')
        thresh_s2: tf.Tensor = tf.boolean_mask(tensor=thresh_s2_all, mask=is_close, name='thresh_s2')
        v: tf.Tensor = tf.divide(s2, thresh_s2, name='v')

        # Row_lengths, for close observations only
        is_close_r = tf.RaggedTensor.from_row_lengths(
            values=is_close, row_lengths=self.row_lengths, name='is_close_r')
        # Compute the row lengths (number of close observations per candidate element)
        row_lengths_close = tf.reduce_sum(tf.cast(is_close_r, tf.int32), axis=1)

        # Assign row_lengths_close
        # self.is_close.assign(is_close)
        self.row_lengths_close.assign(row_lengths_close)

        return s2, v, is_close

    # *******************************************************************************************************
    @tf.function
    def cond_pdf_dist(self, v: tf.Tensor, R_elt: tf.Tensor, thresh_s2_elt: tf.Tensor):
        """
        Calculate conditional probability due to distance
        INPUTS:
            v:                  Distance squared over thresh_s2, in [0, 1]; shape [data_size,]
            R_elt:              Resolution parameter in Cartesian distance; shape [batch_size,]
            thresh_s2_elt:      The threshold in distance squared; shape [batch_size,]
        OUTPUTS:
            q_hit:              Probability density conditional on a hit
            q_miss:             Probability density conditional on a miss
        """
        # Convert R to exponential decay parameter lam
        half_thresh_s2_elt: tf.Tensor = tf.multiply(thresh_s2_elt, 0.5, name='half_thresh_s2_elt')
        R2_elt: tf.Tensor = tf.square(R_elt, name='R2_elt')
        lam_elt: tf.Tensor = tf.divide(half_thresh_s2_elt, R2_elt, name='lam_elt')

        # Shape of close rows
        data_size_close = tf.reduce_sum(self.row_lengths_close)
        data_shape_close: Tuple[tf.int32] = (data_size_close,)
        
        # Upsample lambda 
        lam_rep: tf.Tensor = tf.repeat(input=lam_elt, repeats=self.row_lengths_close, name='lam_rep')
        lam: tf.Tensor = tf.reshape(tensor=lam_rep, shape=data_shape_close, name='lam')

        # Conditional probability based on distance bewteen predicted and observed direction
        # v|hit ~ Expo(lam) truncated to [0, 1]
        lam_v: tf.Tensor = tf.multiply(lam, v, name='lam_v')
        e_minus_lam_v: tf.Tensor = tf.exp(tf.negative(lam_v), name='e_minus_lam_v')
        e_minus_lam: tf.Tensor = tf.exp(tf.negative(lam), name='e_minus_lam')
        q_hit_num: tf.Tensor = tf.multiply(e_minus_lam_v, lam, name='q_hit_num')        
        q_hit_den: tf.Tensor = tf.subtract(1.0, e_minus_lam, name='one_minus_e_minus_lam')
        q_hit: tf.Tensor = tf.divide(q_hit_num, q_hit_den, name='q_hit')

        # Probability of a miss is just 1.0 because v|miss ~ Unif[0, 1]
        q_miss: tf.Tensor = tf.ones_like(q_hit)

        return q_hit, q_miss

    # *******************************************************************************************************
    @tf.function
    def mag_standardize(self, mag_obs_close):
        """Return the mean, standard deviation and z-score of the magnitude by element"""
        # Shape of close rows
        data_size_close = tf.reduce_sum(self.row_lengths_close)
        data_shape_close: Tuple[tf.int32] = (data_size_close,)
        
        # The mean magnitude of all the observations by element
        mag_obs_r: tf.Tensor = tf.RaggedTensor.from_row_lengths(
            values=mag_obs_close, row_lengths=self.row_lengths_close, name='mag_obs_r')
        mag_mean_elt: tf.Tensor = tf.reduce_mean(mag_obs_r, axis=-1)
        # Upsample the mean to the close rows
        mag_mean_rep: tf.Tensor = tf.repeat(input=mag_mean_elt, repeats=self.row_lengths_close, name='mag_mean_rep')
        mag_mean: tf.Tensor = tf.reshape(mag_mean_rep, shape=data_shape_close, name='mag_mean')

        # Compute the variance by element on the close rows
        mag_demean: tf.Tensor = tf.subtract(mag_obs_close, mag_mean, name='mag_demean')
        mag_demean2: tf.Tensor = tf.square(mag_demean, name='mag_demean2')
        mag_demean2_r: tf.Tensor = tf.RaggedTensor.from_row_lengths(
            values=mag_demean2, row_lengths=self.row_lengths_close, name='mag_demean2_r')        
        mag_var_elt: tf.Tensor = tf.reduce_mean(mag_demean2_r, axis=-1, name='mag_var_elt')
        # The estimated standard deviation by element on the close rows
        sigma_hat_elt: tf.Tensor = tf.math.sqrt(mag_var_elt, name='sigma_hat_elt')

        # Upsample the standard deviation
        sigma_hat_rep: tf.Tensor = tf.repeat(input=sigma_hat_elt, repeats=self.row_lengths_close, name='sigma_hat_rep')
        sigma_hat: tf.Tensor = tf.reshape(sigma_hat_rep, shape=data_shape_close, name='sigma_hat')

        # Difference from the mean
        z: tf.Tensor = tf.divide(mag_demean, sigma_hat, name='mag_z')

        return z, mag_mean, sigma_hat

    # *******************************************************************************************************
    @tf.function
    def cond_pdf_mag(self, mag_pred: tf.Tensor, sigma_mag: tf.Tensor, is_close: tf.Tensor):
        """
        Calculate conditional probability due to magnitude
        INPUTS:
            mag_pred:       Predicted magnitude; shape [data_size,]
            sigma_mag_elt:  Standard deviation of magnitude; shape [batch_size,]
            is_close:       Mask indicating which observations are within the threshold distance    
        OUTPUTS:
            q_hit:      Probability density conditional on a hit
            q_miss:     Probability density conditional on a miss
        """
        # Mask the observed and predicted observations, and sigma_mag, to the close rows
        mag_obs_close: tf.Tensor = tf.boolean_mask(tensor=self.mag_obs, mask=is_close, name='mag_obs_close')      
        mag_pred_close: tf.Tensor = tf.boolean_mask(tensor=mag_pred, mask=is_close, name='mag_pred_close')
        sigma_mag_close: tf.Tensor = tf.boolean_mask(tensor=sigma_mag, mask=is_close, name='sigma_mag_close')

        # Difference between predicted and observed magnitude
        mag_diff_close: tf.Tensor = tf.subtract(mag_pred_close, mag_obs_close, name='mag_diff_close')

        # Conditional probability given a hit based on difference between predicted and observed magnitude
        z_hit: tf.Tensor = tf.divide(mag_diff_close, sigma_mag_close, name='z_hit')
        z2_hit: tf.Tensor = tf.square(z_hit, name='z2_hit')
        arg_hit: tf.Tensor = tf.multiply(z2_hit, -0.5, name='arg_hit')
        exp_arg_hit: tf.Tensor = tf.math.exp(arg_hit, name='exp_arg_hit')
        coef_hit: tf.Tensor = tf.divide(self.inv_root_2pi, sigma_mag_close, name='coef_hit')
        q_hit: tf.Tensor = tf.multiply(coef_hit, exp_arg_hit, name='q_hit')

        # Mean and standard deviation of magnitude for PDF conditional on a miss
        z_miss, mag_mean, sigma_hat = self.mag_standardize(mag_obs_close)

        # Conditional probability given a miss based on normal distribution matching the data
        z2_miss: tf.Tensor = tf.square(z_miss, name='z2_miss')
        arg_miss: tf.Tensor = tf.multiply(z2_miss, -0.5, name='arg_miss')
        exp_arg_miss: tf.Tensor = tf.math.exp(arg_miss, name='exp_arg_miss')
        coef_miss: tf.Tensor = tf.divide(self.inv_root_2pi, sigma_hat, name='coef_miss')
        q_miss: tf.Tensor = tf.multiply(coef_miss, exp_arg_miss, name='q_miss')

        # Arrange mag_diff_close as a ragged tensor for computing delta_H required to neutralize it
        delta_H_r: tf.Tensor = tf.RaggedTensor.from_row_lengths(
            values=tf.negative(mag_diff_close), row_lengths=self.row_lengths_close, name='delta_H_r')
        sigma_hat_r: tf.Tensor = tf.RaggedTensor.from_row_lengths(
            values=sigma_hat, row_lengths=self.row_lengths_close, name='sigma_hat_r')

        # Take the mean to get the elementwise delta_H and sigma_hat
        delta_H: tf.Tensor = tf.reduce_mean(delta_H_r, axis=-1)
        sigma_hat: tf.Tensor = tf.reduce_mean(sigma_hat_r, axis=-1)

        # Assign these to variables
        self.delta_H.assign(delta_H)
        self.sigma_hat.assign(sigma_hat)

        # Return the conditional PDFs
        return q_hit, q_miss 

    # *******************************************************************************************************
    # Calculate log likelihood using mixture probability model between hits and misses
    # *******************************************************************************************************

    # *******************************************************************************************************
    @tf.function
    def mixture_model(self, q_hit: tf.Tensor, q_miss: tf.Tensor, num_hits: tf.Tensor):
        """
        Compute log likelihood under a mixture model
        INPUTS:
            q_hit:      Conditional PDF, conditioned on this observation being a hit
            q_miss:     Conditional PDF, conditioned on this observation being a miss
            num_hits:   Number of hits in the mixture model
        """

        # Shape of close rows
        data_size_close = tf.reduce_sum(self.row_lengths_close)
        data_shape_close: Tuple[tf.int32] = (data_size_close,)

        # Row lengths as float
        row_lengths_close_float: tf.Tensor = tf.cast(x=self.row_lengths_close, dtype=dtype)
        # Compute the implied hit rate h from the number of hits and row_lengths_close; shape [batch_size,]
        h_elt: tf.Tensor = tf.divide(num_hits, row_lengths_close_float, name='h_elt')
        # The hit rate must be in [0, 1]
        h_elt = tf.clip_by_value(h_elt, clip_value_min=0.0, clip_value_max=1.0)

        # Upsample h and lambda
        h_rep: tf.Tensor = tf.repeat(input=h_elt, repeats=self.row_lengths_close, name='h_rep')
        h: tf.Tensor = tf.reshape(tensor=h_rep, shape=data_shape_close, name='h')

        # Complement of h
        one_minus_h: tf.Tensor = tf.subtract(1.0, h, name='one_minus_h')

        # Two contributions to PDF under the mixture model: from hits and misses
        p_mix_hit: tf.Tensor = tf.multiply(h, q_hit, name='p_mix_hit')
        p_mix_miss: tf.Tensor = tf.multiply(one_minus_h, q_miss, name='p_mix_miss')
        # PDF according to mixture model and baseline with no hits
        p_mix: tf.Tensor = tf.add(p_mix_hit, p_mix_miss, name='p_mix')
        p_miss: tf.Tensor = q_miss

        # Log of PDF under the mixture model and the baseline 
        log_p_mix: tf.Tensor = tf.math.log(p_mix, name='log_p_mix')
        log_p_miss: tf.Tensor = tf.math.log(p_miss, name='log_p_miss')
        
        # The posterior hit probability is p_hit / p_mix
        p_hit_post: tf.Tensor = tf.divide(p_mix_hit, p_mix)

        # Rearrange to ragged tensors
        log_p_mix_r = tf.RaggedTensor.from_row_lengths(
            values=log_p_mix, row_lengths=self.row_lengths_close, name='log_p_mix_r')
        log_p_miss_r = tf.RaggedTensor.from_row_lengths(
            values=log_p_miss, row_lengths=self.row_lengths_close, name='log_p_miss_r')

        # Compute totals by element
        log_like_mix: tf.Tensor = tf.reduce_sum(log_p_mix_r, axis=-1, name='log_like_mix')
        log_like_miss: tf.Tensor = tf.reduce_sum(log_p_miss_r, axis=-1, name='log_like_miss')

        # Report the differnce in log likelihood for the mixture model over the baseline of all misses
        log_like: tf.Tensor = tf.subtract(log_like_mix, log_like_miss, name='log_like')

        return log_like, p_hit_post

    # *******************************************************************************************************
    @tf.function
    def calc_hits(self, s2: tf.Tensor, p_hit_post: tf.Tensor):
        """
        Compute log likelihood under a mixture model
        INPUTS:
            q_hit:      Conditional PDF, conditioned on this observation being a hit
            s2:         Distance squared; used for counting hits
        """
        # Filter quality hits: only those close enough
        is_near_hit: tf.Tensor = tf.math.less(s2, self.thresh_hit_s2, name='is_near_hit')
        # Filter quality hits: only those with probability above threshold
        # is_high_prob_hit: tf.Tensor = \
        #    tf.math.greater(p_hit_post, self.thresh_hit_prob_post, name='is_high_prob_hit')
        # Quality hits are close enough with high enough posterior probability
        # is_quality_hit: tf.Tensor = tf.math.logical_and(is_near_hit, is_high_prob_hit, name='is_quality_hit')
        is_quality_hit = is_near_hit
        # Count the effective number of quality hits
        p_hit: tf.Tensor = tf.where(condition=is_quality_hit, x=p_hit_post, y=0.0, name='p_hit_quality')

        # Rearrange to ragged tensors
        p_hit_r = tf.RaggedTensor.from_row_lengths(
            values=p_hit, row_lengths=self.row_lengths_close, name='p_hit_quality_r')

        # Compute totals by element
        hits: tf.Tensor = tf.reduce_sum(p_hit_r, axis=1, name='hits')

        return hits

    # *******************************************************************************************************
    @tf.function        
    def calc_score_joint(self, 
             u_pred: tf.Tensor, 
             mag_pred: tf.Tensor,
             num_hits: tf.Tensor, 
             R: tf.Tensor,
             sigma_mag: tf.Tensor) \
             -> Tuple[tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor]:
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            u_pred:     predicted directions with current batch of elements
            mag_pred:   predicted asteroid magnitude
            num_hits:   number of hits in mixture model
            R:          resolution parameter (Cartesian distance)
            sigma_mag:  standard deviation of asteroid magnitude
        """
        # Alias the input R to indicate that it is by element; sigma_mag passed for all rows
        R_elt: tf.Tensor = R

        # Get the threshold by element
        thresh_s2_elt = self.get_thresh_s2()

        # Calculate distance squared s2 and relative distance squared v
        s2, v, is_close = self.calc_dist(u_pred=u_pred, thresh_s2_elt=thresh_s2_elt)

        # Conditional probability q in distance model, conditional on a hit or a miss
        q_hit_d: tf.Tensor
        q_miss_d: tf.Tensor
        q_hit_d, q_miss_d = \
                self.cond_pdf_dist(v=v, R_elt=R_elt, thresh_s2_elt=thresh_s2_elt)

        # Conditional probability q in magnitude model
        q_hit_m: tf.Tensor
        q_miss_m: tf.Tensor
        q_hit_m, q_miss_m = \
            self.cond_pdf_mag(mag_pred=mag_pred, sigma_mag=sigma_mag, is_close=is_close)

        # PDF in the joint model
        q_hit_j: tf.Tensor = tf.multiply(q_hit_d, q_hit_m, name='q_hit')
        q_miss_j: tf.Tensor = tf.multiply(q_miss_d, q_miss_m, name='q_miss')

        # Compute the mixture model for distance
        log_like_d, p_hit_post = self.mixture_model(q_hit=q_hit_d, q_miss=q_miss_d, num_hits=num_hits)

        # Mixture model for magnitude
        log_like_m, _ = self.mixture_model(q_hit=q_hit_m, q_miss=q_miss_m, num_hits=num_hits)

        # Mixture model for joint
        log_like_j, _ = self.mixture_model(q_hit=q_hit_j, q_miss=q_miss_j, num_hits=num_hits)

        # Non-negative part of magnitude log_like and separable log_like
        log_like_m_pos: tf.Tensor = tf.maximum(log_like_m, 0.0)
        log_like_sep: tf.Tensor = tf.add(log_like_d, log_like_m_pos, name='log_like_sep')
        # Better of joint and distance log_like
        log_like_best: tf.Tensor = tf.maximum(log_like_j, log_like_d, 'log_like_j_best')

        # Blended log_like
        log_like_sep_wtd = tf.multiply(log_like_sep, 0.50)
        log_like_best_wtd = tf.multiply(log_like_best, 0.50)
        log_like = tf.add(log_like_sep_wtd, log_like_best_wtd)

        # Count hits
        hits = self.calc_hits(s2=s2, p_hit_post=p_hit_post)

        # Flag indicating whether the magnitude is good
        # mag_is_good = tf.cast(tf.math.greater(log_like_j, log_like_d), dtype)
        mag_is_good = tf.cast(tf.math.greater(log_like_m, 0.0), dtype)
        self.mag_is_good.assign(mag_is_good)

        # Return the log likelihood and hits by element
        return log_like, hits, self.row_lengths_close

    # *******************************************************************************************************
    @tf.function        
    def call(self, 
             u_pred: tf.Tensor, 
             num_hits: tf.Tensor, 
             R: tf.Tensor) \
             -> Tuple[tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor]:
        """
        Score candidate trajectories in current batch based on how well they match observations
        INPUTS:
            u_pred:     predicted directions with current batch of elements
            num_hits:   number of hits in mixture model
            R:          resolution parameter (Cartesian distance)
        """
        # Alias the input R to indicate that it is by element; sigma_mag passed for all rows
        R_elt: tf.Tensor = R

        # Get the threshold by element
        thresh_s2_elt = self.get_thresh_s2()

        # Calculate distance squared s2 and relative distance squared v
        s2, v, is_close = self.calc_dist(u_pred=u_pred, thresh_s2_elt=thresh_s2_elt)

        # Conditional probability q in distance model, conditional on a hit or a miss
        q_hit: tf.Tensor
        q_miss: tf.Tensor
        q_hit, q_miss = \
                self.cond_pdf_dist(v=v, R_elt=R_elt, thresh_s2_elt=thresh_s2_elt)

        # Compute the mixture model
        log_like, p_hit_post = self.mixture_model(q_hit=q_hit, q_miss=q_miss, num_hits=num_hits)

        # Count hits
        hits = self.calc_hits(s2=s2, p_hit_post=p_hit_post)

        # Return the log likelihood and hits by element
        return log_like, hits, self.row_lengths_close

# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
