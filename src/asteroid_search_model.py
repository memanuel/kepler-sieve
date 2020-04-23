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

# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# Utility
import os
import time
from datetime import timedelta

# MSE imports
from asteroid_model import AsteroidDirection, AsteroidMagnitude, make_model_ast_pos
from asteroid_search_layers import CandidateElements, MixtureParameters, TrajectoryScore, R_min_
from asteroid_search_optimizers import make_opt
from candidate_element import elts_np2df
from asteroid_integrate import calc_ast_pos
from candidate_element import perturb_elts
from ztf_element import ztf_elt_hash, make_ztf_batch
from asteroid_search_report import traj_diff
from nearest_asteroid import load_known_ast_pos, nearest_ast_elt_cart, nearest_ast_elt_cov, elt_q_norm
from asteroid_dataframe import calc_ast_data, spline_ast_vec_df
from ztf_ast import load_ztf_nearest_ast
from astro_utils import deg2dist, dist2deg, dist2sec
from tf_utils import tf_quiet, Identity
from tf_astro_utils import tf_dist2deg, tf_dist2sec
from utils import print_header

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

# Ignore irrelevant numpy errors
np.seterr(all='ignore')

# ********************************************************************************************************************* 
# Constants
space_dims: int = 3

# Data type
dtype: tf.dtypes.DType = tf.float32
dtype_np: type = np.float32

# Save directory for candidate elements
save_dir: str = '../data/candidate_elt'

# Load the ZTF asteroid data
ztf_ast = load_ztf_nearest_ast()

# Set plot style variables
mpl.rcParams['figure.figsize'] = [16.0, 10.0]
mpl.rcParams['font.size'] = 16

# Colors for plots
color_mean: str = 'blue'
color_lo: str = 'orange'
color_hi: str = 'green'
color_min: str = 'red'
color_max: str = 'purple'

# ********************************************************************************************************************* 
def candidate_elt_hash(elts: pd.DataFrame, thresh_deg: float) -> int:
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements; and 2 mixture params num_hits, R
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
    OUTPUTS:
        hash_id:    Unique ID for these inputs
    """
    # Columns of the Dataframe to hash
    cols_to_hash: List[str] = ['element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch', 'num_hits', 'R']
    # Tuple of int64; one per orbital element candidate
    hash_df: Tuple = tuple((pd.util.hash_pandas_object(elts[cols_to_hash])).values)
    # Combine the element hash tuple with the threshold
    thresh_deg_int: int = int(thresh_deg*2**48)
    hash_id: int = abs(hash(hash_df + (thresh_deg_int, )))

    return hash_id

# ********************************************************************************************************************* 
# Custom model for Asteroid Search
# ********************************************************************************************************************* 

class AsteroidSearchModel(tf.keras.Model):
    """
    Custom keras model that searches for orbital elements and mixture model parameters
    consistent with a batch of observations generated from the initially guessed candidate elements.
    """

    def __init__(self, 
                 elts: pd.DataFrame, 
                 ztf_elt: pd.DataFrame, 
                 site_name: str='geocenter', 
                 thresh_deg: float = 1.0,
                 optimizer_type: str = 'adam',
                 learning_rate: float = 2.0**-12, 
                 clipnorm: float = 1.0,
                 file_name: str = None,
                 **kwargs):
        """
        INPUTS:
            elts:           DataFrame with initial guess for orbital elements. Columns must include:
                            orbital elements: element_id, a, e, inc, Omega, omega, f, epoch
                            mixture parameters initial guess: num_hits, R, thresh_s
                            magnitude initial guess: H
                            Output of asteroid_elts, perturb_elts or random_elts
            ztf_elt:        DataFrame with ZTF observations within thresh_deg degrees of
                            of the orbits predicted by these elements.
                            Columns must include:
                            ztf_id, element_id, mjd, ux, uy, uz, mag_app
                            Output of make_ztf_batch or load_ztf_batch,
            site_name:      Used for topos adjustment, e.g. 'geocenter' or 'palomar'
            thresh_deg:     Threshold that was used for filtering the ZTF observations.
            optimizer_type: One of 'adam', 'rmsprop', 'adadelta'
            learning_rate:  Initial value of learning rate, default to 2^-15.
            clipnorm:       Initial value of clipnorm for gradient clipping, defaults to 1.0.
        """
        # Initialize tf.keras.Model
        super(AsteroidSearchModel, self).__init__(**kwargs)
        
        # *****************************************************************************************
        # Description of input data 
        # *****************************************************************************************

        # Batch size comes from elts
        self.batch_size: int = elts.shape[0]

        # Shape of the observed trajectories
        self.data_size: int = ztf_elt.shape[0]
        self.traj_shape: Tuple[int, int] = (self.data_size, space_dims)

        # The element_id for the elements in this batch
        self.element_id: keras.backend.constant = keras.backend.constant(value=elts.element_id, dtype=tf.int32)

        # Numpy array and tensor of observation times; flat, shape (data_size,)
        ts_np: np.ndarray = ztf_elt.mjd.values.astype(dtype_np)
        self.ts: keras.backend.constant = keras.backend.constant(value=ts_np, shape=(self.data_size,), dtype=dtype)

        # Get observation count per element
        row_lengths_np: np.ndarray = ztf_elt.element_id.groupby(ztf_elt.element_id).count()
        self.row_lengths: keras.backend.constant = \
            keras.backend.constant(value=row_lengths_np, shape=(self.batch_size,), dtype=tf.int32)
        self.row_lengths_float: keras.backend.constant = \
            keras.backend.constant(value=row_lengths_np, shape=(self.batch_size,), dtype=dtype)

        # Threshold - original, used for data
        self.thresh_deg_data: float = thresh_deg
        # Threshold - live, used for scoring; array of shape [batch_size,]
        self.thresh_deg: np.ndarray = np.full(shape=self.batch_size, fill_value=thresh_deg, dtype=dtype_np)

        # Save the original ztf_elt frame for reuse
        self.ztf_elt: pd.DataFrame = ztf_elt

        # Observed directions; extract from ztf_elt DataFrame
        cols_u_obs: List[str] = ['ux', 'uy', 'uz']
        u_obs_np: np.ndarray = ztf_elt[cols_u_obs].values.astype(dtype_np)

        # Apparent magnitude; extract from ztf_elt DataFrame
        mag_app_np = ztf_elt.mag_app.values.astype(dtype_np)

        # *****************************************************************************************
        # Layers for candidate elements, asteroid direction and score
        # *****************************************************************************************

        # Set of trainable weights with candidate orbital elements; initialize according to elts
        self.candidate_elements: CandidateElements = CandidateElements(elts=elts, name='candidate_elements')

        # Set of trainable weights with candidate mixture parameters
        self.mixture_parameters: MixtureParameters = MixtureParameters(elts=elts, name='mixture_parameters')

        # The predicted direction; shape is [data_size, 3,]
        self.direction: AsteroidDirection = \
            AsteroidDirection(ts_np=ts_np, row_lengths_np=row_lengths_np, site_name=site_name, name='direction')

        # The predicted magnitude; shape is [data_size, ]
        self.magnitude: AsteroidMagnitude = \
               AsteroidMagnitude(ts_np=ts_np, row_lengths_np=row_lengths_np, elts=elts, name='magnitude')

        # Bind the direction layer to this model for legibility
        self.position = self.direction.position

        # Calibration arrays (flat)
        self.cols_q_ast: List[str] = ['qx', 'qy', 'qz']
        self.cols_v_ast: List[str] = ['vx', 'vy', 'vz']
        self.q_ast_np: np.ndarray = ztf_elt[self.cols_q_ast].values.astype(dtype_np)
        self.v_ast_np: np.ndarray = ztf_elt[self.cols_v_ast].values.astype(dtype_np)

        # Run calibration
        self.position.calibrate(elts=elts, q_ast=self.q_ast_np, v_ast=self.v_ast_np)

        # Score layer for these observations
        self.score: TrajectoryScore = \
            TrajectoryScore(row_lengths_np=row_lengths_np, 
                            u_obs_np=u_obs_np, mag_app_np=mag_app_np,
                            thresh_deg=self.thresh_deg, name='score')

        # *****************************************************************************************
        # Variables for adaptive training and training history
        # *****************************************************************************************

        # Save the learning rate on the model object to facilitate adaptive training
        self.learning_rate: float = learning_rate
        self.clipnorm: float = clipnorm

        # Build the selected optimizer with the input learning rate
        self.optimizer_type: str = optimizer_type.lower()
        self.optimizer: keras.optimizers.Optimizer = \
            make_opt(optimizer_type=self.optimizer_type, 
                    learning_rate=self.learning_rate, 
                    clipnorm=self.clipnorm, clipvalue=None)
        # Compile the model with this optimizer instance
        self.recompile()

        # Set the learning rate factors for adaptive training
        self.lr_factor_elt_dn: float = 0.5  # elementwise learning rate

        # Initialize loss history and total training time
        self.training_time: float = 0.0

        # Epoch and episode counters
        self.batches_per_epoch: int = 64
        self.epochs_per_episode: int = 4
        self.samples_per_epoch: int = self.batches_per_epoch * self.batch_size
        self.episode_length: int = 0
        self.current_episode: int = 0
        self.current_epoch: int = 0
        self.current_batch: int = 0

        # Tensor of ones to pass as dummy inputs for evaluating one batch or training one epoch
        self.x_eval: tf.Tensor = tf.ones(shape=self.batch_size, dtype=dtype)
        self.x_trn: tf.Tensor = tf.ones(self.samples_per_epoch, dtype=dtype)

        # Start training timer
        self.t0: float = time.time()
        self.episode_t0: float = time.time()

        # Cached values of log likelihood and loss
        self.log_like_np: np.ndarray = np.zeros(self.batch_size, dtype=dtype_np)
        self.hits_np: np.ndarray = np.zeros(self.batch_size, dtype=dtype_np)
        self.R_sec_np: np.ndarray = np.zeros(self.batch_size, dtype=dtype_np)
        self.thresh_sec_np: np.ndarray = np.zeros(self.batch_size, dtype=dtype_np)
        self.log_like_mean: float = 0.0
        self.hits_mean: float = 0.0
        self.loss: float = 0.0

        # Training mode: one of 'joint', 'element', 'mixture'
        self.training_mode: str = 'joint'

        # Weight factors for training in three modes: joint, element, mixture
        self.weight_joint: tf.Variable = \
                tf.Variable(initial_value=np.ones(self.batch_size, dtype=dtype_np), trainable=False, dtype=dtype)
        self.weight_element: tf.Variable = \
            tf.Variable(initial_value=np.ones(self.batch_size, dtype=dtype_np), trainable=False, dtype=dtype)
        self.weight_mixture: tf.Variable = \
            tf.Variable(initial_value=np.ones(self.batch_size, dtype=dtype_np), trainable=False, dtype=dtype)
        # Mean of active weights; for effective learning rate
        self.active_weight_mean: float = 1.0
        self.effective_learning_rate: float = self.learning_rate * self.active_weight_mean

        # Power of resolution in the denominator of objective function
        self.obj_den_R_power: tf.Variable = \
            tf.Variable(initial_value=0.0, trainable=False, dtype=dtype, name='obj_den_R_power')

        # Power of threshold in the denominator of objective function
        self.obj_den_thresh_power: tf.Variable = \
            tf.Variable(initial_value=0.0, trainable=False, dtype=dtype, name='obj_den_thresh_power')
        # Minimum value of resolution R (Cartesian)
        # self.R_min: keras.backend.constant = \
        #         keras.backend.constant(value=R_min_, dtype=dtype, name='R_min')

        # Early stopping callback
        self.update_early_stop()
        self.callbacks = [self.early_stop]

        # Initialize lists with training history
        # Weights, elements and mixture parameters; for rolling back bad updates
        self.candidate_elements_hist: List[tf.Tensor] = []
        self.mixture_parameters_hist: List[tf.Tensor] = []
        self.magnitude_hist: List[tf.Tensor] = []
        self.score_hist: List[tf.Tensor] = []
      
        # Log likelihoods and losses
        self.log_like_hist: List[np.ndarray] = []
        self.hits_hist: List[np.ndarray] = []
        self.log_like_mean_hist: List[np.ndarray] = []
        self.hits_mean_hist: List[np.ndarray] = []
        self.loss_hist: List[np.ndarray] = []

        # Training counters and times
        self.episode_hist: List[int] = []
        self.epoch_hist: List[int] = []
        self.batch_hist: List[int] = []
        self.training_time_hist: List[float] = []

        # DataFrame with training history
        self.train_hist_summary: pd.DataFrame = pd.DataFrame()
        self.train_hist_elt: pd.DataFrame = pd.DataFrame()

        # Add the first entries with initial weights and losses
        self.save_weights()

        # Threshold for and count of the number of bad training episodes
        self.bad_episode_thresh: float = 0.00
        self.bad_episode_count: int = 0

        # Current fit
        self.elts_fit: pd.DataFrame = self.candidates_df()        

        # Elements of the nearest asteroid
        self.elts_near_ast: pd.DataFrame = pd.DataFrame()

        # Save hash IDs for elts and ztf_elts
        self.elts_hash_id: int = candidate_elt_hash(elts=elts, thresh_deg=self.thresh_deg_data)
        self.ztf_elt_hash_id: int = \
            ztf_elt_hash(elts=elts, ztf=self.ztf_elt, thresh_deg=self.thresh_deg_data, near_ast=False)
        
        # File name for saving training progress
        self.file_name: str = f'candidate_elt_{self.elts_hash_id}.h5' if file_name is None else file_name
        self.file_path: str = os.path.join(save_dir, self.file_name)

        # Initial tuning of magnitude
        self.calc()
        self.tune_mag()

        # Initialize model with frozen elements
        self.freeze_candidate_elements()

        # Fit one time to "prime the pump" on the loss function.
        # Not sure why this seems to be necessary but it is.
        self.fit(self.x_eval, verbose=0)

    @tf.function
    def calc(self, inputs=None):
        """
        Predict directions from current elements and score them.  
        inputs is a dummmy with no affect on the results; it is there to satisfy keras.Model API
        """

        # Extract the candidate elements and mixture parameters; pass dummy inputs to satisfy keras Layer API
        a: tf.Tensor
        e: tf.Tensor
        inc: tf.Tensor
        Omega: tf.Tensor
        omega: tf.Tensor
        f: tf.Tensor
        epoch: tf.Tensor
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=inputs)
        
        # Extract mixture parameters; pass dummy inputs to satisfy keras Layer API
        self.num_hits: tf.Tensor
        self.R: tf.Tensor
        self.num_hits, self.R, self.R_max = self.mixture_parameters(inputs=inputs)
        
        # Stack the current orbital elements.  Shape is [batch_size, 7,]
        orbital_elements: tf.Tensor = keras.backend.stack([a, e, inc, Omega, omega, f, epoch,])

        # Stack mixture model parameters. Shape is [batch_size, 3,]
        mixture_params: tf.Tensor = keras.backend.stack([self.num_hits, self.R, self.R_max])

        # Tensor of predicted directions.  Shape of u_pred is [data_size, 3,]
        # delta_pred is distance from earth to ast.  
        # q_ast also included for use in magnitude calculation. Shape [data_size, 3]
        self.u_pred: tf.Tensor
        self.delta_pred: tf.Tensor
        self.q_ast: tf.Tensor
        self.u_pred, self.delta_pred, self.q_ast, = self.direction(a, e, inc, Omega, omega, f, epoch)

        # Compute the predicted magnitude; shape [data_size,]
        # self.mag_pred, self.sigma_mag = self.magnitude(self.q_ast)
        
        # Compute the log likelihood by element from the predicted direction and mixture model parameters
        # Shape is [elt_batch_size, 3]
        self.log_like: tf.Tensor
        self.hits: tf.Tensor
        self.row_lengths_close: tf.Tensor
        # self.log_like, self.hits, self.row_lengths_close = \
        #     self.score(self.u_pred, mag_pred=self.mag_pred, 
        #                num_hits=self.num_hits, R=self.R, sigma_mag=self.sigma_mag)
        self.log_like, self.hits, self.row_lengths_close = \
            self.score(self.u_pred, num_hits=self.num_hits, R=self.R)
        
        # Add the loss function - the NEGATIVE of the log likelihood weighted by each element's weight
        # Take negative b/c TensorFlow minimizes the loss function, and we want to maximize the log likelihood
        self.elt_weight: tf.Variable = self.get_active_weight()
        # Numerator of objective function is log_like by element, weighted by elt_weight
        obj_num: tf.Tensor = tf.multiply(self.elt_weight, self.log_like)
        
        # Divide by threshold to encourage it getting smaller; objective function to maximize
        thresh_s2: tf.Tensor = self.score.get_thresh_s2()
        thresh_s: tf.Tensor = tf.sqrt(thresh_s2)
        thresh_log1p: tf.Tensor = tf.math.log1p(thresh_s)
        obj_den_1: tf.Tensor = tf.math.pow(x=thresh_log1p, y=self.obj_den_thresh_power, name='obj_den_1')       
        
        # Divide by resolution R to encourage it getting smaller
        # The minimum resolution of 1.0 arc second = 7.77E-7 ~ 2^-20.3
        R_log1p: tf.Tensor = tf.add(tf.math.log1p(self.R), 2.0**-24)
        obj_den_2: tf.Tensor = tf.math.pow(x=R_log1p, y=self.obj_den_R_power, name='obj_den_2')
        
        # Denominator of objective function includes terms for thresh_sec^2 * R_sec^1
        obj_den: tf.Tensor = tf.multiply(obj_den_1, obj_den_2, name='obj_den')
        # obj_den = tf.multiply(obj_den, obj_den_3)

        # Objective function is the quotient of log_like over powers of R and thresh
        obj: tf.Tensor = tf.divide(obj_num, obj_den, name='obj')

        # Take the negative of the objective function as the loss
        loss: tf.Tensor = tf.negative(obj)

        # Stack score outputs. Shape is [batch_size, 4,]
        row_lengths_close_fl = tf.cast(self.row_lengths_close, dtype=dtype)
        score_outputs: tf.Tensor = keras.backend.stack([self.log_like, self.hits, row_lengths_close_fl, loss])

        # Alias outputs and save onto the model
        self.score_outputs = Identity(name='score_outputs')(score_outputs)
        self.orbital_elements = Identity(name='orbital_elements')(orbital_elements)
        self.mixture_params = Identity(name='mixture_parameters')(mixture_params)

        # Wrap outputs
        outputs: Tuple[tf.Tensor, tf.Tensor, tf.Tensor]
        outputs = (self.score_outputs, self.orbital_elements, self.mixture_params)

        return outputs

    # @tf.function
    def call(self, inputs=None):
        """
        Predict directions from current elements and score them.  
        inputs is a dummmy with no affect on the results; it is there to satisfy keras.Model API
        """
        # Calculate the outputs
        score_outputs, orbital_elements, mixture_params = self.calc(inputs=inputs)

        # Extract the loss by candidate element
        loss: tf.Tensor = score_outputs[-1]

        # Filter out nan
        loss_is_nan = tf.math.is_nan(loss)
        # Apply placeholder of +5 (pretty bad score) over NaN until they are rolled back
        loss_clean = tf.where(condition=loss_is_nan, x=5.0, y=loss)

        # Total the clean closs over all the elements in the batch
        loss_sum: tf.Tensor
        loss_sum = tf.reduce_sum(loss_clean, name='loss_sum')

        # Add the blended loss to the model
        self.add_loss(loss_sum)
        
        return loss_sum

    def predict_position(self) -> Tuple[tf.Tensor, tf.Tensor]:
        """Predict position q and velocity v"""
        # Extract the candidate elements and mixture parameters
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=None)

        # Ragged tensor of predicted position and velocity.  Shape is [batch_size, num_obs, 3,]
        q_pred: tf.Tensor
        v_pred: tf.Tensor
        q_pred, v_pred = self.position(a, e, inc, Omega, omega, f, epoch)        

        return (q_pred, v_pred)

    def predict_direction(self) -> Tuple[tf.Tensor, tf.Tensor]:
        """Predict direction u, displacement r and magnitude"""
        # Extract the candidate elements and mixture parameters
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=None)

        # Ragged tensor of predicted directions.  Shape is [batch_size, num_obs, 3,]
        u_pred: tf.Tensor
        delta_pred: tf.Tensor
        u_pred, delta_pred, q_ast = self.direction(a, e, inc, Omega, omega, f, epoch)        

        # Compute the predicted magnitude; shape [data_size,]
        mag_pred, sigma_mag = self.magnitude(q_ast)

        return (u_pred, delta_pred, mag_pred)

    # *********************************************************************************************
    # Methods to calculate outputs, log likelihood, loss
    # *********************************************************************************************

    def calc_score(self) -> tf.Tensor:
        """Calculate a tuple of score outputs; no input required (uses cached x_eval)"""
        outputs = self.calc()
        score_outputs: tf.Tensor = outputs[0]
        return score_outputs
    
    def calc_log_like(self) -> Tuple[tf.Tensor, tf.Tensor]:
        """Calculate the log likelihood as tensor of shape [batch_size,]"""
        score_outputs: tf.Tensor = self.calc_score()
        log_like, hits, num_rows_close, loss = score_outputs
        log_like_mean: tf.Tensor = tf.reduce_mean(log_like).numpy()
        return log_like, log_like_mean

    def calc_loss(self) -> tf.Tensor:
        """Calculate loss function by element with current inputs; shape [batch_size]"""
        score_outputs: tf.Tensor = self.calc_score()
        log_like, hits, loss, num_rows_close = score_outputs
        return loss

    def calc_loss_total(self):
        """Calculate total loss function with current inputs; scalar"""
        loss = self.calc_loss()
        return tf.reduce_sum(loss)

    # def traj_err(self, elts0, elts1):
    #     """Calculate difference in trajectories from two sets of orbital elements"""
    #     return traj_diff(elts0, elts1, self.model_pos)

    # *********************************************************************************************
    def get_orbital_elts(self) \
        -> Tuple[tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor]:
        """Extract the current orbital elements as Numpy arrays"""
        # Extract the candidate elements
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=None)
        return a, e, inc, Omega, omega, f, epoch
  
    # *********************************************************************************************
    def get_mixture_params(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract the current mixture parameters as Numpy arrays"""
        num_hits: tf.Tensor
        R: tf.Tensor
        R_max: tf.Tensor
        num_hits, R, R_max = self.mixture_parameters(inputs=None)
        return num_hits.numpy(), R.numpy(), R_max.numpy()

    def get_thresh_deg(self) -> np.ndarray:
        """Extract the threshold in degrees as a Numpy array"""
        return self.score.get_thresh_deg()

    def get_H(self) -> np.ndarray:
        """Extract the brightness parameter H as a Numpy array"""
        return self.magnitude.get_H().numpy()

    def get_sigma_mag(self) -> np.ndarray:
        """Extract the standard deviation of the magnitude as a Numpy array"""
        return self.magnitude.get_sigma_mag().numpy()

    # *********************************************************************************************
    # Element weights; equivalent to independent learning rates for each element in the batch
    # *********************************************************************************************

    def get_active_weight(self) -> tf.Variable:
        """Get the currently active weight"""
        weight_tbl : Dict[str, tf.Tensor]= {
            'joint': self.weight_joint,
            'element': self.weight_element,
            'mixture': self.weight_mixture,
        }
        weight: tf.Tensor = weight_tbl[self.training_mode]
        return weight

    def set_active_weight(self, weight) -> None:
        """Set the currently active weight"""
        weight_active : Dict[str, tf.Tensor] = {
            'joint': self.weight_joint,
            'element': self.weight_element,
            'mixture': self.weight_mixture,
        }[self.training_mode]
        weight_active.assign(weight)
        # update mean active weight and effective learning rate
        self.active_weight_mean = tf.reduce_mean(weight).numpy()
        self.effective_learning_rate = self.active_weight_mean * self.learning_rate

    def reset_active_weight(self) -> None:
        """Reset the active weight to all 1s"""
        weight_ones: np.ndarray = np.ones(self.batch_size, dtype=dtype_np)
        self.set_active_weight(weight_ones)

    def adjust_active_weight(self, factor: float) -> None:
        """Adjust active weight by a factor"""
        active_weight_old: tf.Tensor = self.get_active_weight()
        active_weight_new: tf.Tensor = tf.multiply(active_weight_old, factor)
        self.set_active_weight(active_weight_new)

    # *********************************************************************************************
    # Methods to change model state in training
    # *********************************************************************************************

    def recompile(self):
        """Recompile this model with its current optimizer"""
        # Note: important not to name this method compile, that breaks relationship with tf.keras.Model
        # tf.keras.Model.compile(self, optimizer=self.optimizer)
        self.compile(optimizer=self.optimizer)

    # *********************************************************************************************
    def set_learning_rate(self, learning_rate: float):
        """Set the learning rate"""
        self.learning_rate = learning_rate
        # Recompile with the new learning rate
        self.optimizer = make_opt(optimizer_type=self.optimizer_type, learning_rate=self.learning_rate, 
                                  clipnorm=self.clipnorm, clipvalue=None)
        self.recompile()
        
    # *********************************************************************************************
    def adjust_learning_rate(self, lr_factor: float, verbose: bool=True):
        """Adjust the learning rate and recompile the model"""
        if verbose:
            print(f'Changing learning rate by factor {lr_factor:8.6f} from '
                  f'{self.learning_rate:8.3e} to {self.learning_rate*lr_factor:8.3e}.')
        learning_rate: float = self.learning_rate * lr_factor
        # Apply the new learning rate
        self.set_learning_rate(learning_rate)

    # *********************************************************************************************
    def set_clipnorm(self, clipnorm):
        """Set the clipnorm parameter for gradient clipping."""
        self.clipnorm = clipnorm
        # Delegate to set_learning_rate to rebuild the optimizer
        self.set_learning_rate(self.learing_rate)
        
    # *********************************************************************************************
    def set_R_deg_max(self, R_deg_max: float):
        """Adjust resolution R_deg to at most R_deg_max"""
        # If R_deg_max was a scalar, promote it to a full numpy array
        if isinstance(R_deg_max, float):
            R_deg_max = np.full(shape=self.batch_size, fill_value=R_deg_max, dtype=dtype_np)            
        # Apply the new R_max value to the mixture parameters layer
        self.mixture_parameters.set_R_deg_max(R_deg_max)

    # *********************************************************************************************
    def set_R_sec_max(self, R_sec_max: float):
        """Adjust R_deg to at most R_sec_max after converting seconds to degrees"""
        R_deg_max = R_sec_max / 3600.0
        self.set_R_deg_max(R_deg_max)

    # *********************************************************************************************
    def set_thresh_deg(self, thresh_deg: np.ndarray):
        """Adjust the threshold for which observations are included in the score function"""
        # If thresh_deg was a scalar, promote it to a full numpy array
        if isinstance(thresh_deg, float) or isinstance(thresh_deg, int):
            thresh_deg = np.full(shape=self.batch_size, fill_value=thresh_deg, dtype=dtype_np)
        # Update the thresh_deg member; this is an array
        self.thresh_deg = thresh_deg
        # Apply this update to the score layer and mixture parameters layer
        self.score.set_thresh_deg(self.thresh_deg)
        self.mixture_parameters.set_thresh_deg(self.thresh_deg)

    # *********************************************************************************************
    def set_thresh_deg_max(self, thresh_deg_max: float):
        """Adjust thresh_deg to at most thresh_deg_max"""
        # If thresh_deg_max was a scalar, promote it to a full numpy array
        if isinstance(thresh_deg_max, float) or isinstance(thresh_deg_max, int):
            thresh_deg_max = np.full(shape=self.batch_size, fill_value=thresh_deg_max, dtype=dtype_np)
        # Update the thresh_deg member; this is an array
        self.thresh_deg = np.minimum(self.thresh_deg, thresh_deg_max)
        # Apply the new max threshold to the score layer
        self.score.set_thresh_deg_max(thresh_deg_max)

    # *********************************************************************************************
    def set_thresh_sec_max(self, thresh_sec_max: float):
        """Adjust thresh_deg to at most thresh_sec_max after converting seconds to degrees"""
        thresh_deg_max = thresh_sec_max / 3600.0
        self.set_thresh_deg_max(thresh_deg_max)

    # *********************************************************************************************
    def set_obj_den_R_power(self, obj_den_R_power: float):
        """Set the power of R in the denominator of the objective function"""
        self.obj_den_R_power.assign(obj_den_R_power)

    # *********************************************************************************************
    def set_obj_den_thresh_power(self, obj_den_thresh_power: float):
        """Set the power of R in the denominator of the objective function"""
        self.obj_den_thresh_power.assign(obj_den_thresh_power)

    # *********************************************************************************************
    def recalibrate(self):
        """Recalibrate Kepler model rebound integration at current orbital elements"""
        # Element IDs as a Numpy array
        element_ids: np.ndarray = self.element_id.numpy()

        # Current candidates in DataFrame format for calc_ast_data
        cols_elt: List[str] = ['element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
        elts = self.candidates_df()[cols_elt]

        # Filter NaN out when training gets messed up; bad settings will be rolled back later
        cols_filter = ['Omega', 'omega', 'f']
        elts[cols_filter] = np.nan_to_num(elts[cols_filter].values, nan=0.0)

        # Get unique times and inverse indices
        mjd: np.ndarray
        unq_idx: np.ndarray
        mjd, unq_idx = np.unique(self.ts.numpy(), return_inverse=True)
        mjd0: float = np.min(mjd) - 1.0
        mjd1: float = np.max(mjd) + 1.0

        # Calculate positions in this date range, sampled daily
        df_ast_daily: pd.DataFrame
        df_earth_daily: pd.DataFrame
        df_sun_daily: pd.DataFrame
        df_ast_daily, df_earth_daily, df_sun_daily = \
            calc_ast_data(elts=elts, mjd0=mjd0, mjd1=mjd1, element_id=element_ids)

        # Spline positions at ztf times
        df_ast, df_earth, df_sun = \
            spline_ast_vec_df(df_ast=df_ast_daily, df_earth=df_earth_daily, df_sun=df_sun_daily, 
                              mjd=mjd, include_elts=False)

        # Iterate over the elements
        element_id: np.int32
        for element_id in element_ids:
            # Mask for observations of this element
            mask_ztf: pd.Series = (self.ztf_elt.element_id == element_id)
            # ztf_i: pd.DataFrame = self.ztf_elt[mask_ztf]
            # Mask for asteroid position of this element
            mask_df: pd.Series = (df_ast.element_id==element_id)
            df_ast_i: pd.DataFrame = df_ast[mask_df]
            # Index into df_ast_i on selected rows
            df_idx: np.ndarray = unq_idx[mask_ztf]
            # Overwrite relevant slices of q_ast and v_ast with the newly integrated q, v
            self.q_ast_np[mask_ztf] = df_ast_i[self.cols_q_ast].iloc[df_idx]
            self.v_ast_np[mask_ztf] = df_ast_i[self.cols_v_ast].iloc[df_idx]            

        # Run calibration with the new q, v
        self.position.calibrate(elts=elts, q_ast=self.q_ast_np, v_ast=self.v_ast_np)

    # *********************************************************************************************
    def tune_mag(self):
        """Tune the magnitude"""
        # Is the magnitude good on each element?
        mag_is_good: tf.Tensor = tf.cast(self.score.mag_is_good, tf.bool)

        # delta_H value to apply on good and bad elements, respectively
        delta_H_bad: tf.Tensor = self.score.delta_H
        delta_H_good: tf.Tensor = tf.zeros_like(delta_H_bad)

        # sigma_mag value to apply on good and bad elements, respectively
        sigma_mag_good: tf.Tensor = self.magnitude.get_sigma_mag()
        sigma_mag_bad: tf.Tensor = self.score.sigma_hat

        # Filter delta_H and sigma_mag between good and bad values
        delta_H: tf.Tensor = tf.where(mag_is_good, delta_H_good, delta_H_good)
        sigma_mag: tf.Tensor = tf.where(mag_is_good, sigma_mag_good, sigma_mag_bad)

        # Apply these updates to magnitude layer
        self.magnitude.shift_H(delta_H)
        self.magnitude.set_sigma_mag(sigma_mag)

    # *********************************************************************************************
    # Freeze and thaw layers relating to orbital elements and mixture parameters
    # *********************************************************************************************

    # *********************************************************************************************
    def freeze_candidate_elements(self):
        """
        Make the candidate orbital elements and magnitude not trainable (frozen);
        only update mixture parameters and score
        """
        # Quit early if already frozen
        if not self.candidate_elements.trainable:
            return
        self.candidate_elements.trainable = False
        self.magnitude.trainable = False
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'mixture'

    def thaw_candidate_elements(self):
        """Make the candidate orbital elements and magnitude trainable (unfrozen)"""
        # Quit early if already thawed
        if self.candidate_elements.trainable:
            return
        self.candidate_elements.trainable = True
        self.magnitude.trainable = True
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'joint' if self.mixture_parameters.trainable else 'element'

    # *********************************************************************************************
    def freeze_mixture_parameters(self):
        """
        Make the mixture parameters and score not trainable (frozen);
        only update candidate oribtal elements and magnitude
        """
        # Quit early if already frozen
        if not self.mixture_parameters.trainable:
            return
        self.mixture_parameters.trainable = False
        self.score.trainable = False
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'element'

    def thaw_mixture_parameters(self):
        """Make the mixture parameters trainable (unfrozen)"""
        # Quit early if already thawed
        if self.mixture_parameters.trainable:
            return
        self.mixture_parameters.trainable = True
        self.score.trainable = True
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'joint' if self.candidate_elements.trainable else 'mixture'

    # *********************************************************************************************
    def thaw_all(self):
        """Make all layers trainable (unfrozen)"""
        # Go through each layer and make it trainable
        self.candidate_elements.trainable = True
        self.mixture_parameters.trainable = True
        self.magnitude.trainable = True
        self.score.trainable = True
        # Recompile the model; update early stopping threshold, and training mode
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'joint'

    # *********************************************************************************************
    # Adaptive training; save weights and history at episode end
    # *********************************************************************************************

    def save_weights(self, is_update: bool=False):
        """Save the current weights, log likelihood and training history"""
        # Generate the outputs
        score_outputs: tf.Tensor
        orbital_elements: tf.Tensor
        mixture_params: tf.Tensor
        score_outputs, orbital_elements, mixture_params = self.calc()
        log_like: tf.Tensor
        hits: tf.Tensor
        # Unpack score_outputs
        log_like, hits, num_rows_close, loss = score_outputs
        # Unpack mixture params
        num_hits, R, R_max = mixture_params = mixture_params

        # The number of hits is at most row_lengths_close
        num_hits_max: tf.Tensor = tf.cast(x=num_rows_close, dtype=dtype)
        num_hits_adj: tf.Tensor = tf.math.minimum(x=num_hits, y=num_hits_max)
        self.mixture_parameters.set_num_hits(num_hits=num_hits_adj)

        # Extract R from mixture_params and convert to arc seconds
        R = mixture_params[1]
        R_sec = dist2sec(R)

        # Extract thresh_deg
        thresh_sec = self.get_thresh_deg() * 3600.0

        # is_update is true when we are updating weights that have already been written
        # this is passed when we need to restore weights that got worse
        is_new: bool = ~is_update
        
        # Calculate mean log likelihood, mean hits, and loss; save them to cache        
        self.log_like_np = log_like.numpy()
        self.hits_np = hits.numpy()
        self.R_sec_np = R_sec
        self.thresh_sec_np = thresh_sec
        self.log_like_mean = tf.reduce_mean(log_like).numpy()
        self.hits_mean = tf.reduce_mean(hits).numpy()
        self.loss = self.calc_loss_total()

        # Write history of weights, elements, and mixture parameters
        if is_new:
            self.candidate_elements_hist.append(self.candidate_elements.get_weights())
            self.mixture_parameters_hist.append(self.mixture_parameters.get_weights())
            self.magnitude_hist.append(self.magnitude.get_weights())
            self.score_hist.append(self.score.get_weights())
        else:
            self.candidate_elements_hist[-1] = self.candidate_elements.get_weights()
            self.mixture_parameters_hist[-1] = self.mixture_parameters.get_weights()
            self.magnitude_hist[-1] = self.magnitude.get_weights()
            self.score_hist[-1] = self.score.get_weights()
        
        # Write history of log likelihood and loss
        if is_new:
            self.log_like_hist.append(self.log_like_np)
            self.hits_hist.append(self.hits_np)
            self.log_like_mean_hist.append(self.log_like_mean)
            self.hits_mean_hist.append(self.hits_mean)
            self.loss_hist.append(self.loss)
        else:
            self.log_like_hist[-1] = self.log_like_np
            self.hits_hist[-1] = self.hits_np
            self.log_like_mean_hist[-1] = self.log_like_mean
            self.hits_mean_hist[-1] = self.hits_mean
            self.loss_hist[-1] = self.loss

        # Write history of training counters and time
        if is_new:
            self.episode_hist.append(self.current_episode)
            self.epoch_hist.append(self.current_epoch)
            self.batch_hist.append(self.current_batch)
            self.training_time_hist.append(self.training_time)

        # Delegate to save_train_hist to write history DataFrames
        self.save_train_hist(score_outputs=score_outputs, orbital_elements=orbital_elements, mixture_params=mixture_params)
        # Deduplicate if we did an update
        self.train_hist_deduplicate()

    # *********************************************************************************************
    def save_train_hist(self, score_outputs, orbital_elements, mixture_params):
        """Save training history, both summary and by element"""

        # Extract score outputs
        log_like = score_outputs[0]
        hits = score_outputs[1]
        num_rows_close = score_outputs[2]

        # Extract mixture parameters
        num_hits = mixture_params[0]
        R = mixture_params[1]
        R_max = mixture_params[2]
        # Transform to degrees or log
        R_deg = dist2deg(R)
        R_sec = dist2sec(R)
        log_R = np.log(R)
        R_deg_max = dist2deg(R_max)

        # Extract thresh_deg
        thresh_deg = self.get_thresh_deg()
        thresh_sec = thresh_deg * 3600.0
        thresh_s = deg2dist(thresh_deg)
        log_thresh = np.log(thresh_s)

        # Extract the brightness and standard deviation of the magnitude
        H = self.get_H()
        sigma_mag = self.get_sigma_mag()

        # Alias candidate elements layer and mixture parameters
        cand = self.candidate_elements
        mixt = self.mixture_parameters

        # DataFrame with detailed training history of this episode by element
        hist_elt_dict = {
            # Key (row number)
            'key': self.current_episode*self.batch_size + np.arange(self.batch_size), 
            # Element number and ID
            'element_num': np.arange(self.batch_size), 
            'element_id': self.element_id.numpy(),
            
            # Training stage: by episode, epoch, batch and time
            'episode': np.full(shape=self.batch_size, fill_value=self.current_episode),
            'epoch': np.full(shape=self.batch_size, fill_value=self.current_epoch),
            'batch': np.full(shape=self.batch_size, fill_value=self.current_batch),
            'training_time': np.full(shape=self.batch_size, fill_value=self.training_time),
            
            # Log likelihood and hits
            'log_like': log_like,
            'hits': hits,
            'num_rows_close': num_rows_close,

            # Orbital elements
            'a': orbital_elements[0].numpy(),
            'e': orbital_elements[1].numpy(),
            'inc': orbital_elements[2].numpy(),
            'Omega': orbital_elements[3].numpy(),
            'omega': orbital_elements[4].numpy(),
            'f': orbital_elements[5].numpy(),
            'epoch_ast': orbital_elements[6].numpy(),
            
            # Mixture parameters
            'num_hits': num_hits,
            'R': R,
            'R_deg': R_deg,
            'R_sec': R_sec,
            'log_R': log_R,
            'minus_log_R': -log_R,
            'R_max': R_max,
            'R_deg_max': R_deg_max,

            # Threshold
            'thresh_deg': thresh_deg,
            'thresh_sec': thresh_sec,
            'thresh_s': thresh_s,
            'minus_log_thresh_s': -np.log(thresh_s),

            # Brightness H and standard deviation of magnitude
            'H': H,
            'sigma_mag': sigma_mag,
            'minus_log_sigma_mag': -np.log(sigma_mag),

            # Control variables - candidate orbital elements
            'a_': cand.a_.numpy(),
            'e_': cand.e_.numpy(),
            'inc_': cand.inc_.numpy(),
            'Omega_': cand.Omega_.numpy(),
            'omega_': cand.omega_.numpy(),
            'f_': cand.f_.numpy(),
            # Control variables - mixture parameters
            'num_hits_': mixt.num_hits_.numpy(),
            'R_': mixt.R_.numpy(),            

            # Weights in three modes
            'weight_joint': self.weight_joint.numpy(),
            'weight_element': self.weight_element.numpy(),
            'weight_mixture': self.weight_mixture.numpy(),
        }
        train_hist_elt_cur = pd.DataFrame(hist_elt_dict, index=hist_elt_dict['key'])
        self.train_hist_elt = pd.concat([self.train_hist_elt, train_hist_elt_cur])

        # DataFrame with summary history for this episode
        hist_sum_dict = {
            # Key; same as the episode
            'key': [self.current_episode],
            # Training stage: by episode, epoch, batch and time
            'episode': [self.current_episode],
            'epoch': [self.current_epoch],
            'batch': [self.current_batch],
            'training_time': [self.training_time],

            # Loss and learning rate
            'loss': np.sum(self.loss),
            'learning_rate': self.learning_rate,

            # Log likelihood summarized over the elements
            'log_like_mean': [np.mean(log_like)],
            'log_like_med': [np.median(log_like)],
            'log_like_std': [np.std(log_like)],
            'log_like_min': [np.min(log_like)],
            'log_like_max': [np.max(log_like)],
            'log_like_q20': [np.quantile(log_like, q=0.2)],
            'log_like_q80': [np.quantile(log_like, q=0.8)],

            # Worst and best element in this batch
            'log_like_argmin': [np.argmin(log_like)],
            'log_like_argmax': [np.argmax(log_like)],

            # Hits summarized over the elements
            'hits_mean': [np.mean(hits)],
            'hits_med': [np.median(hits)],
            'hits_std': [np.std(hits)],
            'hits_min': [np.min(hits)],
            'hits_max': [np.max(hits)],
            'hits_q20': [np.quantile(hits, q=0.2)],
            'hits_q80': [np.quantile(hits, q=0.8)],

            # Summary of mixture parameter num_hits
            'num_hits_mean': [np.mean(num_hits)],
            'num_hits_std': [np.std(num_hits)],
            'num_hits_min': [np.min(num_hits)],
            'num_hits_max': [np.max(num_hits)],
            'num_hits_q20': [np.quantile(num_hits, q=0.2)],
            'num_hits_q80': [np.quantile(num_hits, q=0.8)],

            # Summary of mixture parameter R
            'R_sec_mean': [np.mean(R_sec)],
            'R_sec_std': [np.std(R_sec)],
            'R_sec_min': [np.min(R_sec)],
            'R_sec_max': [np.max(R_sec)],
            'R_sec_q20': [np.quantile(R_sec, q=0.2)],
            'R_sec_q80': [np.quantile(R_sec, q=0.8)],

            # Summary of mixture parameter log(R)
            'log_R_mean': [np.mean(log_R)],
            'log_R_std': [np.std(log_R)],
            'log_R_min': [np.min(log_R)],
            'log_R_max': [np.max(log_R)],
            'log_R_q20': [np.quantile(log_R, q=0.2)],
            'log_R_q80': [np.quantile(log_R, q=0.8)],

            # Summary of threshold_sec
            'thresh_sec_mean' : [np.mean(thresh_sec)],
            'thresh_sec_std' : [np.std(thresh_sec)],
            'thresh_sec_min' : [np.min(thresh_sec)],
            'thresh_sec_max' : [np.max(thresh_sec)],
            'thresh_sec_q20': [np.quantile(thresh_sec, q=0.2)],
            'thresh_sec_q80': [np.quantile(thresh_sec, q=0.8)],

            # Summary of log(threshold)
            'log_thresh_mean': [np.mean(log_thresh)],
            'log_thresh_std' : [np.std(log_thresh)],
            'log_thresh_min' : [np.min(log_thresh)],
            'log_thresh_max' : [np.max(log_thresh)],
            'log_thresh_q20': [np.quantile(log_thresh, q=0.2)],
            'log_thresh_q80': [np.quantile(log_thresh, q=0.8)],

            # Extra data required to rebuild model state
            'candidate_elements_trainable': self.candidate_elements.trainable,
            'mixture_parameters_trainable': self.mixture_parameters.trainable,

        }
        train_hist_summary_cur = pd.DataFrame(hist_sum_dict, index=hist_sum_dict['key'])
        self.train_hist_summary = pd.concat([self.train_hist_summary, train_hist_summary_cur])

    # *********************************************************************************************
    def train_hist_deduplicate(self):
        """Drop duplicate entries from training history"""
        # for hist in [self.train_hist_elt, self.train_hist_summary]:
        self.train_hist_elt = self.train_hist_elt.loc[~self.train_hist_elt.index.duplicated(keep='last')]
        self.train_hist_summary = self.train_hist_summary.loc[~self.train_hist_summary.index.duplicated(keep='last')]

    # *********************************************************************************************
    # Adaptive training; update weights and modify model state at episode end
    # *********************************************************************************************

    def update_weights(self, verbose: bool=True):
        """"Restore the weights for each element to the prior iteration if they got worse"""        
        # If we don't have at least 2 entries on the history lists, terminate early
        n: int = len(self.log_like_hist)
        if n < 2:
            return

        # Tune the magnitude based on the most recent log likelihood
        # self.tune_mag()

        # Get old and new log_like, hits, and loss; to figure out if an element improved
        log_like_old, log_like_new = self.log_like_hist[n-2:n]
        hits_old, hits_new = self.hits_hist[n-2:n]
        loss_old, loss_new = self.loss_hist[n-2:n]
        
        # Get old and new weights; to revert an element if necessary
        candidate_elements_old, candidate_elements_new = self.candidate_elements_hist[n-2:n]
        mixture_parameters_old, mixture_parameters_new = self.mixture_parameters_hist[n-2:n]
        magnitude_old, magnitude_new = self.magnitude_hist[n-2:n]
        score_old, score_new = self.score_hist[n-2:n]

        # Test which elements have gotten worse (usually this should be false)
        is_less_likely = tf.math.less(log_like_new, log_like_old)
        is_less_hits = tf.math.less(hits_new, hits_old)
        is_higher_loss = tf.math.greater(loss_new, loss_old)
        is_nan = tf.math.is_nan(log_like_new)

        # Accumulate different ways an update can be worse than the preceding one
        is_worse = tf.math.logical_or(is_less_likely, is_less_hits)
        is_worse = tf.math.logical_or(is_less_likely, is_higher_loss)
        is_worse = tf.math.logical_or(is_worse, is_nan)        

        # If none of the elements have gotten worse, terminate early
        if not tf.math.reduce_any(is_worse):
            return

        # Calculate the best weights on the trainable layers. Then restore them.
        if self.candidate_elements.trainable:
            candidate_elements_best = tf.where(condition=is_worse, x=candidate_elements_old, y=candidate_elements_new)
            self.candidate_elements.set_weights(candidate_elements_best)
        if self.mixture_parameters.trainable:
            mixture_parameters_best = tf.where(condition=is_worse, x=mixture_parameters_old, y=mixture_parameters_new)
            self.mixture_parameters.set_weights(mixture_parameters_best)
        if self.magnitude.trainable:
            magnitude_best = tf.where(condition=is_worse, x=magnitude_old, y=magnitude_new)
            self.magnitude.set_weights(magnitude_best)
        if self.score.trainable:
            score_best = tf.where(condition=is_worse, x=score_old, y=score_new)
            self.score.set_weights(score_best)

        # If we edited weights in this iteration, we also need to apply the changes to the history tables
        self.save_weights(is_update=True)

        # Reduce the weights on the elements that trained too fast
        weight_old = self.get_active_weight()
        weight_new = tf.where(condition=is_worse, x=weight_old*self.lr_factor_elt_dn, y=weight_old)
        # Median of the proposed new weights before normalization
        # weight_med = tf.reduce_min(tf.nn.top_k(weight_new, self.batch_size//2, sorted=False).values)
        
        # Apply the new active weight 
        self.set_active_weight(weight_new)
        if verbose:
            num_changed = tf.reduce_sum(tf.cast(x=is_worse, dtype=tf.int32)).numpy()
            print(f'Adjusted element weight down on {num_changed} candidate elements. Mean weight = {self.active_weight_mean:6.2e}')
        
        # Change in the mean log likelihood after editing history
        log_like_mean_change = self.log_like_mean_hist[-1] - self.log_like_mean_hist[-2]

        # If the change in the mean log likelihood is below a threhold, increment the bad_episode counter
        if log_like_mean_change < self.bad_episode_thresh:
            self.bad_episode_count += 1
            print(f'Increasing bad_episode_count to {self.bad_episode_count}.')
       
    # *********************************************************************************************
    def update_early_stop(self):
        """Update early stopping monitor"""
        baseline = self.calc_loss_total()
        self.early_stop = tf.keras.callbacks.EarlyStopping(
            monitor='loss', 
            patience=0,
            baseline=baseline,
            min_delta=0.0, 
            restore_best_weights=False)
    
    # *********************************************************************************************
    def episode_end(self, hist, verbose: int = 1):
        """Post-processing after one episode of adaptive training"""
        # Update training counters
        self.episode_length = len(hist.epoch)
        self.current_episode += 1
        self.current_epoch += self.episode_length
        self.current_batch += self.episode_length * self.batches_per_epoch

        # Update timers
        self.episode_time = (time.time() - self.episode_t0)
        self.training_time += self.episode_time
        self.episode_t0 = time.time()

        # Save weights and apply best to each element
        self.save_weights()
        self.update_weights()

        # Count number of good elements (with at least 10 hits)
        hits = np.round(self.hits_np)
        is_good_elt = (hits >= 10.0)
        num_good_elts = np.sum(is_good_elt)

        # Geometric mean of resolution and threshold
        R_sec_mean = np.exp(np.mean(np.log(self.R_sec_np)))
        R_sec_good = np.exp(np.mean(np.log(self.R_sec_np[is_good_elt])))
        R_sec_bad = np.exp(np.mean(np.log(self.R_sec_np[~is_good_elt])))
        thr_sec_mean = np.exp(np.mean(np.log(self.thresh_sec_np)))
        thr_sec_good = np.exp(np.mean(np.log(self.thresh_sec_np[is_good_elt])))
        thr_sec_bad = np.exp(np.mean(np.log(self.thresh_sec_np[~is_good_elt])))

        # Mean log_like on all and good elements
        log_like_mean = np.mean(self.log_like_np)
        log_like_good = np.mean(self.log_like_np[is_good_elt])
        log_like_bad = np.mean(self.log_like_np[~is_good_elt])
        hits_mean = np.mean(self.hits_np)
        hits_good = np.mean(self.hits_np[is_good_elt])
        hits_bad = np.mean(self.hits_np[~is_good_elt])

        # Update early_stop
        self.update_early_stop()

        # Recalibrate the model if the orbital elements are not frozen
        if self.candidate_elements.trainable:
            self.recalibrate()
            # Also need to save weights; loss is going to get worse immediately after a recalibration
            self.save_weights()

        # Status message
        if verbose > 0:
            print(f'                    \  All Elts : Bad Elts : Good Elts ({int(num_good_elts)})')
            print(f'Geom Mean Resolution:  {R_sec_mean:8.2f} : {R_sec_bad:8.2f} : {R_sec_good:8.2f} arc seconds')
            print(f'Geom Mean Threshold :  {thr_sec_mean:8.2f} : {thr_sec_bad:8.2f} : {thr_sec_good:8.2f} arc seconds')
            print(f'Mean Log Likelihood :  {log_like_mean:8.2f} : {log_like_bad:8.2f} : {log_like_good:8.2f}')
            print(f'Mean Hits           :  {hits_mean:8.2f} : {hits_bad:8.2f} : {hits_good:8.2f}')
            print(f'Good Elements       :  {num_good_elts:8.2f}')

    # *********************************************************************************************
    # Main adaptive search routine; called by external consumers
    # *********************************************************************************************

    # *********************************************************************************************
    def train_one_episode(self, verbose=1):
        """One episode of training"""
        # Status update for this episode
        print(f'\nTraining episode {self.current_episode}: Epoch {self.current_epoch:4}, Batch {self.current_batch:6}')
        print(f'effective_learning_rate={self.effective_learning_rate:8.3e}, training_time {self.training_time:0.0f} sec.')

        # Scheduled adjustment to thresh_deg if applicable
        if self.thresh_deg_schedule:
            episode_batch = self.current_batch - self.episode_batch_start
            episode_t = episode_batch / float(self.episode_max_batches)
            thresh_factor = np.exp(episode_t * self.thresh_deg_log_range)
            thresh_deg_max = thresh_factor * self.thresh_deg_start
            thresh_deg_mean = np.mean(thresh_deg_max)
            thresh_factor_mean = np.mean(thresh_factor)
            print(f'Updating thresh_deg. episode_batch {episode_batch} / {self.episode_max_batches} ({episode_t*100:4.1f}%); '
                  f'thresh_deg_max = {thresh_deg_mean:8.6f} (factor {thresh_factor_mean:8.6f} from start).')
            # Apply the new thresh_deg_max to the model
            self.set_thresh_deg_max(thresh_deg_max)
            # Save weights to account for jump in loss function after resetting thresh_deg
            self.save_weights()

        # Train for another episode
        hist = self.fit(x=self.x_trn, batch_size=self.batch_size, epochs=self.current_epoch+self.epochs_per_episode, 
                        steps_per_epoch=self.batches_per_epoch, initial_epoch=self.current_epoch,
                        callbacks=self.callbacks, shuffle=False, verbose=verbose)

        # Episode end processing; includes LR adjustment
        self.episode_end(hist)

    # *********************************************************************************************
    def search_adaptive(self, 
                        max_batches: int = 1024, 
                        batches_per_epoch: int = 64, 
                        epochs_per_episode: int = 4,
                        thresh_deg_start: float = None,
                        thresh_deg_end: float = None,
                        max_bad_episodes: int = 3,
                        learning_rate: Optional[float] = None,
                        min_learning_rate: Optional[float] = None,
                        reset_active_weight: bool = False,
                        save_at_end: bool = False,
                        verbose: int = 1):
        """
        Run asteroid search model adaptively.  
        Start with a high learning rate, gradually reduce it if early stopping triggered.
        INPUTS:
            max_batches: The maximum number of batches processed in the adaptive training.
            batches_per_epoch: The number of batches processed in one epoch of training
            epochs_per_episode: The number of epochs that comprise one episode of training
            max_bad_episodes: The number of episodes with poor progress before terminating training early
            learning_rate:   If specified, overwrite the prior learning_rate with this one
            min_learning_rate: Minimum for the learning rate; terminate early if LR drops below this
            verbose: Integer controlling verbosity level (passed on to tf.keras.model.fit)
        """
        # Save epochs_per_episode
        self.epochs_per_episode = epochs_per_episode

        # Apply learning_rate input if it was specified
        original_learning_rate = self.learning_rate
        if learning_rate is not None and learning_rate != self.learning_rate:
            log2_lr = np.log2(learning_rate)
            print(f'Applying learning_rate {learning_rate:6.2e} (2.0^{log2_lr:5.1f}) for adaptive training.')
            self.set_learning_rate(learning_rate)

        # Reset the active weights to 1 if requested
        if reset_active_weight:
            self.reset_active_weight()
        
        # If min_learning_rate was not set, default it to ratio of the current learning rate
        min_learning_rate = self.learning_rate / 128.0 if min_learning_rate is None else min_learning_rate

        # Update batches_per_epoch
        self.batches_per_epoch = batches_per_epoch

        # Early stopping callback
        self.update_early_stop()

        # Define one epoch as a number of batches        
        self.samples_per_epoch = self.batches_per_epoch * self.batch_size
        # Maximum number of episodes is twice the expected number if all episodes are full
        max_episodes = self.current_episode + (max_batches - self.current_batch)*2 // (batches_per_epoch * epochs_per_episode)

        # Reset the bad episode counter
        self.bad_episode_count = 0

        # Schedule thresh_deg if specified
        if (thresh_deg_start is not None) and (thresh_deg_end is not None):
            self.thresh_deg_schedule = True
            self.episode_batch_start = self.current_batch
            self.episode_max_batches = max_batches - self.episode_batch_start
            self.thresh_deg_start = thresh_deg_start
            self.thresh_deg_log_range = np.log(thresh_deg_end / self.thresh_deg_start)
        else:
            self.thresh_deg_schedule = False

        # Dummy training inputs
        self.x_trn = tf.ones(self.samples_per_epoch, dtype=dtype)

        # Save the weights; don't want to reset back to old weights b/c e.g. thresh_deg changed
        self.save_weights()

        # Continue training until max_epochs or max_episodes have elapsed, or learning_rate has dropped too low
        while (self.current_batch < max_batches) and \
              (self.current_episode < max_episodes) and \
              (self.bad_episode_count < max_bad_episodes) and \
              (min_learning_rate < self.effective_learning_rate):
              self.train_one_episode(verbose=verbose)

        # Report cause of training termination
        if self.current_batch >= max_batches:
            print_header(f'Terminating: Completed {self.current_batch} batches.')
        elif self.current_episode >= max_episodes:
            print_header(f'Terminating: Completed {self.current_episode} episodes.')
        elif self.bad_episode_count >= max_bad_episodes:
            print_header(f'Terminating: Had {self.bad_episode_count} bad episodes.')
        elif self.effective_learning_rate <= min_learning_rate:
            print_header(f'Terminating: Effective Learning Rate '
                         f'{self.effective_learning_rate:8.3e} <= minimum {min_learning_rate:8.3e}.')

        # Save training progress to disk
        if save_at_end:
            self.save_state(verbose=True)

        # Restore original learning_rate if it was altered
        if self.learning_rate != original_learning_rate:
            print(f'Restoring original learning rate {original_learning_rate:6.2e}.')
            self.set_learning_rate(original_learning_rate)

    # *********************************************************************************************
    def sieve_round(self, 
                    round: int,
                    num_batches: int, 
                    batches_per_epoch: int,
                    epochs_per_episode: int,
                    training_mode: str,
                    learning_rate: float, 
                    min_learning_rate: float,
                    reset_active_weight: bool=False,
                    obj_den_R_power: Optional[float]=None,
                    obj_den_thresh_power: Optional[float]=None,
                    R_sec_max: Optional[float]=None,
                    thresh_sec_max: Optional[float]=None,
                    charts: bool=False):
        """One round of sieving"""
        # Status update
        log2_lr: int = int(np.round(np.log(learning_rate) / np.log(2.0)))
        msg: str = f'Round {round}: {num_batches} batches @ LR 2^{log2_lr} in {training_mode} mode'
        msg_suffix: str = '.' if thresh_sec_max is None else f'; thresh_sec_max = {thresh_sec_max}'
        print_header(f'{msg}{msg_suffix}')

        # Thaw or freeze elements as requested
        if training_mode == 'mixture':
            self.set_obj_den_R_power(obj_den_R_power or 1.0)
            self.set_obj_den_thresh_power(obj_den_thresh_power or 1.0)
            self.save_state()
            self.thaw_mixture_parameters()
            self.freeze_candidate_elements()
        elif training_mode == 'joint':
            self.set_obj_den_R_power(obj_den_R_power or 0.0)
            self.set_obj_den_thresh_power(obj_den_thresh_power or 0.0)
            self.save_state()
            self.thaw_all()
        elif training_mode == 'element':            
            self.thaw_candidate_elements()
            self.freeze_mixture_parameters()
        else:
            raise ValueError('Bad training_mode.  Must be one of mixture, joint, element')

        # Set R_sec_max if it was specified
        if R_sec_max is not None:
            self.set_R_sec_max(R_sec_max)

        # Set thresh_sec_max if it was specified
        if thresh_sec_max is not None:
            self.set_thresh_sec_max(thresh_sec_max)

        # Delegate to search_adaptive
        max_batches = self.current_batch + num_batches
        self.search_adaptive(
            max_batches=max_batches,
            batches_per_epoch=batches_per_epoch,
            epochs_per_episode=epochs_per_episode,
            learning_rate=learning_rate,
            min_learning_rate=min_learning_rate,
            reset_active_weight=reset_active_weight)

        # Report progress
        self.report()
        if charts:
            self.plot_bar('log_like', sorted=False)
            self.plot_bar('hits', sorted=False)
        self.save_state()

    # *********************************************************************************************
    def sieve(self, base_size: int=64, nearest_ast: bool=True, charts: bool=False):
        """
        One shot method for automated search
        INPUTS:
            nearest_ast: Run comparison to nearest known asteroids at end
        """
        # Episode sizes
        num_batches_mixture = base_size * 8 # 512
        num_batches_joint = base_size * 32  # 2048

        # Set batches_per_epoch parameter
        batches_per_epoch: int = 64
        epochs_per_episode: int = 4

        # Set learning rates
        learning_rate_mixture = 2.0**-12
        learning_rate_joint = 2.0**-16
        
        # Minimum learning rate
        min_learning_rate_mixture = 2.0**-20
        min_learning_rate_joint = 2.0**-24

        # Set all layers thawed
        self.thaw_all()

        # Schedule of max thresholds and resolution
        sched = pd.DataFrame()
        sched['thresh_sec'] = np.array([7200, 5400, 3600, 2400], dtype=dtype_np)
        sched['R_sec'] = np.maximum(sched.thresh_sec * 0.50, 50.0)
        # Lenght of the schedule
        schedule_len = sched.shape[0]

        for i in range(schedule_len):
            # The round number, R_sec_max and thresh_sec_max
            round = 2*i
            R_sec_max = sched.R_sec[i]
            thresh_sec_max = sched.thresh_sec[i]

            # Reset the active weight every on mixture
            reset_active_weight_mixture = True
            reset_active_weight_joint = False
            
            # Mixture parameters
            self.sieve_round(round=round+1, 
                             num_batches=num_batches_mixture, 
                             batches_per_epoch=batches_per_epoch,
                             epochs_per_episode=epochs_per_episode,
                             training_mode='mixture', 
                             learning_rate=learning_rate_mixture, 
                             min_learning_rate=min_learning_rate_mixture,
                             reset_active_weight=reset_active_weight_mixture,
                             R_sec_max=R_sec_max,
                             thresh_sec_max=thresh_sec_max,
                             charts=charts)

            # Joint (mixture and candidate elements)
            self.sieve_round(round=round+2, 
                             num_batches=num_batches_joint, 
                             batches_per_epoch=batches_per_epoch,
                             epochs_per_episode=epochs_per_episode,
                             training_mode='joint',
                             learning_rate=learning_rate_joint, 
                             min_learning_rate=min_learning_rate_joint,
                             reset_active_weight=reset_active_weight_joint,
                             charts=charts)
        
        # Powers for denominator adjustment during fine tuning at the end
        fine_tuning = pd.DataFrame()
        fine_tuning['obj_den_R_power_mixture'] = np.array([2.0, 3.0, 4.0])
        fine_tuning['obj_den_thresh_power_mixture'] = np.array([4.0, 8.0, 16.0])
        fine_tuning['obj_den_R_power_joint'] = np.array([1.0, 2.0, 3.0])
        fine_tuning['obj_den_thresh_power_joint'] = np.array([2.0, 4.0, 8.0])
        fine_tuning['lr_adj_mixture'] = np.array([0.5, 0.5, 0.25])
        fine_tuning['lr_adj_joint'] = np.array([0.5, 0.5, 0.25])

        for i in range(fine_tuning.shape[0]):            
            # Fine tuning at the end
            obj_den_R_power = fine_tuning.obj_den_R_power_mixture.values[i]
            obj_den_thresh_power = fine_tuning.obj_den_R_power_mixture.values[i]
            lr_adj = fine_tuning.lr_adj_mixture.values[i]
            self.sieve_round(round=2*(schedule_len+i)+1,
                            num_batches=num_batches_mixture, 
                            batches_per_epoch=batches_per_epoch,
                            epochs_per_episode=epochs_per_episode,
                            training_mode='mixture', 
                            learning_rate=learning_rate_mixture*lr_adj, 
                            min_learning_rate=min_learning_rate_mixture*lr_adj,
                            reset_active_weight=True,
                            obj_den_R_power=obj_den_R_power,
                            obj_den_thresh_power=obj_den_thresh_power,
                            charts=charts)

            obj_den_R_power = fine_tuning.obj_den_R_power_joint.values[i]
            obj_den_thresh_power = fine_tuning.obj_den_R_power_joint.values[i]
            lr_adj = fine_tuning.lr_adj_mixture.values[i]
            self.sieve_round(round=2*(schedule_len+i)+2, 
                            num_batches=num_batches_joint, 
                            batches_per_epoch=batches_per_epoch,
                            epochs_per_episode=epochs_per_episode,
                            training_mode='joint',
                            learning_rate=learning_rate_joint*lr_adj,
                            min_learning_rate=min_learning_rate_joint*lr_adj,
                            reset_active_weight=True,
                            obj_den_R_power=obj_den_R_power,
                            obj_den_thresh_power=obj_den_thresh_power,
                            charts=charts)

        # Full round of charts
        self.plot_bar('log_like', sorted=False)
        self.plot_bar('hits', sorted=False)
        self.plot_bar('minus_log_R', sorted=False)
        self.plot_hist('log_like')
        self.plot_hist('hits')

        # Compare to nearest known asteroid if requested
        if nearest_ast:
            _ = self.nearest_ast()

            # Plot position error vs. nearest known asteroid
            self.plot_q_error(plot_type='cart', is_log=True, use_near_ast_dist=True)
            # Plot covariance error
            self.plot_q_error(plot_type='cov', is_log=True, use_near_ast_dist=True)

            # Plot error in orbital elements
            self.plot_elt_error_bar(elt_name='a', is_log=True)
            self.plot_elt_error_bar(elt_name='e', is_log=True)
            self.plot_elt_error_bar(elt_name='inc', is_log=True)
            self.plot_elt_error_bar(elt_name='Omega', is_log=True)
            self.plot_elt_error_bar(elt_name='omega', is_log=True)
            self.plot_elt_error_bar(elt_name='f', is_log=True)

    # *********************************************************************************************
    # Output candidate elements as DataFrame; save and load model state (including training history)
    # *********************************************************************************************

    # *********************************************************************************************
    def candidates_df(self):
        """The current candidates as a DataFrame."""
        # Generate the outputs
        score_outputs, orbital_elements, mixture_params = self.calc()

        # Unpack score_outputs
        log_like, hits, num_rows_close, loss = score_outputs
        # Extract thresh_deg from score layer
        thresh_deg = self.get_thresh_deg()
        # Extract the brightness H and standard deviation of magnitude from magnitude layer
        H = self.get_H()
        sigma_mag = self.get_sigma_mag()

        # Build DataFrame of orbital elements
        elts = elts_np2df(orbital_elements.numpy().T)
        
        # Add column with the element_id
        elts.insert(loc=0, column='element_id', value=self.element_id.numpy())
       
        # Add columns for the mixture parameters
        elts['num_hits'] = mixture_params[0].numpy()
        elts['R'] = mixture_params[1].numpy()
        elts['R_deg'] = dist2deg(elts.R)
        elts['R_sec'] = 3600.0 * elts.R_deg
        elts['R_max'] = mixture_params[2].numpy()
        elts['R_deg_max'] = dist2deg(elts.R_max)
        elts['thresh_s'] = deg2dist(thresh_deg)
        elts['thresh_deg'] = thresh_deg
        elts['thresh_sec'] = 3600.0 * elts.thresh_deg

        # Add columns with log likelihood and hits
        elts['log_like'] = log_like.numpy()
        elts['hits'] = hits.numpy()
        elts['num_rows_close'] = num_rows_close.numpy()

        # The brightness parameter H and standard deviation of magnitude
        elts['H'] = H
        elts['sigma_mag'] = sigma_mag

        # Add columns with the weights in the three modes
        elts['weight_joint'] = self.weight_joint.numpy()
        elts['weight_element'] = self.weight_element.numpy()
        elts['weight_mixture'] = self.weight_mixture.numpy()

        # Save this onto the model as the attribute elts_fit
        self.elts_fit = elts

        return elts

    # *********************************************************************************************
    def save_state(self, verbose: bool=False):
        """Save model state to disk: candidate elements and training history"""
        # Generate DataFrame with current candidates
        elts = self.candidates_df()
        # Status message
        if verbose:
            print(f'Saving candidate elements DataFrame in {self.file_path}.')
        
        # Save to disk
        elts.to_hdf(self.file_path, key='elts', mode='w')

        # Also save training history (summary and by element) after deduplicating
        self.train_hist_deduplicate()
        self.train_hist_summary.to_hdf(self.file_path, key='train_hist_summary', mode='a')
        self.train_hist_elt.to_hdf(self.file_path, key='train_hist_elt', mode='a')
        # Save elts_fit and elts_near_ast
        self.elts_fit.to_hdf(self.file_path, key='elts_fit', mode='a')
        self.elts_near_ast.to_hdf(self.file_path, key='elts_near_ast', mode='a')

    # *********************************************************************************************
    def set_elements(self, elts: pd.DataFrame):
        """Set training state """
        # Restore weights on candidate_elements layer of this model
        self.candidate_elements.load(elts=elts)
        # Restore weights on mixture_parameters layer of this model
        self.mixture_parameters.load(elts=elts)
        # Restore weights on magnitude layer
        self.magnitude.load(elts=elts)
        # Set the threshold in degrees
        self.score.set_thresh_deg(thresh_deg=self.thresh_deg)

        # Restore element weight in the three modes
        self.weight_joint.assign(elts['weight_joint'].values)
        self.weight_element.assign(elts['weight_element'].values)
        self.weight_mixture.assign(elts['weight_mixture'].values)

        # Alias history
        hist = self.train_hist_summary

        # Restore training counters
        self.current_episode = hist.episode.values[-1]
        self.current_epoch = hist.epoch.values[-1]
        self.current_batch = hist.batch.values[-1]
        self.training_time = hist.training_time.values[-1]

        # Restore learning rate and threshold members
        self.learning_rate = hist.learning_rate.values[-1]
        self.thresh_deg = elts['thresh_deg'].values

        # Recompile model
        self.recompile()
        # Update early stopping callback
        self.update_early_stop()

        # Check whether any layers are frozen
        if not self.train_hist_summary.candidate_elements_trainable.values[-1]:
            self.freeze_candidate_elements()
        if not self.train_hist_summary.mixture_parameters_trainable.values[-1]:
            self.freeze_mixture_parameters()

        # Run the save_weights method so the current weights are available to restore a bad training episode
        self.save_weights()

    # *********************************************************************************************
    def load(self, verbose: bool=True):
        """Load model state from disk: candidate elements and training history"""
        try:
            # Load DataFrames for candidate elements and auxiliary data
            elts = pd.read_hdf(self.file_path, key='elts')
            self.train_hist_summary = pd.read_hdf(self.file_path, key='train_hist_summary')
            self.train_hist_elt = pd.read_hdf(self.file_path, key='train_hist_elt')
            # Set the threshold for scoring
            self.thresh_deg = elts['thresh_deg']
            # Load nearest asteroid information; empty DataFrame when not yet calculated
            self.elts_fit = pd.read_hdf(self.file_path, key='elts_fit')
            self.elts_near_ast = pd.read_hdf(self.file_path, key='elts_near_ast')
        except FileNotFoundError:
            if verbose:
                print(f'Unable to find {self.file_path}.')
            return

        # Set candidate elements and mixture paramters
        self.set_elements(elts)

        # Status message
        if verbose:
            print(f'Loaded candidate elements and training history from {self.file_path}.')

    # *********************************************************************************************
    def revert_training(self, episode: int):
        """Revert training to the specified epoch"""
        # Mask for this epoch on training history
        hist = self.train_hist_elt
        mask = (hist.episode == episode)

        # Required columns for training state
        cols = self.candidates_df().columns
        # Candidate  elements used to restore training state
        elts = hist[cols][mask]
        # Rename the column epoch_ast to epoch; there is a name clash between
        # epoch of training and epoch in astronomy!
        elts['epoch'] = hist['epoch_ast'][mask]

        # Mask for history up to and including this episode
        mask = (hist.episode <= episode)
        self.train_hist_elt = hist[mask]

        # Summary history
        summary = self.train_hist_summary
        mask = (summary.episode <= episode)
        self.train_hist_summary = summary[mask]

        # Apply these elements with the shortened history
        self.set_elements(elts)

    # *********************************************************************************************
    def calc_ztf_hits(self, thresh_sec: float = 10.0):
        """Summarize the Hits in a DataFrame"""
        # The elements
        elts_fit = self.candidates_df()
        # The ztf_elt data used to build this model
        ztf_elt = self.ztf_elt.copy()
        
        # Unique ZTF identifiers
        ztf_id_unq = np.unique(ztf_elt.ztf_id)        
        # Slice of relevant ZTF observations
        ztf_slice = ztf_ast.loc[ztf_id_unq]

        # Build the hits withing the threshold
        thresh_deg: float = thresh_sec / 3600.0
        ztf_hits = make_ztf_batch(elts=elts_fit, ztf=ztf_slice, thresh_deg=thresh_deg)
        # Filter to relevant columns
        cols = \
            ['element_id', 'ztf_id', 'ObjectID', 'CandidateID', 'TimeStampID', 
            'mjd', 'ra', 'dec', 'mag_app', 'ux', 'uy', 'uz',
            'qx', 'qy', 'qz', 'vx', 'vy', 'vz', 
            'elt_ux', 'elt_uy', 'elt_uz', 'elt_r', 's_sec']
        ztf_hits = ztf_hits[cols]

        # Predicted outputs; use this for the predicted magnitude
        u_pred, delta_pred, mag_pred = self.predict_direction()
        # Add predicted magnitude to ztf_elt
        ztf_elt['mag_pred'] = mag_pred
        # Align with mag_pred with ztf_hits
        ztf_elt.set_index(keys=['ztf_id', 'element_id'], drop=True, inplace=True)
        keys = list(zip(ztf_hits.ztf_id, ztf_hits.element_id))
        mag_pred_hit = ztf_elt.loc[keys, 'mag_pred'].values
        # loc = ztf_hits.columns.get_loc('s_sec')
        # ztf_hits.insert(loc=loc, column='mag_pred', value=mag_pred_hit)
        ztf_hits['mag_pred'] = mag_pred_hit
        ztf_hits['mag_diff'] = np.abs(ztf_hits.mag_app - ztf_hits.mag_pred)

        # Display columns to visualize this without getting overwhelmed
        # cols = \
            # ['ztf_id', 'element_id', 'mjd',
            # 'ra', 'dec', 'mag_app', 'ux', 'uy', 'uz',
            # 'elt_ux', 'elt_uy', 'elt_uz', 
            # 's_sec']

        return ztf_hits

    # *********************************************************************************************
    # Plot training results: bar or time series charts
    # *********************************************************************************************

    # *********************************************************************************************
    def plot_bar(self, att_name: str = 'log_like', sorted: bool=False, episode=None):
        """
        Bar chart for one of the attributes monitored during training
        INPUTS:
            att_name: Name of the attribute to plot, e.g. 'log_like', 'hits', 'R_deg'
            episode:  The episode as of which to plot.  Default to the last episode.
        """
        # Alias the training history
        hist = self.train_hist_elt

        # Default epoch to the last one
        if episode is None:
            episode = np.max(hist.episode)

        # Filter to the requested epoch only
        mask = (hist.episode == episode)
        hist = hist[mask]

        # The values for this attribute
        values = hist[att_name].values

        # Sort the values in descending order if sorting was requested
        values_plot = np.sort(values)[::-1] if sorted else values

        # Mean of the plotted values
        values_mean = np.mean(values)

        # Score attributed formatted for inclusion in chart title or labels
        att_title_tbl = {
            'log_like': 'Log Likelihood',
            'hits': 'Hits',
            'R_deg': 'Resolution (degrees)',
        }
        att_title = att_title_tbl.get(att_name, att_name)
        
        # The x-label varies based on sorted input
        xlabel = 'Rank' if sorted else 'Element Num'

        # Bar plot of log likelihood after last episode
        fig, ax = plt.subplots()
        ax.set_title(f'{att_title} by Element (Episode {episode})')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'{att_title}')
        ax.bar(x=hist.element_num, height=values_plot, color='blue', label='elt')
        ax.axhline(y=values_mean, color='red', label=f'mean')
        # ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/{att_name}_bar.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    def plot_hist(self, att_name: str = 'log_like', x_axis: str = 'batch'):
        """
        Learning curve chart (progress vs. episode) for one of the score attributes log_like or hits.
        INPUTS:
            score_att: Name of the attribute to plot.  One of 'log_like', 'hits', 'R_deg'
            x_axis:    x-axis for plot; one of 'batch', 'epoch', 'episode', 'time'
        """
        # Alias the training history
        hist = self.train_hist_summary

        # Extract score_mean, score_std, score_min, score_max from hist
        score_mean = hist[f'{att_name}_mean'].values
        score_std = hist[f'{att_name}_std'].values
        score_min = hist[f'{att_name}_min'].values
        score_max = hist[f'{att_name}_max'].values
        # score_lo = score_mean - score_std
        # score_hi = score_mean + score_std
        score_lo = hist[f'{att_name}_q20']
        score_hi = hist[f'{att_name}_q80']

        # Only tabulate argmin for log likelihood
        min_elt = hist.log_like_argmin.values[-1]
        max_elt = hist.log_like_argmax.values[-1]
        min_label = f'min ({min_elt})' if att_name=='log_like' else 'min'
        max_label = f'max ({max_elt})' if att_name=='log_like' else 'max'

        # Score attributed formatted for inclusion in chart title or labels
        att_title_tbl = {
            'log_like': 'Log Likelihood',
            'hits': 'Hits',
            'R_deg': 'Resolution (degrees)'
        }
        att_title = att_title_tbl.get(att_name, att_name)

        # x-axis for plot
        if x_axis == 'batch':
            x_label = 'Batch Trained'
            x_plot = hist.batch
        elif x_axis == 'epoch':
            x_label = 'Epoch Trained'
            x_plot = hist.epoch
        elif x_axis == 'episode':
            x_label = 'Episode Trained'
            x_plot = hist.episode
        elif x_axis == 'time':
            x_label = 'Training Time (sec)'
            x_plot = hist.training_time
        else:
            raise ValueError('Bad x_axis input; must be one of batch, epoch, episode, time.')

        # Plot total log likelihood over training
        fig, ax = plt.subplots()
        ax.set_title(f'Training Progress: {att_title} by Element')
        ax.set_xlabel(x_label)
        ax.set_ylabel(att_title)
        # Plot mean +/- 1 SD
        ax.plot(x_plot, score_mean, color=color_mean, label='Mean')
        # Plot quantiles 20 and 80
        ax.plot(x_plot, score_lo, color=color_lo, label='Quantile 20')
        ax.plot(x_plot, score_hi, color=color_hi, label='Quantile 80')
        # Plot min and max
        ax.plot(x_plot, score_min, color=color_min, label=min_label)
        ax.plot(x_plot, score_max, color=color_max, label=max_label)
        # Legend etc
        ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/{att_name}_hist.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    # Find asteroid with nearest orbital elements; calculate position error vs. known elements
    # *********************************************************************************************

    # *********************************************************************************************
    def nearest_ast(self, search_type='cart'):
        """
        Search for the asteroid that is nearest in orbital trajectory to the candidate elements
        INPUTS:
            search_type: Search method for nearest element; one of 'cart', 'cov'
        """
        # The current candidate elements, i.e. after fitting process has been run
        self.elts_fit = self.candidates_df()

        # Drop irrelevant columns
        self.elts_fit.drop(columns=['weight_joint', 'weight_element', 'weight_mixture'], inplace=True)

        # Search for nearest asteroids to the fitted elements
        if search_type == 'cart':
            self.elts_near_ast = nearest_ast_elt_cart(self.elts_fit)
            # If we ran a Cartesian search, add an extra column with the Q_norm
            q_norm = elt_q_norm(elts=self.elts_fit, ast_num=self.elts_fit.nearest_ast_num)
            self.elts_fit['nearest_ast_q_norm'] = q_norm
            # self.elts_near_ast['nearest_ast_q_norm'] = q_norm
            loc = self.elts_near_ast.columns.get_loc('nearest_ast_dist')+1
            self.elts_near_ast.insert(loc=loc, column='nearest_ast_q_norm', value=q_norm)            
        elif search_type == 'cov':
            self.elts_near_ast = nearest_ast_elt_cov(self.elts_fit)
        else:
            raise ValueError(f'Bad search_type {search_type}, must be one of \'cart\', \'cov\'.')

        # Return the inputs and nearest asteroids
        return self.elts_fit, self.elts_near_ast

    # *********************************************************************************************
    def calc_error(self, elts_true):
        """
        Calculate error vs. known orbital elements
        INPUTS:
            elts_true: DataFrame of true orbital elements
        OUTPUTS:
            hist_err:  DataFrame with errors by element over the training episode
            q_err:     Position error of current orbital elements
        """
        # Create a copy of the training history that will be augmented with error columns
        hist_err = self.train_hist_elt.copy()
        # The columns we're computing errors for: six orbital elements
        cols_elt = ['a', 'e', 'inc', 'Omega', 'omega', 'f']
        # Iterate through the elements; calulate error and log error
        for col in cols_elt:
            # col_true = f'{col}_true'
            col_err = f'{col}_err'
            col_log_err = f'{col}_log_err'
            # hist[col_true] = elts_true.loc[hist.element_num, col].values
            hist_err[col_err] = np.abs(hist_err[col] - elts_true.loc[hist_err.element_num, col].values)
            hist_err[col_log_err] = np.log(hist_err[col_err])

        # Unpack true orbital elements
        a = elts_true.a.values
        e = elts_true.e.values
        inc = elts_true.inc.values
        Omega = elts_true.Omega.values
        omega = elts_true.omega.values
        f = elts_true.f.values
        epoch = elts_true.epoch.values

        # Calculate position of the candidate elements
        q_pred, v_pred = self.predict_position()

        # Calculate position of true elements
        q_true, v_true = self.position(a, e, inc, Omega, omega, f, epoch)

        # Position error - flat; shape [data_size, 3,]
        q_err_flat = q_pred - q_true
        # Norm of position error l flat; shape [data_size,]
        q_err_norm_flat = tf.linalg.norm(q_err_flat, axis=-1)
        # Position error as a ragged tensor; shape [batch_size, num_obs,]
        q_err_r = tf.RaggedTensor.from_row_lengths(values=q_err_norm_flat, row_lengths=self.row_lengths, name='q_err_r')
        # Mean position error by element; shape [batch_size,]
        q_err = tf.reduce_mean(q_err_r, axis=-1, name='q_err')

        return hist_err, q_err

    # *********************************************************************************************
    # Plot error vs. known orbital elements or selected element
    # *********************************************************************************************

    # *********************************************************************************************
    def plot_q_error(self, 
                     elts_true=None, 
                     plot_type: str = 'cart',
                     is_log: bool=False, 
                     use_near_ast_dist: bool=True):
        """
        Bar chart of position error by orbital element.
        INPUTS:
            elts_true: DataFrame of true orbital elements
            plot_type: One of 'cart' or 'cov' for Cartesian Distance in AU or dimenionless covariance metric
            is_log:    Whether to plot log(q_err) instead of q_err
            use_near_ast_dist:
                    use_ast_dist=True means plot the distance to the nearest asteroid
                    saved when near_ast_dist was called.
                    These are based on 20 years of data sampled every 3 months.
                    use_near_ast_dist = False means to regenerate two new
                    trajectories at times matching the observation times.
        """
        # If elts_true wasn't specified, it defaults to nearest asteroid elements
        elts_true = self.elts_near_ast

        # Calculate distance error with selected method: live or from near_ast
        if use_near_ast_dist:
            q_err = self.elts_near_ast.nearest_ast_dist
            calc_method = '20 Years of Monthly Dates'
        else:
            hist_err, q_err = self.calc_error(elts_true)
            calc_method = 'Observation Dates'

        # Compute the q_norm error
        q_norm_err = self.elts_near_ast.nearest_ast_q_norm

        # Choose between Cartesian distance and q_norm
        err_tbl = {
            'cart': q_err,
            'cov': q_norm_err,
        }
        err = err_tbl[plot_type]

        # Get element numbers and mean error over the whole batch
        element_num = np.arange(self.batch_size, dtype=np.int32)
        mean_err = np.exp(np.mean(np.log(err+2**-40))) if is_log else np.mean(err)

        # Type of error
        is_q_norm: bool = (plot_type == 'cov')

        # Reference errors for plot
        ref_err = 0.01 if is_q_norm else 0.001

        # Heights for plot
        height = 1.0 / err if is_log else err
        mean_height = 1.0 / mean_err if is_log else mean_err
        ref_height = 1.0 / ref_err if is_log else ref_err

        # Chart title
        error_tag = 'Covariance Error' if is_q_norm else 'Mean Position Error in AU'
        title_tag_log = 'Precision of Covariance' if is_q_norm else 'Mean Precision of Position in AU'
        title_tag = title_tag_log if is_log else error_tag
        title = f'{title_tag} vs. Elements on {calc_method}'

        # Caption for y-axis
        ylabel_prefix = 'Inverse ' if is_log else ''
        # ylabel_suffix = ' (log scaled)' if is_log else ''
        ylabel = f'{ylabel_prefix}{error_tag}'

        # Plot selected error
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel('Element')
        ax.set_ylabel(ylabel)
        ax.bar(x=element_num, height=height, label='elt', color='blue')
        ax.axhline(y=mean_height, label='mean', color='red')
        ax.axhline(y=ref_height, label=f'err={ref_err}', color='black')
        if is_log:
            ax.set_yscale('log', basey=2)
        ax.grid()
        ax.legend()
        # fig.savefig('../figs/training/q_error.png', bbox_inches='tight')
        plt.show()
        return fig, ax
        
    # *********************************************************************************************
    def plot_elt_error_bar(self,
                           elt_name: str, 
                           elts_true=None,
                           is_log: bool=False):
        """Bar chart of error in a specified element"""
        # If elts_true wasn't specified, it defaults to nearest asteroid elements
        elts_true = self.elts_near_ast

        # Get element numbers and mean error over the whole batch
        element_num = np.arange(self.batch_size, dtype=np.int32)

        # Fitted value of the selected orbital element
        sel_elt_fit = self.elts_fit[elt_name].values

        # True value of selected orbital element
        sel_elt_true = elts_true[elt_name].values

        # Error in selected element
        is_angle = elt_name in ['inc', 'Omega', 'omega', 'f']
        if is_angle:
            diff_deg = np.rad2deg(sel_elt_fit - sel_elt_true)
            # Subtract out full circles to get error in the interval [-180, 180)
            n = (diff_deg + 180.0) // 360.0
            # The error in degrees is the absolute value; in [0, 180)
            elt_err = np.abs(diff_deg - n * 360.0)
        else:
            elt_err = np.abs(sel_elt_fit - sel_elt_true)

        # Error to plot
        values_plot = 1.0 / elt_err if is_log else elt_err
        
        # Mean
        mean_plot = np.exp(np.mean(np.log(values_plot))) if is_log else np.mean(elt_err)

        # Table with reference values
        ref_err_tbl = {
            'a': 0.001,
            'e': 0.0001,
            'inc': 0.1,
            'Omega': 0.1,
            'omega': 0.1,
            'f': 0.1,
        }
        ref_err = ref_err_tbl[elt_name]
        ref_plot = 1.0 / ref_err if is_log else ref_err

        # Chart titles
        title_prefix = 'Precision' if is_log else 'Error'
        title_suffix = ' (degrees)' if is_angle else ''
        title = f'{title_prefix} in Orbital Element {elt_name}{title_suffix}'
        err_label = 'Inverse Error in' if is_log else 'Absolute Error in'

        # Bar plot of error in selected orbital element
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel('Candidate Element')
        ax.set_ylabel(f'{err_label} {elt_name}{title_suffix}')
        ax.bar(x=element_num, height=values_plot, color='blue', label='elt')
        ax.axhline(y=mean_plot, color='red', label='mean')
        ax.axhline(y=ref_plot, color='black', label=f'error= {ref_err}')
        if is_log:
            ax.set_yscale('log', basey=2)
        ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/{att_name}_bar.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    def plot_elt_error_hist(self, elt_name: str, elts_true=None, is_log: bool=True, elt_num: int=None):
        """
        Plot (log) error of selected orbital element vs. known elements
        INPUTS:
            elt_name:  Name of orbital element; one of 'a', 'e', 'inc', 'Omega', 'omega', 'f'
            elts_true: DataFrame of true orbital elements; defaults to self.near_ast
            is_log:    Flag; whether to plot log(err) (True) or err (False)
            elt_num:   Specific element_num to plot; plot just this one rather than mean, std, etc.
        """
        # If elts_true wasn't specified, it defaults to nearest asteroid elements
        if elts_true is None:
            elts_true = self.elts_near_ast

        # Calculate error
        hist_err, q_err = self.calc_error(elts_true)
        
        # Column names with error
        col_err = f'{elt_name}_err'
        col_log_err = f'{elt_name}_log_err'
        col_plot = col_log_err if is_log else col_err
        title = f'Log Error of {elt_name}' if is_log else f'Error of {elt_name}'
        ylabel = f'Log Error' if is_log else f'Error'

        # Selected error; reshape to size (episode_count, batch_size)
        err = hist_err[col_plot].values.reshape((-1, self.batch_size))

        # Style of plot: summary or selected element
        is_summary: bool = elt_num is None

        # The errors to plot
        if is_summary:
            err_mean = np.mean(err, axis=1)
            err_std = np.std(err, axis=1)
            err_min = np.min(err, axis=1)
            err_max = np.max(err, axis=1)
            err_lo = np.quantile(err, q=0.20, axis=1)
            err_hi = np.quantile(err, q=0.80, axis=1)
        else:
            err_selected = err[:, elt_num]

        # Plot selected error
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel('Episode')
        ax.set_ylabel(ylabel)
        if is_summary:
            ax.plot(err_mean, color=color_mean, label='mean')
            # ax.plot(err_mean - err_std, color=color_lo, label='mean -1 SD')
            # ax.plot(err_mean + err_std, color=color_hi, label='mean +1 SD')
            ax.plot(err_lo, color=color_lo, label='quantile 20')
            ax.plot(err_hi, color=color_hi, label='quantile 80')
            ax.plot(err_min, color=color_min, label='min')
            ax.plot(err_max, color=color_max, label='max')
        else:
            element_num = np.arange(self.batch_size, dtype=np.int32)
            ax.plot(err_selected, color='blue', label=f'element_num {element_num}')
        ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/log_like_bar.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    def plot_control(self, element_num):
        """Plot control variables"""
        # Get control variables for selected element_num
        hist = self.train_hist_elt
        mask = (hist.element_num == element_num)
        hist = hist[mask]

        # Plot control variables
        fig, ax = plt.subplots()
        ax.set_title(f'Control Variables for Element {element_num}')
        ax.set_xlabel('Episode')
        ax.set_ylabel('Control Variable ([0, 1] Scale)')
        # Plot mixture parameters
        ax.plot(hist.episode, hist.num_hits_, label='num_hits')
        ax.plot(hist.episode, hist.R_, label='R')
        # Plot elements
        ax.plot(hist.episode, hist.a_, label='a')
        ax.plot(hist.episode, hist.e_, label='e')
        ax.plot(hist.episode, hist.inc_, label='inc')
        ax.plot(hist.episode, hist.Omega_, label='Omega')
        ax.plot(hist.episode, hist.omega_, label='omega')
        ax.plot(hist.episode, hist.f_, label='f')
        # Legend, scale etc
        ax.set_ylim([0, 1])
        ax.grid()
        # Put legend to the right of the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # fig.savefig('../figs/training/control_variables.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    def plot_element_att(self, att_name: str, element_num: int):
        """
        Plot an attribute for a specified candidate element
        INPUTS:
            att_name:    Name of the attribute to plot.  One of:
                         'log_like', 'hits', 'R_deg', 'R_sec'
            element_num: Which element in the batch to plot
        """
        # Get history for selected element_num
        hist = self.train_hist_elt
        mask = (hist.element_num == element_num)
        hist = hist[mask]

        # Extract param_mean, param_std, param_min, param_max from hist
        att_values = hist[att_name].values

        # Attribute formatted for inclusion in chart title or labels
        att_title = {
            'log_like': 'Log Likelihood',
            'hits':     'Hits',
            'R_deg':    'Resolution (Degrees)'
        }[att_name]

        # Plot the named attribute
        fig, ax = plt.subplots()
        ax.set_title(f'{att_title} for Element {element_num}')
        ax.set_xlabel('Episode')
        ax.set_ylabel(f'{att_title}')
        # Plot mixture parameters
        ax.plot(hist.episode, att_values, label='f{att_name}')
        ax.grid()
        # fig.savefig('../figs/training/element_{att_name}.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    #  Diagnostic
    # *********************************************************************************************

    def report(self):
        """Run model and report a few summary outputs to console"""
        # Run model on current candidates
        score_outputs, candidate_elements, mixture_params = self.calc()

        # Unpack score_outputs
        log_like, hits, num_rows_close, loss = score_outputs

        # Count number of good elements (with at least 10 hits)
        hits = np.round(hits.numpy())
        is_good_elt = (hits >= 10.0)
        num_good_elts = np.sum(is_good_elt)

        # Number of quality elements
        num_good_elts = np.sum(is_good_elt)

        # Summarize hits
        hits_good = np.mean(hits[is_good_elt])       
        hits_bad = np.mean(hits[~is_good_elt])
        hits_mean = np.mean(hits)
        hits_med = np.median(hits)
        hits_geo = np.exp(np.mean(np.log(hits+1.0)))-1.0
        # hits_std = np.std(hits)
        hits_min = np.min(hits)
        hits_max = np.max(hits)

        # Summarize log likelihood
        log_like = log_like.numpy()
        log_like_good = np.mean(log_like[is_good_elt])
        log_like_bad = np.mean(log_like[~is_good_elt])
        log_like_mean = np.mean(log_like)
        log_like_med = np.median(log_like)
        log_like_geo = np.exp(np.mean(np.log(np.maximum(log_like, 1.0))))
        # log_like_std = np.std(log_like)
        log_like_min = np.min(log_like)
        log_like_max = np.max(log_like)

        # Summarize resolution
        R = mixture_params[1].numpy()
        R_sec = dist2sec(R)
        R_sec_good = np.mean(R_sec[is_good_elt])
        R_sec_bad = np.mean(R_sec[~is_good_elt])
        R_sec_mean = np.mean(R_sec)
        R_sec_med = np.median(R_sec)
        R_sec_geo = np.exp(np.mean(np.log(R_sec+1.0)))-1.0
        # R_sec_std = np.std(R_sec)
        R_sec_min = np.min(R_sec)
        R_sec_max = np.max(R_sec)

        # Threshold in arc seconds
        thresh_deg = self.get_thresh_deg()
        thresh_sec = thresh_deg * 3600.0
        thresh_sec_good = np.mean(thresh_sec[is_good_elt])
        thresh_sec_bad = np.mean(thresh_sec[~is_good_elt])
        thresh_sec_mean = np.mean(thresh_sec)
        thresh_sec_med = np.median(thresh_sec)
        thresh_sec_geo = np.exp(np.mean(np.log(thresh_sec+1.0)))
        # thresh_sec_std = np.std(thresh_sec)
        thresh_sec_min = np.min(thresh_sec)
        thresh_sec_max = np.max(thresh_sec)

        # Report on log likelihood and resolution
        print(f'\nGood elements (hits >= 10): {num_good_elts:6.2f}\n')
        print(f'         \\  log_like :  hits  :    R_sec : thresh_sec')
        print(f'Mean Good: {log_like_good:8.2f}  : {hits_good:6.2f} : {R_sec_good:8.2f} : {thresh_sec_good:8.2f}')
        print(f'Mean Bad : {log_like_bad:8.2f}  : {hits_bad:6.2f} : {R_sec_bad:8.2f} : {thresh_sec_bad:8.2f}')
        print(f'Mean     : {log_like_mean:8.2f}  : {hits_mean:6.2f} : {R_sec_mean:8.2f} : {thresh_sec_mean:8.2f}')
        print(f'Median   : {log_like_med:8.2f}  : {hits_med:6.2f} : {R_sec_med:8.2f} : {thresh_sec_med:8.2f}')
        print(f'GeoMean  : {log_like_geo:8.2f}  : {hits_geo:6.2f} : {R_sec_geo:8.2f} : {thresh_sec_geo:8.2f}')
        # print(f'Std      : {log_like_std:8.2f}  : {hits_std:6.2f} : {R_sec_std:8.2f} : {thresh_sec_std:8.2f}')
        print(f'Min      : {log_like_min:8.2f}  : {hits_min:6.2f} : {R_sec_min:8.2f} : {thresh_sec_min:8.2f}')
        print(f'Max      : {log_like_max:8.2f}  : {hits_max:6.2f} : {R_sec_max:8.2f} : {thresh_sec_max:8.2f}')

        print(f'Trained for {self.current_batch} batches over {self.current_epoch} epochs '
              f'and {self.current_episode} episodes (elapsed time {int(np.round(self.training_time))} seconds).')

    def review_members(self):
        """Print diagnostic review of members to console"""
        preview_size = 5

        # Input data description
        print(f'batch_size: {self.batch_size}')
        print(f'\element_id: shape {self.element_id.shape}')
        print(self.element_id[0:preview_size])
        print(f'\nrow_lengths: shape {self.row_lengths.shape}')
        print((self.row_lengths[0:preview_size]))
        print(f'\ndata_size: {self.data_size}')
        print(f'traj_shape: {self.traj_shape}')

        # Learning rate and training configuration
        print(f'\nlearning_rate: {self.learning_rate:6.2e}')
        print(f'clipnorm: {self.clipnorm}')
        print(f'episode_length: {self.episode_length}')
        print(f'\ncurrent_episode: {self.current_episode}')
        print(f'current_epoch: {self.current_epoch}')
        print(f'current_batch: {self.current_batch}')

        # Weights history
        print(f'\ncandidate_elements_hist: list of tensors shape '
              f'[{len(self.candidate_elements_hist)}, {len(self.candidate_elements_hist[0])}, {len(self.candidate_elements_hist[0][0])}]')
        print([x for x in self.candidate_elements_hist[0][0][0:preview_size]])
        print(f'\nmixture_parameters_hist: list of tensors shape '
              f'[{len(self.mixture_parameters_hist)}, {len(self.mixture_parameters_hist[0])}, {len(self.mixture_parameters_hist[0][0])}]')
        print([x for x in self.mixture_parameters_hist[0][0][0:preview_size]])

        # Serialization
        print(f'\nelts_hash_id: {self.elts_hash_id}')
        print(f'ztf_elt_hash_id: {self.ztf_elt_hash_id}')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
