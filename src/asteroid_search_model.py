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
import time
from datetime import timedelta

# MSE imports
from asteroid_model import AsteroidDirection, AsteroidMagnitude, make_model_ast_pos
from asteroid_search_layers import CandidateElements, MixtureParameters, TrajectoryScore
from candidate_element import elts_np2df
from asteroid_integrate import calc_ast_pos
from candidate_element import perturb_elts
from ztf_element import ztf_elt_hash
from asteroid_search_report import traj_diff
from nearest_asteroid import nearest_ast_elt_cart, nearest_ast_elt_cov, elt_q_norm
from asteroid_dataframe import calc_ast_data, spline_ast_vec_df
from astro_utils import deg2dist, dist2deg, dist2sec
from tf_utils import tf_quiet, Identity
from utils import print_header

# Typing
from typing import List, Tuple, Dict, Optional, Union

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

# Save directory for candidate elements
save_dir = '../data/candidate_elt'

# Set plot style variables
mpl.rcParams['figure.figsize'] = [16.0, 10.0]
mpl.rcParams['font.size'] = 16

# Colors for plots
color_mean = 'blue'
color_lo = 'orange'
color_hi = 'green'
color_min = 'red'
color_max = 'purple'

# ********************************************************************************************************************* 
def make_opt_adam(learning_rate: float, 
                  clipnorm: float=1.0, 
                  clipvalue:Optional[float]=None) \
                  -> keras.optimizers.Optimizer:
    """
    Build Adam optimizer for training 
    Default settings are:
    learning_rate = 1.0E-3
    clipnorm = None
    clipvalue = None
    These are changed based on trial and error.
    Other arguments left at default settings
    """

    # Settings for other arguments; leave at defaults
    beta_1 = 0.900          # default 0.900
    beta_2 = 0.999          # default 0.999
    epsilon = 1.0E-7        # default 1.0E-7
    amsgrad = False         # default False

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args = {
        'learning_rate': learning_rate,
        'beta_1': beta_1,
        'beta_2': beta_2,
        'epsilon': epsilon,
        'amsgrad': amsgrad,
    }
    # Add clipnorm if it was set
    if clipnorm is not None:
        opt_args['clipnorm'] = clipnorm
    # Add clipvalue if it was set
    if clipvalue is not None:
        opt_args['clipvalue'] = clipvalue

    # Build the optimizer
    opt = keras.optimizers.Adam(**opt_args)
    return opt
    
# ********************************************************************************************************************* 
def make_opt_rmsprop(learning_rate: float, 
                     clipnorm: float=1.0, 
                     clipvalue: Optional[float]=None) \
                     -> keras.optimizers.Optimizer:
    """
    Build RMSprop optimizer for training 
    Default settings are:
    learning_rate = 1.0E-3
    rho = 0.90
    momentum = 0.0
    epsilon = 1.0E-7
    centered = False
    clipnorm = None
    clipvalue = None
    """

    # Settings for other arguments; leave at defaults
    rho = 0.900             # default 0.900
    momentum = 0.000        # default 0.000
    epsilon = 2.0**-23      # default 1.0E-7; nearest power of 2
    centered = False

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args = {
        'learning_rate': learning_rate,
        'rho': rho,
        'momentum': momentum,
        'epsilon': epsilon,
        'centered': centered,
    }
    # Add clipnorm if it was set
    if clipnorm is not None:
        opt_args['clipnorm'] = clipnorm
    # Add clipvalue if it was set
    if clipvalue is not None:
        opt_args['clipvalue'] = clipvalue

    # Build the optimizer
    opt = keras.optimizers.RMSprop(**opt_args)
    return opt

# ********************************************************************************************************************* 
def make_opt_adadelta(learning_rate: float, 
                      clipnorm: float=1.0, 
                      clipvalue: Optional[float]=None) \
                      -> keras.optimizers.Optimizer:
    """
    Build Adadelta optimizer for training 
    Default settings are:
    learning_rate = 1.0E-3
    rho = 0.950
    epsilon = 1.0E-7
    clipnorm = None
    clipvalue = None
    """

    # Settings for other arguments; leave at defaults
    rho = 0.950             # default 0.950
    epsilon = 2.0**-23      # default 1.0E-7; nearest power of 2

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args = {
        'learning_rate': learning_rate,
        'rho': rho,
        'epsilon': epsilon,
    }
    # Add clipnorm if it was set
    if clipnorm is not None:
        opt_args['clipnorm'] = clipnorm
    # Add clipvalue if it was set
    if clipvalue is not None:
        opt_args['clipvalue'] = clipvalue

    # Build the optimizer
    opt = keras.optimizers.Adadelta(**opt_args)
    return opt

# ********************************************************************************************************************* 
def make_opt(optimizer_type: str, 
             learning_rate: float, 
             clipnorm: float=1.0, 
             clipvalue: Optional[float]=None) \
             -> keras.optimizers.Optimizer:
    """
    Create an instance of the specified optimizer.
    INPUTS:
        learning_rate:  learning_rate for the optimizer
        optimizer_type: one of 'adam', 'rmsprop', 'adadelta'
        clipnorm:       gradient clipping by norm of gradient vector
        clipvalue:      gradient clipping by element of gradient vector
    """
    # Table of factory functions keyed by optimizer_type
    optimizer_func = {
        'adam': make_opt_adam,
        'rmsprop': make_opt_rmsprop,
        'adadelta': make_opt_adadelta,
    }
    # The factor function for the selected optimizer_type
    optimizer_func = optimizer_func[optimizer_type]
    # Instantiate this optimizer with selected input parameters
    return optimizer_func(learning_rate=learning_rate, clipnorm=clipnorm, clipvalue=clipvalue)

# ********************************************************************************************************************* 
def candidate_elt_hash(elts: pd.DataFrame, thresh_deg: float):
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements; and 2 mixture params num_hits, R
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
    OUTPUTS:
        hash_id:    Unique ID for these inputs
    """
    # Columns of the Dataframe to hash
    cols_to_hash = ['element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch', 'num_hits', 'R']
    # Tuple of int64; one per orbital element candidate
    hash_df = tuple((pd.util.hash_pandas_object(elts[cols_to_hash])).values)
    # Combine the element hash tuple with the threshold
    thresh_deg_int = int(thresh_deg*2**48)
    hash_id = abs(hash(hash_df + (thresh_deg_int, )))

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
                 learning_rate: float = 2.0**-13, 
                 clipnorm: float = 1.0,
                 **kwargs):
        """
        INPUTS:
            elts:           DataFrame with initial guess for orbital elements.
                            Columns: element_id, a, e, inc, Omega, omega, f, epoch
                            Output of asteroid_elts, perturb_elts or random_elts
            ztf_elt:        DataFrame with ZTF observations within thresh_deg degrees of
                            of the orbits predicted by these elements.
                            Output of make_ztf_batch or load_ztf_batch
            site_name:      Used for topos adjustment, e.g. 'geocenter' or 'palomar'
            thresh_deg:     Threshold used for filtering the ZTF observations
            optimizer_type: One of 'adam', 'rmsprop', 'adadelta'
            learning_rate:  Initial value of learning rate
            clipnorm:       Initial value of clipnorm for gradient clipping
        """
        # Initialize tf.keras.Model
        super(AsteroidSearchModel, self).__init__(**kwargs)
        
        # *****************************************************************************************
        # Description of input data 
        # *****************************************************************************************

        # Batch size comes from elts
        self.batch_size = elts.shape[0]

        # Shape of the observed trajectories
        self.data_size: int = ztf_elt.shape[0]
        self.traj_shape = (self.data_size, space_dims)

        # The element_id for the elements in this batch
        self.element_id = keras.backend.constant(value=elts.element_id, dtype=tf.int32)

        # Numpy array and tensor of observation times; flat, shape (data_size,)
        ts_np = ztf_elt.mjd.values.astype(dtype_np)
        self.ts = keras.backend.constant(value=ts_np, shape=(self.data_size,), dtype=dtype)

        # Get observation count per element
        row_lengths_np = ztf_elt.element_id.groupby(ztf_elt.element_id).count()
        self.row_lengths = keras.backend.constant(value=row_lengths_np, shape=(self.batch_size,), dtype=tf.int32)
        self.row_lengths_float = tf.cast(x=self.row_lengths, dtype=dtype)

        # Threshold - original, used for data
        self.thresh_deg_data = thresh_deg
        # Threshold - live, used for scoring; array of shape [batch_size,]
        self.thresh_deg = np.full(shape=self.batch_size, fill_value=thresh_deg, dtype=dtype_np)

        # Save the original ztf_elt frame for reuse
        self.ztf_elt = ztf_elt

        # Observed directions; extract from ztf_elt DataFrame
        cols_u_obs = ['ux', 'uy', 'uz']
        u_obs_np = ztf_elt[cols_u_obs].values.astype(dtype_np)

        # Apparent magnitude; extract from ztf_elt DataFrame
        cols_mag_app = ['mag_app']
        mag_app_np = ztf_elt[cols_mag_app].values.astype(dtype_np)

        # *****************************************************************************************
        # Layers for candidate elements, asteroid direction and score
        # *****************************************************************************************

        # Set of trainable weights with candidate orbital elements; initialize according to elts
        self.candidate_elements = CandidateElements(elts=elts, name='candidate_elements')

        # Set of trainable weights with candidate mixture parameters
        # self.mixture_parameters = MixtureParameters(elts=elts, thresh_deg=self.thresh_deg, name='mixture_parameters')
        self.mixture_parameters = MixtureParameters(elts=elts, name='mixture_parameters')

        # The predicted direction; shape is [data_size, 3,]
        self.direction = AsteroidDirection(ts_np=ts_np, row_lengths_np=row_lengths_np, 
                                           site_name=site_name, name='direction')

        # The predicted magnitude; shape is [data_size, ]
        self.magnitude = AsteroidMagnitude(ts_np=ts_np, row_lengths_np=row_lengths_np, elts=elts, name='magnitude')

        # Bind the direction layer to this model for legibility
        self.position = self.direction.position

        # Calibration arrays (flat)
        self.cols_q_ast = ['qx', 'qy', 'qz']
        self.cols_v_ast = ['vx', 'vy', 'vz']
        self.q_ast = ztf_elt[self.cols_q_ast].values.astype(dtype_np)
        self.v_ast = ztf_elt[self.cols_v_ast].values.astype(dtype_np)

        # Run calibration
        self.position.calibrate(elts=elts, q_ast=self.q_ast, v_ast=self.v_ast)

        # Score layer for these observations
        self.score = TrajectoryScore(row_lengths_np=row_lengths_np, 
                                     u_obs_np=u_obs_np, mag_app_np=mag_app_np,
                                     thresh_deg=self.thresh_deg, name='score')

        # Position model - used for comparing trajectories of fitted and known orbital elements
        self.model_pos: keras.Model = make_model_ast_pos(ts_np=ts_np, row_lengths_np=row_lengths_np)

        # *****************************************************************************************
        # Variables for adaptive training and training history
        # *****************************************************************************************

        # Save the learning rate on the model object to facilitate adaptive training
        self.learning_rate: float = learning_rate
        self.clipnorm: float = clipnorm

        # Build the selected optimizer with the input learning rate
        self.optimizer_type = optimizer_type.lower()
        self.optimizer = make_opt(optimizer_type=self.optimizer_type, 
                                  learning_rate=self.learning_rate, 
                                  clipnorm=self.clipnorm, clipvalue=None)
        # Compile the model with this optimizer instance
        self.recompile()

        # Set the learning rate factors for adaptive training
        self.lr_factor_dn: float = 0.5      # global learning rate; currently not used
        self.lr_factor_elt_dn: float = 0.5  # elementwise learning rate

        # Initialize loss history and total training time
        self.training_time: float = 0.0

        # Epoch and episode counters
        self.batches_per_epoch: int = 100
        self.epochs_per_episode: int = 5
        self.samples_per_epoch = self.batches_per_epoch * self.batch_size
        self.episode_length: int = 0
        self.current_episode: int = 0
        self.current_epoch: int = 0
        self.current_batch: int = 0

        # Tensor of ones to pass as dummy inputs for evaluating one batch or training one epoch
        self.x_eval = tf.ones(shape=self.batch_size, dtype=dtype)
        self.x_trn = tf.ones(self.samples_per_epoch, dtype=dtype)

        # Start training timer
        self.t0: float = time.time()
        self.episode_t0: float = time.time()

        # Cached values of log likelihood and loss
        self.log_like: np.ndarray = np.zeros(self.batch_size, dtype=dtype_np)
        self.hits: np.ndarray = np.zeros(self.batch_size, dtype=dtype_np)
        self.log_like_mean: float = 0.0
        self.hits_mean: float = 0.0
        self.loss: float = 0.0

        # Training mode: one of 'joint', 'element', 'mixture'
        self.training_mode = 'joint'

        # Weight factors for training in three modes: joint, element, mixture
        self.weight_joint = tf.Variable(initial_value=np.ones(self.batch_size, dtype=dtype_np), trainable=False, dtype=dtype)
        self.weight_element = tf.Variable(initial_value=np.ones(self.batch_size, dtype=dtype_np), trainable=False, dtype=dtype)
        self.weight_mixture = tf.Variable(initial_value=np.ones(self.batch_size, dtype=dtype_np), trainable=False, dtype=dtype)
        # Mean of active weights; for effective learning rate
        self.active_weight_mean: float = 1.0
        self.effective_learning_rate: float = self.learning_rate * self.active_weight_mean

        # Early stopping callback
        self.update_early_stop()
        self.callbacks = [self.early_stop]

        # Initialize lists with training history
        # Weights, elements and mixture parameters
        self.candidate_elements_hist = []
        self.mixture_parameters_hist = []

        # Log likelihoods and losses
        self.log_like_hist = []
        self.hits_hist = []
        self.log_like_mean_hist = []
        self.hits_mean_hist = []
        self.loss_hist = []

        # Training counters and times
        self.episode_hist = []
        self.epoch_hist = []
        self.batch_hist = []
        self.training_time_hist = []

        # DataFrame with training history
        self.train_hist_summary = pd.DataFrame()
        self.train_hist_elt = pd.DataFrame()

        # Add the first entries with initial weights and losses
        self.save_weights()

        # Threshold for and count of the number of bad training episodes
        self.bad_episode_thresh: float = 0.00
        self.bad_episode_count: int = 0

        # Elements in search for nearest asteroid: fit (inputs) and nearest asteroid (outputs)
        self.elts_fit = None
        self.elts_near_ast = None

        # Save hash IDs for elts and ztf_elts
        self.elts_hash_id: int = candidate_elt_hash(elts=elts, thresh_deg=self.thresh_deg_data)
        self.ztf_elt_hash_id: int = ztf_elt_hash(elts=elts, ztf=self.ztf_elt, thresh_deg=self.thresh_deg_data, near_ast=False)
        # File name for saving training progress
        self.file_path = f'{save_dir}/candidate_elt_{self.elts_hash_id}.h5'

        # Initialize model with frozen score layer
        self.freeze_score()

    # @tf.function
    def call(self, inputs=None):
        """
        Predict directions from current elements and score them.  
        inputs is a dummmy with no affect on the results; it is there to satisfy keras.Model API
        """

        # Extract the candidate elements and mixture parameters; pass dummy inputs to satisfy keras Layer API
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=inputs)
        
        # Extract mixture parameters; pass dummy inputs to satisfy keras Layer API
        num_hits, R = self.mixture_parameters(inputs=inputs)
        
        # Stack the current orbital elements.  Shape is [batch_size, 7,]
        orbital_elements = keras.backend.stack([a, e, inc, Omega, omega, f, epoch,])

        # Stack mixture model parameters. Shape is [batch_size, 2,]
        mixture_parameters = keras.backend.stack([num_hits, R,])

        # Tensor of predicted directions.  Shape of u_pred is [data_size, 3,]
        # delta_pred is distance from earth to ast.  
        # q_ast also included for use in magnitude calculation. Shape [data_size, 3]
        u_pred, delta_pred, q_ast, = self.direction(a, e, inc, Omega, omega, f, epoch)        
        
        # Compute the log likelihood by element from the predicted direction and mixture model parameters
        # Shape is [elt_batch_size, 3]
        log_like, log_like_wtd, hits, row_lengths_close = self.score(u_pred, num_hits=num_hits, R=R)
        # Get the mean log like
        # log_like_mean = tf.divide(log_like, tf.cast(x=row_lengths_close, dtype=dtype))
        # Re-scale to the original number of rows; this way doesn't shrink as resolution changes
        # log_like_scaled = tf.multiply(log_like_mean, self.row_lengths_float)
        
        # Stack score outputs. Shape is [batch_size, 2,]
        score_outputs = keras.backend.stack([log_like, hits])

        # Add the loss function - the NEGATIVE of the log likelihood weighted by each element's weight
        elt_weight = self.get_active_weight()
        # mean log_like by element, weighted by elt_weight
        loss_by_elt_log_like = tf.multiply(elt_weight, log_like)
        loss_by_elt_log_like_wtd = tf.multiply(elt_weight , log_like_wtd)
        # Take negative b/c TensorFlow minimizes the loss function, and we want to maximize the log likelihood
        loss_log_like = -tf.reduce_sum(loss_by_elt_log_like)
        loss_log_like_wtd = -tf.reduce_sum(loss_by_elt_log_like_wtd)
        # Create a combined loss with two terms
        loss_factor_log_like = 0.90
        loss_factor_log_like_wtd = 0.10
        # Add these losses to the model
        self.add_loss(loss_log_like*loss_factor_log_like)
        self.add_loss(loss_log_like_wtd*loss_factor_log_like_wtd)

        # Wrap outputs
        outputs = (score_outputs, orbital_elements, mixture_parameters)
        
        return outputs

    def predict_position(self):
        """Predict position q and velocity v"""
        # Extract the candidate elements and mixture parameters
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=None)

        # Ragged tensor of predicted position and velocity.  Shape is [batch_size, num_obs, 3,]
        q_pred, v_pred = self.position(a, e, inc, Omega, omega, f, epoch)        

        return (q_pred, v_pred)

    def predict_direction(self):
        """Predict direction u and displacement r"""
        # Extract the candidate elements and mixture parameters
        a, e, inc, Omega, omega, f, epoch, = self.candidate_elements(inputs=None)

        # Ragged tensor of predicted directions.  Shape is [batch_size, num_obs, 3,]
        u_pred, r_pred, q_ast, = self.direction(a, e, inc, Omega, omega, f, epoch)        

        return (u_pred, r_pred)

    # *********************************************************************************************
    # Methods to calculate outputs, log likelihood, loss
    # *********************************************************************************************

    def calc_outputs(self):
        """Calculate the outputs; no input required (uses cached x_eval)"""
        return self(self.x_eval)
    
    def calc_log_like(self):
        """Calculate the log likelihood as tensor of shape [batch_size,]"""
        score_outputs, orbital_elements, mixture_params = self.calc_outputs()
        log_like, hits = score_outputs
        log_like_mean = tf.reduce_mean(log_like).numpy()
        return log_like, log_like_mean

    def calc_loss(self):
        """Calculate loss function with current inputs"""
        return self.evaluate(self.x_eval, verbose=0)

    def traj_err(self, elts0, elts1):
        """Calculate difference in trajectories from two sets of orbital elements"""
        return traj_diff(elts0, elts1, self.model_pos)

    def get_orbital_elts(self):
        """Extract the current orbital elements as Numpy arrays"""
        pass

    
    def get_mixture_params(self) -> Tuple[np.ndarray, np.ndarray]:
        """Extract the current mixture parameters as Numpy arrays"""
        num_hits, R = self.mixture_parameters(inputs=None)
        return num_hits.numpy(), R.numpy()

    def get_thresh_deg(self) -> np.ndarray:
        """Extract the threshold in degrees as a Numpy array"""
        return self.score.get_thresh_deg()

    def get_H(self) -> np.ndarray:
        """Extract the brightness parameter H as a Numpy array"""
        return self.magnitude.get_H().numpy()

    # *********************************************************************************************
    # Element weights; equivalent to independent learning rates for each element in the batch
    # *********************************************************************************************

    def get_active_weight(self):
        """Get the currently active weight"""
        weight = {
            'joint': self.weight_joint,
            'element': self.weight_element,
            'mixture': self.weight_mixture,
        }[self.training_mode]
        return weight

    def set_active_weight(self, weight):
        """Set the currently active weight"""
        weight_active = {
            'joint': self.weight_joint,
            'element': self.weight_element,
            'mixture': self.weight_mixture,
        }[self.training_mode]
        weight_active.assign(weight)
        # update mean active weight and effective learning rate
        self.active_weight_mean = tf.reduce_mean(weight).numpy()
        self.effective_learning_rate = self.active_weight_mean * self.learning_rate

    def reset_active_weight(self):
        """Reset the active weight to all 1s"""
        weight_ones = np.ones(self.batch_size, dtype=dtype_np)
        self.set_active_weight(weight_ones)

    # *********************************************************************************************
    # Methods to change model state in training
    # *********************************************************************************************

    def recompile(self):
        """Recompile this model with its current optimizer"""
        # Note: important not to name this method compile, that breaks relationship with tf.keras.Model
        # tf.keras.Model.compile(self, optimizer=self.optimizer)
        self.compile(optimizer=self.optimizer)

    def set_learning_rate(self, learning_rate):
        """Set the learning rate"""
        self.learning_rate = learning_rate
        
    def set_clipnorm(self, clipnorm):
        """Set the clipnorm parameter for gradient clipping"""
        self.clipnorm = clipnorm
        
    def adjust_learning_rate(self, lr_factor, verbose: bool=True):
        """Adjust the learning rate and recompile the model"""
        if verbose:
            print(f'Changing learning rate by factor {lr_factor:8.6f} from '
                  f'{self.learning_rate:8.3e} to {self.learning_rate*lr_factor:8.3e}.')
        self.learning_rate = self.learning_rate * lr_factor
        # Recompile with the new learning rate
        self.optimizer = make_opt(optimizer_type=self.optimizer_type, learning_rate=self.learning_rate, 
                                  clipnorm=self.clipnorm, clipvalue=None)
        self.recompile()

    def set_R_deg_max(self, R_deg_max: float):
        """Adjust resolution R_deg to at most R_deg_max"""
        # If R_deg_max was a scalar, promote it to a full numpy array
        if isinstance(R_deg_max, float):
            R_deg_max = np.full(shape=self.batch_size, fill_value=R_deg_max, dtype=dtype_np)            
        # Apply the new R_max value to the mixture parameters layer
        self.mixture_parameters.set_R_deg_max(R_deg_max)

    def set_thresh_deg(self, thresh_deg: np.ndarray):
        """Adjust the threshold for which observations are included in the score function"""
        # If thresh_deg was a scalar, promote it to a full numpy array
        if isinstance(thresh_deg, float):
            thresh_deg = np.full(shape=self.batch_size, fill_value=thresh_deg, dtype=dtype_np)
        # Update the thresh_deg member; this is an array
        self.thresh_deg = thresh_deg
        # Apply this update to the score layer and mixture parameters layer
        self.score.set_thresh_deg(self.thresh_deg)
        self.mixture_parameters.set_thresh_deg(self.thresh_deg)

    def set_thresh_deg_max(self, thresh_deg_max: float):
        """Adjust thresh_deg to at most thresh_deg_max"""
        # If thresh_deg_max was a scalar, promote it to a full numpy array
        if isinstance(thresh_deg_max, float):
            thresh_deg_max = np.full(shape=self.batch_size, fill_value=thresh_deg_max, dtype=dtype_np)
        # Update the thresh_deg member; this is an array
        self.thresh_deg = np.minimum(self.thresh_deg, thresh_deg_max)
        # Apply the new max threshold to the score layer
        self.score.set_thresh_deg_max(thresh_deg_max)
        # Apply the updated value to the mixture parameters layer
        self.mixture_parameters.set_thresh_deg(self.thresh_deg)

    def recalibrate(self):
        """Recalibrate Kepler model rebound integration at current orbital elements"""
        # Element IDs as a Numpy array
        element_ids = self.element_id.numpy()

        # Current candidates in DataFrame format for calc_ast_data
        cols_elt = ['element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch']
        elts = self.candidates_df()[cols_elt]

        # Get unique times and inverse indices
        mjd, unq_idx = np.unique(self.ts.numpy(), return_inverse=True)
        mjd0 = np.min(mjd) - 1.0
        mjd1 = np.max(mjd) + 1.0

        # Calculate positions in this date range, sampled daily
        df_ast_daily, df_earth_daily, df_sun_daily = calc_ast_data(elts=elts, mjd0=mjd0, mjd1=mjd1, element_id=element_ids)

        # Spline positions at ztf times
        df_ast, df_earth, df_sun = spline_ast_vec_df(df_ast=df_ast_daily, df_earth=df_earth_daily, df_sun=df_sun_daily, 
                                                    mjd=mjd, include_elts=False)

        # Iterate over the elements
        for element_id in element_ids:
            # Mask for observations of this element
            mask_ztf = (self.ztf_elt.element_id == element_id)
            ztf_i = self.ztf_elt[mask_ztf]
            # Mask for asteroid position of this element
            mask_df = (df_ast.element_id==element_id)
            df_ast_i = df_ast[mask_df]
            # Index into df_ast_i on selected rows
            df_idx = unq_idx[mask_ztf]
            # Overwrite relevant slices of q_ast and v_ast with the newly integrated q, v
            self.q_ast[mask_ztf] = df_ast_i[self.cols_q_ast].iloc[df_idx]
            self.v_ast[mask_ztf] = df_ast_i[self.cols_v_ast].iloc[df_idx]            

        # Run calibration with the new q, v
        self.position.calibrate(elts=elts, q_ast=self.q_ast, v_ast=self.v_ast)

    # *********************************************************************************************
    # Freeze and thaw layers relating to orbital elements and mixture parameters
    # *********************************************************************************************

    def freeze_candidate_elements(self):
        """Make the candidate orbital elements not trainable; only update mixture parameters"""
        self.candidate_elements.trainable = False
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'mixture'

    def thaw_candidate_elements(self):
        """Make the candidate orbital elements trainable (unfrozen)"""
        self.candidate_elements.trainable = True
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'joint' if self.mixture_parameters.trainable else 'element'

    def freeze_mixture_parameters(self):
        """Make the mixture parameters not trainable; only update candidate oribtal elements"""
        self.mixture_parameters.trainable = False
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'element'

    def thaw_mixture_parameters(self):
        """Make the mixture parameters trainable (unfrozen)"""
        self.mixture_parameters.trainable = True
        self.recompile()
        self.update_early_stop()
        self.training_mode = 'joint' if self.candidate_elements.trainable else 'mixture'

    def freeze_magnitude(self):
        """Make the magnitude layer not trainable"""
        self.magnitude.trainable = False
        self.recompile()
        self.update_early_stop()

    def thaw_magnitude(self):
        """Make the magnitude layer trainable (unfrozen)"""
        self.magnitude.trainable = True
        self.recompile()
        self.update_early_stop()

    def freeze_score(self):
        """Make the score layer (e.g. thresh_deg) not trainable"""
        self.score.trainable = False
        self.recompile()
        self.update_early_stop()

    def thaw_score(self):
        """Make the score layer (e.g. thresh_deg) trainable (unfrozen)"""
        self.score.trainable = True
        self.recompile()
        self.update_early_stop()

    # *********************************************************************************************
    # Adaptive training; save weights and history at episode end
    # *********************************************************************************************

    def save_weights(self, is_update: bool=False):
        """Save the current weights, log likelihood and training history"""
        # Generate the outputs
        score_outputs, orbital_elements, mixture_params = self.calc_outputs()
        log_like, hits = score_outputs

        # is_update is true when we are updating weights that have already been written
        # this is passed when we need to restore weights that got worse
        is_new: bool = ~is_update
        
        # Calculate mean log likelihood, mean hits, and loss; save them to cache        
        self.log_like = log_like.numpy()
        self.hits = hits.numpy()
        self.log_like_mean = tf.reduce_mean(log_like).numpy()
        self.hits_mean = tf.reduce_mean(hits).numpy()
        self.loss = self.calc_loss()

        # Write history of weights, elements, and mixture parameters
        if is_new:
            self.candidate_elements_hist.append(self.candidate_elements.get_weights())
            self.mixture_parameters_hist.append(self.mixture_parameters.get_weights())
        else:
            self.candidate_elements_hist[-1] = self.candidate_elements.get_weights()
            self.mixture_parameters_hist[-1] = self.mixture_parameters.get_weights()
        
        # Write history of log likelihood and loss
        if is_new:
            self.log_like_hist.append(self.log_like)
            self.hits_hist.append(self.hits)
            self.log_like_mean_hist.append(self.log_like_mean)
            self.hits_mean_hist.append(self.hits_mean)
            self.loss_hist.append(self.loss)
        else:
            self.log_like_hist[-1] = self.log_like
            self.hits_hist[-1] = self.hits
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

    def save_train_hist(self, score_outputs, orbital_elements, mixture_params):
        """Save training history, both summary and by element"""

        # Extract score outputs
        log_like = score_outputs[0]
        hits = score_outputs[1]

        # Extract mixture parameters
        num_hits = mixture_params[0]
        R = mixture_params[1]
        R_deg = dist2deg(R)
        R_sec = dist2sec(R)
        log_R = np.log(R)

        # Extract thresh_deg
        thresh_deg = self.get_thresh_deg()

        # Extract the brightness
        H = self.get_H()

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

            # Orbital elements
            'a': orbital_elements[0].numpy(),
            'e': orbital_elements[1].numpy(),
            'inc': orbital_elements[2].numpy(),
            'Omega': orbital_elements[3].numpy(),
            'omega': orbital_elements[4].numpy(),
            'f': orbital_elements[5].numpy(),
            
            # Mixture parameters
            'num_hits': num_hits,
            'R': R,
            'R_deg': R_deg,
            'R_sec': R_sec,
            'log_R': log_R,

            # Threshold
            'thresh_deg': thresh_deg,
            'thresh_s': deg2dist(thresh_deg),

            # Brightness H
            'H': H,

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
        # train_hist_elt_cur.set_index('key', inplace=True)
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
            'loss': self.loss,
            'learning_rate': self.learning_rate,

            # Log likelihood summarized over the elements
            'log_like_mean': [np.mean(log_like)],
            'log_like_med': [np.median(log_like)],
            'log_like_std': [np.std(log_like)],
            'log_like_min': [np.min(log_like)],
            'log_like_max': [np.max(log_like)],

            # Worst and best element in this batch
            'log_like_argmin': [np.argmin(log_like)],
            'log_like_argmax': [np.argmax(log_like)],

            # Hits summarized over the elements
            'hits_mean': [np.mean(hits)],
            'hits_med': [np.median(hits)],
            'hits_std': [np.std(hits)],
            'hits_min': [np.min(hits)],
            'hits_max': [np.max(hits)],

            # Summary of mixture parameter num_hits
            'num_hits_mean': [np.mean(num_hits)],
            'num_hits_std': [np.std(num_hits)],
            'num_hits_min': [np.min(num_hits)],
            'num_hits_max': [np.max(num_hits)],

            # Summary of mixture parameter R
            'R_deg_mean': [np.mean(R_deg)],
            'R_deg_std': [np.std(R_deg)],
            'R_deg_min': [np.min(R_deg)],
            'R_deg_max': [np.max(R_deg)],

            # Summary of mixture parameter R
            'log_R_mean': [np.mean(log_R)],
            'log_R_std': [np.std(log_R)],
            'log_R_min': [np.min(log_R)],
            'log_R_max': [np.max(log_R)],

            # Summary of threshold_deg
            'thresh_deg_mean' : [np.mean(thresh_deg)],
            'thresh_deg_std' : [np.std(thresh_deg)],
            'thresh_deg_min' : [np.min(thresh_deg)],
            'thresh_deg_max' : [np.max(thresh_deg)],

            # Extra data required to rebuild model state
            'candidate_elements_trainable': self.candidate_elements.trainable,
            'mixture_parameters_trainable': self.mixture_parameters.trainable,

        }
        train_hist_summary_cur = pd.DataFrame(hist_sum_dict, index=hist_sum_dict['key'])
        self.train_hist_summary = pd.concat([self.train_hist_summary, train_hist_summary_cur])

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

        # Get old and new log_like, candidate_elements and mixture_parameters
        log_like_old, log_like_new = self.log_like_hist[n-2:n]
        log_like_mean_old, log_like_mean_new = self.log_like_mean_hist[n-2:n]
        candidate_elements_old, candidate_elements_new = self.candidate_elements_hist[n-2:n]
        mixture_parameters_old, mixture_parameters_new = self.mixture_parameters_hist[n-2:n]

        # Test which elements have gotten worse (usually this should be false)
        is_worse = tf.math.less(log_like_new, log_like_old)

        # If none of the elements have gotten worse, terminate early
        if not tf.math.reduce_any(is_worse):
            return

        # If we get here, at least one candidate got worse.  Want to restore it.
        # Calculate the best weights on the candidate elements and mixture parameters
        log_like_best = tf.math.maximum(x=log_like_old, y=log_like_new)
        candidate_elements_best = tf.where(condition=is_worse, x=candidate_elements_old, y=candidate_elements_new)
        mixture_parameters_best = tf.where(condition=is_worse, x=mixture_parameters_old, y=mixture_parameters_new)

        # Apply the best weights to the trainable layers if they are not frozen
        if self.candidate_elements.trainable:
            self.candidate_elements.set_weights(candidate_elements_best)
        if self.mixture_parameters.trainable:
            self.mixture_parameters.set_weights(mixture_parameters_best)

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
       
    def update_early_stop(self):
        """Update early stopping monitor"""
        self.early_stop = tf.keras.callbacks.EarlyStopping(
            monitor='loss', 
            patience=0, 
            baseline=self.calc_loss(), 
            min_delta=0.0, 
            restore_best_weights=False)
    
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

        # Geometric mean of resolution
        R_deg_geomean = dist2deg(np.exp(self.train_hist_summary.log_R_mean.values[-1]))

        # Update early_stop
        self.update_early_stop()

        # Recalibrate the model if the orbital elements are not frozen
        if self.candidate_elements.trainable:
            self.recalibrate()
            # Also need to save weights; loss is going to get worse immediately after a recalibration
            self.save_weights()

        # Status message
        if verbose > 0:
            # print(f'Epoch {self.current_epoch:4}. Elapsed time = {self.episode_time:0.0f} sec')
            # print(f'Total Log Likelihood: {log_like_total:8.2f}')
            print(f'Geom Mean Resolution: {R_deg_geomean:8.6f} degrees ({R_deg_geomean*3600:6.1f} arc seconds)')
            print(f'Mean Hits          :  {self.hits_mean:8.2f}')
            print(f'Mean Log Likelihood:  {self.log_like_mean:8.2f}')

    # *********************************************************************************************
    # Main adaptive search routine; called by external consumers
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

    def search_adaptive(self, 
                        max_batches: int = 1000, 
                        batches_per_epoch: int = 100, 
                        epochs_per_episode: int = 5,
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
        if learning_rate is not None and learning_rate != self.learning_rate:
            self.learning_rate = learning_rate
            self.recompile()

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
        max_episodes = (max_batches*2) // (batches_per_epoch * epochs_per_episode)

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
        elif self.learning_rate <= self.effective_learning_rate:
            print_header(f'Terminating: Effective Learning Rate '
                         f'{self.effective_learning_rate:8.3e} <= minimum {min_learning_rate:8.3e}.')

        # Save training progress to disk
        if save_at_end:
            self.save_state(verbose=True)

    def sieving_round(self, thresh_deg_max: float, batches_elt: int=5000, batches_mix: int=1000):
        """
        One round of sieving
        INPUTS:
            thresh_deg_max: Maximum of threshold parameter
            batches_elt:    Number of batches to train orbital elements
            batches_mix:    Number of batches to train just mixture
        """

        # Adaptive search parameters
        batches_per_epoch = 100
        epochs_per_episode = 5
        learning_rate = 2.0**-15
        reset_active_weight = True
        save_at_end = False
        verbose = 1

        # Set threshold
        self.set_thresh_deg_max(thresh_deg_max)

        # Train orbital elements and mixture parameters
        self.thaw_candidate_elements()
        self.thaw_mixture_parameters()
        self.search_adaptive(
            max_batches=self.current_batch+batches_elt, 
            batches_per_epoch=batches_per_epoch,
            epochs_per_episode=epochs_per_episode,
            learning_rate=learning_rate,
            reset_active_weight=reset_active_weight,
            save_at_end=save_at_end,
            verbose=verbose)

        # Train with frozen orbital elements
        self.freeze_candidate_elements()
        self.search_adaptive(
            max_batches=self.current_batch+batches_mix, 
            batches_per_epoch=batches_per_epoch,
            epochs_per_episode=epochs_per_episode,
            learning_rate=learning_rate,
            reset_active_weight=reset_active_weight,
            save_at_end=save_at_end,
            verbose=verbose)

        # Report progress and save
        self.report()
        self.save_state()

    def sieve(self):
        """Run a predetermined sequence of adaptive searches"""

        # Training rounds on schedule
        training_schedule = [
            (2.00, 0, 2000,),
            (2.00, 5000, 1000,),
            (1.75, 5000, 1000,),
            (1.50, 5000, 1000,),
            (1.25, 5000, 1000,),
            (1.00, 5000, 1000,),
            (0.80, 5000, 1000,),
            (0.70, 5000, 1000,),
            (0.60, 5000, 1000,),
            (0.50, 5000, 1000,),
        ]

        # Iterate over training schedule
        for thresh_deg_max, batches_elt, batches_mix in training_schedule:
            self.sieving_round(thresh_deg_max=thresh_deg_max, batches_elt=batches_elt, batches_mix=batches_mix)

    # *********************************************************************************************
    # Output candidate elements as DataFrame; save and load model state (including training history)
    # *********************************************************************************************

    def candidates_df(self):
        """The current candidates as a DataFrame."""
        # Generate the outputs
        score_outputs, orbital_elements, mixture_params = self.calc_outputs()
        # Unpack score_outputs
        log_like, hits = score_outputs
        # Extract thresh_deg from score layer
        thresh_deg = self.get_thresh_deg()
        # Extract the brightness H from magnitude layer
        H = self.get_H()

        # Build DataFrame of orbital elements
        elts = elts_np2df(orbital_elements.numpy().T)
        
        # Add column with the element_id
        elts.insert(loc=0, column='element_id', value=self.element_id.numpy())
       
        # Add columns for the mixture parameters
        elts['num_hits'] = mixture_params[0].numpy()
        elts['R'] = mixture_params[1].numpy()
        elts['R_deg'] = dist2deg(elts.R)
        elts['R_sec'] = 3600.0 * elts.R_deg
        elts['thresh_s'] = deg2dist(thresh_deg)
        elts['thresh_deg'] = thresh_deg
        elts['thresh_sec'] = 3600.0 * elts.thresh_deg

        # Add columns with log likelihood and hits
        elts['log_like'] = log_like.numpy()
        elts['hits'] = hits.numpy()

        # The brightness parameter H
        elts['H'] = H.numpy()

        # Add columns with the weights in the three modes
        elts['weight_joint'] = self.weight_joint.numpy()
        elts['weight_element'] = self.weight_element.numpy()
        elts['weight_mixture'] = self.weight_mixture.numpy()

        return elts

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
        # Save elts_fit and elts_near_ast if populated
        if self.elts_fit is not None:
            self.elts_fit.to_hdf(self.file_path, key='elts_fit', mode='a')
        if self.elts_near_ast is not None:
            self.elts_near_ast.to_hdf(self.file_path, key='elts_near_ast', mode='a')

    def load(self, verbose: bool=False):
        """Load model state from disk: candidate elements and training history"""
        try:
            # Load DataFrames for candidate elements and auxiliary data
            elts = pd.read_hdf(self.file_path, key='elts')
            self.train_hist_summary = pd.read_hdf(self.file_path, key='train_hist_summary')
            self.train_hist_elt = pd.read_hdf(self.file_path, key='train_hist_elt')
            # Set the threshold for scoring
            self.thresh_deg = elts['thresh_deg']
            # Load nearest asteroid information if available
            try:
                self.elts_fit = pd.read_hdf(self.file_path, key='elts_fit')
                self.elts_near_ast = pd.read_hdf(self.file_path, key='elts_near_ast')
            except KeyError:
                pass
        except FileNotFoundError:
            if verbose:
                print(f'Unable to find {self.file_path}.')
            return

        # Status message
        if verbose:
            print(f'Loaded candidate elements and training history from {self.file_path}.')

        # Regenerate candidate_elements layer of this model
        self.candidate_elements = CandidateElements(elts=elts, name='candidate_elements')
        # Regenerate mixture_parameters layer of this model
        self.mixture_parameters = MixtureParameters(elts=elts, name='mixture_parameters')
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
    # Plot training results: bar or time series charts
    # *********************************************************************************************

    def plot_bar(self, att_name: str = 'log_like', sorted: bool = True, episode=None):
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
        ax.bar(x=hist.element_num, height=values_plot, color='blue')
        # ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/{att_name}_bar.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    def plot_hist(self, att_name: str = 'log_like'):
        """
        Learning curve chart (progress vs. episode) for one of the score attributes log_like or hits.
        INPUTS:
            score_att: Name of the attribute to plot.  One of 'log_like', 'hits', 'R_deg'
        """
        # Alias the training history
        hist = self.train_hist_summary

        # Extract score_mean, score_std, score_min, score_max from hist
        score_mean = hist[f'{att_name}_mean'].values
        score_std = hist[f'{att_name}_std'].values
        score_min = hist[f'{att_name}_min'].values
        score_max = hist[f'{att_name}_max'].values
        score_lo = score_mean - score_std
        score_hi = score_mean + score_std

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

        # Plot total log likelihood over training
        fig, ax = plt.subplots()
        ax.set_title(f'Training Progress: {att_title} by Element')
        ax.set_xlabel('Batch Trained')
        ax.set_ylabel(f'{att_title}')
        # Plot mean +/- 1 SD
        ax.plot(hist.batch, score_mean, color=color_mean, label='Mean')
        ax.plot(hist.batch, score_lo, color=color_lo, label='Mean -1 SD')
        ax.plot(hist.batch, score_hi, color=color_hi, label='Mean +1 SD')
        # Plot min and max
        ax.plot(hist.batch, score_min, color=color_min, label=min_label)
        ax.plot(hist.batch, score_max, color=color_max, label=max_label)
        # Legend etc
        ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/{att_name}_hist.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    # *********************************************************************************************
    # Find asteroid with nearest orbital elements; calculate position error vs. known elements
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

        # Table of functions to find nearest asteroid
        nearest_ast_elt_tbl = {
            'cart' : nearest_ast_elt_cart,
            'cov' : nearest_ast_elt_cov,
        }
        nearest_ast_elt_func = nearest_ast_elt_tbl[search_type]

        # Search for nearest asteroids to the fitted elements
        self.elts_near_ast = nearest_ast_elt_func(self.elts_fit)

        # If we ran a Cartesian search, add an extra column with the Q_norm
        if search_type == 'cart':
            q_norm = elt_q_norm(elts=self.elts_fit, ast_num=self.elts_fit.nearest_ast_num)
            self.elts_fit['nearest_ast_q_norm'] = q_norm

        # Return the inputs and nearest asteroids
        return self.elts_fit, self.elts_near_ast

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

    def plot_elt_error(self, elt_name: str, elts_true=None, is_log: bool=True, elt_num: int=None):
        """
        Plot (log) error of selected orbital element vs. known elements
        INPUTS:
            elt_name:  Name of orbital element; one of 'a', 'e', 'inc', 'Omega', 'omega', 'f'
            elts_true: DataFrame of true orbital elements; defaults to self.near_ast
            is_log:    Flag; whether to plot log(err) (True) or err (False)
            elt_num:   Specific element_num to plot; plot just this one rather than mean, std, etc.
        """
        # If elts_true wasn't specified, it defaults to nearest asteroid elements
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
        else:
            err_selected = err[:, elt_num]

        # Plot selected error
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel('Episode')
        ax.set_ylabel(ylabel)
        if is_summary:
            ax.plot(err_mean, color=color_mean, label='mean')
            ax.plot(err_mean - err_std, color=color_lo, label='mean -1 SD')
            ax.plot(err_mean + err_std, color=color_hi, label='mean +1 SD')
            ax.plot(err_min, color=color_min, label='min')
            ax.plot(err_max, color=color_max, label='max')
        else:
            ax.plot(err_selected, color='blue', label=f'element_num {element_num}')
        ax.legend()
        ax.grid()
        # fig.savefig('../figs/training/log_like_bar.png', bbox_inches='tight')
        plt.show()
        return fig, ax

    def plot_q_error(self, elts_true=None, is_log: bool=False, use_near_ast_dist: bool=True):
        """
        Bar chart of position error by orbital element.
        INPUTS:
            elts_true: DataFrame of true orbital elements
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
            calc_method = 'Observation Dates'
        else:
            hist_err, q_err = self.calc_error(elts_true)
            calc_method = 'Near Ast. Dates'
        
        # Get element numbers and mean error over the whole batch
        element_num = np.arange(self.batch_size, dtype=np.int32)
        mean_err = np.exp(np.mean(np.log(q_err+2**-40))) if is_log else np.mean(q_err)

        # Plot selected error
        fig, ax = plt.subplots()
        ax.set_title(f'Mean Position Error vs. Elements on {calc_method}')
        ax.set_xlabel('Element')
        ax.set_ylabel('Mean Position Error in AU')
        ax.bar(x=element_num, height=q_err, label='elt', color='blue')
        ax.axhline(y=mean_err, label='mean', color='red')
        if is_log:
            ax.set_yscale('log', basey=2)
        ax.grid()
        ax.legend()
        # fig.savefig('../figs/training/q_error.png', bbox_inches='tight')
        plt.show()
        return fig, ax
        

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
        score_outputs, candidate_elements, mixture_parameters = self.calc_outputs()
        # Unpacke score_outputs
        log_like, hits = score_outputs

        # Summarize log likelihood
        log_like = log_like.numpy()
        log_like_mean = np.mean(log_like)
        log_like_std = np.std(log_like)
        log_like_min = np.min(log_like)
        log_like_max = np.max(log_like)

        # Summarize hits
        hits = hits.numpy()
        hits_mean = np.mean(hits)
        hits_std = np.std(hits)
        hits_min = np.min(hits)
        hits_max = np.max(hits)

        # Summarize resolution
        # num_hits = mixture_parameters[0].numpy()
        R_deg = dist2deg(mixture_parameters[1].numpy())
        R_deg_mean = np.mean(R_deg)
        R_deg_std = np.std(R_deg)
        R_deg_min = np.min(R_deg)
        R_deg_max = np.max(R_deg)

        # Threshold in degrees
        thresh_deg = self.get_thresh_deg()
        thresh_deg_mean = np.mean(thresh_deg)
        thresh_deg_std = np.std(thresh_deg)
        thresh_deg_min = np.min(thresh_deg)
        thresh_deg_max = np.max(thresh_deg)

        # Report on log likelihood and resolution
        print(f'     \  log_like :  hits  :  R_deg    :    R_sec : thresh_deg')
        print(f'Mean : {log_like_mean:8.2f}  : {hits_mean:6.2f} :  {R_deg_mean:8.6f} : {R_deg_mean*3600:8.2f} : {thresh_deg_mean:8.6f}')
        print(f'Std  : {log_like_std:8.2f}  : {hits_std:6.2f} :  {R_deg_std:8.6f} : {R_deg_std*3600:8.2f} : {thresh_deg_std:8.6f}')
        print(f'Min  : {log_like_min:8.2f}  : {hits_min:6.2f} :  {R_deg_min:8.6f} : {R_deg_min*3600:8.2f} : {thresh_deg_min:8.6f}')
        print(f'Max  : {log_like_max:8.2f}  : {hits_max:6.2f} :  {R_deg_max:8.6f} : {R_deg_max*3600:8.2f} : {thresh_deg_max:8.6f}')
        print(f'Trained for {self.current_batch} batches over {self.current_epoch} epochs and {self.current_episode} episodes.')

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
