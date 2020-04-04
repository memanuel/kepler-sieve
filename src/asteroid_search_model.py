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

# MSE imports
from asteroid_model import AsteroidDirection, make_model_ast_pos
from candidate_element import elts_np2df
from asteroid_search_layers import CandidateElements, TrajectoryScore
from asteroid_integrate import calc_ast_pos
from candidate_element import perturb_elts
from ztf_element import ztf_elt_hash
from asteroid_search_report import traj_diff
from astro_utils import deg2dist, dist2deg, dist2sec
from tf_utils import tf_quiet, Identity
from utils import print_header

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

# Save directory for candidate elements
# savedir = '../data/candidate_elt'

# ********************************************************************************************************************* 
def make_adam_opt(learning_rate=1.0E-4, clipnorm=1.0, clipvalue=None):
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
def candidate_elt_hash(elts: pd.DataFrame, thresh_deg: float):
    """
    Load or generate a ZTF batch with all ZTF observations within a threshold of the given elements
    INPUTS:
        elts:       Dataframe including element_id; 6 orbital elements; and 2 mixture params h, R
        thresh_deg: Threshold in degrees; only observations this close to elements are returned
    OUTPUTS:
        hash_id:    Unique ID for these inputs
    """
    # Columns of the Dataframe to hash
    cols_to_hash = ['element_id', 'a', 'e', 'inc', 'Omega', 'omega', 'f', 'epoch', 'h', 'R']
    # Tuple of int64; one per orbital element candidate
    hash_df = tuple((pd.util.hash_pandas_object(elts[cols_to_hash])).values)
    # Combine the element hash tuple with the threshold
    thresh_int = int(thresh_deg*2**48)
    hash_id = hash(hash_df + (thresh_int, ))

    return hash_id
    
# ********************************************************************************************************************* 
# Custom model for Asteroid Search
# ********************************************************************************************************************* 

class AsteroidSearchModel(tf.keras.Model):
    """
    Custom keras model that searches for orbital elements and mixture model parameters
    consistent with a batch of observations generated from the initially guessed candidate elements.
    """

    def __init__(self, elts: pd.DataFrame, ztf_elt: pd.DataFrame, 
                 site_name: str='geocenter', thresh_deg: float = 1.0, 
                 learning_rate: float = 2.0**-13, clipnorm: float = 1.0,
                 **kwargs):
        """
        INPUTS:
            elts:       DataFrame with initial guess for orbital elements.
                        Columns: element_id, a, e, inc, Omega, omega, f, epoch
                        Output of asteroid_elts, perturb_elts or random_elts
            ztf_elt:    DataFrame with ZTF observations within thresh_deg degrees of
                        of the orbits predicted by these elements.
                        Output of make_ztf_batch or load_ztf_batch
            site_name:  Used for topos adjustment, e.g. 'geocenter' or 'palomar'
            h:          Initial value of hit probability in mixture model
            R_deg:      Initial value of resolution parameter (in degrees) in mixture model
            learning_rate: Initial value of learning rate
            clipnorm:   Initila value of clipnorm for gradient clipping
        """
        # Initialize tf.keras.Model
        super(AsteroidSearchModel, self).__init__(**kwargs)
        
        # ***************************************************
        # Description of input data 
        # ***************************************************

        # Batch size comes from elts
        self.batch_size = elts.shape[0]

        # The element_id for the elements in this batch
        self.elts_element_id = keras.backend.constant(value=elts.element_id, dtype=tf.int32)

        # Numpy array and tensor of observation times; flat, shape (data_size,)
        ts_np = ztf_elt.mjd.values.astype(dtype_np)

        # Get observation count per element
        row_lengths_np = ztf_elt.element_id.groupby(ztf_elt.element_id).count()
        self.row_lengths = keras.backend.constant(value=row_lengths_np, shape=(self.batch_size,), dtype=tf.int32)

        # Shape of the observed trajectories
        self.data_size: int = ztf_elt.shape[0]
        self.traj_shape = (self.data_size, space_dims)

        # Observed directions; extract from ztf_elt DataFrame
        cols_u_obs = ['ux', 'uy', 'uz']
        u_obs_np = ztf_elt[cols_u_obs].values.astype(dtype_np)

        # ***************************************************
        # Layers for candidate elements, asteroid direction and score
        # ***************************************************

        # Set of trainable weights with candidate orbital elements; initialize according to elts
        self.candidate_elements = CandidateElements(elts=elts, name='candidate_elements')

        # The predicted direction; shape is [data_size, 3,]
        self.direction = AsteroidDirection(ts_np=ts_np, row_lengths_np=row_lengths_np, 
                                           site_name=site_name, name='direction')

        # Calibration arrays (flat)
        cols_q_ast = ['qx', 'qy', 'qz']
        cols_v_ast = ['vx', 'vy', 'vz']
        q_ast = ztf_elt[cols_q_ast].values.astype(dtype_np)
        v_ast = ztf_elt[cols_v_ast].values.astype(dtype_np)

        # Run calibration
        self.direction.q_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)

        # Score layer for these observations
        self.score = TrajectoryScore(row_lengths_np=row_lengths_np, u_obs_np=u_obs_np,
                                     thresh_deg=thresh_deg, name='score')

        # Position model - used for comparing trajectories of fitted and known orbital elements
        self.model_pos: keras.Model = make_model_ast_pos(ts_np=ts_np, row_lengths_np=row_lengths_np)

        # ***************************************************
        # Variables for adaptive training and training history
        # ***************************************************

        # Save the learning rate on the model object to facilitate adaptive training
        self.learning_rate: float = learning_rate
        self.clipnorm: float = clipnorm

        # Initialize loss history and total training time
        self.total_training_time: float = 0.0

        # Compile the model with its learning rate and clipnorm
        self.recompile()

        # Tensor of ones to pass as dummy inputs for evaluating one batch
        self.x_eval = tf.ones(shape=self.batch_size, dtype=dtype)

        # Epoch and episode counters
        self.batches_per_epoch: int = 100
        self.episode_length: int = 0
        self.current_episode: int = 0
        self.current_epoch: int = 0
        self.current_batch: int = 0

        # Start training timer
        self.t0: float = time.time()
        self.episode_t0: float = time.time()

        # Cached values of total log likelihood and loss
        self.total_log_like: float = 0.0
        self.loss: float = 0.0

        # Initialize lists with training history

        # Weights, elements and mixture parameters
        self.weights_hist = []

        # Log likelihoods and losses
        self.total_log_like_hist = []
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

        # Save hash IDs for elts and ztf_elts
        self.elts_hash_id: int = candidate_elt_hash(elts=elts, thresh_deg=thresh_deg)
        self.ztf_elt_hash_id: int = ztf_elt_hash(elts=elts, thresh_deg=thresh_deg, near_ast=False)
        # File name for saving training progress
        self.file_path = f'../data/candidate_elt/candidate_elt_{self.elts_hash_id}.h5'
    
    # @tf.function
    def call(self, inputs=None):
        """
        Predict directions from current elements and score them.  
        inputs is a dummmy with no affect on the results; it is there to satisfy keras.Model API
        """

        # Extract the candidate elements and mixture parameters; pass dummy inputs to satisfy keras Layer API
        a, e, inc, Omega, omega, f, epoch, h, lam, R, = self.candidate_elements(inputs=inputs)
        
        # Stack the current orbital elements.  Shape is [batch_size, 7,]
        orbital_elements = keras.backend.stack([a, e, inc, Omega, omega, f, epoch,])

        # Stack mixture model parameters. Shape is [batch_size, 3,]
        mixture_params = keras.backend.stack([h, lam, R])

        # Tensor of predicted directions.  Shape is [data_size, 3,]
        u_pred, r_pred = self.direction(a, e, inc, Omega, omega, f, epoch)        
        
        # Compute the log likelihood by element from the predicted direction and mixture model parameters
        # Shape is [elt_batch_size,]
        log_like = self.score(u_pred, h=h, lam=lam)
        
        # Add the loss function - the NEGATIVE of the log likelihood
        # (Take negative b/c TensorFlow minimizes the loss function, and we want to maximize the log likelihood)
        self.add_loss(-tf.reduce_sum(log_like))
        
        # Wrap outputs
        outputs = (log_like, orbital_elements, mixture_params)
        
        return outputs

    # *************************************************************************
    # Methods to calculate outputs, log likelihood, loss
    # *************************************************************************

    def calc_outputs(self):
        """Calculate the outputs; no input required (uses cached x_eval)"""
        return self(self.x_eval)
    
    def calc_log_like(self):
        """Calculate the log likelihood as tensor of shape [batch_size,]"""
        log_like, orbital_elements, mixture_params = self.calc_outputs()
        total_log_like = tf.reduce_sum(log_like).numpy()
        return log_like, total_log_like

    def calc_loss(self):
        """Calculate loss function with current inputs"""
        return self.evaluate(self.x_eval, verbose=0)

    def traj_err(self, elts0, elts1):
        """Calculate difference in trajectories from two sets of orbital elements"""
        return traj_diff(elts_ast, elts0, elts1)

    # *************************************************************************
    # Methods to change model state in training
    # *************************************************************************

    def recompile(self):
        """Recompile this model with its current learning rate"""
        # Note: important not to name this method compile, that breaks relationship with tf.keras.Model
        optimizer = make_adam_opt(learning_rate=self.learning_rate, clipnorm=self.clipnorm)
        tf.keras.Model.compile(self, optimizer=optimizer)

    def set_learning_rate(self, learning_rate):
        """Set the learning rate"""
        self.learning_rate = learning_rate
        
    def set_clipnorm(self, clipnorm):
        """Set the clipnorm parameter for gradient clipping"""
        self.clipnorm = clipnorm
        
    def adjust_learning_rate(self, lr_factor, verbose: int =1):
        """Adjust the learning rate and recompile the model"""
        if verbose:
            print(f'Changing learning rate by factor {lr_factor:8.6f} from '
                  f'{self.learning_rate:8.3e} to {self.learning_rate*lr_factor:8.3e}.')
        self.learning_rate = self.learning_rate * lr_factor
        # Recompile with the new learning rate
        self.recompile()

    # *************************************************************************
    # Adaptive training; save weights and history at episode end
    # *************************************************************************

    def save_weights(self):
        """Save the current weights, log likelihood and training history"""
        # Generate the outputs
        log_like, orbital_elements, mixture_params = self.calc_outputs()
        # Caclulat total log likelihood and loss; save them to cache
        self.total_log_like = tf.reduce_sum(log_like).numpy()
        self.loss = self.calc_loss()

        # Write history of weights, elements, and mixture parameters
        self.weights_hist.append(self.candidate_elements.get_weights())
        
        # Write history of log likelihood and loss
        self.total_log_like_hist.append(self.total_log_like)
        self.loss_hist.append(self.loss)

        # Write history of training counters and time
        self.episode_hist.append(self.current_episode)
        self.epoch_hist.append(self.current_epoch)
        self.batch_hist.append(self.current_batch)
        self.training_time_hist.append(self.total_training_time)

        # Delegate to save_train_hist to write history DataFrames
        self.save_train_hist(log_like=log_like, orbital_elements=orbital_elements, mixture_params=mixture_params)

    def save_train_hist(self, log_like, orbital_elements, mixture_params):
        """Save training history, both summary and by element"""
        # DataFrame with summary history for this episode
        hist_sum_dict = {
            # Training stage: by episode, epoch, batch and time
            'episode': self.current_episode,
            'epoch': self.current_epoch,
            'batch': self.current_batch,
            'training_time': self.total_training_time,
            # Summary of mixture parameters
            'h_mean': np.mean(mixture_params[0]),
            'h_std': np.std(mixture_params[0]),
            'R_mean': np.mean(mixture_params[1]),
            'R_std': np.std(mixture_params[1]),
            # Loss and learning rate
            'loss': self.loss,
            'learning_rate': self.learning_rate,
            # Log likelihood summarized over the elements
            'log_like_total': self.total_log_like,
            'log_like_mean': np.mean(log_like),
            'log_like_med': np.median(log_like),
            'log_like_std': np.std(log_like),
            'log_like_min': np.min(log_like),
            'log_like_max': np.max(log_like),
            'log_like_argmin': np.argmin(log_like),
        }
        train_hist_summary_cur = pd.DataFrame(hist_sum_dict, index=[self.current_episode])
        self.train_hist_summary = pd.concat([self.train_hist_summary, train_hist_summary_cur])

        # DataFrame with detailed training history of this episode by element
        hist_elt_dict = {
            # Element number and ID
            'element_num': np.arange(self.batch_size), 
            'element_id': self.elts_element_id.numpy(),
            # Training stage: by episode, epoch, batch and time
            'episode': np.full(shape=self.batch_size, fill_value=self.current_episode),
            'epoch': np.full(shape=self.batch_size, fill_value=self.current_epoch),
            'batch': np.full(shape=self.batch_size, fill_value=self.current_batch),
            'training_time': np.full(shape=self.batch_size, fill_value=self.total_training_time),
            # Orbital elements
            'a': orbital_elements[0],
            'e': orbital_elements[1],
            'inc': orbital_elements[2],
            'Omega': orbital_elements[3],
            'omega': orbital_elements[4],
            'f': orbital_elements[5],
            # Mixture parameters
            'h': mixture_params[0],
            'lam': mixture_params[1],
            'R': mixture_params[2],
            'R_deg': dist2deg(mixture_params[2]),
            'R_sec': dist2sec(mixture_params[2]),
            # Log likelihood
            'log_like': log_like
        }
        train_hist_elt_cur = pd.DataFrame(hist_elt_dict)
        self.train_hist_elt = pd.concat([self.train_hist_elt, train_hist_elt_cur])

    # *************************************************************************
    # Adaptive training; update weights and modify model state at episode end
    # *************************************************************************

    def update_weights(self):
        """"Restore the weights for each element to the prior iteration if they got worse"""
        # Quit early if we don't have at least 2 entries on the lists of weights and log likes
        n = len(self.weights_hist)
        if n < 2:
            return

        # # Get weights, log likelihoods and loss from the last 2 iterations
        # weights_old, weights_new = self.weights_list[n-2:n]
        # log_like_old, log_like_new = self.log_likes[n-2:n]
        # loss_old, loss_new = self.loss_hist[n-2:n]

        # # Test whether each element has improved
        # is_better = tf.math.greater(log_like_new, log_like_old)

        # # Calculate the best weights and associated log likelihoods
        # best_weights = tf.where(condition=is_better, x=weights_new, y=weights_old)
        # best_log_like = tf.where(condition=is_better, x=log_like_new, y=log_like_old)

        # # Apply the best weights to the model
        # self.candidate_elements.set_weights(best_weights)
        
        # # Overwrite the weights that got worse so they will be used at most once in an update
        # weights_old2 = self.weights_list[n-3] if n > 2 else weights_old
        # weights_new_update = tf.where(condition=is_better, x=weights_new, y=weights_old2)
        # self.weights[-1] = weights_new_update

        # # Update the log_likes and loss_hist lists to reflect these updates
        # self.log_likes[-1] = best_log_like
        # self.loss_hist[-1] = -tf.reduce_sum(best_log_like)
        
        # # If the overall loss has gotten worse, reduce the learning rate
        # if loss_new > loss_old:
        #     self.adjust_learning_rate(self.lr_factor_dn, verbose=True)
       
    def update_early_stop(self):
        """Update early stopping monitor"""
        self.early_stop = tf.keras.callbacks.EarlyStopping(
            monitor='loss', 
            patience=0, 
            baseline=self.calc_loss(), 
            min_delta=1.0, 
            restore_best_weights=True)
    
    def episode_end(self, hist, verbose: int = 1):
        """Post-processing after one episode of adaptive training"""
        # Update training counters
        self.episode_length = len(hist.epoch)
        self.current_episode += 1
        self.current_epoch += self.episode_length
        self.current_batch += self.episode_length * self.batches_per_epoch

        # Update timers
        self.episode_time = (time.time() - self.episode_t0)
        self.total_training_time += self.episode_time
        self.episode_t0 = time.time()

        # Save weights and apply best to each element
        self.save_weights()
        self.update_weights()

        # Current total log likelihood
        total_log_like = self.total_log_like_hist[-1]

        # Update early_stop
        self.update_early_stop()

        # Status message
        if verbose > 0:
            # print(f'Epoch {self.current_epoch:4}. Elapsed time = {self.episode_time:0.0f} sec')
            print(f'Total Log Likelihood: {total_log_like:8.2f}')        

    # *********************************************************************************************
    # Main adaptive search routine; called by external consumers
    # *********************************************************************************************
        
    def search_adaptive(self, 
                        max_batches: int = 10000, 
                        batches_per_epoch: int = 100, 
                        epochs_per_episode: int = 5,
                        min_learning_rate: float = 2.0**-20,
                        regenerate: bool=False,
                        verbose: int = 1):
        """
        Run asteroid search model adaptively.  
        Start with a high learning rate, gradually reduce it if early stopping triggered.
        INPUTS:
            max_batches: The maximum number of batches processed in the adaptive training.
            batches_per_epoch: The number of batches processed in one epoch of training
            epochs_per_episode: The number of epochs that comprise one episode of training
            min_learning_rate: Minimum for the learning rate; terminate early if LR drops below this
            verbose: Integer controlling verbosity level (passed on to tf.keras.model.fit)
        """
        # Update batches_per_epoch
        self.batches_per_epoch = batches_per_epoch

        # Try to load candidates
        if not regenerate:
            self.load_candidates(verbose=True)

        # Early stopping callback
        self.update_early_stop()
        callbacks = [self.early_stop]

        # Define one epoch as a number of batches        
        samples_per_epoch: int = batches_per_epoch * self.batch_size
        max_epochs: int = max_batches // batches_per_epoch + self.current_epoch
        x_trn = tf.ones(samples_per_epoch, dtype=dtype)

        # Set the learning rate factor
        self.lr_factor_dn: float = 0.5
        # self.lr_factor_up: float = 1.0

        # Run first episode of training
        hist = self.fit(x=x_trn, batch_size=self.batch_size, epochs=epochs_per_episode, 
                        steps_per_epoch=batches_per_epoch, initial_epoch=self.current_epoch,
                        callbacks=callbacks, shuffle=False, verbose=verbose)
        # Episode end processing; includes LR adjustment
        self.episode_end(hist)

        # Continue training until max_epochs have elapsed or learning_rate has dropped too low
        while self.current_epoch < max_epochs and min_learning_rate < self.learning_rate:
            # Status update for this episode
            print(f'\nTraining episode {self.current_episode}: Epoch {self.current_epoch:4}')
            print(f'learning_rate={self.learning_rate:8.3e}, total training time {self.total_training_time:0.0f} sec.')
            # Train for another episode
            hist = self.fit(x=x_trn, batch_size=self.batch_size, epochs=self.current_epoch+epochs_per_episode, 
                            steps_per_epoch=batches_per_epoch, initial_epoch=self.current_epoch,
                            callbacks=callbacks,  shuffle=False, verbose=verbose)
            # Episode end processing; includes LR adjustment
            self.episode_end(hist)

        # Report cause of training termination
        if self.current_epoch == max_epochs:
            print_header(f'Terminating: Reached final epoch {max_epochs}.')
        if self.learning_rate <= min_learning_rate:
            print_header(f'Terminating: Learning rate {self.learning_rate:8.3e} <= minimum {min_learning_rate:8.3e}.')

        # Save training progress to disk
        self.save_state(verbose=True)

    # *********************************************************************************************
    # Output candidate elements as DataFrame; save and load model state (including training history)
    # *********************************************************************************************

    def candidates_df(self):
        """The current candidates as a DataFrame."""
        # Generate the outputs
        log_like, orbital_elements, mixture_params = self.calc_outputs()
        
        # Build DataFrame of orbital elements
        elts = elts_np2df(orbital_elements.numpy().T)
        # print(f'elts.shape={elts.shape}')
        # print(f'self.elts_element_id.numpy().shape = {self.elts_element_id.numpy().shape}')
        
        # Add column with the element_id
        elts.insert(loc=0, column='element_id', value=self.elts_element_id.numpy())
        
        # Add columns for the mixture parameters
        elts['h'] = mixture_params[0].numpy()
        elts['lam'] = mixture_params[1].numpy()
        elts['R'] = mixture_params[2].numpy()
        elts['R_deg'] = dist2deg(elts.R)
        elts['R_sec'] = 3600.0 * elts.R_deg

        # Add column with log likelihood
        elts['log_like'] = log_like.numpy()

        return elts

    def save_state(self, verbose: bool=False):
        """Save model state to disk: candidate elements and training history"""
        # Generate DataFrame with current candidates
        candidate_elts = self.candidates_df()       
        # Status message
        if verbose:
            print(f'Saving candidate elements DataFrame in {self.file_path}.')
        
        # Save to disk
        candidate_elts.to_hdf(self.file_path, key='candidate_elts', mode='w')

        # Also save training history (summary and by element)
        self.train_hist_summary.to_hdf(self.file_path, key='train_hist_summary', mode='a')
        self.train_hist_elt.to_hdf(self.file_path, key='train_hist_elt', mode='a')

    def load(self, verbose: bool=False):
        """Load model state from disk: candidate elements and training history"""
        try:
            # Load DataFrames for candidate elements and auxiliary data
            candidate_elts = pd.read_hdf(self.file_path, key='candidate_elts')
            self.train_hist_summary = pd.read_hdf(self.file_path, key='train_hist_summary')
            self.train_hist_elt = pd.read_hdf(self.file_path, key='train_hist_elt')
        except FileNotFound:
            if verbose:
                print(f'Unable to find {self.file_path}.')
            return

        # Status message
        if verbose:
            print(f'Loaded candidate elements and training history from {self.file_path}.')

        # Regenerate candidate_elements layer of this model
        self.candidate_elements = CandidateElements(elts=candidate_elts, name='candidate_elements')

        # Restore learning rate
        self.learning_rate = self.train_hist_summary.learning_rate.values[-1]

    # *********************************************************************************************
    # Review model and training
    # *********************************************************************************************

    def plot_log_likes(self):
        pass

    def review_members(self):
        """Print diagnostic review of members to console"""
        preview_size = 5

        # Input data description
        print(f'batch_size: {self.batch_size}')
        print(f'\nelts_element_id: shape {self.elts_element_id.shape}')
        print(self.elts_element_id[0:preview_size])
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
        print(f'\nweights_list: list of tensors shape [{len(self.weights_list)}, {len(self.weights_list[0])}, {len(self.weights_list[0][0])}]')
        print([x for x in self.weights_list[0][0][0:preview_size]])

        # Serialization
        print(f'\nelts_hash_id: {self.elts_hash_id}')
        print(f'ztf_elt_hash_id: {self.ztf_elt_hash_id}')

# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
