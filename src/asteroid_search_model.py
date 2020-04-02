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
from asteroid_search_layers import CandidateElements, TrajectoryScore
from asteroid_integrate import calc_ast_pos
from search_score_functions import score_mean, score_var, score_mean_var
from candidate_element import perturb_elts
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
# Custom model for Asteroid Search
# ********************************************************************************************************************* 

class AsteroidSearchModel(tf.keras.Model):
    def __init__(self, elts: pd.DataFrame, ztf_elt: pd.DataFrame, site_name: str='geocenter', 
                 thresh_deg: float = 1.0, h: float = 0.01, R_deg: float = 1.0, 
                 learning_rate: float = 1.0E-4, clipnorm: float = 1.0,
                 **kwargs):
        """
        Functional API model for scoring elements
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
            R_deg:      Initial value of resolution parameter (in degrees) in mixture model
            learning_rate: Initial value of learning rate
            clipnorm:   Initila value of clipnorm for gradient clipping
        """
        # Initialize tf.keras.Model
        super(AsteroidSearchModel, self).__init__(**kwargs)
        
        # Batch size comes from elts
        self.batch_size = elts.shape[0]

        # Numpy array and tensor of observation times; flat, shape (data_size,)
        ts_np = ztf_elt.mjd.values.astype(dtype_np)

        # Get observation count per element
        row_lengths_np = ztf_elt.element_id.groupby(ztf_elt.element_id).count()
        self.row_lengths = keras.backend.constant(value=row_lengths_np, shape=(self.batch_size,), dtype=tf.int32)

        # Shape of the observed trajectories
        self.data_size = ztf_elt.shape[0]
        self.traj_shape = (self.data_size, space_dims)

        # Observed directions; extract from ztf_elt DataFrame
        cols_u_obs = ['ux', 'uy', 'uz']
        u_obs_np = ztf_elt[cols_u_obs].values.astype(dtype_np)

        # Convert resolution from degrees to Cartesian distance
        R_s = deg2dist(R_deg)

        # Set of trainable weights with candidate orbital elements; initialize according to elts
        self.candidate_elements = CandidateElements(elts=elts, h=h, R=R_s, name='candidate_elements')

        # Stack the current orbital elements; shape is [elt_batch_size, 7,]
        self.orbital_elements = tf.keras.layers.Concatenate(axis=-1, name='orbital_elements') 

        # Stack mixture model parameters
        self.mixture_params = tf.keras.layers.Concatenate(axis=-1, name='mixture_params')

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

        # Save the learning rate on the model object to facilitate adaptive training
        self.learning_rate = learning_rate
        self.clipnorm = clipnorm

        # Initialize loss history and total training time
        self.loss_hist: list = []
        self.total_training_time: float = 0.0

        # Compile the model with its learning rate and clipnorm
        self.recompile()

        # Epoch and episode counters
        self.training_episode: int = 1
        self.episode_length: int = 0
        self.current_epoch: int = 0

        # Initialize best_loss and best_weights
        self.x_eval = tf.ones(shape=self.batch_size, dtype=dtype)
        self.best_loss = self.current_loss()
        self.best_weights = self.candidate_elements.get_weights()        
        
    def call(self, inputs=None):
        # Extract the candidate elements and mixture parameters; pass dummy inputs to satisfy keras Layer API
        a, e, inc, Omega, omega, f, epoch, h, lam, R, = self.candidate_elements(inputs=inputs)
        
        # Stack the current orbital elements
        orbital_elements = self.orbital_elements([a, e, inc, Omega, omega, f, epoch,])
        # Stack mixture model parameters
        mixture_params = self.mixture_params([h, lam, R,])

        # Tensor of predicted directions.  Shape is [data_size, 3,]
        u_pred, r_pred = self.direction(a, e, inc, Omega, omega, f, epoch)        
        
        # Compute the log likelihood by element from the predicted direction and mixture model parameters
        # Shape is [elt_batch_size,]
        log_like = self.score(u_pred, h=h, lam=lam)
        
        # Add the loss function - the NEGATIVE of the log likelihood
        # (Take negative b/c TensorFlow minimizes the loss function)
        self.add_loss(-tf.reduce_sum(log_like))
        
        # Wrap outputs
        outputs = (log_like, orbital_elements, mixture_params)
        
        return outputs

    def current_loss(self):
        return self.evaluate(self.x_eval, verbose=0)

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
        """Reduce the learning rate and recompile the model"""
        if verbose:
            print(f'Changing learning rate by factor {lr_factor:8.6f} from '
                  f'{self.learning_rate:8.3e} to {self.learning_rate*lr_factor:8.3e}.')
        self.learning_rate = self.learning_rate * lr_factor
        # Recompile with the new learning rate
        self.recompile()

    def restore_best_weights(self):
        """Restore the best weights"""
        self.candidate_elements.set_weights(self.best_weights)
    
    def adaptive_episode_end(self, hist, verbose: int = 1):
        """Post-processing after one episode of adaptive training"""
        # Update training counters
        self.episode_length = len(hist.epoch)
        self.current_epoch += self.episode_length
        self.training_episode += 1
        self.episode_time = (time.time() - self.t0)

        # Update best loss and weights
        current_loss: float = self.current_loss()
        if verbose > 1:
            print(f'Updating best_loss & best_weights at end of episode {self.training_episode-1}')
        if current_loss < self.best_loss:
            if verbose > 1:
                print(f'New best_loss = {current_loss:8.2f}.  Old best_loss was {self.best_loss:8.2f}')
            self.best_loss = current_loss
            self.best_weights = self.candidate_elements.get_weights()
        else:
            if verbose > 0:
                print(f'Manually restoring best weights to recover best_loss = {self.best_loss:8.2f}')
            self.restore_best_weights()
            current_loss = self.current_loss()

        # Log likelihood is negative of the loss
        log_like = -current_loss        

        # Update early_stop
        early_stop = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=0, baseline=self.best_loss, 
                                                      restore_best_weights=True)

        # Update loss history and training time
        self.loss_hist += hist.history['loss']
        self.total_training_time += self.episode_time

        # Status message
        if verbose > 0:
            # print(f'Epoch {self.current_epoch:4}. Elapsed time = {self.episode_time:0.0f} sec')
            print(f'Log Likelihood: {log_like:8.2f}')        
        
    def search_adaptive(self, 
                        max_batches: int = 10000, 
                        batches_per_epoch: int = 100, 
                        epochs_per_episode: int = 5,
                        verbose: int = 1):
        """
        Run asteroid search model adaptively.  
        Start with a high learning rate, gradually reduce it if early stopping triggered
        """

        # Start timer
        self.t0 = time.time()

        # Early stopping callback
        early_stop = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=0, restore_best_weights=True)
        callbacks = [early_stop]

        # Define one epoch as a number of batches        
        samples_per_epoch: int = batches_per_epoch * self.batch_size
        max_epochs: int = max_batches // batches_per_epoch
        x_trn = tf.ones(samples_per_epoch, dtype=dtype)

        # Set the learning rate factor
        lr_factor_up: float = 1.0
        lr_factor_dn: float = 0.5

        # Run first episode of training
        hist = self.fit(x=x_trn, batch_size=self.batch_size, epochs=epochs_per_episode, 
                        steps_per_epoch=batches_per_epoch, 
                        callbacks=callbacks, shuffle=False, verbose=verbose)            
        self.adaptive_episode_end(hist)

        # Continue training until max_epochs have elapsed
        while self.current_epoch < max_epochs:        
            # If the last training ran without early stopping, increase the learning rate
            if self.episode_length == epochs_per_episode:
                if lr_factor_up != 1.0:
                    self.adjust_learning_rate(lr_factor_up, verbose=False)
            # If the last training hit early stopping, decrease the learning rate
            else:
                self.adjust_learning_rate(lr_factor_dn, verbose=True)
            # Train for another episode
            print(f'\nTraining episode {self.training_episode}: Epoch {self.current_epoch:4}')
            print(f'learning_rate={self.learning_rate:8.3e}, total training time {self.total_training_time:0.0f} sec.')
            hist = self.fit(x=x_trn, batch_size=self.batch_size, epochs=self.current_epoch+epochs_per_episode, 
                            steps_per_epoch=batches_per_epoch, initial_epoch=self.current_epoch,
                            callbacks=callbacks,  shuffle=False, verbose=verbose)        
            self.adaptive_episode_end(hist)        
# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
