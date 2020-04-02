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
# Functional API model
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_model_asteroid_search(elts: pd.DataFrame,
                               ztf_elt: pd.DataFrame,
                               site_name: str='geocenter',
                               thresh_deg: float = 1.0,
                               h: float = 0.01,
                               R_deg: float = 1.0):
    """
    Make functional API model for scoring elements
    INPUTS:
        elts:       DataFrame with initial guess for candidate orbital elements.
                    Columns: element_id, a, e, inc, Omega, omega, f, epoch
                    Output of orbital_element_batch, perturb_elts or random_elts
        ztf_elt:    DataFrame with ZTF observations within thresh_deg degrees of
                    of the orbits predicted by these elements.
                    Output of make_ztf_batch or load_ztf_batch
        site_name:  Used for topos adjustment, e.g. 'geocenter' or 'palomar'
        h:          Initial value of hit probability in mixture model
        lam:        Initial value of exponential decay parameter in mixture model
        R_deg:      Initial value of resolution parameter (in degrees) in mixture model
    """

    # Element batch size comes from elts
    elt_batch_size = elts.shape[0]

    # Dummy input; this does absolutely nothing, but otherwise keras pukes
    # Did I mention that sometimes I HATE keras?
    x = keras.Input(shape=(), batch_size=elt_batch_size, name='x')

    # Numpy array and tensor of observation times; flat, shape (data_size,)
    ts_np = ztf_elt.mjd.values.astype(dtype_np)

    # Get observation count per element
    row_lengths_np = ztf_elt.element_id.groupby(ztf_elt.element_id).count()

    # Shape of the observed trajectories
    data_size = ztf_elt.shape[0]
    traj_shape = (data_size, space_dims)

    # Observed directions; extract from ztf_elt DataFrame
    cols_u_obs = ['ux', 'uy', 'uz']
    u_obs_np = ztf_elt[cols_u_obs].values.astype(dtype_np)

    # Convert resolution from degrees to Cartesian distance
    R_s = deg2dist(R_deg)

    # Set of trainable weights with candidate orbital elements; initialize according to elts
    elements_layer = CandidateElements(elts=elts, h=h, R=R_s, name='candidates')
    
    # Extract the candidate elements and mixture parameters; pass dummy inputs to satisfy keras Layer API
    a, e, inc, Omega, omega, f, epoch, h, lam, R, = elements_layer(inputs=x)
    
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
    R = Identity(name='R')(R)

    # Stack the current orbital elements; shape is [elt_batch_size, 7,]
    elts_tf = tf.keras.layers.Concatenate(axis=-1, name='elts') ([a, e, inc, Omega, omega, f, epoch,])

    # Stack mixture model parameters; shape [elt_batch_size, 3,]
    mixture = tf.keras.layers.Concatenate(axis=-1, name='mixture') ([h, lam, R,])

    # The predicted direction; shape is [data_size, 3,]
    direction_layer = AsteroidDirection(ts_np=ts_np, row_lengths_np=row_lengths_np, 
                                        site_name=site_name, name='direction')

    # Calibration arrays (flat)
    cols_q_ast = ['qx', 'qy', 'qz']
    cols_v_ast = ['vx', 'vy', 'vz']
    q_ast = ztf_elt[cols_q_ast].values.astype(dtype_np)
    v_ast = ztf_elt[cols_v_ast].values.astype(dtype_np)

    # Run calibration
    direction_layer.q_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)

    # Tensor of predicted directions.  Shape is [data_size, 3,]
    u_pred, r_pred = direction_layer(a, e, inc, Omega, omega, f, epoch)

    # Score layer for these observations
    score_layer = TrajectoryScore(row_lengths_np=row_lengths_np, u_obs_np=u_obs_np,
                                  thresh_deg=thresh_deg, name='score')

    # Compute the log likelihood by element from the predicted direction and mixture model parameters
    # Shape is [elt_batch_size,]
    log_like = score_layer(u_pred, h=h, lam=lam)

    # Wrap inputs and outputs
    inputs = (x,)
    outputs = (log_like, elts_tf, mixture)

    # Create model with functional API
    model = keras.Model(inputs=inputs, outputs=outputs)

    # Bind the custom layers to model
    model.elements = elements_layer
    model.direction = direction_layer
    model.position = model.direction.q_layer
    model.score = score_layer

    # Save the batch size to the model
    model.elt_batch_size = elt_batch_size
  
    # Save the learning rate on the model object to facilitate adaptive training
    model.learning_rate = 1.0E-4
    model.clipnorm = 1.0
    # model.clipvalue = None

    # Add the loss function - the NEGATIVE of the log likelihood
    # (Take negative b/c TensorFlow minimizes the loss function)
    model.add_loss(-tf.reduce_sum(log_like))

    # Compile the model
    compile_search_model(model)

    # Initialize best_loss and best_weights
    x_eval = tf.ones(elt_batch_size, dtype=dtype)
    model.best_loss = model.evaluate(x_eval, verbose=0)
    model.best_weights = model.elements.get_weights()

    return model

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
def compile_search_model(model):
    """Recompile this model with the current learning rate"""
    optimizer = make_adam_opt(learning_rate=model.learning_rate, clipnorm=model.clipnorm)
    model.compile(optimizer=optimizer)

# ********************************************************************************************************************* 
def adjust_learning_rate(model, lr_factor, verbose=False):
    """Reduce the learning rate and recompile the model"""
    if verbose:
        print(f'Changing learning rate by factor {lr_factor:8.6f} from '
              f'{model.learning_rate:8.3e} to {model.learning_rate*lr_factor:8.3e}.')
    model.learning_rate = model.learning_rate * lr_factor
    compile_search_model(model)

# ********************************************************************************************************************* 
def ast_search_adaptive(model, learning_rate=None, clipnorm=None,
                        batch_size: int = 64, max_epochs: int = 100):
    """
    Run asteroid search model adaptively.  
    Start with a high learning rate, gradually reduce it if early stopping triggered
    """
    
    # Start timer
    t0 = time.time()
    
    # Set the learning rate and clipnorm to the model if they were specified
    if learning_rate is not None:
        model.learning_rate = learning_rate
    if clipnorm is not None:
        model.clipnorm = clipnorm
    
    # Recompile the model
    compile_search_model(model)
    
    # Early stopping callback
    early_stop = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=0, restore_best_weights=True)
    callbacks = [early_stop]

    # Define one epoch as a number of batches
    steps_per_epoch: int = 200
    # steps_per_epoch: int = 10
    samples_per_epoch: int = steps_per_epoch * batch_size
    x_trn: np.ndarray = np.ones(samples_per_epoch)

    # Set number of epochs for one episode of training
    epochs_per_episode: int = 5 # min(10, max_epochs)
    training_episode: int = 1

    # Epoch and episode counters
    episode_length: int = 0
    model.current_epoch: int = 0
    model.best_loss: float = np.inf
    model.best_weights: np.ndarray

    # Set the learning rate factor
    lr_factor_up: float = 1.0 # 2.0**0.125
    lr_factor_dn: float = 0.5
        
    # Verbose flag for training
    verbose = 1
    
    def after_episode():
        nonlocal episode_length, training_episode
        # Update training counters
        episode_length = hist.epoch[-1] + 1
        model.current_epoch += episode_length
        training_episode += 1
        elapsed_time = (time.time() - t0)

        # Update best loss and weights
        current_loss: float = hist.history['loss'][-1]
        print(f'Updating best_loss at end of episode {training_episode}')
        if current_loss < model.best_loss:
            print(f'New best_loss = {current_loss:8.2f}.  Old best_loss was {model.best_loss:8.2f}')
            model.best_loss = current_loss
            model.best_weights = model.elements.get_weights()
        else:
            print(f'Manually restoring best weights to recover best_loss = {model.best_loss:8.2f}')
            model.elements.set_weights(model.best_weights)
            current_loss = model.best_loss

        # Log likelihood is negative of the loss
        log_like = -current_loss        
        
        # Update early_stop
        early_stop = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=0, baseline=model.best_loss, restore_best_weights=True)

        # Status message
        print(f'Epoch {model.current_epoch:4}. Elapsed time = {elapsed_time:0.0f} sec')
        print(f'Log Likelihood: {log_like:8.2f}')

    # Run first episode of training
    hist = model.fit(x=x_trn, batch_size=batch_size, epochs=epochs_per_episode, steps_per_epoch=steps_per_epoch, 
                     callbacks=callbacks, shuffle=False, verbose=verbose)            
    after_episode()
    
    # Continue training until max_epochs have elapsed
    while model.current_epoch < max_epochs:        
        # If the last training ran without early stopping, increase the learning rate
        if episode_length == epochs_per_episode:
            # adjust_learning_rate(model, lr_factor_up, verbose=False)
            pass
        # If the last training hit early stopping, decrease the learning rate
        else:
            adjust_learning_rate(model, lr_factor_dn, verbose=False)
        # Train for another episode
        print(f'\nTraining episode {training_episode}:')
        print(f'Epoch {model.current_epoch:4}, learning_rate={model.learning_rate:8.3e}, clipnorm={model.clipnorm:6.3f}.')
        hist = model.fit(x=x_trn, batch_size=batch_size, epochs=epochs_per_episode, steps_per_epoch=steps_per_epoch, 
                         callbacks=callbacks, shuffle=False, verbose=verbose)        
        after_episode()
        
# ********************************************************************************************************************* 
if __name__ == '__main__':
    pass
