"""
Harvard IACS Masters Thesis
asteroid_search_model.py: Tensorflow layers and models used to search for asteroid orbital elements.

Michael S. Emanuel
Thu Oct 17 15:24:10 2019
"""

# Core
import numpy as np

# Tensorflow / ML
import tensorflow as tf

# Utility
from datetime import timedelta

# MSE imports
from tf_utils import tf_quiet

# Typing
from typing import List, Tuple, Dict, Optional, Union

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Run TF quietly
tf_quiet()

# ********************************************************************************************************************* 
def make_opt_adam(learning_rate: float, 
                  clipnorm: float=1.0, 
                  clipvalue:Optional[float]=None) \
                  -> keras.optimizers.Optimizer:
    """
    Build Adam optimizer for training
    INPUTS:
        learning_rate: The learning rate for the optimizer
        clipnorm:      Parameter for gradient clipping by norm.
        clipvalue:     Parameter for gradient clipping by value; probably better to use clipnorm.
    Default settings in tensorflow are:
    learning_rate = 1.0E-3
    clipnorm = None
    clipvalue = None
    These are changed based on trial and error.
    Other arguments left at default settings
    """

    # Settings for other arguments; leave at defaults
    beta_1: float = 0.900          # default 0.900
    beta_2: float = 0.999          # default 0.999
    epsilon: float = 1.0E-7        # default 1.0E-7
    amsgrad: float = False         # default False

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args: Dict[str, float] = {
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
    opt: keras.optimizers.Optimizer = keras.optimizers.Adam(**opt_args)
    return opt
    
# ********************************************************************************************************************* 
def make_opt_sgd(learning_rate: float, 
                 clipnorm: float=1.0, 
                 clipvalue: Optional[float]=None) \
                 -> keras.optimizers.Optimizer:
    """
    Build SGD optimizer for training 
    INPUTS:
        learning_rate: The learning rate for the optimizer. Recommended default is 2^-15.
        clipnorm:      Parameter for gradient clipping by norm.
        clipvalue:     Parameter for gradient clipping by value; probably better to use clipnorm.
    Default settings in tensorflow are:
    learning_rate = 1.0E-2
    momentum = 0.0
    nestorov = False
    clipnorm = None
    clipvalue = None
    """

    # Settings for other arguments; leave at defaults
    momentum: float = 0.000        # default 0.000
    nesterov: bool = False

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args: Dict = {
        'learning_rate': learning_rate,
        'momentum': momentum,
        'nesterov': nesterov,
    }
    # Add clipnorm if it was set
    if clipnorm is not None:
        opt_args['clipnorm'] = clipnorm
    # Add clipvalue if it was set
    if clipvalue is not None:
        opt_args['clipvalue'] = clipvalue

    # Build the optimizer
    opt: keras.optimizers.Optimizer = keras.optimizers.SGD(**opt_args)
    return opt

# ********************************************************************************************************************* 
def make_opt_rmsprop(learning_rate: float, 
                     clipnorm: float=1.0, 
                     clipvalue: Optional[float]=None) \
                     -> keras.optimizers.Optimizer:
    """
    Build RMSprop optimizer for training 
    INPUTS:
        learning_rate: The learning rate for the optimizer. Recommended default is 2^-15.
        clipnorm:      Parameter for gradient clipping by norm.
        clipvalue:     Parameter for gradient clipping by value; probably better to use clipnorm.
    Default settings in tensorflow are:
    learning_rate = 1.0E-3
    rho = 0.90
    momentum = 0.0
    epsilon = 1.0E-7
    centered = False
    clipnorm = None
    clipvalue = None
    """

    # Settings for other arguments; leave at defaults
    rho: float = 0.900             # default 0.900
    momentum: float = 0.000        # default 0.000
    epsilon: float = 2.0**-23      # default 1.0E-7; nearest power of 2
    centered: bool = False

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args: Dict = {
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
    opt: keras.optimizers.Optimizer = keras.optimizers.RMSprop(**opt_args)
    return opt

# ********************************************************************************************************************* 
def make_opt_adadelta(learning_rate: float, 
                      clipnorm: float=1.0, 
                      clipvalue: Optional[float]=None) \
                      -> keras.optimizers.Optimizer:
    """
    Build Adadelta optimizer for training 
    INPUTS:
        learning_rate: The learning rate for the optimizer. Recommended default is 2^-15.
        clipnorm:      Parameter for gradient clipping by norm.
        clipvalue:     Parameter for gradient clipping by value; probably better to use clipnorm.
    Default settings in tensorflow are:
    learning_rate = 1.0E-3
    rho = 0.950
    epsilon = 1.0E-7
    clipnorm = None
    clipvalue = None
    """

    # Settings for other arguments; leave at defaults
    rho: float = 0.950             # default 0.950
    epsilon: float = 2.0**-23      # default 1.0E-7; nearest power of 2

    # Optimizer arguments - wrap as dict.  Point is to permit clip=None
    opt_args: Dict = {
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
        learning_rate: The learning rate for the optimizer. Recommended default is 2^-15.
        optimizer_type: one of 'adam', 'sgd', 'rmsprop', 'adadelta'
        clipnorm:       gradient clipping by norm of gradient vector
        clipvalue:      gradient clipping by element of gradient vector
    """
    # Table of factory functions keyed by optimizer_type
    optimizer_func: Dict = {
        'adam': make_opt_adam,
        'sgd': make_opt_sgd,
        'rmsprop': make_opt_rmsprop,
        'adadelta': make_opt_adadelta,
    }
    # The factor function for the selected optimizer_type
    optimizer_func = optimizer_func[optimizer_type]
    # Instantiate this optimizer with selected input parameters
    return optimizer_func(learning_rate=learning_rate, clipnorm=clipnorm, clipvalue=clipvalue)
