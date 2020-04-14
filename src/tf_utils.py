"""
Harvard IACS Masters Thesis
Tensorflow Utilites

Michael S. Emanuel
Wed Jun 26 20:24:42 2019
"""

# Tensorflow
import tensorflow as tf
from silence_tensorflow import silence_tensorflow

# Core
import numpy as np

# Utility
import matplotlib.pyplot as plt
import time
import datetime

# *************************************************************************************************
# Aliases
keras = tf.keras

# Data type (for degrees / radians etc.)
dtype = tf.float32

# ********************************************************************************************************************* 
def tf_quiet():
    """Silence excessive TensorFlow warnings and status messages"""
    silence_tensorflow()

# *************************************************************************************************
def gpu_grow_memory(verbose : bool=False):
    """Set TensorFlow to grow memory of GPUs rather than grabbing it all at once."""
    # Get available GPUs
    gpus = tf.config.list_physical_devices('GPU')
    # List GPUs if requested
    if verbose:
        print(f'Found {len(gpus)} GPUs.  Setting memory growth = True.')
    # Set all selected GPUs to 
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)

# ********************************************************************************************************************* 
def get_gpu_device(gpu_num: int):
    """Select a specific GPU"""
    # Name of this device
    device_name = f':/gpu:{gpu_num}'
    # Select the desired GPU and return it; use it in a with tf.device context
    gpu_device = tf.device(device_name)
    return gpu_device

# ********************************************************************************************************************* 
def plot_loss_hist(hist,  model_name, key='loss', baseline=None):
    """Plot loss vs. wall time"""
    # Extract loss and wall time arrays
    loss = hist[key]
    time = hist['time']
    
    # Plot loss vs. wall time
    fig, ax = plt.subplots(figsize=[16,9])
    ax.set_title(f'Loss vs. Wall Time for {model_name}')
    ax.set_xlabel('Wall Time (Seconds)')
    ax.set_ylabel(f'{key}')
    ax.plot(time, loss, color='blue', label=key)
    if baseline is not None:
        ax.axhline(baseline, color='red', label='baseline')
    ax.set_yscale('log')
    ax.grid()
    ax.legend()

    return fig, ax

# ********************************************************************************************************************* 
# Angles between radians and degrees
coef_rad2deg = tf.constant(value=(np.pi / 180.0), dtype=dtype)
coef_deg2rad = tf.constant(value=(180.0 / np.pi), dtype=dtype)

# ********************************************************************************************************************* 
@tf.function
def tf_deg2rad(deg):  
    return tf.multiply(deg, coef_deg2rad)

# ********************************************************************************************************************* 
@tf.function
def tf_rad2deg(rad):
    return tf.multiply(rad, coef_rad2deg)
    
# ********************************************************************************************************************* 
@tf.function
def tf_dist2rad(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to radians in [0, pi]"""
    half_dist = tf.multiply(dist, 0.5)
    half_x_rad = tf.asin(half_dist)
    x_rad = tf.multiply(half_x_rad, 2.0)
    return x_rad

# ********************************************************************************************************************* 
@tf.function
def tf_rad2dist(x_rad):
    """Convert a distance on unit sphere from radians in [0, pi] to cartesian distance in [0, 2]"""
    half_x_rad = tf.multiply(x_rad, 0.5)
    sin_half_x = tf.sin(half_x_rad)
    dist = tf.multiply(sin_half_x, 2.0)
    return dist

# ********************************************************************************************************************* 
@tf.function
def tf_dist2deg(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to degrees in [0, 180]"""
    x_rad = tf_dist2rad(dist)
    return tf_rad2deg(x_rad)

# ********************************************************************************************************************* 
@tf.function
def tf_deg2dist(x_deg):
    """Convert a distance on unit sphere from degrees in [0, 180] to cartesian distance in [0, 2]"""
    x_rad = tf_deg2rad(x_deg)
    return tf_rad2dist(x_rad)

# ********************************************************************************************************************* 
@tf.function
def tf_dist2sec(dist):
    """Convert a cartesian distance on unit sphere in [0, 2] to arc seconds in [0, 180*3600]"""
    x_rad = tf_dist2rad(dist)
    x_deg = tf_rad2deg(x_rad)
    x_sec = tf.multiply(x_deg, 3600.0)
    return x_sec

# *************************************************************************************************
# Custom Callbacks
# *************************************************************************************************

# *************************************************************************************************
# https://stackoverflow.com/questions/43178668/record-the-computation-time-for-each-epoch-in-keras-during-model-fit 
class TimeHistory(keras.callbacks.Callback):
    """Save the wall time after every epoch"""
    def on_train_begin(self, logs={}):
        self.times = []
        self.train_time_start = time.time()

    def on_epoch_end(self, batch, logs={}):
        self.times.append(time.time() - self.train_time_start)

# *************************************************************************************************
class EpochLoss(tf.keras.callbacks.Callback):
    """Log the loss every N epochs"""
    def __init__(self, interval=10, newline=False):
        super(EpochLoss, self).__init__()
        self.interval = interval
        self.train_time_start = time.time()
        self.log_prefix = '\n' if newline else ''

    def log_to_screen(self, epoch, logs):
        loss = logs['loss']
        elapsed = time.time() - self.train_time_start
        elapsed_str = str(datetime.timedelta(seconds=np.round(elapsed)))
        print(f'{self.log_prefix}Epoch {epoch:04}; loss {loss:5.2e}; elapsed {elapsed_str}') 

    def on_epoch_end(self, epoch, logs=None):
        epoch = epoch+1
        if (epoch % self.interval == 0) or (epoch == 1):
            self.log_to_screen(epoch, logs)

# ********************************************************************************************************************* 
# Custom Layers
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
class Identity(keras.layers.Activation):
    """Identity layer for labeling outputs in models"""
    def __init__(self, **kwargs):
        super(Identity, self).__init__(activation = tf.identity, **kwargs)
        
    def get_config(self):
        return dict()

# ********************************************************************************************************************* 
class Divide(keras.layers.Layer):
    """Division"""
    def call(self, inputs):
        # Unpack inputs
        x, y = inputs
        
        # Return quotient x / y
        return tf.divide(x, y)
    
    def get_config(self):
        return dict()
