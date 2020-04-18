"""
Harvard IACS Masters Thesis
Tensorflow Utilites

Michael S. Emanuel
Fri Apr 17 2020 20:37
"""

# Core
import numpy as np

# Tensorflow
import tensorflow as tf

# ********************************************************************************************************************* 
# Data type
dtype = tf.float32

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
