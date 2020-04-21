"""
Harvard IACS Masters Thesis
asteroid_search.py: Search for orbital elements of asteroids given observational data.

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
import argparse

# Local imports
# from asteroid_search_model import make_model_asteroid_search
# from ztf_data import load_ztf_easy_batch, make_ztf_batch, report_ztf_score
from asteroid_search_model import AsteroidSearchModel
from ztf_element import make_ztf_batch, report_ztf_score
from asteroid_data import make_ztf_dataset
from asteroid_element import load_ast_elt
from asteroid_integrate import calc_ast_pos
from asteroid_search_report import report_model, report_training_progress
from candidate_element import perturb_elts
from utils import print_header
from tf_utils import tf_quiet, gpu_grow_memory, get_gpu_device

# Typing
from typing import Dict

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Tensorflow config
# Run TF quietly
tf_quiet()

# Configure TensorFlow to use GPU memory variably; this is done in asteroid_model
# gpu_grow_memory(verbose=True)

# ********************************************************************************************************************* 
# Constants
ast_elt = load_ast_elt()

# ********************************************************************************************************************* 
if __name__ == '__main__':    
    # Process command line arguments
    parser = argparse.ArgumentParser(description='Search for asteroids.')
    parser.add_argument('-thresh_deg', nargs='?', metavar='thresh_deg', type=float, default=1.0,
                        help='threshold in degrees; observations further from candidate elements than this are ignored')
    parser.add_argument('-R_deg', nargs='?', metavar='R_deg', type=float, default=0.2,
                        help='initial value of resolution factor, R, in degrees')
    parser.add_argument('-elt_batch_size', nargs='?', metavar='elt_batch_size', type=int, default=64,
                        help='the number of candidate orbital elements per batch')
    parser.add_argument('-epochs', nargs='?', metavar='elt_batch_size', type=int, default=10,
                        help='the number of epochs to train')
    parser.add_argument('-gpu_num', nargs='?', metavar='gpu_num', type=int, default=1,
                        help='the GPU to use; defaults to 1 to avoid clash with Jupyter sessions')
    parser.add_argument('--test', default=True, action='store_true',
                        help='test on easy batch')
    args = parser.parse_args()

    # Alias arguments
    R_deg = args.R_deg
    thresh_deg = args.thresh_deg
    R_is_trainable = args.R_is_trainable
    elt_batch_size = args.elt_batch_size
    epochs = args.epochs
    cycles_per_epoch = args.cycles_per_epoch
    gpu_num = args.gpu_num

    # The selected gpu device
    gpu_device = get_gpu_device(gpu_num)

    # Set time_batch_size to None always to use all available snapshots (with threshold this is way better)
    time_batch_size = None

    # Set use_calibration; should be True unless doing quick code testing
    use_calibration = True

    # If run in test mode, test on the easy batch
    if args.test:
        with gpu_device:
            print(f'Running test_easy_batch with gpu {gpu_num}.')