"""
Harvard IACS Masters Thesis
Solar Asteroid Model: Predict the movement of a test particle (e.g. asteroid) in the solar system
using the Kepler approximation with the sun as a fixed central attractor.

Michael S. Emanuel
Sun Oct 13 11:56:50 2019
"""

# Core
import tensorflow as tf
import pandas as pd
import numpy as np

# Astronomy
import astropy
from astropy.units import au, day, year

# Local imports
from orbital_element import MeanToTrueAnomaly, TrueToMeanAnomaly
from asteroid_data import get_earth_pos, get_sun_pos_vel
from candidate_element import orbital_element_batch
from ra_dec import calc_topos
from tf_utils import tf_quiet, gpu_grow_memory, Identity

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Run TF quietly
tf_quiet()

# Configure TensorFlow to use GPU memory variably
gpu_grow_memory(verbose=True)

# ********************************************************************************************************************* 
# Constants

# The gravitational constant in ('day', 'AU', 'Msun') coordinates
# sim = rebound.Simulation()
# sim.units = ('day', 'AU', 'Msun')
# G_ = sim.G
# Hard code G
G_ = 0.00029591220828559104
# The gravitational field strength mu = G * (m0 + m1)
# For massless asteroids orbiting the sun with units Msun, m0=1.0, m1=0.0, and mu = G
mu = tf.constant(G_)

# Speed of light; express this in AU / day
light_speed_au_day = astropy.constants.c.to(au / day).value

# Number of spatial dimensions
space_dims = 3

# Data type
dtype = tf.float32

# ********************************************************************************************************************* 
# Custom Layers
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
class ElementToPosition(keras.layers.Layer):
    def __init__(self, **kwargs):
        super(ElementToPosition, self).__init__(**kwargs)

    @tf.function
    def call(self, inputs):
        """Compute position from orbital elements (a, e, inc, Omega, omega, f)"""
        # Unpack inputs
        # a: semimajor axis
        # e: eccentricity
        # inc: inclination
        # Omega: longitude of ascending node
        # omega: argument of pericenter
        # f: true anomaly
        a, e, inc, Omega, omega, f = inputs

        # Shape of input
        shape = a.shape
        
        # sine and cosine of the angles inc, Omega, omega, and f
        ci = keras.layers.Activation(activation=tf.cos, name='cos_inc')(inc)
        si = keras.layers.Activation(activation=tf.sin, name='sin_inc')(inc)
        cO = keras.layers.Activation(activation=tf.cos, name='cos_Omega')(Omega)
        sO = keras.layers.Activation(activation=tf.sin, name='sin_Omega')(Omega)
        co = keras.layers.Activation(activation=tf.cos, name='cos_omega')(omega)
        so = keras.layers.Activation(activation=tf.sin, name='sin_omega')(omega)
        cf = keras.layers.Activation(activation=tf.cos, name='cos_f')(f)
        sf = keras.layers.Activation(activation=tf.sin, name='sin_f')(f)

        # Distance from center
        e2 = keras.layers.Activation(activation=tf.square, name='e2')(e)
        one = tf.broadcast_to(1.0, shape)
        one_minus_e2 = tf.subtract(one, e2, name='one_minus_e2')
        e_cos_f = tf.multiply(e, cf, name='e_cos_f')
        one_plus_e_cos_f = tf.add(one, e_cos_f, name='one_plus_e_cos_f')
        a_x_one_minus_e2 = tf.multiply(a, one_minus_e2, name='a_x_one_minus_e2')
        r = tf.divide(a_x_one_minus_e2, one_plus_e_cos_f, name='r')
        
        # Position
        cocf = tf.multiply(co ,cf, name='cocf')
        sosf = tf.multiply(so, sf, name='sosf')
        cocf_sosf = tf.subtract(cocf, sosf, name='cocf_sosf')

        socf = tf.multiply(so, cf, name='socf')
        cosf = tf.multiply(co, sf, name='cosf')
        socf_cosf = tf.add(socf, cosf, name='socf_cosf')

        cO_x_cocf_sosf = tf.multiply(cO, cocf_sosf, name='cO_x_cocf_sosf')
        sO_x_socf_cosf = tf.multiply(sO, socf_cosf, name = 'sO_x_socf_cosf')
        sO_x_socf_cosf_x_ci = tf.multiply(sO_x_socf_cosf, ci, name='sO_x_socf_cosf_x_ci')       
        sO_x_cocf_sosf = tf.multiply(sO, cocf_sosf, name='sO_x_cocf_sosf')
        cO_x_socf_cosf = tf.multiply(cO, socf_cosf, name='cO_x_socf_cosf')
        cO_x_socf_cosf_x_ci = tf.multiply(cO_x_socf_cosf, ci, name='cO_x_socf_cosf_x_ci')

        # Direction components
        ux = tf.subtract(cO_x_cocf_sosf, sO_x_socf_cosf_x_ci, name='ux')
        uy = tf.add(sO_x_cocf_sosf, cO_x_socf_cosf_x_ci, name='uy')
        uz = tf.multiply(socf_cosf, si, name='socf_cosf_x_si')

        # Position components; units are in AU
        qx = tf.multiply(r, ux, name='qx')
        qy = tf.multiply(r, uy, name='qy')
        qz = tf.multiply(r, uz, name='qz')

        # Assemble the position vector
        q = keras.layers.concatenate(inputs=[qx, qy, qz], axis=-1, name='q')

        # Calculate the velocity
        # Current speed
        v0 = tf.sqrt(mu / a / one_minus_e2)
        # The term e+cf appears three times
        epcf = tf.add(e, cf)
        # The term cocO appears twice
        cocO = tf.multiply(co, cO)
        # The term cosO appears twice
        cosO = tf.multiply(co, sO)
        # The term so*sO appears twice
        sosO = tf.multiply(so, sO)
        # The terms socO appears twice
        socO = tf.multiply(so, cO)
        # Simplified expression for velocity with substitutions
        vx = v0*(epcf*(-ci*cosO - socO) - sf*(cocO - ci*sosO))
        vy = v0*(epcf*(ci*cocO - sosO)  - sf*(cosO + ci*socO))
        vz = v0*(epcf*co*si - sf*si*so)

        # Assemble the velocity vector; units are in AU/year
        v = keras.layers.concatenate(inputs=[vx, vy, vz], axis=-1, name='v')

        return q, v

    def get_config(self):
        return dict()

# ********************************************************************************************************************* 
class AsteroidPosition(keras.layers.Layer):
    """
    Compute orbit positions for asteroids in the solar system from the initial orbital elements with the Kepler model.
    Inputs for the model are 6 orbital elements, the epoch, and the desired times for position outputs.
    Outputs of the model are the position of the asteroid in the BARYCENTRIC frame.
    """
    def __init__(self, ts_np: np.ndarray, row_lengths_np: np.ndarray, **kwargs):
        """
        INPUTS:
            ts_np: Numpy array of time snapshots at which to simulate the position.
                   Must be passed as a Numpy array for compatibility with get_sun_pos_vel()
                   Shape [data_size, 1]; flattened to match ztf_elt Dataframe
            row_lengths_np: Number of observations for each element; shape [elt_batch_size]
                            elt_batch_size is the number of elements simulated, e.g. 64.   
        """
        super(AsteroidPosition, self).__init__(**kwargs)

        # Configuration for serialization
        self.cfg = {
            'ts_np': ts_np,
            'row_lengths_np': row_lengths_np
        }

        # Save ts and row_lenghts as a keras constants
        self.ts = keras.backend.constant(value=ts_np, shape=ts_np.shape, dtype=dtype)
        self.row_lengths = keras.backend.constant(value=row_lengths_np, shape=row_lengths_np.shape, dtype=tf.int32)

        # Infer elt_batch_size from shape of row_lengths
        self.elt_batch_size = keras.backend.constant(value=self.row_lengths.shape[0], dtype=tf.int32)

        # Save the data size
        self.data_size = keras.backend.constant(value=tf.reduce_sum(row_lengths_np), dtype=tf.int32)
        
        # Take a one time snapshot of the sun's position and velocity at these times
        q_sun_np, v_sun_np = get_sun_pos_vel(ts_np)

        # Convert q_sun and v_sun into keras constants
        traj_shape = (self.data_size, space_dims)
        self.q_sun = keras.backend.constant(q_sun_np, dtype=dtype, shape=traj_shape, name='q_sun')
        self.v_sun = keras.backend.constant(v_sun_np, dtype=dtype, shape=traj_shape, name='v_sun')

        # The adjustments dq and dv to correct the Kepler approximation to match the numerical integration
        self.dq = tf.Variable(initial_value=np.zeros(traj_shape), dtype=dtype, trainable=False, name='dq')
        self.dv = tf.Variable(initial_value=np.zeros(traj_shape), dtype=dtype, trainable=False, name='dv')

    def update_dq_dv(self, dq, dv):
        """Update the value of dq and dv"""
        self.dq.assign(dq)
        self.dv.assign(dv)

    def calibrate(self, elts: pd.DataFrame, q_ast: np.ndarray, v_ast: np.ndarray):
        """Calibrate this model by setting dq to recover q_ast"""
        # Unpack elements
        a = elts['a']
        e = elts['e']
        inc = elts['inc']
        Omega = elts['Omega']
        omega = elts['omega']
        f = elts['f']
        epoch = elts['epoch']

        # Zero out calibration and predict with these elements
        self.update_dq_dv(self.dq*0.0, self.dv*0.0)
        q_pred, v_pred = self.call(a, e, inc, Omega, omega, f, epoch)

        # We expected to get the numerically integrated barycentric coordinates of the asteroids
        # Compute the offset and apply it to the model
        dq = q_ast - q_pred
        dv = v_ast - v_pred
        self.update_dq_dv(dq, dv)

    @tf.function
    def call(self, a, e, inc, Omega, omega, f, epoch):
        """
        Simulate the orbital trajectories.
        Snapshot times t shared by all the input elements.  
        The inputs orbital elements and reference epoch should all have size [data_size,].
        That is, inputs are flat, not ragged, e.g. 92000 entries for 64 candidate elements.
        Outputs are the barycentric position and velocity: q, v.
        The elements generate a heliocentric calculation, and the constants q_sun, v_sun are added.
        These are also flat, with shape [data_size, 3,]
        """
        # Alias row_lengths for legibility
        row_lengths = self.row_lengths

        # Times, elements all have shape (data_size, 1,)
        elt_shape = (self.data_size, 1,)

        # Time relative to epoch; call this t
        epoch_t = tf.reshape(tf.repeat(epoch, row_lengths), elt_shape)
        t = keras.layers.subtract([self.ts, epoch_t], name='t')

        # Compute eccentric anomaly E from f and e
        M = TrueToMeanAnomaly(name='TrueToMeanAnomaly')([f, e])
        
        # Compute mean motion N from mu and a
        a3 = tf.math.pow(a, 3, name='a3')
        mu_over_a3 = tf.divide(mu, a3, name='mu_over_a3')
        N = tf.sqrt(mu_over_a3, name='N')

        # ******************************************************************
        # Predict orbital elements over time

        # Repeat the constant orbital elements to be vectors of shape elt_shape = (data_size, 1)
        a_t = tf.reshape(tensor=tf.repeat(a, row_lengths), shape=elt_shape, name='a_t')
        e_t = tf.reshape(tensor=tf.repeat(e, row_lengths), shape=elt_shape, name='e_t')
        inc_t = tf.reshape(tensor=tf.repeat(inc, row_lengths), shape=elt_shape, name='inc_t')
        Omega_t = tf.reshape(tensor=tf.repeat(Omega, row_lengths), shape=elt_shape, name='Omega_t')
        omega_t = tf.reshape(tensor=tf.repeat(omega, row_lengths), shape=elt_shape, name='omega_t')
        
        # Repeat initial mean anomaly M0 and mean motion N0 to match shape of outputs
        M0_t = tf.reshape(tensor=tf.repeat(M, row_lengths), shape=elt_shape, name='M0_t')
        N0_t = tf.reshape(tensor=tf.repeat(N, row_lengths), shape=elt_shape, name='N0_t')
        # Compute the mean anomaly M(t) as a function of time
        N_mult_t = keras.layers.multiply(inputs=[N0_t, t], name='N_mult_t')
        M_t = keras.layers.add(inputs=[M0_t, N_mult_t], name='M_t')
    
        # Compute the true anomaly from the mean anomly and eccentricity
        f_t = MeanToTrueAnomaly(name='mean_to_true_anomaly')([M_t, e_t])
    
        # Wrap orbital elements into one tuple of inputs for layer converting to cartesian coordinates
        elt_t = (a_t, e_t, inc_t, Omega_t, omega_t, f_t,)
        
        # Convert orbital elements to heliocentric cartesian coordinates
        q_helio, v_helio = ElementToPosition(name='qv_helio')(elt_t)

        # Add solar position and velocity to get q, v in barycentric coordinates
        # Also add the optional correction factors dq, dv from calibration
        q = keras.layers.add(inputs=[q_helio, self.q_sun, self.dq], name='q_flat')
        v = keras.layers.add(inputs=[v_helio, self.v_sun, self.dv], name='v_flat')

        return q, v

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class AsteroidDirection(keras.layers.Layer):
    """
    Layer to compute the direction from earth to asteroid.
    """
    def __init__(self, ts_np: np.ndarray, row_lengths_np: np.ndarray, site_name: str, **kwargs):
        """
        INPUTS:
            ts: Ragged tensor of time snapshots at which to simulate the position.
                Shape [elt_batch_size, (traj_size),]; each element has a different number of time snaps
            row_lengths: Number of observations for each element; shape [elt_batch_size]
            site_name: name of the observatory site, used for topos adjustment, e.g. 'geocenter' or 'palomar'            
        """
        super(AsteroidDirection, self).__init__(**kwargs)

        # Configuration for serialization
        self.cfg = {
            'ts_np': ts_np,
            'row_lengths_np': row_lengths_np,
            'site_name': site_name,
        }

        # Save ts and row_lenghts as a keras constants
        self.ts = keras.backend.constant(value=ts_np, shape=ts_np.shape, dtype=dtype)
        self.row_lengths = keras.backend.constant(value=row_lengths_np, shape=row_lengths_np.shape, dtype=tf.int32)

        # Infer elt_batch_size from shape of row_lengths
        self.elt_batch_size = keras.backend.constant(value=self.row_lengths.shape[0], dtype=tf.int32)

        # Save the data size
        self.data_size = keras.backend.constant(value=tf.reduce_sum(self.row_lengths), dtype=tf.int32)

        # Shape of trajectories is flat: (data_size, 3,)
        traj_shape = (self.data_size, space_dims)

        # Build layer to compute positions
        self.q_layer = AsteroidPosition(ts_np=ts_np, row_lengths_np=row_lengths_np, name='q_ast')
        
        # Take a one time snapshot of the earth's position at these times in barycentric coordinates
        q_earth_np = get_earth_pos(ts_np)

        # Take a one time snapshot of the topos adjustment; displacement from geocenter to selected observatory
        dq_topos_ap = calc_topos(obstime_mjd=ts_np, site_name=site_name)
        # Convert dq_topos to a numpy array with units of au
        dq_topos_np = dq_topos_ap.to(au).value

        # Position of the observatory in barycentric frame as a Keras constant
        q_obs_np = q_earth_np + dq_topos_np
        self.q_obs = keras.backend.constant(value=q_obs_np, shape=traj_shape, dtype=dtype, name='q_obs')

    def calibrate(self, elts: pd.DataFrame, q_ast: np.ndarray, v_ast: np.ndarray):
        """Calibrate this model by calibrating the underlying position layer"""
        # Calibrate the position model
        self.q_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)

    @tf.function
    def call(self, a, e, inc, Omega, omega, f, epoch):
        """
        Simulate direction from observatory to asteroid with these orbital elements.
        Snapshot times t shared by all the input elements.  
        The inputs orbital elements and reference epoch should all have size [data_size,].
        That is, inputs are flat, not ragged, e.g. 92000 entries for 64 candidate elements.
        Outputs are the direction and distance: u, r.
        These are also flat, with u.shape = [data_size, 3,] r.shape = [data_size, 1]
        """
        # Calculate position and velocity of the asteroid
        q_ast, v_ast = self.q_layer(a, e, inc, Omega, omega, f, epoch)

        # Relative displacement from observatory to asteroid; instantaneous, before light time adjustment
        # q_rel_inst = tf.subtract(q_ast, self.q_obs, name='q_rel_inst')
        q_rel_inst = keras.layers.subtract(inputs=[q_ast, self.q_obs], name='q_rel_inst')

        # Distance between earth and asteroid, before light time adjustment
        r_inst = tf.norm(q_rel_inst, axis=-1, keepdims=True, name='r_earth_inst')

        # Light time in days from asteroid to earth in days (time units is days)
        light_time = tf.divide(r_inst, light_speed_au_day)
        
        # Adjusted relative position, accounting for light time; simulation velocity units are AU /day
        dq_lt = tf.multiply(v_ast, light_time)
        q_rel = tf.subtract(q_rel_inst, dq_lt)
        # Adjusted distance to earth, accounting for light time
        r = tf.norm(q_rel, axis=-1, keepdims=True, name='r')

        # Convert q_rel and r to ragged tensors
        # q_rel = tf.RaggedTensor.from_row_lengths(values=q_rel, row_lengths=self.row_lengths, name='q_rel_r')
        # r_r = tf.RaggedTensor.from_row_lengths(values=r, row_lengths=self.row_lengths, name='r_r')

        # Direction from earth to asteroid as unit vectors u = (ux, uy, uz)    
        u = tf.divide(q_rel, r, name='u')

        return u, r
    
    def get_config(self):
        return self.cfg
    
# ********************************************************************************************************************* 
# Functional API Models
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_model_ast_pos(ts_np: np.ndarray, row_lengths_np: np.ndarray) -> keras.Model:
    """
    Compute orbit positions for asteroids in the solar system from
    the initial orbital elements with the Kepler model.
    Factory function that returns a functional model.
    Inputs for the model are 6 orbital elements, the epoch, and the desired times for position outputs.
    Outputs of the model are the position of the asteroid relative to the sun.
    INPUTS:
        ts_np:  Numpy array of time snapshots at which to simulate the position.
                Shape [data_size, 1]; flattened to match ztf_elt Dataframe
        row_lengths_np: Number of observations for each element; shape [elt_batch_size]
    OUTPUTS:
        model: A tf.keras model to predict asteroid position and velocity (q, v,)
    """
    # Infer elt_batch_size from row_lengths
    elt_batch_size: int = row_lengths_np.shape[0]

    # Inputs: 6 orbital elements; epoch;
    a = keras.Input(shape=(), batch_size=elt_batch_size, name='a')
    e = keras.Input(shape=(), batch_size=elt_batch_size, name='e')
    inc = keras.Input(shape=(), batch_size=elt_batch_size, name='inc')
    Omega = keras.Input(shape=(), batch_size=elt_batch_size, name='Omega')
    omega = keras.Input(shape=(), batch_size=elt_batch_size, name='omega')
    f = keras.Input(shape=(), batch_size=elt_batch_size, name='f')
    epoch = keras.Input(shape=(), batch_size=elt_batch_size, name='epoch')

    # Wrap these up into one tuple of inputs for the model
    inputs = (a, e, inc, Omega, omega, f, epoch)
    
    # Build asteroid position layer
    ast_pos_layer = AsteroidPosition(ts_np=ts_np, row_lengths_np=row_lengths_np, name='ast_pos_layer')

    # Define output tensors q, v by applying the position layer to the input tensors with the elements
    q_flat, v_flat = ast_pos_layer(a, e, inc, Omega, omega, f, epoch)

    # Convert q, v to ragged tensors matching the element batch
    # q = tf.RaggedTensor.from_row_lengths(values=q_flat, row_lengths=row_lengths, name='q')
    # v = tf.RaggedTensor.from_row_lengths(values=v_flat, row_lengths=row_lengths, name='v')
    ragged_map_func = lambda x : tf.RaggedTensor.from_row_lengths(values=x, row_lengths=row_lengths_np)
    q = tf.keras.layers.Lambda(function=ragged_map_func, name='q')(q_flat)
    v = tf.keras.layers.Lambda(function=ragged_map_func, name='v')(v_flat)
    
    # Wrap up the outputs
    outputs = (q, v,)

    # Wrap this into a model
    model = keras.Model(inputs=inputs, outputs=outputs, name='model_asteroid_pos')
    
    # Bind the asteroid position layer
    model.ast_pos_layer = ast_pos_layer
    
    return model


# ********************************************************************************************************************* 
def make_model_ast_dir(ts_np: np.ndarray, row_lengths_np: np.ndarray, site_name: str) -> keras.Model:
    """
    Compute direction from earth to asteroids in the solar system from
    the initial orbital elements with the Kepler model.
    Factory function that returns a functional model.
    Inputs for the model are 6 orbital elements, the epoch, and the desired times for position outputs.
    Outputs of the model are the unit vector (direction) pointing from earth to the asteroid
    INPUTS:
        ts_np:  Numpy array of time snapshots at which to simulate the position.
                Shape [data_size, 1]; flattened to match ztf_elt Dataframe
        row_lengths_np: Number of observations for each element; shape [elt_batch_size]
        site_name: name of the observatory site for topos adjustment, e.g. 'geocenter' or 'palomar'
    """
    # Convert ts and row_lengths to keras constants
    # ts = keras.backend.constant(value=ts, shape=ts.shape, dtype=dtype)
    # row_lengths = keras.backend.constant(value=row_lengths_np, shape=row_lengths_np.shape, dtype=tf.int32)

    # Infer elt_batch_size from row_lengths
    elt_batch_size: int = row_lengths_np.shape[0]

    # Inputs: 6 orbital elements; epoch
    a = keras.Input(shape=(), batch_size=elt_batch_size, name='a')
    e = keras.Input(shape=(), batch_size=elt_batch_size, name='e')
    inc = keras.Input(shape=(), batch_size=elt_batch_size, name='inc')
    Omega = keras.Input(shape=(), batch_size=elt_batch_size, name='Omega')
    omega = keras.Input(shape=(), batch_size=elt_batch_size, name='omega')
    f = keras.Input(shape=(), batch_size=elt_batch_size, name='f')
    epoch = keras.Input(shape=(), batch_size=elt_batch_size, name='epoch')

    # Wrap these up into one tuple of inputs for the model
    inputs = (a, e, inc, Omega, omega, f, epoch)
    
    # Build asteroid direction layer
    ast_dir_layer = AsteroidDirection(ts_np=ts_np, row_lengths_np=row_lengths_np, 
                                      site_name=site_name, name='ast_dir_layer')

    # Define output tensors u, r by applying the position layer to the input tensors with the elements 
    u_flat, r_flat = ast_dir_layer(a, e, inc, Omega, omega, f, epoch)

    # Convert u, r to ragged tensors matching the element batch
    # u = tf.RaggedTensor.from_row_lengths(values=u_flat, row_lengths=row_lengths, name='u')
    # r = tf.RaggedTensor.from_row_lengths(values=r_flat, row_lengths=row_lengths, name='r')
    ragged_map_func = lambda x : tf.RaggedTensor.from_row_lengths(values=x, row_lengths=row_lengths_np)
    u = tf.keras.layers.Lambda(function=ragged_map_func, name='u')(u_flat)
    r = tf.keras.layers.Lambda(function=ragged_map_func, name='r')(r_flat)
    
    # Wrap the outputs
    outputs = (u, r)
    
    # Wrap this into a model
    model = keras.Model(inputs=inputs, outputs=outputs, name='model_asteroid_dir')
    
    # Bind the asteroid direction layer and aasteroid position layer
    model.ast_dir_layer = ast_dir_layer
    model.ast_pos_layer = ast_dir_layer.q_layer
    
    return model
