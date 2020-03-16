"""
Harvard IACS Masters Thesis
Solar Asteroid Model: Predict the movement of a test particle (e.g. asteroid) in the solar system
using the Kepler approximation with the sun as a fixed central attractor.

Michael S. Emanuel
Sun Oct 13 11:56:50 2019
"""

# Core
import tensorflow as tf
import numpy as np
from silence_tensorflow import silence_tensorflow

# Astronomy
import astropy
from astropy.units import au, day

# Local imports
from orbital_element import MeanToTrueAnomaly, TrueToMeanAnomaly
from asteroid_data import make_dataset_ast_pos, make_dataset_ast_dir, get_earth_pos, get_sun_pos_vel
from asteroid_data import orbital_element_batch
from ra_dec import calc_topos
from tf_utils import gpu_grow_memory, Identity

# ********************************************************************************************************************* 
# Aliases
keras = tf.keras

# ********************************************************************************************************************* 
# Turn off massive amount of superfluous tensorflow warnings and status messages
silence_tensorflow()
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

    # @tf.function
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

        # See module OrbitalElements for original version that includes velocity as well
        # This is pared down for speed
        
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

        # Position components
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

        # Assemble the velocity vector
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
    def __init__(self, ts, batch_size: int, **kwargs):
        """
        INPUTS:
            ts: fixed tensor of time snapshots at which to simulate the position
            batch_size: the number of elements to simulate at a time, e.g. 64; not to be confused with traj_size!
        """
        super(AsteroidPosition, self).__init__(**kwargs)

        # Configuration for serialization
        self.cfg = {
            'ts': ts,
            'batch_size': batch_size,
        }

        # Get trajectory size from ts
        self.ts = ts
        self.traj_size = ts.shape[0]
        self.batch_size = batch_size

        # Take a one time snapshot of the sun's position and velocity at these times
        q_sun_np, v_sun_np = get_sun_pos_vel(ts)
        traj_size = ts.shape[0]
        q_sun_np = q_sun_np.reshape(1, traj_size, space_dims)
        v_sun_np = v_sun_np.reshape(1, traj_size, space_dims)

        # Convert q_sun and v_sun into keras constants
        self.q_sun = keras.backend.constant(q_sun_np, dtype=dtype, shape=q_sun_np.shape, name='q_sun')
        self.v_sun = keras.backend.constant(v_sun_np, dtype=dtype, shape=v_sun_np.shape, name='v_sun')

        # Reshape ts to (batch_size, traj_size, 1)
        target_shape = (-1, 1)

        # print(f'ts.shape = {ts.shape}')
        # First repeat ts batch_size times; now size is (traj_size, batch_size, 1)
        t_rep = keras.layers.RepeatVector(n=batch_size, name='ts_rep')(keras.backend.reshape(ts, target_shape))
        # print(f't_rep.shape = {t_rep.shape}')
        # Transpose axes to make shape (batch_size, traj_size, 1)
        self.t_vec = tf.transpose(t_rep, perm=(1,0,2))
        # print(f't_vec.shape = {self.t_vec.shape}')

        # The adjustment dq to correct the Kepler approximation to match the numerical integration
        # self.dq = keras.backend.variable(np.zeros((batch_size, self.traj_size, space_dims,)), dtype=tf.float32, name = 'dq')
        self.dq = tf.Variable(initial_value=np.zeros((batch_size, self.traj_size, space_dims,)), 
                              dtype=dtype, trainable=False, name='dq')

        # The adjustment vq to correct the Kepler approximation to match the numerical integration
        self.dv = tf.Variable(initial_value=np.zeros((batch_size, self.traj_size, space_dims,)), 
                              dtype=dtype, trainable=False, name='dv')

    def update_dq_dv(self, dq, dv):
        """Update the value of dq and dv"""
        self.dq.assign(dq)
        self.dv.assign(dv)

    def calibrate(self, elts, q_ast, v_ast):
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

    def call(self, a, e, inc, Omega, omega, f, epoch):
        """
        Simulate the orbital trajectories.  
        Snapshot times t shared by all the input elements.  
        The inputs orbital elements and reference epoch should all have size (batch_size,).
        Output is the barycentric position; the elements generate a heliocentric calculation, and the
        constant with the sun's position is added to it.
        """
        # Alias traj_size, batch_size for legibility
        traj_size = self.traj_size

        # Reshape epoch to (batch_size, traj_size, 1)
        target_shape = (-1, 1)
        epoch_vec = keras.layers.RepeatVector(n=traj_size, name='epoch_vec')(keras.backend.reshape(epoch, target_shape))
        
        # Subtract epoch from t_vec; now it is relative to the epoch
        t = keras.layers.subtract([self.t_vec, epoch_vec], name='t')        

        # Compute eccentric anomaly E from f and e
        M = TrueToMeanAnomaly(name='TrueToMeanAnomaly')([f, e])
        
        # Compute mean motion N from mu and a
        a3 = tf.math.pow(a, 3, name='a3')
        mu_over_a3 = tf.divide(mu, a3, name='mu_over_a3')
        N = tf.sqrt(mu_over_a3, name='N')

        # Reshape t to (batch_size, traj_size, 1)
        target_shape = (-1, 1)
        # ******************************************************************
        # Predict orbital elements over time
        
        # Repeat the constant orbital elements to be vectors of shape (batch_size, traj_size, 1)
        target_shape = (-1, 1)
        a_t = keras.layers.RepeatVector(n=traj_size, name='a_t')(keras.backend.reshape(a, target_shape))
        e_t = keras.layers.RepeatVector(n=traj_size, name='e_t')(keras.backend.reshape(e, target_shape))
        inc_t = keras.layers.RepeatVector(n=traj_size, name='inc_t')(keras.backend.reshape(inc, target_shape))
        Omega_t = keras.layers.RepeatVector(n=traj_size, name='Omega_t')(keras.backend.reshape(Omega, target_shape))
        omega_t = keras.layers.RepeatVector(n=traj_size, name='omega_t')(keras.backend.reshape(omega, target_shape))
        
        # Repeat initial mean anomaly M0 and mean motion N0 to match shape of outputs
        M0_t = keras.layers.RepeatVector(n=traj_size, name='M0_t')(keras.backend.reshape(M, target_shape))
        N0_t = keras.layers.RepeatVector(n=traj_size, name='N0_t')(keras.backend.reshape(N, target_shape))
        # Compute the mean anomaly M(t) as a function of time
        N_mult_t = keras.layers.multiply(inputs=[N0_t, t])
        M_t = keras.layers.add(inputs=[M0_t, N_mult_t])
    
        # Compute the true anomaly from the mean anomly and eccentricity
        f_t = MeanToTrueAnomaly(name='mean_to_true_anomaly')([M_t, e_t])
    
        # Wrap orbital elements into one tuple of inputs for layer converting to cartesian coordinates
        elt_t = (a_t, e_t, inc_t, Omega_t, omega_t, f_t,)
        
        # Convert orbital elements to heliocentric cartesian coordinates
        q_helio, v_helio = ElementToPosition(name='q_helio')(elt_t)
        # Add solar position to get q in barycentric coordinates
        q_bary = tf.add(q_helio, self.q_sun)
        # Add solar velocity get v in barycentric coordinates
        v_bary = tf.add(v_helio, self.v_sun)
   
        # The estimated barycentric position and velocity includes the optional correction factor dq
        q = tf.add(q_bary, self.dq)
        v = tf.add(v_bary, self.dv)

        return q, v

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class AsteroidDirection(keras.layers.Layer):
    """
    Layer to compute the direction from earth to asteroid.
    """
    def __init__(self, ts, site_name: str, batch_size: int, **kwargs):
        """
        INPUTS:
            ts: fixed tensor of time snapshots at which to simulate the position
            site_name: name of the observatory site, used for topos adjustment, e.g. 'geocenter' or 'palomar'
            batch_size: the number of elements to simulate at a time, e.g. 64; not to be confused with traj_size!
        """
        super(AsteroidDirection, self).__init__(**kwargs)

        # Configuration for serialization
        self.cfg = {
            'ts': ts,
            'site_name': site_name,
            'batch_size': batch_size,
        }
        # Save the times
        self.ts = ts
        
        # Build layer to compute positions
        self.q_layer = AsteroidPosition(ts=ts, batch_size=batch_size, name='q_ast')
        
        # Take a one time snapshot of the earth's position at these times in barycentric coordinates
        q_earth_np = get_earth_pos(ts)
        traj_size = ts.shape[0]
        q_earth_np = q_earth_np.reshape(1, traj_size, space_dims)
        # self.q_earth = keras.backend.constant(q_earth_np, dtype=tf.float32, shape=q_earth_np.shape, name='q_earth')

        # Take a one time snapshot of the topos adjustment; displacement from geocenter to selected observatory
        dq_topos_ap = calc_topos(obstime_mjd=ts, site_name=site_name)
        # Convert dq_topos to a numpy array with units of au
        dq_topos_np = dq_topos_ap.to(au).value

        # Position of the observatory in barycentric frame as a Keras constant
        q_obs_np = q_earth_np + dq_topos_np
        self.q_obs = keras.backend.constant(q_obs_np, dtype=dtype, shape=q_obs_np.shape, name='q_obs')

    def calibrate(self, elts, q_ast, v_ast):
        """Calibrate this model by calibrating the underlying position layer"""
        # Calibrate the position model
        self.q_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)

    def call(self, a, e, inc, Omega, omega, f, epoch):
        # Calculate position and velocity of the asteroid
        q_ast, v_ast = self.q_layer(a, e, inc, Omega, omega, f, epoch)

        # Relative displacement from observatory to asteroid; instantaneous, before light time adjustment
        q_rel_inst = tf.subtract(q_ast, self.q_obs, name='q_rel_inst')
        # Distance between earth and asteroid, before light time adjustment
        r_inst = tf.norm(q_rel_inst, axis=-1, keepdims=True, name='r_earth_inst')

        # Light time in days from asteroid to earth in days (time units is days)
        light_time = tf.divide(r_inst, light_speed_au_day)
        # Adjusted relative position, accounting for light time; simulation velocity units are AU /day
        dq_lt = tf.multiply(v_ast, light_time)
        q_rel = tf.subtract(q_rel_inst, dq_lt)
        # Adjusted distance to earth, accounting for light time
        r = tf.norm(q_rel, axis=-1, keepdims=True, name='r')
        # Direction from earth to asteroid as unit vectors u = (ux, uy, uz)    
        u = tf.divide(q_rel, r)

        return u, r
    
    def get_config(self):
        return self.cfg
    
# ********************************************************************************************************************* 
# Functional API Models
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def make_model_ast_pos(ts: tf.Tensor, batch_size:int =64) -> keras.Model:
    """
    Compute orbit positions for asteroids in the solar system from
    the initial orbital elements with the Kepler model.
    Factory function that returns a functional model.
    Inputs for the model are 6 orbital elements, the epoch, and the desired times for position outputs.
    Outputs of the model are the position of the asteroid relative to the sun.
    INPUTS;
        ts: times to evaluate asteroid position in heliocentric coordinates
        batch_size: defaults to None for variable batch size
    """
    # Inputs: 6 orbital elements; epoch;
    a = keras.Input(shape=(), batch_size=batch_size, name='a')
    e = keras.Input(shape=(), batch_size=batch_size, name='e')
    inc = keras.Input(shape=(), batch_size=batch_size, name='inc')
    Omega = keras.Input(shape=(), batch_size=batch_size, name='Omega')
    omega = keras.Input(shape=(), batch_size=batch_size, name='omega')
    f = keras.Input(shape=(), batch_size=batch_size, name='f')
    epoch = keras.Input(shape=(), batch_size=batch_size, name='epoch')

    # Wrap these up into one tuple of inputs for the model
    inputs = (a, e, inc, Omega, omega, f, epoch)
    
    # Output times are a constant
    ts = keras.backend.constant(ts, name='ts')

    # Call asteroid position layer
    ast_pos_layer = AsteroidPosition(ts, batch_size, name='ast_pos_layer')
    q, v = ast_pos_layer(a, e, inc, Omega, omega, f, epoch)
    # Alias outputs
    q = Identity(name='q')(q)
    v = Identity(name='v')(v)
    
    # Wrap up the outputs
    outputs = (q, v,)

    # Wrap this into a model
    model = keras.Model(inputs=inputs, outputs=outputs, name='model_asteroid_pos')
    
    # Bind the asteroid position layer
    model.ast_pos_layer = ast_pos_layer
    
    return model


# ********************************************************************************************************************* 
def make_model_ast_dir(ts: tf.Tensor, site_name: str = 'geocenter', batch_size:int =64) -> keras.Model:
    """
    Compute direction from earth to asteroids in the solar system from
    the initial orbital elements with the Kepler model.
    Factory function that returns a functional model.
    Inputs for the model are 6 orbital elements, the epoch, and the desired times for position outputs.
    Outputs of the model are the unit vector (direction) pointing from earth to the asteroid
    INPUTS;
        ts: times to evaluate asteroid direction from earth
        site_name: name of the observatory site for topos adjustment, e.g. 'geocenter' or 'palomar'
        batch_size: defaults to None for variable batch size
    """
    # Inputs: 6 orbital elements; epoch
    a = keras.Input(shape=(), batch_size=batch_size, name='a')
    e = keras.Input(shape=(), batch_size=batch_size, name='e')
    inc = keras.Input(shape=(), batch_size=batch_size, name='inc')
    Omega = keras.Input(shape=(), batch_size=batch_size, name='Omega')
    omega = keras.Input(shape=(), batch_size=batch_size, name='omega')
    f = keras.Input(shape=(), batch_size=batch_size, name='f')
    epoch = keras.Input(shape=(), batch_size=batch_size, name='epoch')

    # Wrap these up into one tuple of inputs for the model
    inputs = (a, e, inc, Omega, omega, f, epoch)
    
    # Output times are a constant
    ts = keras.backend.constant(ts, name='ts')

    # All the work done in a single layer
    # u = AsteroidDirection(ts, batch_size, name='u')(a, e, inc, Omega, omega, f, epoch)
    ast_dir_layer = AsteroidDirection(ts=ts, site_name=site_name, batch_size=batch_size, name='ast_dir_layer')
    u, r = ast_dir_layer(a, e, inc, Omega, omega, f, epoch)

    # Alias outputs
    u = Identity(name='u')(u)
    r = Identity(name='r')(r)

    # Wrap the outputs
    outputs = (u, r)
    
    # Wrap this into a model
    model = keras.Model(inputs=inputs, outputs=outputs, name='model_asteroid_dir')
    
    # Bind the asteroid direction layer and aasteroid position layer
    model.ast_dir_layer = ast_dir_layer
    model.ast_pos_layer = ast_dir_layer.q_layer
    
    return model
