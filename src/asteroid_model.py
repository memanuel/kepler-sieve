"""
Harvard IACS Masters Thesis
Solar Asteroid Model: Predict the movement of a test particle (e.g. asteroid) in the solar system
using the Kepler approximation with the sun as a fixed central attractor.

Michael S. Emanuel
Sun Oct 13 11:56:50 2019
"""

# Library imports
import tensorflow as tf
import numpy as np
import time
from silence_tensorflow import silence_tensorflow

# Local imports
# from tf_utils import Identity
from orbital_element import MeanToTrueAnomaly, TrueToMeanAnomaly
from asteroid_data import make_dataset_ast_pos, make_dataset_ast_dir, get_earth_pos, get_sun_pos_vel
from asteroid_data import orbital_element_batch
from asteroid_integrate import calc_ast_pos
from  tf_utils import gpu_grow_memory, Identity
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
   
        # The estimated barycentric position includes the optional correction factor dq
        q = tf.add(q_bary, self.dq)
        v = tf.add(v_bary, self.dv)

        return q, v

    def get_config(self):
        return self.cfg

# ********************************************************************************************************************* 
class DirectionUnitVector(keras.layers.Layer):
    """
    Layer to compute the direction from object 1 (e.g. earth) to object 2 (e.g. asteroid)
    """
    
    def __init__(self, **kwargs):
        super(DirectionUnitVector, self).__init__(**kwargs)

    # don't declare this tf.function because it breaks when using it with q_earth
    # still not entirely sure how tf.function works ...
    # @tf.function
    def call(self, q1, q2):
        # Relative displacement from earth to asteroid
        q_rel = tf.subtract(q2, q1, name='q_rel')
        # Distance between objects
        r = tf.norm(q_rel, axis=-1, keepdims=True, name='r')
        # Unit vector pointing from object 1 to object 2
        u = tf.divide(q_rel, r, name='q_rel_over_r')
        return u
    
    def get_config(self):
        return dict()       


## ********************************************************************************************************************* 
class AsteroidDirection(keras.layers.Layer):
    """
    Layer to compute the direction from earth to asteroid.
    """
    def __init__(self, ts, batch_size: int, **kwargs):
        """
        INPUTS:
            ts: fixed tensor of time snapshots at which to simulate the position
            batch_size: the number of elements to simulate at a time, e.g. 64; not to be confused with traj_size!
        """
        super(AsteroidDirection, self).__init__(**kwargs)

        # Configuration for serialization
        self.cfg = {
            'ts': ts,
            'batch_size': batch_size,
        }
        self.ts = ts
        
        # Build layer to compute positions
        self.q_layer = AsteroidPosition(ts=ts, batch_size=batch_size, name='q_ast')
        
        # Take a one time snapshot of the earth's position at these times; done in heliocentric coordinates
        q_earth_np = get_earth_pos(ts)
        traj_size = ts.shape[0]
        q_earth_np = q_earth_np.reshape(1, traj_size, space_dims)
        self.q_earth = keras.backend.constant(q_earth_np, dtype=tf.float32, shape=q_earth_np.shape, name='q_earth')
        
    def calibrate(self, elts, q_ast, v_ast):
        """Calibrate this model by calibrating the underlying position layer"""
        # Calibrate the position model
        self.q_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)

    def call(self, a, e, inc, Omega, omega, f, epoch):
        # Calculate position
        q_ast = self.q_layer(a, e, inc, Omega, omega, f, epoch)

        # Unit displacement vector (direction) from earth to asteroid
        u = DirectionUnitVector(name='u')(self.q_earth, q_ast)
        
        return u
    
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
def make_model_ast_dir(ts: tf.Tensor, batch_size:int =64) -> keras.Model:
    """
    Compute direction from earth to asteroids in the solar system from
    the initial orbital elements with the Kepler model.
    Factory function that returns a functional model.
    Inputs for the model are 6 orbital elements, the epoch, and the desired times for position outputs.
    Outputs of the model are the unit vector (direction) pointing from earth to the asteroid
    INPUTS;
        ts: times to evaluate asteroid direction from earth
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
    ast_dir_layer = AsteroidDirection(ts, batch_size, name='u')
    u = ast_dir_layer(a, e, inc, Omega, omega, f, epoch)

    # Wrap the outputs
    outputs = (u,)
    
    # Wrap this into a model
    model = keras.Model(inputs=inputs, outputs=outputs, name='model_asteroid_dir')
    
    # Bind the asteroid direction layer and aasteroid position layer
    model.ast_dir_layer = ast_dir_layer
    model.ast_pos_layer = ast_dir_layer.q_layer
    
    return model

# ********************************************************************************************************************* 
# Tests
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def test_ast_pos_layer():
    """Test custom layer for asteroid positions"""
    ast_pos_layer = AsteroidPosition(batch_size=64)
    ts = np.arange(51544, 58744, dtype=np.float32)
    a,e,inc,Omega,omega,f, epoch = orbital_element_batch(1).values()
    q_ast, v_ast = ast_pos_layer(ts,a,e,inc,Omega,omega,f,epoch)
    return q_ast, v_ast

# ********************************************************************************************************************* 
def test_ast_pos() -> bool:
    """Test asteroid position model"""
    # Load data for the first 1000 asteroids
    ds: tf.data.Dataset = make_dataset_ast_pos(n0=0, num_files=1, include_vel=True)
    # Get reference times
    batch_in, batch_out = list(ds.take(1))[0]
    ts = batch_in['ts'][0]
    batch_size = 64
    # Create the model to predict asteroid trajectories
    model: keras.Model = make_model_ast_pos(ts=ts, batch_size=batch_size)
    # Compile with MSE (mean squared error) loss
    model.compile(loss='MSE')
    # Evaluate this model
    mse_qv, mse_q, mse_v = model.evaluate(ds)
    # RMS errors for q and v
    rmse_q: float = np.sqrt(mse_q)
    rmse_v: float = np.sqrt(mse_v)
    # Threshold for passing
    thresh_q: float = 0.125
    thresh_v: float = 0.125
    isOK_mse: bool = (rmse_q < thresh_q) and (rmse_v < thresh_v)
    # Report results
    msg: str = 'PASS' if isOK_mse else 'FAIL'
    print(f'Root MSE for asteroid model on first 1000 asteroids:')
    print(f'q in AU     = {rmse_q:8.6f}')
    print(f'v in AU/day = {rmse_v:8.6f}')
    print(f'***** {msg} *****')

    # Evaluate q, v on first batch before adjustments
    q_true = batch_out['q']
    v_true = batch_out['v']
    q_pred, v_pred = model.predict(batch_in)
    err_q_pre = np.mean(np.linalg.norm(q_pred - q_true, axis=2))
    err_v_pre = np.mean(np.linalg.norm(v_pred - v_true, axis=2))

    # Assemble orbital elements
    elts = {
        'a': batch_in['a'].numpy(),
        'e': batch_in['e'].numpy(),
        'inc': batch_in['inc'].numpy(),
        'Omega': batch_in['Omega'].numpy(),
        'omega': batch_in['omega'].numpy(),
        'f': batch_in['f'].numpy(),
        'epoch': batch_in['epoch'].numpy(),
        }
    # Epoch must be constant in a batch
    epoch = batch_in['epoch'][0].numpy()
    
    # Compute numerical orbit for calibration
    t0 = time.time()
    q_ast, q_earth, v_ast = calc_ast_pos(elts=elts, epoch=epoch, ts=ts)
    t1 = time.time()
    calc_time = t1 - t0
    print(f'\nNumerical integration of orbits took {calc_time:5.3f} seconds.')
    
    # Calibrate the position model
    model.ast_pos_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)
    # Evaluate error after calibration
    q_pred, v_pred = model.predict(batch_in)
    err_q_post = np.mean(np.linalg.norm(q_pred - q_true, axis=2))
    err_v_post = np.mean(np.linalg.norm(v_pred - v_true, axis=2))

    # Relative errors
    err_q_pre_rel = err_q_pre / np.mean(np.linalg.norm(q_true, axis=2))
    err_v_pre_rel = err_v_pre / np.mean(np.linalg.norm(v_true, axis=2))
    err_q_post_rel = err_q_post / np.mean(np.linalg.norm(q_true, axis=2))
    err_v_post_rel = err_v_post / np.mean(np.linalg.norm(v_true, axis=2))

    # Report results for position after calibration
    thresh_q = 2.0E-6
    isOK_q: bool = (err_q_post < thresh_q)
    msg = 'PASS' if isOK_q else 'FAIL'
    print(f'\nMean position error on first batch of 64 asteroids in AU:')
    print(f'Before calibration: {err_q_pre:5.3e} ({err_q_pre_rel:5.3e} relative)')
    print(f'After calibration:  {err_q_post:5.3e} ({err_q_post_rel:5.3e} relative)')
    print(f'***** {msg} *****')

    # Report results for velocity after calibration
    thresh_v = 1.0E-8
    isOK_v: bool = (err_v_post < thresh_v)
    msg = 'PASS' if isOK_v else 'FAIL'
    print(f'\nMean velocity error on first batch of 64 asteroids in AU/day:')
    print(f'Before calibration: {err_v_pre:5.3e} ({err_v_pre_rel:5.3e} relative)')
    print(f'After calibration:  {err_v_post:5.3e} ({err_v_post_rel:5.3e} relative)')
    print(f'***** {msg} *****')

    return (isOK_mse and isOK_q and isOK_v)

# ********************************************************************************************************************* 
def test_ast_dir() -> bool:
    """Test the asteroid direction model"""
    # Load data for the first 1000 asteroids
    ds: tf.data.Dataset = make_dataset_ast_dir(0, 1)
    # Get reference times
    batch_in, batch_out = list(ds.take(1))[0]
    ts = batch_in['ts'][0]
    # Create the model to predict asteroid trajectories
    model: keras.Model = make_model_ast_dir(ts=ts)
    # Set number of steps
    num_ast: int = 1000
    batch_size: int = 64
    steps = num_ast // batch_size
    # Compile with MSE (mean squared error) loss
    model.compile(loss='MSE')
    # Evaluate this model
    mse: float = model.evaluate(ds, steps=steps)
    rmse: float = np.sqrt(mse)
    # Convert error from unit vector to angle
    rmse_rad = 2.0 * np.arcsin(rmse / 2.0)
    rmse_deg = np.rad2deg(rmse_rad)
    rmse_sec = rmse_deg * 3600
    # Threshold for passing
    thresh: float = 2.5
    isOK_1: bool = (rmse_deg < thresh)
    
    # Report results
    msg: str = 'PASS' if isOK_1 else 'FAIL'
    print(f'MSE for asteroid model on first 1000 asteroids = {mse:8.6f}')
    print(f'Angle error = {rmse_rad:5.3e} rad / {rmse_deg:8.6f} degrees / {rmse_sec:6.2f} arc seconds')
    print(f'***** {msg} *****')

    # Evaluate on first batch before adjustments
    u1_true = batch_out['u']
    u1_pred = model.predict(batch_in)
    # Error in degrees
    err1_pre = np.rad2deg(np.mean(np.linalg.norm(u1_pred - u1_true, axis=2)))
    
    # Assemble orbital elements for calibration
    elts = {
        'a': batch_in['a'].numpy(),
        'e': batch_in['e'].numpy(),
        'inc': batch_in['inc'].numpy(),
        'Omega': batch_in['Omega'].numpy(),
        'omega': batch_in['omega'].numpy(),
        'f': batch_in['f'].numpy(),
        'epoch': batch_in['epoch'].numpy(),
        }
    # Epoch must be common in a batch to allow for numerical integration
    epoch = elts['epoch'][0]
    
    # Compute correction factor dq; also time this operation
    t0 = time.time()
    q_ast, q_earth, v_ast = calc_ast_pos(elts=elts, epoch=epoch, ts=ts)
    t1 = time.time()
    calc_time = t1 - t0
    print(f'\nNumerical integration of orbits took {calc_time:5.3f} seconds.')
    
    # Evaluate error after correction
    model.ast_dir_layer.calibrate(elts=elts, q_ast=q_ast, v_ast=v_ast)
    u1_pred = model.predict(batch_in)
    err1_post = np.rad2deg(np.mean(np.linalg.norm(u1_pred - u1_true, axis=2)))
    err1_post_sec = err1_post * 3600
    
    # Report results
    thresh = 1.0
    isOK_2: bool = (err1_post_sec < thresh)
    msg = 'PASS' if isOK_2 else 'FAIL'
    print(f'Mean angle error on first batch of 64 asteroids in degrees:')
    print(f'Before calibration: {err1_pre:5.3e}')
    print(f'After calibration:  {err1_post:5.3e}  ({err1_post_sec:5.3f} arc-seconds)')
    print(f'***** {msg} *****')
    
    return (isOK_1 and isOK_2)

# ********************************************************************************************************************* 
def main():
    test_ast_pos()
    # test_ast_dir()
    
# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()

