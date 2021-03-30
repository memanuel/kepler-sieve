"""
Utilities for integrating Rebound simulations.

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Core
import numpy as np
import pandas as pd
import sqlalchemy

# Astronomy
import rebound

# Utilities
from datetime import datetime, timedelta

# UI
# from tqdm.auto import tqdm as tqdm_auto

# MSE imports
from astro_utils import mjd_to_datetime
from orbital_element import OrbitalElement_aeiOofM as OrbitalElement # a, e, inc, Omega, omega, f, M

# Typing
from typing import List, Tuple, Dict, Set, Optional

# ********************************************************************************************************************* 
def integrate_mjds(sim_epoch: rebound.Simulation, mjds: np.array, save_elements: bool, progbar: bool):
    """
    Perform an integration in rebound to a specific set of MJDs and save to Numpy arrays.
    INPUTS:
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjds:           array of MJDs on which to write simulation outputs
        save_elements:  flag indicating whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        Numpy arrays with output of integration
        body_ids:       N integer body IDs of the bodies that were integrated
        body_names:     N body names
        q:              MxNx3 array of positions in AU
        v:              MxNx3 array of velocities in AU / day
    """
    # Look up the epoch from the base simulation
    epoch: float = float(sim_epoch.epoch)
    # print(f'epoch={epoch}')

    # Look up body IDs and names
    body_ids: np.ndarray = sim_epoch.body_ids
    body_names: np.ndarray = sim_epoch.body_names

    # Time is measured as an offset from the epoch of the base simulation
    t_epoch: float = 0.0
    
    # Create copies of the simulation to integrate forward and backward
    sim_fwd: rebound.Simulation = sim_epoch.copy()
    sim_back: rebound.Simulation = sim_epoch.copy()

    # Set the time counter on both simulation copies to the epoch time
    sim_fwd.t = t_epoch
    sim_back.t = t_epoch
    # Flip sign of dt on sim_back
    sim_back.dt *= -1.0

    # Set the times for snapshots in both directions;
    ts: np.array = mjds - epoch

    # Get the index placement for the left vs. right sides.
    # Left side will be integrated backwards from epoch, right side integrated forward from epoch
    idx: int = np.searchsorted(ts, t_epoch, side='left')
    ts_fwd: np.array = ts[idx:]
    ts_back: np.array = ts[:idx][::-1]

    # Number of snapshots
    M_back: int = len(ts_back)
    M_fwd: int = len(ts_fwd)
    M: int = M_back + M_fwd
    # Number of particles
    N: int = sim_epoch.N

    # Initialize arrays for the position and velocity
    shape_qv: Tuple[int] = (M, N, 3)
    q: np.array = np.zeros(shape_qv, dtype=np.float64)
    v: np.array = np.zeros(shape_qv, dtype=np.float64)
    
    # Initialize arrays for orbital elements if applicable  
    if save_elements:
        # Arrays for a, e, inc, Omega, omega, f, M
        shape_elt: Tuple[int] = (M, N)
        orb_a: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)
        orb_e: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)
        orb_inc: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)
        orb_Omega: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)
        orb_omega: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)
        orb_f: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)
        orb_M: np.array = np.full(shape_elt, fill_value=np.nan, dtype=np.float64)

        # Wrap these into a named tuple
        elts: OrbitalElement = \
            OrbitalElement(a=orb_a, e=orb_e, inc=orb_inc,
                           Omega=orb_Omega, omega=orb_omega, f=orb_f, M=orb_M)
    else:
        elts = None

    # Subfunction: process one row of the loop    
    def process_row(sim: rebound.Simulation, t: float, row: int):
        # Integrate to the current time step with an exact finish time
        sim.integrate(t, exact_finish_time=1)
        
        # Serialize the position and velocity
        # While the rebound interface includes methods for the accelerations, these do NOT appear to be
        # supported at this time! They just return 0.0, even after calculating the orbit.
        sim.serialize_particle_data(xyz=q[row])
        sim.serialize_particle_data(vxvyvz=v[row])

        # Save the orbital elements if applicable
        if save_elements:
            # Compute the orbital elements vs. the sun as primary
            primary = sim.particles['Sun']
            orbits = sim.calculate_orbits(primary=primary, jacobi_masses=False)
            # Compute the elements of the Moon using the Earth as primary
            orb_moon = sim.particles['Moon'].calculate_orbit(primary=sim.particles['Earth'])
            # Save the elements on this date as an Nx7 array
            # This approach only traverses the (slow) Python list orbits one time
            elt_array = np.array([[orb.a, orb.e, orb.inc, orb.Omega, orb.omega, orb.f, orb.M] \
                                  for orb in orbits])
            # Handle the special entry for the moon
            orb = orb_moon
            elt_array[3, :] = np.array([orb.a, orb.e, orb.inc, orb.Omega, orb.omega, orb.f, orb.M])

            # Save the elements into the current row of the named orbital elements arrays
            # The LHS row mask starts at 1 b/c orbital elements are not computed for the first object (Sun)
            orb_a[row, 1:] = elt_array[:, 0]
            orb_e[row, 1:] = elt_array[:, 1]
            orb_inc[row, 1:] = elt_array[:, 2]
            orb_Omega[row, 1:] = elt_array[:, 3]
            orb_omega[row, 1:] = elt_array[:, 4]
            orb_f[row, 1:] = elt_array[:, 5]
            orb_M[row, 1:] = elt_array[:, 6]

    # Integrate the simulation forward in time
    idx_fwd: List[Tuple[int, float]] = list(enumerate(ts_fwd))
    if progbar:
        idx_fwd = tqdm_auto(idx_fwd)
    for i, t in idx_fwd:
        # Row index for position data
        row: int = M_back + i
        # Process this row
        process_row(sim=sim_fwd, t=t, row=row)
        
    # Integrate the simulation backward in time
    idx_back: List[Tuple[int, float]] = list(enumerate(ts_back))
    if progbar:
        idx_back = tqdm_auto(idx_back)
    for i, t in idx_back:
        # Row index for position data
        row: int = M_back - (i+1)
        # Process this row
        process_row(sim=sim_back, t=t, row=row)

    return body_ids, body_names, q, v, elts

# ********************************************************************************************************************* 
def integrate_numpy(sim_epoch: rebound.Simulation, 
                    mjd0: int, 
                    mjd1: int,
                    interval_p: int,
                    interval_q: int,
                    save_elements: bool,
                    progbar: bool):
    """
    Perform an integration in rebound and save to Numpy arrays.
    INPUTS:
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjd0:           the earliest MJD to simulate back to
        mjd1:           the latest MJD to simulate forward to
        interval_p:     Numerator of the time interval between steps in days
        interval_q:     Denominator of the time interval between steps in days
                        Time steps are saved to output array at multiples of tau = (interval_p / interval_q)
                        The internal time step of the simulation is adaptive
        save_elements:  flag indicating whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        Numpy arrays with output of integration
        body_ids:       N integer body IDs of the bodies that were integrated
        body_names:     N body names
        epochs:         M times as of which the integration was saved; MJDs
        q:              MxNx3 array of positions in AU
        v:              MxNx3 array of velocities in AU / day
    """
    # Look up the epoch from the base simulation
    epoch: int = sim_epoch.epoch

    # Convert epoch, start and end times relative to a base date of the simulation start
    # This way, time is indexed from t0=0 to t1 = (dt1-dt0)
    t0: float = 0.0
    t1: float = mjd1 - mjd0
    t_epoch: float = float(epoch) - mjd0
    
    # Set the times for snapshots in both directions;
    step_count: int = ((mjd1 - mjd0)*interval_q) // interval_p
    time_step: np.float64 = np.float64(interval_p) / np.float64(interval_q)
    ts: np.array = np.arange(step_count+1) * time_step

    # The mjds corresponding to the times in ts
    epochs: np.array = ts + (epoch - t_epoch)
    # dt0: datetime = mjd_to_datetime(mjd0)
    # epochs_dt: np.array = np.array([dt0 + timedelta(t) for t in ts])

    # Delegate to integrate_mjds
    body_ids, body_names, q, v, elts = \
        integrate_mjds(sim_epoch=sim_epoch, mjds=epochs, save_elements=save_elements, progbar=progbar)

    return body_ids, body_names, epochs, q, v, elts

# ********************************************************************************************************************* 
def integration_np2df(body_ids: np.array, body_names: np.array, epochs: np.array, 
                      q: np.array, v: np.array, elts: Optional[np.array]):
    """
    Arrange Numpy arrays with integration output into a Pandas DataFrame with one row per observation.\
    INPUTS:
        body_ids:       N integer body IDs of the bodies that were integrated
        body_names:     N body names
        epochs:         M times as of which the integration was saved; MJDs
        q:              MxNx3 array of positions in AU
        v:              MxNx3 array of velocities in AU / day
        elts:           MxNx7 array of orbital elements; None when not included
    OUTPUTS:
        df: DataFrame with columns
        BodyID:       N integer body IDs of the bodies that were integrated
        BodyName:     N body names
        mjd:          M times as of which the integration was saved; MJDs
        qx, qy, qx:   Positions in AU in the BME
        vx, vy, vz:   Velocities in AU / day
        a, e, inc, Omega, omega, f, M: orbital elements when included
    """
    # Array sizes
    M: int = epochs.shape[0]
    N: int = body_ids.shape[0]

    # The time stamps
    mjd = epochs.repeat(N)
    TimeID = np.rint(mjd*24*60).astype(np.int32)

    # The ID and name of each body
    BodyID = np.tile(body_ids, M)
    BodyName = np.tile(body_names, M)
    
    # The positions
    qx = q[:,:,0].flatten()
    qy = q[:,:,1].flatten()
    qz = q[:,:,2].flatten()

    # The velocities
    vx = v[:,:,0].flatten()
    vy = v[:,:,1].flatten()
    vz = v[:,:,2].flatten()    

    # Wrap into a dictionary
    data_dict = {
        'TimeID': TimeID,
        'BodyID': BodyID,
        # 'BodyName': BodyName,
        'mjd': mjd,
        'qx': qx,
        'qy': qy,
        'qz': qz,
        'vx': vx,
        'vy': vy,
        'vz': vz,
    }
    
    # Add the orbital elements if they were included
    if elts is not None:
        data_dict['a'] = elts.a.flatten()
        data_dict['e'] = elts.e.flatten()
        data_dict['inc'] = elts.inc.flatten()
        data_dict['Omega'] = elts.Omega.flatten()
        data_dict['omega'] = elts.omega.flatten()
        data_dict['f'] = elts.f.flatten()
        data_dict['M'] = elts.M.flatten()

    # Convert to a DataFrame
    df = pd.DataFrame(data_dict)
    return df

# ********************************************************************************************************************* 
def integrate_df(sim_epoch: rebound.Simulation, 
                 mjd0: int, 
                 mjd1: int, 
                 interval_p: int,
                 interval_q: int,
                 save_elements: bool,
                 progbar: bool) -> None:
    """
    Perform an integration in rebound and save to Pandas DataFrame.
    INPUTS:
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjd0:           the earliest MJD to simulate back to
        mjd1:           the latest MJD to simulate forward to
        interval_p:     Numerator of the time interval between steps in days
        interval_q:     Denominator of the time interval between steps in days
        save_elements:  flag indiciting whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        df: DataFrame with columns
        BodyID:       N integer body IDs of the bodies that were integrated
        BodyName:     N body names
        mjd:          M times as of which the integration was saved; MJDs
        qx, qy, qx:   Positions in AU in the BME
        vx, vy, vz:   Velocities in AU / day
    """
    # Delegate to integrate_numpy
    body_ids, body_names, epochs, q, v, elts = \
            integrate_numpy(sim_epoch=sim_epoch, mjd0=mjd0, mjd1=mjd1, 
            interval_p=interval_p, interval_q=interval_q, save_elements=save_elements, progbar=progbar)
    # Delegate to integration_np2df
    df = integration_np2df(body_ids=body_ids, body_names=body_names, epochs=epochs, q=q, v=v, elts=elts)

    return df

# ********************************************************************************************************************* 
def integrate_mjds_df(sim_epoch: rebound.Simulation, mjds: np.array, save_elements: bool, progbar: bool) -> None:
    """
    Perform an integration in rebound to a specific set of MJDs and save to Pandas DataFrame.
    INPUTS:
        sim_epoch:      rebound simulation object as of the epoch time; to be integrated in both directions
        mjds:           array of MJDs on which to write simulation outputs
        save_elements:  flag indicating whether to save orbital elements
        progbar:        flag - whether to display a progress bar
    OUTPUTS:
        df: DataFrame with columns
        BodyID:       N integer body IDs of the bodies that were integrated
        BodyName:     N body names
        mjd:          M times as of which the integration was saved; MJDs
        qx, qy, qx:   Positions in AU in the BME
        vx, vy, vz:   Velocities in AU / day
    """
    # Delegate to integrate_numpy
    body_ids, body_names, q, v, elts = \
            integrate_mjds(sim_epoch=sim_epoch, mjds=mjds, save_elements=save_elements, progbar=progbar)
    # Delegate to integration_np2df
    df = integration_np2df(body_ids=body_ids, body_names=body_names, epochs=mjds, q=q, v=v, elts=elts)

    return df
