"""
Michael S. Emanuel
Fri Oct 25 19:20:56 2019
"""
# Library imports
import tensorflow as tf
import numpy as np

# Local imports
from observation_data import random_direction
from tf_utils import tf_quiet

# ********************************************************************************************************************* 
# Functions for scoring predicted asteroid directions vs. observations
# ********************************************************************************************************************* 

# ********************************************************************************************************************* 
def score_mean(A: tf.Tensor, thresh: float=2.0):
    """
    Expected value of the score function 
    f(s) = exp(-1/2 A s^2)
    integrated on the region where s^2 < thresh^2
    s is the Euclidean distance between the predicted and observed direction; 
    s^2 = dx^2 + dy^2 + dz^2 with dx = u_obs_x - u_pred_x, etc.
    thesh is the maximum Cartesian distance considered before observations are cut off.
    If we use cylindrical coordinates (z, phi), with r^2 +z^2 = 1, x^2 + y^2 = r^2, and
    x = r cos(phi), y = r sin(phi), z = z
    and set the true direction at (0, 0, 1), we can get the EV by integrating over the sphere.
    s is a right triangle with sides (1-z) and r, regardless of phi, so
    s^2 = (1-z)^2 + r^2 = 1 -2z + z^2 + r^2 = 1 + (r^2+z^2) - 2z = 2 - 2z = 2(1-z)
    The integrand is therefore
    f(z) = exp(a(z-1))
    This can be integrated symbolically to F(z) = exp(a(z-1))/a
    Need to divide by 2 for the expected value on interval z=[-1, +1]
    With a threshold, the integral is over z ranging from 1-thresh^2 / 2 to 1
    Note that when thresh=2, z_min = 1 - 2^2/2 = -1 as usual.
    """
    # The expected value is (1 - exp(-0.5 * A * thresh^2)) / 2A
    # return (1.0 - np.exp(-2.0*A)) / (2*A)
    # Save original data type and cast to double
    try:
        dtype = A.dtype
    except:
        dtype = tf.float32
    A = tf.cast(A, tf.float64)

    # Do calculation in double
    # expm1_m2a = tf.math.expm1(minus_two_A, name='expm1_m2a')   
    # Coefficient that multiplies A
    coefficient = -0.5 * thresh**2
    # Argument of the exponential function
    arg = tf.multiply(coefficient, A)
    # Use the built-in exp minus one function for better precision
    # ans = (1 - exp(arg)) / (2A) = (exp(arg)-a) / (-2A) = expm1(arg) / (-2A)
    expm1_arg = tf.math.expm1(arg, name='expm1_arg')
    # The denominator in fraction with expm1_arg
    minus_two_A = tf.multiply(A, -2.0, name='minus_two_A')
    # Take the mean over the batch
    mean = tf.divide(expm1_arg, minus_two_A)
    return tf.cast(mean, dtype)

# ********************************************************************************************************************* 
def score_var(A: tf.Tensor, thresh: float = 2.0):
    """
    Variance of the score function 
    f(s) = exp(-1/2 A s^2)
    integrated on the region where s^2 < thresh^2
    This can be integrated symbolically.
    Because the integrand f_a(z) = exp(a(z-1)), f(z)^2 = f_2a(z)
    If mu(a, thresh) = E_T[f(z)] with the given threshold, then
    E[f(z)^2] = mu(2a, thresh) so 
    Var[f(z)] = E[f(z)^2] - E[f(z)]^2 = mu(2a, thresh) - mu(a, thresh)**2
    """
    # Expected value E[f(z)^2]
    e_f2 = score_mean(2.0*A, thresh=thresh)
    # Expected value of E[f(z)]^2
    ef_2 = tf.square(score_mean(A, thresh=thresh))
    # The variance is the difference E[f(z)^2] - E[f(z)]^2
    var = e_f2 - ef_2

    return var

# ********************************************************************************************************************* 
def score_mean_var(A: tf.Tensor, thresh: float = 2.0):
    """
    Mean and variance of the score function in one call for efficiency
    f(s) = exp(-1/2 A s^2)
    integrated on the region where s^2 < thresh^2
    This can be integrated symbolically.
    Because the integrand f_a(z) = exp(a(z-1)), f(z)^2 = f_2a(z)
    If mu(a, thresh) = E_T[f(z)] with the given threshold, then
    E[f(z)^2] = mu(2a, thresh) so 
    Var[f(z)] = E[f(z)^2] - E[f(z)]^2 = mu(2a, thresh) - mu(a, thresh)**2
    """
    # The mean
    mean = score_mean(A, thresh=thresh)
    # Expected value E[f(z)^2]
    e_f2 = score_mean(2.0*A, thresh=thresh)
    # Expected value of E[f(z)]^2
    ef_2 = tf.square(mean)
    # The variance is the difference E[f(z)^2] - E[f(z)]^2
    var = e_f2 - ef_2

    return mean, var

# ********************************************************************************************************************* 
def score_std(A: tf.Tensor, thresh: float = 2.0):
    """
    Standard deviation of the score function 
    f(s) = exp(-1/2 A s^2)
    integrated on the region where s^2 < thresh^2
    This can be integrated symbolically.
    """
    # The variance is (A*Coth(A)) * E[f]^2
    # The std deviation is sqrt(A*Coth(A)) * E[f]
    # B = (A / tf.tanh(A)) - 1.0
    # Save original data type and cast to double
    
    # Delegate to score_var to compute the variance
    var = score_var(A, thresh=thresh)
    # Standard deviation is square root of variance
    std = tf.sqrt(var)
    return std

# ********************************************************************************************************************* 
def test_score_moments_one(A: float, thresh: float, mean_exp: float, var_exp: float) -> bool:
    """Test the mean and variance for the score function with known values"""
    # Calculate mean and variance
    mean_calc, var_calc = score_mean_var(A, thresh=thresh)

    # Description of this test
    description = f'A = {A:10.6f}, thresh={thresh:10.6f}'

    # Test the mean
    # mean_calc = score_mean(A)
    isOK: bool = np.isclose(mean_exp, mean_calc, atol=1.0E-7)     
    if not isOK:
        print(f'Failed on {description}, expected mean {mean_exp}, got mean {mean_calc}.')

    # Test the variance
    # var_calc = score_var(A)
    isOK = isOK and np.isclose(var_exp, var_calc, atol=1.0E-7)
    if not isOK:
        print(f'Failed on {description}, expected var {var_exp}, got var {var_calc}.')

    # Message if we passed
    if isOK:
        print(f'Passed on {description}.  mean = {mean_calc:8.6f}, var = {var_calc:8.6f}.')

    return isOK

# ********************************************************************************************************************* 
def test_score_moments():
    """Test the mean and variance for the score function with known values"""
    # When A = 1.0, mean = 0.432332, var = 0.0585098
    A = np.array(1.0)
    thresh = 2.0
    mean_exp = 0.432332
    var_exp = 0.0585098
    isOK: bool = test_score_moments_one(A, thresh=thresh, mean_exp=mean_exp, var_exp=var_exp)
    
    # When A = 32.828063500117445 (10 degrees), mean = 0.0152309, var = 0.00738346
    res = np.deg2rad(10.0)
    A = np.array(1.0 / res**2)
    thresh = 2.0
    mean_exp = 0.015230870989335427
    var_exp = 0.00738346
    isOK = isOK and test_score_moments_one(A, thresh=thresh, mean_exp=mean_exp, var_exp=var_exp)
    
    # When A = 32.828063500117445 (10 degrees), thresh = 0.17453292519943295 (10 degrees)
    # mean = 0.005992880760175801, var = 0.004777958720806503
    res = np.deg2rad(10.0)
    A = np.array(1.0 / res**2)
    thresh = np.deg2rad(10.0)
    mean_exp = 0.005992880760175801
    var_exp = 0.004777958720806503
    isOK = isOK and test_score_moments_one(A, thresh=thresh, mean_exp=mean_exp, var_exp=var_exp)
    
    # Report results
    msg: str = 'PASS' if isOK else 'FAIL'
    print(f'\nTest score moments:\n***** {msg} *****')    
    return isOK

# ********************************************************************************************************************* 
def test_score_num(N: int, R_deg: float):
    """
    Numerical test of score_mean and score_var functions
    INPUTS:
        N: number of points to test
        R_deg: resolution factor in degrees
    """
    # Convert resolution factor to A
    R: float = np.deg2rad(R_deg)
    A: float = 1.0 / R**2
    # Calculated mean and variance
    mean_calc: float = score_mean(A)
    std_calc: float = score_std(A)

    # Numerical samples
    directions: np.array = random_direction(N)
    # Difference with north pole    
    north_pole = np.array([0.0, 0.0, 1.0])
    epsilon = directions- north_pole
    epsilon2 = np.sum(epsilon*epsilon, axis=-1)

    # Simulated mean and stdev
    scores = np.exp(-0.5 * A * epsilon2)
    mean_num: float = np.mean(scores)
    mean_error_rel = np.abs(mean_num-mean_calc) / mean_calc
    std_num: float = np.std(scores)
    std_error_rel = np.abs(std_num-std_calc) / std_calc

    # Test the mean
    rtol: float = 2.0 / N**0.4
    isOK: bool = np.isclose(mean_calc, mean_num, rtol=rtol)
    print(f'\nTest score function with N={N}, R={R_deg} degrees:')
    print(f'Mean:')
    print(f'calc      = {mean_calc:10.8f}')
    print(f'num       = {mean_num:10.8f}')
    print(f'error_rel = {mean_error_rel:10.8f}')
    
    # Test the standard deviation
    rtol = 2.0 / N**0.4
    isOK = isOK and np.isclose(mean_calc, mean_num, rtol=rtol)
    print(f'\nStandard Deviation:')
    print(f'calc      = {std_calc:10.8f}')
    print(f'num       = {std_num:10.8f}')
    print(f'error_rel = {std_error_rel:10.8f}')
    
    # Summary results    
    msg: str = 'PASS' if isOK else 'FAIL'
    print(f'***** {msg} *****')
    return isOK

# ********************************************************************************************************************* 
def main():
    # Run TensorFlow quietly
    tf_quiet()

    # Test the score function
    test_score_num(10000, 10.0)
    test_score_moments()

# ********************************************************************************************************************* 
if __name__ == '__main__':
    main()
