"""
Utilities for working with orbital elements.

Michael S. Emanuel
2021-01-29
"""

# Core
import numpy as np
# Utility
from collections import namedtuple

# ********************************************************************************************************************* 
# Named tuple data type for orbital elements with just the six 
# a, e, inc, Omega, omega, f
OrbitalElement_aeiOof = namedtuple('OrbitalElement', 'a e inc Omega omega f')

# Named tuple data type for orbital elements that includes the mean anomaly M
# a, e, inc, Omega, omega, f, M
OrbitalElement_aeiOofM = namedtuple('OrbitalElement', 'a e inc Omega omega f M')
