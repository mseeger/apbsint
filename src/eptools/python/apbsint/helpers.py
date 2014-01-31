"""
helpers
=======

Collection of internal helper functions and classes, used by different
modules.

"""

import numpy as np

__all__ = ['check_vecsize', 'maxreldiff', 'Struct']

def check_vecsize(v,n=None):
    """
    Check whether 'v' is a 1D numpy array. If 'n' is given, also check
    whether its length is equal to that.
    """
    return (isinstance(v,np.ndarray) and v.ndim == 1 and
            (n is None or v.shape[0] == n))

def maxreldiff(a,b):
    """
    Returns maximum relative difference (component-wise) over all entries of
    numpy arrays 'a', 'b' (same dimensionalities).
    """
    return (np.abs(a-b)/np.maximum(np.maximum(np.abs(a),np.abs(b)),1e-8)).max()

# Provides type for structure without class attributes, can have any number
# and type of object attributes.
class Struct:
    pass
