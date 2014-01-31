# -------------------------------------------------------------------
# APBTEST_EXT
# -------------------------------------------------------------------
# Cython code wrapping external C++ functions.
# Collect functions which are used in test code for the C++ code
# (allows to write test code in Python)
# Author: Matthias Seeger
# -------------------------------------------------------------------

import cython
import numpy as np
cimport numpy as np

# Declarations: Static methods as functions

cdef extern from "src/eptools/potentials/SpecfunServices.h" namespace "SpecfunServices":
    double logCdfNormal(double z)
    double derivLogCdfNormal(double z)

# Cython functions

@cython.boundscheck(False)
@cython.wraparound(False)
def specfun_logcdfnormal(np.ndarray[np.double_t,ndim=1] z not None,
                         np.ndarray[np.double_t,ndim=1] res not None):
    cdef int i, sz
    # Ensure that input/output arguments are contiguous
    sz = z.shape[0]
    if not (res.flags.c_contiguous and res.shape[0]==sz):
        # HIER: Define own exception, say EptwrapError
        raise TypeError('RES must be contiguous array of size of Z')
    z = np.ascontiguousarray(z)
    # Loop: Compute function values
    for i in range(sz):
        res[i] = logCdfNormal(z[i])

@cython.boundscheck(False)
@cython.wraparound(False)
def specfun_derivlogcdfnormal(np.ndarray[np.double_t,ndim=1] z not None,
                              np.ndarray[np.double_t,ndim=1] res not None):
    cdef int i, sz
    # Ensure that input/output arguments are contiguous
    sz = z.shape[0]
    if not (res.flags.c_contiguous and res.shape[0]==sz):
        # HIER: Define own exception, say EptwrapError
        raise TypeError('RES must be contiguous array of size of Z')
    z = np.ascontiguousarray(z)
    # Loop: Compute function values
    for i in range(sz):
        res[i] = derivLogCdfNormal(z[i])
