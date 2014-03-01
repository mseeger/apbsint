# Included by ptannotate_ext_*.pyx (normal or workaround)

import cython

# Potential annotations are passed to/from C++ code as void* (they are
# wrapped in Python classes to keep the objects alive). We cast them to
# 'uintptr_t' when going through Python (using the Capsules API seems
# very complicated).
from libc.stdint cimport uintptr_t

# Cython classes, wrapping external classes.
# They are derived from the abstract PotentialAnnotation, whose 'getptr'
# method returns a void* (as 'uintptr_t') to the C++ object.

# TODO:
# - Proper treatment of exceptions thrown by the C++ classes. Right now,
#   we rely on the automatic 'except +' mechanism of Cython (creates a
#   'RuntimeError')

cdef class PotentialAnnotation:
    def getptr(self):
        return <uintptr_t>NULL
