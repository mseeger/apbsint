# -------------------------------------------------------------------
# CPTANNOTATE_EXT
# -------------------------------------------------------------------
# Header file for ptannotate_ext extension module.
# External class declarations.
# Author: Matthias Seeger
# -------------------------------------------------------------------

# External classes. These are wrapped as potential annotation types in
# ptannotate_ext.

#cdef extern from "src/eptools/XXX/YYY.h":
#    cdef cppclass YYY:
#        YYY(...) except +

IF INCLUDE_WORKAROUND:
   include "cptannotate_ext_workaround.pxi"
