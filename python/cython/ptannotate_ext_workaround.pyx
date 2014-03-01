# -------------------------------------------------------------------
# PTANNOTATE_EXT
# -------------------------------------------------------------------
# Wrapper classes for potential annotation types.
# In the C++ code wrapped by eptools_ext, potential objects can be
# annotated by additional objects. The corresponding C++ class
# are wrapped by Python classes defined here, in order to keep them
# persistent.
# Author: Matthias Seeger
# -------------------------------------------------------------------

DEF INCLUDE_WORKAROUND = True

# Declarations
IF INCLUDE_WORKAROUND:
    # Workaround: Additional declarations
    import cptannotate_ext_workaround as pax

# Definitions
include "ptannotate_ext.pxi"
IF INCLUDE_WORKAROUND:
    # Workaround: Additional definitions
    include "ptannotate_ext_workaround.pxi"
