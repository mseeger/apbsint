/* -------------------------------------------------------------------
 * EPTWRAP_DEBUG_CASTANNOBJ
 *
 * Debug function to test the cast
 *   np.uint64 -> void* -> QuadratureServices
 * If the pointer is non-zero, we call 'debug_method'.
 *
 * Input:
 * - ANNOBJ: void*
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_debug_castannobj.h"
#include "src/eptools/potentials/quad/QuadratureServices.h"

void eptwrap_debug_castannobj(void* annobj,W_ERRORARGS)
{
  try {
    if (annobj==0)
      W_RETERROR(1,"ANNOBJ is NULL");
    ((QuadratureServices*) annobj)->debug_method();
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
