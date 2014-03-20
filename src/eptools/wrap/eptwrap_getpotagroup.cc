/* -------------------------------------------------------------------
 * EPTWRAP_GETPOTAGROUP
 *
 * Returns argument group ID for potential with ID PID. Argument
 * group IDs are defined in 'EPScalarPotential::atypeXXX'. All
 * potentials within one group have the same input and return
 * arguments (number, types).
 *
 * Input:
 * - PID:    Potential ID
 *
 * Return:
 * - AGID:   Argument group ID, or -1 if PID is not valid ID
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_getpotagroup.h"
#include "src/eptools/potentials/EPPotentialFactory.h"

void eptwrap_getpotagroup(int ain,int aout,int pid,int* agid,W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain!=1)
      W_RETERROR(2,"Need 1 input argument");
    if (aout!=1)
      W_RETERROR(2,"Need 1 return argument");
    if (EPPotentialFactory::isValidID(pid))
      *agid = EPPotentialFactory::getArgumentGroup(pid);
    else
      *agid = -1;
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
