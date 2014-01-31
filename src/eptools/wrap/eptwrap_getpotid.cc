/* -------------------------------------------------------------------
 * EPTWRAP_GETPOTID
 *
 * Potential names <--> IDs maintained in 'EPPotentialNamedFactory'
 *
 * Input:
 * - NAME:    Potential name
 *
 * Return:
 * - PID:     Potential ID, or -1 if NAME is not potential name
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_getpotid.h"
#include "src/eptools/potentials/EPPotentialNamedFactory.h"

void eptwrap_getpotid(int ain,int aout,char* name,int* pid,W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain!=1)
      W_RETERROR(2,"Need 1 input argument");
    if (aout!=1)
      W_RETERROR(2,"Need 1 return argument");
    *pid = EPPotentialNamedFactory::getID4Name(name);
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
