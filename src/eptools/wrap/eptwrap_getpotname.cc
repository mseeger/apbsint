/* -------------------------------------------------------------------
 * EPTWRAP_GETPOTNAME
 *
 * Potential names <--> IDs maintained in 'EPPotentialNamedFactory'
 *
 * Input:
 * - PID:     Potential ID
 *
 * Return:
 * - NAME:    Potential name, or "" if PID is not valid ID
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_getpotname.h"
#include "src/eptools/potentials/EPPotentialNamedFactory.h"

const char* eptwrap_getpotname_emptyStr = "";

void eptwrap_getpotname(int ain,int aout,int pid,char** name,W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain!=1)
      W_RETERROR(2,"Need 1 input argument");
    if (aout!=1)
      W_RETERROR(2,"Need 1 return argument");
    if (!EPPotentialFactory::isValidID(pid))
      *name = (char*) eptwrap_getpotname_emptyStr;
    else
      *name = (char*) EPPotentialNamedFactory::getName4ID(pid).c_str();
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
