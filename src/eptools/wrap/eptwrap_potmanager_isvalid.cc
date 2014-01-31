/* -------------------------------------------------------------------
 * EPTWRAP_POTMANAGER_ISVALID
 *
 * Potential manager defined by POTIDS, NUMPOT, PARVEC, PARSHRD. See
 * 'PotManagerFactory' comments for full details.
 * Here, the validity of this representation is checked. If an error
 * is detected, an error string is returned, containing the
 * coordinate (block and position within block) where things are
 * wrong. Otherwise, the return string is empty.
 *
 * Use POSOFF == 1 if scripting language uses 1-based indexing (e.g.,
 * Matlab). POSOFF is added to coordinates stated in RETSTR.
 *
 * Input:
 * - POTIDS:  Potential manager representation [int32 array]
 * - NUMPOT:  " [int32 array]
 * - PARVEC:  " [double array]
 * - PARSHRD: " [int32 array]
 * - POSOFF:  Optional. Offset added to positions stated in RETSTR.
 *            Def.: 0
 *
 * Return:
 * - RETSTR:  S.a. Empty if no error
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_potmanager_isvalid.h"
#include "src/eptools/potentials/PotManagerFactory.h"

const char* eptwrap_potmanager_isvalid_emptyStr = "";

/*
 * We use 'errstr' for '*retstr'
 */
void eptwrap_potmanager_isvalid(int ain,int aout,W_IARRAY(potids),
				W_IARRAY(numpot),W_DARRAY(parvec),
				W_IARRAY(parshrd),int posoff,char** retstr,
				W_ERRORARGS)
{
  ArrayHandle<int> potidsA,numpotA,parshrdA;
  ArrayHandle<double> parvecA;

  try {
    /* Read arguments */
    if (ain==4)
      posoff=0;
    else if (ain!=5)
      W_RETERROR(2,"Wrong number of input arguments");
    if (aout!=1)
      W_RETERROR(2,"Need 1 return argument");
    W_CHKSIZE(numpot,npotids,"NUMPOT");
    W_MASKARRAY(potids);
    W_MASKARRAY(numpot);
    W_MASKARRAY(parvec);
    W_MASKARRAY(parshrd);
    *retstr = (char*) eptwrap_potmanager_isvalid_emptyStr;
    try {
      PotManagerFactory::checkRepres(potidsA,numpotA,parvecA,parshrdA,posoff);
    } catch (InvalidParameterException ex) {
      /* Use 'errstr' as buffer for RETSTR */
      strcpy(errstr,ex.msg());
      *retstr = errstr;
    }
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
