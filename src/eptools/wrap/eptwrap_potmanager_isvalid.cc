/* -------------------------------------------------------------------
 * EPTWRAP_POTMANAGER_ISVALID
 *
 * Potential manager defined by POTIDS, NUMPOT, PARVEC, PARSHRD,
 * ANNOBJ. See 'PotManagerFactory' comments for full details.
 * Here, the validity of this representation is checked. If an error
 * is detected, an error string is returned, containing the
 * coordinate (block and position within block) where things are
 * wrong. Otherwise, the return string is empty.
 *
 * Use POSOFF == 1 if scripting language uses 1-based indexing (e.g.,
 * Matlab). POSOFF is added to coordinates stated in RETSTR.
 * For TAUIND, see 'PotManagerFactory::checkRepres' and
 * 'FactorizedEPRepresentation' comments.
 *
 * Input:
 * - POTIDS:  Potential manager representation [int32 array]
 * - NUMPOT:  " [int32 array]
 * - PARVEC:  " [double array]
 * - PARSHRD: " [int32 array]
 * - ANNOBJ:  " [void* array]
 * - POSOFF:  Optional. Offset added to positions stated in RETSTR.
 *            Def.: 0
 * - TAUIND:  Must be given iff the PM contains bivariate precision
 *            potentials
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
				W_IARRAY(parshrd),W_ARRAY(annobj,void*),
				int posoff,W_IARRAY(tauind),char** retstr,
				W_ERRORARGS)
{
  ArrayHandle<int> potidsA,numpotA,parshrdA,tauindA;
  ArrayHandle<double> parvecA;
  ArrayHandle<void*> annobjA;

  try {
    /* Read arguments */
    if (ain==5)
      posoff=0;
    else if (ain!=6 && ain!=7)
      W_RETERROR(2,"Wrong number of input arguments");
    if (aout!=1)
      W_RETERROR(2,"Need 1 return argument");
    W_CHKSIZE(numpot,npotids,"NUMPOT");
    W_CHKSIZE(annobj,npotids,"ANNOBJ");
    W_MASKARRAY(potids);
    W_MASKARRAY(numpot);
    W_MASKARRAY(parvec);
    W_MASKARRAY(parshrd);
    W_MASKARRAY(annobj);
    if (ain==7)
      W_MASKARRAY(tauind);
    *retstr = (char*) eptwrap_potmanager_isvalid_emptyStr;
    try {
      PotManagerFactory::checkRepres(potidsA,numpotA,parvecA,parshrdA,annobjA,
				     posoff,tauindA);
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
