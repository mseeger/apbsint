/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_SINGLE
 *
 * Local EP update for a single potential t(s)
 *
 * Two ways to pass potential type and parameters. (1) Potential type PID
 * either ID or Name. PARS is parameter vector for potential.
 * (2) Potential manager by POTIDS, NUMPOT, PARVEC, PARSHRD, see
 * EPTOOLS_EPUPDATE_PARALLEL, then potential selected by PIND (0-floor).
 *
 * Input:
 * - (1):
 *   - PID:     Potential ID or Name (if string)
 *   - PARS:    Parameter vector
 * - Or (2):
 *   - POTIDS:  Potential manager representation [int32 array]
 *   - NUMPOT:  " [int32 array]
 *   - PARVEC:  " [double array]
 *   - PARSHRD: " [int32 array]
 *   - PIND:    Potential number
 * - CMU:     Cavity mean
 * - CRHO:    Cavity variance
 *
 * Return:
 * - RSTAT:   Return status (1: Success, 0: Failure)
 * - ALPHA:   Value for alpha
 * - NU:      Value for nu
 * - LOGZ:    Value for log Z (optional)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_epupdate_single.h"
#include "src/eptools/potentials/EPPotentialNamedFactory.h"
#include "src/eptools/potentials/PotentialManager.h"

void eptwrap_epupdate_single1(int ain,int aout,int pid,W_DARRAY(pars),
			      double cmu,double crho,int* rstat,double* alpha,
			      double* nu,double* logz,W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain!=4)
      W_RETERROR(2,"Need 4 input arguments");
    if (aout==3)
      logz=0;
    else if (aout==4) {
      if (logz==0)
	W_RETERROR(2,"LOGZ must be given");
    } else
      W_RETERROR(2,"Wrong number of return arguments");
    if (!EPPotentialFactory::isValidID(pid))
      W_RETERROR(1,"PID: Invalid potential ID");
    /* Create potential object, do update */
    Handle<EPScalarPotential> epPot;
    try {
      epPot.changeRep(EPPotentialFactory::create(pid,pars));
    } catch (...) {
      W_RETERROR(1,"Cannot create potential object");
    }
    *rstat = epPot->compMoments(cmu,crho,*alpha,*nu,logz);
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}

void eptwrap_epupdate_single2(int ain,int aout,char* pname,W_DARRAY(pars),
			      double cmu,double crho,int* rstat,double* alpha,
			      double* nu,double* logz,W_ERRORARGS)
{
  eptwrap_epupdate_single1(ain,aout,EPPotentialNamedFactory::getID4Name(pname),
			   W_ARR(pars),cmu,crho,rstat,alpha,nu,logz,W_ERRARGS);
}

void eptwrap_epupdate_single3(int ain,int aout,W_IARRAY(potids),
			      W_IARRAY(numpot),W_DARRAY(parvec),
			      W_IARRAY(parshrd),int pind,double cmu,
			      double crho,int* rstat,double* alpha,double* nu,
			      double* logz,W_ERRORARGS)
{
  Handle<PotentialManager> potMan;

  try {
    /* Read arguments */
    if (ain!=7)
      W_RETERROR(2,"Need 7 input arguments");
    if (aout==3)
      logz=0;
    else if (aout==4) {
      if (logz==0)
	W_RETERROR(2,"LOGZ must be given");
    } else
      W_RETERROR(2,"Wrong number of return arguments");
    createPotentialManager(W_ARR(potids),W_ARR(numpot),W_ARR(parvec),
			   W_ARR(parshrd),potMan,W_ERRARGS);
    if (pind<0 || pind>=potMan->size())
      W_RETERROR(1,"PIND out of range");
    /* Do update */
    *rstat = potMan->getPot(pind).compMoments(cmu,crho,*alpha,*nu,logz);
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
