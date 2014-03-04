/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_SINGLE_BVPREC
 *
 * Local EP update for a single potential t(s,tau). Same as
 * EPTWRAP_EPUPDATE_SINGLE, but for bivariate potentials with
 * precision parameter.
 * Potential must be in group 'atypeBivarPrec'.
 *
 * Input:
 * - (1):
 *   - PID:     Potential ID or Name (if string)
 *   - PARS:    Parameter vector
 *   - ANNOBJ:  Annotation (0: None)
 * - Or (2):
 *   - POTIDS:  Potential manager representation [int32 array]
 *   - NUMPOT:  " [int32 array]
 *   - PARVEC:  " [double array]
 *   - PARSHRD: " [int32 array]
 *   - ANNOBJ:  " [void* array]
 *   - PIND:    Potential number
 * - CMU:     Cavity mean
 * - CRHO:    Cavity variance
 * - CA:      Cavity a param.
 * - CC:      Cavity c param.
 *
 * Return:
 * - RSTAT:   Return status (1: Success, 0: Failure)
 * - ALPHA:   Value for alpha
 * - NU:      Value for nu
 * - HATA:    Value for a_hat
 * - HATC:    Value for c_hat
 * - LOGZ:    Value for log Z (optional)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_epupdate_single_bvprec.h"
#include "src/eptools/potentials/EPPotentialNamedFactory.h"
#include "src/eptools/potentials/PotentialManager.h"

void eptwrap_epupdate_single_bvprec1(int ain,int aout,int pid,W_DARRAY(pars),
				     void* annobj,double cmu,double crho,
				     double ca,double cc,int* rstat,
				     double* alpha,double* nu,double* hata,
				     double* hatc,double* logz,W_ERRORARGS)
{
  Handle<EPScalarPotential> epPot;
  double inp[4],ret[4];

  try {
    /* Read arguments */
    if (ain!=7)
      W_RETERROR(2,"Need 7 input arguments");
    if (aout==5)
      logz=0;
    else if (aout==6) {
      if (logz==0)
	W_RETERROR(2,"LOGZ must be given");
    } else
      W_RETERROR(2,"Wrong number of return arguments");
    if (!EPPotentialFactory::isValidID(pid))
      W_RETERROR(1,"PID: Invalid potential ID");
    /* Create potential object, do update */
    try {
      epPot.changeRep(EPPotentialFactory::create(pid,pars,annobj));
    } catch (...) {
      W_RETERROR(1,"Cannot create potential object");
    }
    if (epPot->getArgumentGroup()!=EPScalarPotential::atypeBivarPrec)
      W_RETERROR(1,"Potential must be in group 'atypeBivarPrec'");
    inp[0]=cmu; inp[1]=crho; inp[2]=ca; inp[3]=cc;
    *rstat = epPot->compMoments(inp,ret,logz);
    *alpha=ret[0]; *nu=ret[1]; *hata=ret[2]; *hatc=ret[3];
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}

void eptwrap_epupdate_single_bvprec2(int ain,int aout,char* pname,
				     W_DARRAY(pars),void* annobj,double cmu,
				     double crho,double ca,double cc,int* rstat,
				     double* alpha,double* nu,double* hata,
				     double* hatc,double* logz,W_ERRORARGS)
{
  eptwrap_epupdate_single_bvprec1(ain,aout,
				  EPPotentialNamedFactory::getID4Name(pname),
				  W_ARR(pars),annobj,cmu,crho,ca,cc,rstat,
				  alpha,nu,hata,hatc,logz,W_ERRARGS);
}

void eptwrap_epupdate_single_bvprec3(int ain,int aout,W_IARRAY(potids),
				     W_IARRAY(numpot),W_DARRAY(parvec),
				     W_IARRAY(parshrd),W_ARRAY(annobj,void*),
				     int pind,double cmu,double crho,double ca,
				     double cc,int* rstat,double* alpha,
				     double* nu,double* hata,double* hatc,
				     double* logz,W_ERRORARGS)
{
  Handle<PotentialManager> potMan;
  double inp[4],ret[4];

  try {
    /* Read arguments */
    if (ain!=10)
      W_RETERROR(2,"Need 10 input arguments");
    if (aout==5)
      logz=0;
    else if (aout==6) {
      if (logz==0)
	W_RETERROR(2,"LOGZ must be given");
    } else
      W_RETERROR(2,"Wrong number of return arguments");
    createPotentialManager(W_ARR(potids),W_ARR(numpot),W_ARR(parvec),
			   W_ARR(parshrd),W_ARR(annobj),potMan,W_ERRARGS);
    if (pind<0 || pind>=potMan->size())
      W_RETERROR(1,"PIND out of range");
    if (potMan->getPot(pind).getArgumentGroup()!=
	EPScalarPotential::atypeBivarPrec)
      W_RETERROR(1,"Potential must be in group 'atypeBivarPrec'");
    /* Do update */
    inp[0]=cmu; inp[1]=crho; inp[2]=ca; inp[3]=cc;
    *rstat = potMan->getPot(pind).compMoments(inp,ret,logz);
    *alpha=ret[0]; *nu=ret[1]; *hata=ret[2]; *hatc=ret[3];
     W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
