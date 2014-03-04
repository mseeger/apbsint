/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_PARALLEL_BVPREC
 *
 * Local EP updates (in parallel) for all potentials
 * t_j(s_j,tau_{k(j)}) of potential manager.
 *
 * Same as EPTWRAP_EPUPDATE_PARALLEL, but for bivariate potentials
 * with precision parameter. Additional inputs are CA, CC (Gamma
 * parameters of cavity marginals over tau variables), additional
 * returns are HATA, HATC (Gamma parameters of updated marginals).
 * All potentials must be in the argument group 'atypeBivarPrec'.
 *
 * Input:
 * - POTIDS:  Potential manager representation [int32 array]
 * - NUMPOT:  " [int array]
 * - PARVEC:  " [double array]
 * - PARSHRD: " [int32 array]
 * - ANNOBJ:  " [void* array]
 * - CMU:     Vector cavity means
 * - CRHO:    Vector cavity variances
 * - CA:      Vector cavity a parameters
 * - CC:      Vector cavity c parameters
 * - UPDIND:  S.a. Optional [int32 array]
 *
 * Return:
 * - RSTAT:   Vector of return stati (1: Success, 0: Failure)
 *            [int32 array]
 * - ALPHA:   Vector of alpha values
 * - NU:      Vector of nu values
 * - HATA:    Vector of a_hat values
 * - HATC:    Vector of c_hat values
 * - LOGZ:    Vector of log Z values (optional)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_epupdate_parallel_bvprec.h"
#include "src/eptools/potentials/PotentialManager.h"

void eptwrap_epupdate_parallel_bvprec(int ain,int aout,W_IARRAY(potids),
				      W_IARRAY(numpot),W_DARRAY(parvec),
				      W_IARRAY(parshrd),W_ARRAY(annobj,void*),
				      W_DARRAY(cmu),W_DARRAY(crho),
				      W_DARRAY(ca),W_DARRAY(cc),
				      W_IARRAY(updind),W_IARRAY(rstat),
				      W_DARRAY(alpha),W_DARRAY(nu),
				      W_DARRAY(hata),W_DARRAY(hatc),
				      W_DARRAY(logz),W_ERRORARGS)
{
  int i,j,totsz;
  double temp;
  Handle<PotentialManager> potMan;
  double inp[4],ret[4];

  try {
    /* Read arguments */
    if (ain<9 || ain>10)
      W_RETERROR(2,"Wrong number of input arguments");
    if (aout<5 || aout>6)
      W_RETERROR(2,"Wrong number of return arguments");
    /* Create potential manager */
    createPotentialManager(W_ARR(potids),W_ARR(numpot),W_ARR(parvec),
			   W_ARR(parshrd),W_ARR(annobj),potMan,W_ERRARGS);
    //cout << "Wrap: Done createPotentialManager" << endl; // DEBUG!
    totsz=ncmu;
    W_CHKSIZE(crho,totsz,"CRHO");
    W_CHKSIZE(ca,totsz,"CA");
    W_CHKSIZE(cc,totsz,"CC");
    if (potMan->numArgumentGroup(EPScalarPotential::atypeBivarPrec)!=totsz)
      W_RETERROR(1,"All potentials must be in group 'atypeBivarPrec'");
    if (ain>9) {
      if (updind==0)
	W_RETERROR(2,"UPDIND missing");
      W_CHKSIZE(updind,totsz,"UPDIND");
      Interval<int> ivM(0,potMan->size()-1,IntVal::ivClosed,IntVal::ivClosed);
      if (ivM.check(updind,nupdind)!=0)
	W_RETERROR(1,"UPDIND: Entries out of range");
    } else {
      updind=0;
      if (potMan->size()!=totsz)
	W_RETERROR(1,"CMU, potential manager: Different sizes");
    }
    /* Return arguments */
    W_CHKSIZE(rstat,totsz,"RSTAT");
    W_CHKSIZE(alpha,totsz,"ALPHA");
    W_CHKSIZE(nu,totsz,"NU");
    W_CHKSIZE(hata,totsz,"HATA");
    W_CHKSIZE(hatc,totsz,"HATC");
    if (aout>5)
      W_CHKSIZE(logz,totsz,"LOGZ");
    else
      logz=0;

    /* Main loop over all potentials */
    for (i=0; i<totsz; i++) {
      j=(updind==0)?i:updind[i];
      //cout << "i=" << i << endl;
      inp[0]=cmu[i]; inp[1]=crho[i]; inp[2]=ca[i]; inp[3]=cc[i];
      rstat[i] = potMan->getPot(j).compMoments(inp,ret,&temp);
      alpha[i]=ret[0]; nu[i]=ret[1]; hata[i]=ret[2]; hatc[i]=ret[3];
      if (rstat[i] && aout>5)
	logz[i]=temp;
    }
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
