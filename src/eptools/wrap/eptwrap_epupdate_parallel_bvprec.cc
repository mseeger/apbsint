/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_PARALLEL_BVPREC
 *
 * Generalization of EPTWRAP_EPUPDATE_PARALLEL to potential managers
 * which include bivariate potentials t_j(s_j,tau_k(j)), where tau_k
 * are precision parameters (argument group 'atypeBivarPrec').
 *
 * The PM may contain standard potentials t_j(s_j) as well (group
 * 'atypeUnivariate'). The precision potentials must come last.
 * Additional inputs are CA, CC (Gamma parameters of cavity marginals
 * over tau variables), additional returns are HATA, HATC (Gamma
 * parameters of updated marginals). These are flat vectors (size:
 * number of precision potentials).
 *
 * UPDIND (optional) plays the same role as in
 * EPTWRAP_EPUPDATE_PARALLEL. It must only select standard
 * potentials (group 'atypeUnivariate'): we always update on all
 * precision potentials. In fact, if there are m_st standard, then
 * m_prec precision potentials (so m = m_st + m_prec), we use the
 * effective subselection index
 *   [UPDIND; m_st:(m-1)] (size sz(UPDIND) + m_prec).
 * If UPDIND is not given, substitute 0:(m_st-1).
 * The variables CMU, CRHO, RSTAT, ALPHA, NU, LOGZ follow this
 * ordering (if j is l-th entry in effective index, then
 * alpha_j == ALPHA(l), etc.).
 * The variables CA, CC, HATA, HATC are of size m_prec (so that
 * hat{a}_j == HATA(j - m_st), etc.).
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
  int i,j,totsz,numBVPrec,numSt,thresSt;
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
    if (ain<=9 && totsz!=potMan->size())
      W_RETERROR(1,"CMU: Wrong size");
    numBVPrec=potMan->numArgumentGroup(EPScalarPotential::atypeBivarPrec);
    if (numBVPrec==0)
      W_RETERROR(1,"Potential manager must contain precision parameter potentials");
    numSt=potMan->numArgumentGroup(EPScalarPotential::atypeUnivariate);
    if (numSt+numBVPrec!=potMan->size())
      W_RETERROR(1,"Potentials of unsupported argument group are present");
    W_CHKSIZE(crho,totsz,"CRHO");
    W_CHKSIZE(ca,numBVPrec,"CA");
    W_CHKSIZE(cc,numBVPrec,"CC");
    if (ain>9) {
      if (updind==0)
	W_RETERROR(2,"UPDIND missing");
      W_CHKSIZE(updind,totsz-numBVPrec,"UPDIND");
      Interval<int> ivM(0,potMan->size()-1,IntVal::ivClosed,IntVal::ivClosed);
      if (ivM.check(updind,numSt)!=0)
	W_RETERROR(1,"UPDIND: Entries out of range");
    } else
      updind=0;
    /* Return arguments */
    W_CHKSIZE(rstat,totsz,"RSTAT");
    W_CHKSIZE(alpha,totsz,"ALPHA");
    W_CHKSIZE(nu,totsz,"NU");
    W_CHKSIZE(hata,numBVPrec,"HATA");
    W_CHKSIZE(hatc,numBVPrec,"HATC");
    if (aout>5)
      W_CHKSIZE(logz,totsz,"LOGZ");
    else
      logz=0;

    /* Main loop over all potentials */
    thresSt=(updind==0)?numSt:nupdind; // Number standard pots to upd. on
    for (i=0; i<totsz; i++) {
      inp[0]=cmu[i]; inp[1]=crho[i];
      if (i<thresSt)
	j=(updind==0)?i:updind[i];
      else {
	// Precision potential: 'updind' does not apply
	j=i-thresSt;
	inp[2]=ca[j]; inp[3]=cc[j];
	j+=numSt;
      }
      rstat[i] = potMan->getPot(j).compMoments(inp,ret,&temp);
      alpha[i]=ret[0]; nu[i]=ret[1];
      if (i>=thresSt) {
	hata[j-numSt]=ret[2]; hatc[j-numSt]=ret[3];
      }
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
