/* -------------------------------------------------------------------
 * EPTOOLS_FACT_SEQUPDATES
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * EP with factorized Gaussian backbone. Run a number of updates on
 * potentials, one of the after (sequential updating).
 * Operates on ([I]: Input, [I/O]: Input/output, some vectors are
 * overwritten):
 * - Potential manager [I]. PM_POTIDS, PM_NUMPOT, PM_PARVEC, PM_PARSHRD,
 *   see EPTOOLS_EPUPDATE_PARALLEL comments.
 *   NOTE: Need full potential manager, even if updates only done on
 *   subset of potentials.
 * - Representation [I/O]. Structure and content of coupling factor
 *   B [I], EP (message) parameters [I/O]
 * - Marginals on variables [I/O]. MARGPI, MARGBETA
 * - Support data structure for selective damping mechanism [I/O].
 *   SD_NUMVALID, SD_TOPIND, SD_TOPVAL
 *
 * There are M potentials (factors), N variables. We run EP updates on
 * potentials in UPDJIND, one after the other. Messages and marginals
 * are factorized Gaussians, given by natural parameters (pi, beta).
 * An update modifies marginals MARGPI, MARGBETA and EP (message)
 * parameters RP_PI, RP_BETA. If DAMPFACT>0, the update is damped. There
 * may also be selective damping (see below).
 * Updates can fail for various reasons. They are either skipped, or
 * selective damping is applied:
 * - Cavity marginal undefined: If pi < PIMINTHRES/2
 * - New marginal undefined: If pi < PIMINTHRES/2
 * RSTAT is return status for each update. Codes defined in
 * 'FactorizedEPDriver':
 * - 0 [updSuccess]: Update successful
 * - 1 [updCavityInvalid]: Cavity marginal undefined. Update skipped
 * - 2 [updNumericalError]: Local EP update raises error. Update skipped
 * - 3 [updMarginalsInvalid]: New marginals undefined. Update skipped
 * - 4 [updCavCondSkipped]: Selective damping requires skipping
 * DELTA is relative change in moments for each non-skipped update, or 0
 * for skipped ones. The entry is the maximum relative difference for
 * means and stddevs (before and after update).
 *
 * Representation:
 * Consists of RP_ROWIND, RP_COLIND, RP_BVALS, RP_PI, RP_BETA. Details in
 * 'FactorizedEPRepresentation' comments. Internal representation
 * automatically compiled by Matlab code, see EPT.BFACT_INTREPRES.
 * RP_PI, RP_BETA (EP parameters) are I/O, their content is overwritten.
 *
 * Selective damping (optional):
 * SD_NUMVALID, SD_TOPIND, SD_TOPVAL, SD_SUBIND, SD_SUBEXCL. Details in
 * technical report, 'FactorizedEPDriver', 'FactEPMaximumPiValues'
 * comments. Idea is to ensure that for all EP parameters and marginals:
 * pi >= PIMINTHRES. This is a precondition (not checked). If the
 * condition is violated after an update, we apply the
 * minimum amount of extra damping (may be on top of DAMPFACT). In the
 * extreme case, the update is skipped (RSTAT: updCavCondSkipped).
 * Effective damping factor used for each update can be returned in
 * SD_DAMPFACT.
 * Requires 'FactEPMaximumPiValues' data structure, represented by SD_XXX
 * variables (I/O, content is overwritten).
 * SD_NUPD, SD_NREC return statistics about this datastructure (number of
 * update calls and block recomputations).
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - UPDJIND:     Update on these potentials, in order
 * - PM_POTIDS:   Potential manager [int32 array]
 * - PM_NUMPOT:   " [int32 array]
 * - PM_PARVEC:   " [double array]
 * - PM_PARSHRD:  " [int32 array]
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array; I/O]
 * - RP_BETA:     " [double array; I/O]
 * - MARGPI:      Variable marginals [I/O]
 * - MARGBETA:    " [I/O]
 * - PIMINTHRES:  See above. Positive
 * - DAMPFACT:    Damping factor, in [0,1). Optional, def. is 0
 * - SD_NUMVALID: Selective damping. Optional [int32 array; I/O]
 * - SD_TOPIND:   " [int32 array; I/O]
 * - SD_TOPVAL:   " [double array; I/O]
 * - SD_SUBIND    " [int32 array]
 * - SD_SUBEXCL   ". Def.: false
 *
 * Return:
 * - RSTAT:       Return stati for each update. Optional
 * - DELTA:       See above. Optional
 * - SD_DAMPFACT: See above. Optional, only if selective damping
 * - SD_NUPD:     "
 * - SD_NREC:     "
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_fact_sequpdates.h"

char errMsg[512];

/* Main function EPTOOLS_FACT_SEQUPDATES */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n,m,argidx;
  double piminthres,dampfact=0.0,sd_subexcl=0;
  int* updjind,*pm_potids,*pm_numpot,*pm_parshrd,*rp_rowind,*rp_colind;
  double* pm_parvec,*rp_bvals,*rp_pi,*rp_beta,*margpi,*margbeta;
  int nupdjind,npm_potids,npm_numpot,npm_parshrd,nrp_rowind,nrp_colind,
    npm_parvec,nrp_bvals,nrp_pi,nrp_beta,nmargpi,nmargbeta;
  int* sd_numvalid=0,*sd_topind=0,*sd_subind=0;
  double* sd_topval=0;
  int nsd_numvalid=0,nsd_topind=0,nsd_subind=0,nsd_topval=0;
  int* rstat=0;
  double* delta=0,*sd_dampfact=0;
  int nrstat,ndelta,nsd_dampfact;
  int sd_nupd,sd_nrec;
  int errcode;
  char* errstr;

  errstr = errMsg;
  /* Read arguments */
  if (nrhs<15)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs>5)
    mexErrMsgTxt("Too many return arguments");
  argidx = -1; /* ++argidx in macros */
  M_GETISCAL(n,"N");
  M_GETISCAL(m,"M");
  M_GETIARRAY(updjind,"UPDJIND");
  M_GETIARRAY(pm_potids,"PM_POTIDS");
  M_GETIARRAY(pm_numpot,"PM_NUMPOT");
  M_GETDARRAY(pm_parvec,"PM_PARVEC");
  M_GETIARRAY(pm_parshrd,"PM_PARSHRD");
  M_GETIARRAY(rp_rowind,"RP_ROWIND");
  M_GETIARRAY(rp_colind,"RP_COLIND");
  M_GETDARRAY(rp_bvals,"RP_BVALS");
  M_GETDARRAY(rp_pi,"RP_PI");
  M_GETDARRAY(rp_beta,"RP_BETA");
  M_GETDARRAY(margpi,"MARGPI");
  M_GETDARRAY(margbeta,"MARGBETA");
  M_GETDSCAL(piminthres,"PIMINTHRES");
  if (nrhs>15) {
    M_GETDSCAL(dampfact,"DAMPFACT");
    if (nrhs>16) {
      /* Selective damping */
      if (nrhs<19)
	mexErrMsgTxt("Need all SD_XXX or none");
      M_GETIARRAY(sd_numvalid,"SD_NUMVALID");
      M_GETIARRAY(sd_topind,"SD_TOPIND");
      M_GETDARRAY(sd_topval,"SD_TOPVAL");
      if (nrhs>19) {
	M_GETIARRAY(sd_subind,"SD_SUBIND");
	if (nrhs>20)
	  M_GETISCAL(sd_subexcl,"SD_SUBEXCL");
      }
    }
  }
  /* Create return arguments */
  argidx = -1; /* ++argidx in macro */
  if (nlhs>0) {
    nrstat=ndelta=nsd_dampfact=nupdjind;
    M_MAKEIARRAY(rstat);
    if (nlhs>1) {
      M_MAKEDARRAY(delta);
      if (nlhs>2)
	M_MAKEDARRAY(sd_dampfact);
    }
  }
  /* Call C++ wrapper, deal with error */
  /*sprintf(errstr,"nupdjind=%d. Call wrapper",nupdjind);
    printMsgStdout(errstr);*/
  eptwrap_fact_sequpdates(std::min(nrhs,21),nlhs,n,m,M_ARR(updjind),
			  M_ARR(pm_potids),M_ARR(pm_numpot),M_ARR(pm_parvec),
			  M_ARR(pm_parshrd),M_ARR(rp_rowind),M_ARR(rp_colind),
			  M_ARR(rp_bvals),M_ARR(rp_pi),M_ARR(rp_beta),
			  M_ARR(margpi),M_ARR(margbeta),piminthres,dampfact,
			  M_ARR(sd_numvalid),M_ARR(sd_topind),
			  M_ARR(sd_topval),M_ARR(sd_subind),sd_subexcl,
			  M_ARR(rstat),M_ARR(delta),M_ARR(sd_dampfact),
			  &sd_nupd,&sd_nrec,&errcode,errstr);
  /*printMsgStdout("Exit from wrapper");*/
  if (errcode!=0)
    mexErrMsgTxt(errstr);
  if (nlhs>3) {
    argidx = 2; /* ++argidx in macro */
    M_SETISCAL(sd_nupd);
    if (nlhs>4)
      M_SETISCAL(sd_nrec);
  }
}
