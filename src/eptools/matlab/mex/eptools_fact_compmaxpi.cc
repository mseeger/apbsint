/* -------------------------------------------------------------------
 * EPTOOLS_FACT_COMPMAXPI
 *
 * EP with factorized Gaussian backbone.
 * Computes top-K values in 'FactEPMaximumPiValues' data structure from
 * scratch ('FactEPMaximumPiValues::recompute').
 * This data structure is used for selective damping, see
 * EPTOOLS_FACT_SEQUPDATES.
 *
 * If SD_SUBIND is given, it is a subset of 0:(M-1), sorted in
 * ascending order. See 'FactEPMaximumPiValues', fields 'subInd' and
 * 'subExcl'
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array]
 * - RP_BETA:     " [double array]
 * - SD_K:        Value K (must be >1)
 * - SD_SUBIND    See above. Optional [int32 array]
 * - SD_SUBEXCL   ". Def.: false
 *
 * Return:
 * - SD_NUMVALID: Max pi data structure [int32 array]
 * - SD_TOPIND:   " [int32 array]
 * - SD_TOPVAL:   " [double array]
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_fact_compmaxpi.h"

char errMsg[512];

/* Main function EPTOOLS_FACT_COMPMAXPI */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int argidx,n,m,sd_k,sd_subexcl=0;
  int* rp_rowind,*rp_colind,*sd_subind=0;
  double* rp_bvals,*rp_pi,*rp_beta;
  int nrp_rowind,nrp_colind,nsd_subind,nrp_bvals,nrp_pi,nrp_beta;
  int* sd_numvalid,*sd_topind;
  double* sd_topval;
  int nsd_numvalid,nsd_topind,nsd_topval;
  int errcode;
  char* errstr;

  errstr = errMsg;
  /* Read arguments */
  if (nrhs<8)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=3)
    mexErrMsgTxt("Need 3 return arguments");
  argidx = -1;
  M_GETISCAL(n,"N");
  M_GETISCAL(m,"M");
  M_GETIARRAY(rp_rowind,"RP_ROWIND");
  M_GETIARRAY(rp_colind,"RP_COLIND");
  M_GETDARRAY(rp_bvals,"RP_BVALS");
  M_GETDARRAY(rp_pi,"RP_PI");
  M_GETDARRAY(rp_beta,"RP_BETA");
  M_GETISCAL(sd_k,"SD_K");
  if (nrhs>8) {
    M_GETIARRAY(sd_subind,"SD_SUBIND");
    if (nrhs>9)
      M_GETISCAL(sd_subexcl,"SD_SUBEXCL");
  }
  /* Return arguments */
  argidx = -1;
  nsd_numvalid = n;
  M_MAKEIARRAY(sd_numvalid);
  nsd_topind = nsd_topval = n*(sd_k+1);
  M_MAKEIARRAY(sd_topind);
  M_MAKEDARRAY(sd_topval);
  eptwrap_fact_compmaxpi(std::min(nrhs,10),nlhs,n,m,M_ARR(rp_rowind),
			 M_ARR(rp_colind),M_ARR(rp_bvals),M_ARR(rp_pi),
			 M_ARR(rp_beta),sd_k,M_ARR(sd_subind),sd_subexcl,
			 M_ARR(sd_numvalid),M_ARR(sd_topind),M_ARR(sd_topval),
			 &errcode,errstr);
  if (errcode!=0)
    mexErrMsgTxt(errstr);
}
