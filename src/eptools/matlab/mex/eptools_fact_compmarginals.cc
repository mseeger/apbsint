/* -------------------------------------------------------------------
 * EPTOOLS_FACT_COMPMARGINALS
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * EP with factorized Gaussian backbone.
 * Compute marginals on variables from EP (message) parameters, overwrite
 * MARGPI, MARGBETA.
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array]
 * - RP_BETA:     " [double array]
 * - MARGPI:      Marginal pi parameters written here
 * - MARGBETA:    Marginal beta parameters written here
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_fact_compmarginals.h"

char errMsg[512];

/* Main function EPTOOLS_FACT_COMPMARGINALS */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n,m,argidx;
  int* rp_rowind,*rp_colind;
  double* rp_bvals,*rp_pi,*rp_beta,*margpi,*margbeta;
  int nrp_rowind,nrp_colind,nrp_bvals,nrp_pi,nrp_beta,nmargpi,nmargbeta;
  int errcode;
  char* errstr;

  errstr = errMsg;
  /* Read arguments */
  if (nrhs<9)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs>0)
    mexErrMsgTxt("Too many return arguments");
  argidx = -1;
  M_GETISCAL(n,"N");
  M_GETISCAL(m,"M");
  M_GETIARRAY(rp_rowind,"RP_ROWIND");
  M_GETIARRAY(rp_colind,"RP_COLIND");
  M_GETDARRAY(rp_bvals,"RP_BVALS");
  M_GETDARRAY(rp_pi,"RP_PI");
  M_GETDARRAY(rp_beta,"RP_BETA");
  M_GETDARRAY(margpi,"MARGPI");
  M_GETDARRAY(margbeta,"MARGBETA");
  eptwrap_fact_compmarginals(9,0,n,m,M_ARR(rp_rowind),M_ARR(rp_colind),
			     M_ARR(rp_bvals),M_ARR(rp_pi),M_ARR(rp_beta),
			     M_ARR(margpi),M_ARR(margbeta),&errcode,errstr);
  if (errcode!=0)
    mexErrMsgTxt(errstr);
}
