/* -------------------------------------------------------------------
 * EPTOOLS_EPUPDATE_SINGLE
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
#include "matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_epupdate_single.h"

/* Main function EPTOOLS_EPUPDATE_SINGLE */

char errMsg[512];

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int argidx,rstat;
  double temp,cmu,crho,nu,alpha,logz;
  bool usePMan=(nrhs>=7);
  int errcode;
  char* errstr;

  errstr = errMsg;
  /* Read arguments */
  if (!usePMan && nrhs!=4)
    mexErrMsgTxt("Wrong number of input arguments");
  if (nlhs<3 || nlhs>4)
    mexErrMsgTxt("Wrong number of return arguments");
  if (!usePMan) {
    /* (1) Potential by PID, PARS */
    double* pars;
    int npars;
    argidx = 0; /* Skip 1st argument for now */
    M_GETDARRAY(pars,"PARS");
    M_GETDSCAL(cmu,"CMU");
    M_GETDSCAL(crho,"CRHO");
    if (mxIsChar(prhs[0])) {
      // PID is string: Name
      char* pidStr = (char*) getString(prhs[0],"PID");
      eptwrap_epupdate_single2(4,nlhs,pidStr,M_ARR(pars),cmu,crho,&rstat,
			       &alpha,&nu,&logz,&errcode,errstr);
      mxFree((void*) pidStr);
    } else {
      // PID is number: ID
      int pid = getScalInt(prhs[0],"PID");
      eptwrap_epupdate_single1(4,nlhs,pid,M_ARR(pars),cmu,crho,&rstat,
			       &alpha,&nu,&logz,&errcode,errstr);
    }
  } else {
    /* (2) Potential manager to identify potential type and parameters */
    int pind;
    int* potids,*numpot,*parshrd;
    double* parvec;
    int npotids,nnumpot,nparshrd,nupdind,nparvec;
    argidx = -1;
    M_GETIARRAY(potids,"POTIDS");
    M_GETIARRAY(numpot,"NUMPOT");
    M_GETDARRAY(parvec,"PARVEC");
    M_GETIARRAY(parshrd,"PARSHRD");
    M_GETISCAL(pind,"PIND");
    M_GETDSCAL(cmu,"CMU");
    M_GETDSCAL(crho,"CRHO");
    eptwrap_epupdate_single3(7,nlhs,M_ARR(potids),M_ARR(numpot),
			     M_ARR(parvec),M_ARR(parshrd),pind,cmu,crho,
			     &rstat,&alpha,&nu,&logz,&errcode,errstr);
  }
  if (errcode!=0)
    mexErrMsgTxt(errstr);
  /* Write return arguments */
  argidx = -1;
  M_SETISCAL(rstat);
  M_SETDSCAL(alpha);
  M_SETDSCAL(nu);
  if (nlhs>3)
    M_SETDSCAL(logz);
}
