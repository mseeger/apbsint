/* -------------------------------------------------------------------
 * EPTOOLS_GETPOTID
 *
 * Potential names <--> IDs maintained in 'EPPotentialNamedFactory'
 *
 * Input:
 * - NAME:    Potential name
 *
 * Return:
 * - PID:     Potential ID, or -1 if NAME is not potential name
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_getpotid.h"

char errMsg[512];

/* Main function EPTOOLS_GETPOTID */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int argidx,pid;
  char* name;
  int errcode;
  char* errstr;

  errstr = errMsg;
  /* Read arguments */
  if (nrhs<1)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=1)
    mexErrMsgTxt("Need 1 return argument");
  name = (char*) getString(prhs[0],"NAME");
  eptwrap_getpotid(1,1,name,&pid,&errcode,errstr);
  mxFree((void*) name);
  if (errcode!=0)
    mexErrMsgTxt(errstr);
  argidx = -1;
  M_SETISCAL(pid);
}
