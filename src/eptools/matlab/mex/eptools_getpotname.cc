/* -------------------------------------------------------------------
 * EPTOOLS_GETPOTNAME
 *
 * Potential names <--> IDs maintained in 'EPPotentialNamedFactory'
 *
 * Input:
 * - PID:     Potential ID
 *
 * Return:
 * - NAME:    Potential name, or "" if PID is not valid ID
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_getpotname.h"

char errMsg[512];

/* Main function EPTOOLS_GETPOTNAME */

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
  argidx = -1;
  M_GETISCAL(pid,"PID");
  eptwrap_getpotname(1,1,pid,&name,&errcode,errstr);
  if (errcode!=0)
    mexErrMsgTxt(errstr);
  plhs[0] = mxCreateString(name);
}
