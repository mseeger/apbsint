/* -------------------------------------------------------------------
 * EPTOOLS_POTMANAGER_ISVALID
 *
 * Potential manager defined by POTIDS, NUMPOT, PARVEC, PARSHRD. See
 * 'PotManagerFactory' comments for full details.
 * Here, the validity of this representation is checked. If an error
 * is detected, an error string is returned, containing the
 * coordinate (block and position within block) where things are
 * wrong. Otherwise, the return string is empty.
 *
 * Use POSOFF == 1 if scripting language uses 1-based indexing (e.g.,
 * Matlab). POSOFF is added to coordinates stated in RETSTR.
 *
 * Input:
 * - POTIDS:  Potential manager representation [int32 array]
 * - NUMPOT:  " [int32 array]
 * - PARVEC:  " [double array]
 * - PARSHRD: " [int32 array]
 * - POSOFF:  Optional. Offset added to positions stated in RETSTR.
 *            Def.: 0
 *
 * Return:
 * - RETSTR:  S.a. Empty if no error
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_potmanager_isvalid.h"

char errMsg[512];

/* Main function EPTOOLS_POTMANAGER_ISVALID */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int argidx,posoff=0;
  int* potids,*numpot,*parshrd;
  double* parvec;
  int npotids,nnumpot,nparshrd,nparvec;
  char* retstr;
  int errcode;
  char* errstr;

  errstr=errMsg;
  /* Read arguments */
  if (nrhs<4)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=1)
    mexErrMsgTxt("Need one return argument");
  argidx = -1;
  M_GETIARRAY(potids,"POTIDS");
  M_GETIARRAY(numpot,"NUMPOT");
  M_GETDARRAY(parvec,"PARVEC");
  M_GETIARRAY(parshrd,"PARSHRD");
  if (nrhs>4)
    M_GETISCAL(posoff,"POSOFF");
  /* Call C++ wrapper, deal with error */
  eptwrap_potmanager_isvalid(4,1,M_ARR(potids),M_ARR(numpot),M_ARR(parvec),
			     M_ARR(parshrd),posoff,&retstr,&errcode,errstr);
  if (errcode!=0)
    mexErrMsgTxt(errstr);
  /* Return argument */
  plhs[0] = mxCreateString(retstr);
}
