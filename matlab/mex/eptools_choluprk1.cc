/* -------------------------------------------------------------------
 * EPTOOLS_CHOLUPRK1
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * If
 *   A = L*L', A_ = A + v*v' = L_ L_', A, L n-by-n,
 * where L is lower triangular, the method computes L_ from L. This
 * is called Cholesky rank one update. L or L' (upper triangular) can
 * be passed, only the relevant triagle is accessed. L (or L') is
 * passed in L, v in VEC.
 * NOTE: For the present implementation, the method is more efficient
 * when a lower triangular matrix is used.
 *
 * Dragging along:
 * If Z (r-by-n) is given, so must be the r-vector y. In this case,
 * we overwrite Z by Z_, where
 *   Z_ L_' = Z L' + y v'.
 *
 * Working array:
 * The method uses n Givens rotations, param. by angles c_k, s_k. The
 * change L -> L' is specified by these 2*n numbers. They are written
 * into CVEC, SVEC.
 * Requires a working vector of size >= max(n,r), passed in WORKV.
 * NOTE: The same vector can be passed for VEC and WORKV, in which case
 * VEC is overwritten in an undefined way. If r > n, VEC can be of size
 * r, containing v in the first n components.
 *
 * Input:
 * - L:     Factor L (or L'), overwritten by L_ (or L_'). Must be
 *          lower (upper) triangular, str. code UPLO
 * - VEC:   Vector v. Can have size >n, only first n elem. are used
 * - CVEC:  Vector [n]. c_k ret. here
 * - SVEC:  Vector [n]. s_k ret. here
 * - WORKV: Working vector of size max(n,r). Can be same as VEC
 * - Z:     Dragging along matrix [r-by-n]. Optional
 * - Y:     Dragging along vector [r]. Iff Z is given
 *
 * Return:
 * - STAT:  0 (OK), 1 (Numerical error)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "matlab/mex/mex_helper.h"
#include "src/eptools/wrap/eptwrap_choluprk1.h"
#include "blas.h"

char errMsg[512];

/* Main function EPTOOLS_CHOLUPRK1 */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int argidx,stat;
  fst_matrix lmat,zmat;
  double* vec,*cvec,*svec,*workv,*yvec=0;
  int nvec,ncvec,nsvec,nworkv,nyvec=0;
  int errcode;
  char* errstr;

  errstr = errMsg;
  /* Read arguments */
  if (nrhs<5)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs>1)
    mexErrMsgTxt("Too many return arguments");
  argidx = -1;
  parseBLASMatrix(prhs[++argidx],"L",&lmat,-1,-1);
  M_GETDARRAY(vec,"VEC");
  M_GETDARRAY(cvec,"CVEC");
  M_GETDARRAY(svec,"SVEC");
  M_GETDARRAY(workv,"WORKV");
  if (nrhs>5) {
    if (nrhs<7) mexErrMsgTxt("Need both Z, Y or none");
    parseBLASMatrix(prhs[++argidx],"Z",&zmat,-1,lmat.n);
    M_GETDARRAY(yvec,"Y");
  }
  /* Call C++ wrapper, deal with error */
  eptwrap_choluprk1((nrhs<=7)?nrhs:7,1,&lmat,M_ARR(vec),M_ARR(cvec),
		    M_ARR(svec),M_ARR(workv),&zmat,M_ARR(yvec),&stat,
		    &BLASFUNC(dcopy),&BLASFUNC(drotg),&BLASFUNC(drot),
		    &errcode,errstr);
  if (errcode!=0)
    mexErrMsgTxt(errstr);
  if (nlhs==1) {
    argidx = -1;
    M_SETISCAL(stat);
  }
}
