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

#include "src/eptools/wrap/eptools_helper_basic.h"
#include "src/eptools/wrap/eptwrap_choluprk1.h"

/*
 * The method is adapted from LINPACK dchud. We did the following
 * modifications:
 * - Using BLAS drot in order to avoid any explicit O(n^2) loops
 * - Keeping diag(L_) positive, by flipping angles c_k, s_k whenever
 *   a negative element pops up
 * See the TR
 *   M. Seeger
 *   Low Rank Updates for the Cholesky Decomposition
 */

void eptwrap_choluprk1(int ain,int aout,fst_matrix* lmat,W_DARRAY(vvec),
		       W_DARRAY(cvec),W_DARRAY(svec),W_DARRAY(wkvec),
		       fst_matrix* zmat,W_DARRAY(yvec),int* stat,
		       dcopy_type f_dcopy,drotg_type f_drotg,drot_type f_drot,
		       W_ERRORARGS)

{
  blasint_t i,n,r=0,stp,sz,ione=1;
  double temp;
  bool islower;
  double* tbuff;

  /* Read arguments */
  if (ain<5 || ain>7)
    W_RETERROR(2,"Wrong number of input arguments");
  if (aout!=1)
    W_RETERROR(2,"Need one return argument");
  if ((n=lmat->n)!=lmat->m || lmat->n==0 ||
      (!(islower=(UPLO(lmat->strcode)=='L')) && UPLO(lmat->strcode)!='U'))
    W_RETERROR(1,"L: Wrong size or structure code");
  if (nvvec!=n)
    W_RETERROR(1,"VEC: Wrong size");
  if (ncvec!=n || nsvec!=n)
    W_RETERROR(1,"CVEC, SVEC: Wrong size");
  if (ain>5) {
    if (ain<7)
      W_RETERROR(1,"Need both Z, Y or none");
    r=zmat->m;
    if (zmat->n!=n || r==0)
      W_RETERROR(1,"Z: Wrong size");
    W_CHKSIZE(yvec,r,"Y");
  } else {
    zmat=0; yvec=0; nyvec=0;
  }
  if (nwkvec<n || nwkvec<r)
    W_RETERROR(1,"WORKV: Wrong size");
  *stat=0; /* OK so far */

  /* Generate Givens rotations, update L */
  f_dcopy(&n,vvec,&ione,wkvec,&ione);
  stp = islower?1:lmat->stride;
  for (i=0,sz=n,tbuff=lmat->buff; i<n-1; i++) {
    /* drotg(a,b,c,s): J = [c s; -s c], s.t. J [a; b] = [r; 0]
       a overwritten by r, b by some other information (NOT 0!) */
    if (*tbuff==0.0 && wkvec[i]==0.0) {
      *stat=1; break;
    }
    f_drotg(tbuff,wkvec+i,cvec+i,svec+i);
    /* Do not want negative elements on factor diagonal */
    if ((temp=*tbuff)<0.0) {
      *tbuff=-temp; cvec[i]=-cvec[i]; svec[i]=-svec[i];
    } else if (temp==0.0) {
      *stat=1; break;
    }
    /* drot(x,y,c,s): J = [c s; -s c]. [x_i; y_i] overwritten by
       J [x_i; y_i], for all i
       BAD: Slower for upper triangular! */
    sz--;
    f_drot(&sz,tbuff+stp,&stp,wkvec+(i+1),&ione,cvec+i,svec+i);
    tbuff+=(lmat->stride+1);
  }
  if (*stat==0 && (*tbuff!=0.0 || wkvec[n-1]!=0.0)) {
    f_drotg(tbuff,wkvec+(n-1),cvec+i,svec+i);
    if ((temp=*tbuff)<0.0) {
      *tbuff=-temp; cvec[i]=-cvec[i]; svec[i]=-svec[i];
    } else if (temp==0.0)
      *stat=1;
  } else
    *stat=1;

  /* Dragging along */
  if (r>0 && *stat==0) {
    f_dcopy(&r,yvec,&ione,wkvec,&ione);
    for (i=0; i<n; i++)
      f_drot(&r,zmat->buff+(i*zmat->stride),&ione,wkvec,&ione,cvec+i,
	     svec+i);
  }

  W_RETOK;
}
