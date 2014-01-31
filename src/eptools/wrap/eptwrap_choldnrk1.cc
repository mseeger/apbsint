/* -------------------------------------------------------------------
 * EPTOOLS_CHOLDNRK1
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * If
 *   A = L*L', A_ = A - v*v' = L_ L_', A, L n-by-n,
 * where L is lower triangular, the method computes L_ from L. This
 * is called Cholesky rank one downdate. L or L' (upper triangular) can
 * be passed, only the relevant triagle is accessed. L (or L') is
 * passed in L, v in VEC.
 * We require p = L\v. If ISP==true, VEC contains p rather than v.
 * Otherwise, p is computed locally, stored in WORKV.
 * NOTE: For the present implementation, the method is more efficient
 * when a lower triangular matrix is used.
 *
 * Dragging along:
 * If Z (r-by-n) is given, so must be the r-vector y. In this case,
 * we overwrite Z by Z_, where
 *   Z_ L_' = Z L' - y v'.
 *
 * Working array:
 * Requires a working vector of size >= max(n,r), passed in WORKV.
 * NOTE: The same vector can be passed for VEC and WORKV, in which case
 * VEC is overwritten in an undefined way. If r > n, VEC can be of size
 * r, containing v (or p) in the first n components.
 * The method uses n Givens rotations, param. by angles c_k, s_k. The
 * change L -> L' is specified by these 2*n numbers. They are written
 * into CVEC, SVEC.
 * NOTE: The original routine LINPACK dchdd can produce negative
 * values in DIAG(L), something we want to avoid here. If this happens,
 * the corr. column of L_ is flipped. In the present implementation, this
 * is not reported back, so the change L -> L' cannot always be
 * reconstructed from CVEC, SVEC alone.
 *
 * Input:
 * - L:     Factor L (or L'), overwritten by L_ (or L_'). Must be
 *          lower (upper) triangular, str. code UPLO
 * - VEC:   Vector v. Can have size >n, only first n elem. are used
 * - CVEC:  Vector [n]. c_k ret. here
 * - SVEC:  Vector [n]. s_k ret. here
 * - WORKV: Working vector of size max(n,r). Can be same as VEC
 * - ISP:   S.a. Def.: false
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
#include "src/eptools/wrap/eptwrap_choldnrk1.h"
#include <stdlib.h>
#include <math.h>

/*
 * The method is adapted from LINPACK dchdd. We did the following
 * modifications:
 * - Using BLAS drot in order to avoid any explicit O(n^2) loops
 * - Keeping diag(L_) positive, by flipping columns of L_ whenever a
 *   negative element pops up
 * See the TR
 *   M. Seeger
 *   Low Rank Updates for the Cholesky Decomposition
 */

void eptwrap_choldnrk1(int ain,int aout,fst_matrix* lmat,W_DARRAY(vvec),
		       W_DARRAY(cvec),W_DARRAY(svec),W_DARRAY(wkvec),int isp,
		       fst_matrix* zmat,W_DARRAY(yvec),int* stat,
		       dcopy_type f_dcopy,dtrsv_type f_dtrsv,ddot_type f_ddot,
		       drotg_type f_drotg,drot_type f_drot,dscal_type f_dscal,
		       daxpy_type f_daxpy,W_ERRORARGS)

{
  blasint_t i,n,r=0,stp,sz,ione=1,nxi,npos;
  double qs,cval,sval,c1,c2;
  bool islower;
  double* tbuff,*zcol;
  const char* diag="N";
  char trans[2];
  int* flind=0;

  /* Read arguments */
  if (ain<5 || ain>8)
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
  if (ain<7) {
    zmat=0; yvec=0; nyvec=0;
    if (ain<6)
      isp=0;
  }
  if (ain>6) {
    if (ain<8)
      W_RETERROR(1,"Need both Z, Y or none");
    r=zmat->m;
    if (zmat->n!=n || r==0)
      W_RETERROR(1,"Z: Wrong size");
    W_CHKSIZE(yvec,r,"Y");
  }
  if (nwkvec<n || nwkvec<r)
    W_RETERROR(1,"WORKV: Wrong size");
  *stat=0; /* OK so far */

  /* Compute p (if not given)
   * NOTE: This requires dtrsv, which may not be given
   */
  f_dcopy(&n,vvec,&ione,wkvec,&ione);
  if (!isp) {
    if (f_dtrsv==0)
      W_RETERROR(2,"Internal error: Need BLAS dtrsv");
    trans[1]=0;
    trans[0]=islower?'N':'T';
    i=(blasint_t) lmat->stride;
    f_dtrsv(&UPLO(lmat->strcode),trans,(char*) diag,&n,lmat->buff,&i,wkvec,
	    &ione);
  }
  /* Generate Givens rotations */
  qs = 1.0 - f_ddot(&n,wkvec,&ione,wkvec,&ione);
  if (qs<=0.0)
    *stat=1;
  else {
    qs = sqrt(qs);
    for (i=n-1; i>=0; i--) {
      f_drotg(&qs,wkvec+i,cvec+i,svec+i);
      /* 'qs' must remain positive */
      if (qs<0.0) {
	qs=-qs; cvec[i]=-cvec[i]; svec[i]=-svec[i];
      }
    }
    /* NOTE: 'qs' should be 1 now */

    /* Update L. If there are any flips of L_ cols, we alloc. 'flind' are
       store their pos. there */
    fillVec(wkvec,n,0.0);
    stp = islower?1:lmat->stride;
    for (i=n-1,sz=0,tbuff=lmat->buff+((n-1)*(lmat->stride+1)); i>=0; i--) {
      /* BAD: Slower for upper triangular! */
      sz++;
      if (*tbuff<=0.0) {
	*stat=1; break;
      }
      f_drot(&sz,wkvec+i,&ione,tbuff,&stp,cvec+i,svec+i);
      /* Do not want negative elements on diagonal */
      if (*tbuff<0.0) {
	if (flind==0) {
	  /* Does this ever happen?
	     Allocate 'flind'. Size n, to make sure. Grows from the right */
	  flind = (int*) malloc(n*sizeof(int));
	  npos = n;
	}
	flind[--npos]=i;
	qs=-1.0;
	f_dscal(&sz,&qs,tbuff,&stp);
      } else if (*tbuff==0.0) {
	*stat=1; break;
      }
      tbuff-=(lmat->stride+1);
    }
    /* NOTE: Should have v in 'wkvec' now */

    /* Dragging along */
    if (r>0 && *stat==0) {
      f_dcopy(&r,yvec,&ione,wkvec,&ione);
      nxi=(flind!=0)?flind[npos]:-1;
      for (i=0; i<n; i++) {
	zcol=zmat->buff+(i*zmat->stride);
	cval=cvec[i]; sval=svec[i];
	qs=-sval;
	f_daxpy(&r,&qs,wkvec,&ione,zcol,&ione);
	if (nxi==i) {
	  if (++npos<n) nxi=flind[npos];
	  c1=-1.0/cval; c2=sval;
	} else {
	  c1=1.0/cval; c2=-sval;
	}
	f_dscal(&r,&c1,zcol,&ione);
	if (i<n-1) {
	  f_dscal(&r,&cval,wkvec,&ione);
	  f_daxpy(&r,&c2,zcol,&ione,wkvec,&ione);
	}
      }
    }
  }

  if (flind!=0)
    free((void*) flind);
  W_RETOK;
}
