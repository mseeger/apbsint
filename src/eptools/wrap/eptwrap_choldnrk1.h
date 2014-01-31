/* -------------------------------------------------------------------
 * EPTWRAP_CHOLDNRK1
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_CHOLDNRK1_H
#define EPTWRAP_CHOLDNRK1_H

#include "src/eptools/wrap/eptools_helper_macros.h"
#include "src/eptools/wrap/matrix_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_choldnrk1(int ain,int aout,fst_matrix* lmat,W_DARRAY(vvec),
			 W_DARRAY(cvec),W_DARRAY(svec),W_DARRAY(wkvec),int isp,
			 fst_matrix* zmat,W_DARRAY(yvec),int* stat,
			 dcopy_type f_dcopy,dtrsv_type f_dtrsv,
			 ddot_type f_ddot,drotg_type f_drotg,drot_type f_drot,
			 dscal_type f_dscal,daxpy_type f_daxpy,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
