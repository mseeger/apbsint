/* -------------------------------------------------------------------
 * EPTWRAP_CHOLUPRK1
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_CHOLUPRK1_H
#define EPTWRAP_CHOLUPRK1_H

#include "src/eptools/wrap/eptools_helper_macros.h"
#include "src/eptools/wrap/matrix_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_choluprk1(int ain,int aout,fst_matrix* lmat,W_DARRAY(vvec),
			 W_DARRAY(cvec),W_DARRAY(svec),W_DARRAY(wkvec),
			 fst_matrix* zmat,W_DARRAY(yvec),int* stat,
			 dcopy_type f_dcopy,drotg_type f_drotg,drot_type f_drot,
			 W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
