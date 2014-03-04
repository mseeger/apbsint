/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_PARALLEL_BVPREC
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_EPUPDATE_PARALLEL_BVPREC_H
#define EPTWRAP_EPUPDATE_PARALLEL_BVPREC_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_epupdate_parallel_bvprec(int ain,int aout,W_IARRAY(potids),
					W_IARRAY(numpot),W_DARRAY(parvec),
					W_IARRAY(parshrd),W_ARRAY(annobj,void*),
					W_DARRAY(cmu),W_DARRAY(crho),
					W_DARRAY(ca),W_DARRAY(cc),
					W_IARRAY(updind),W_IARRAY(rstat),
					W_DARRAY(alpha),W_DARRAY(nu),
					W_DARRAY(hata),W_DARRAY(hatc),
					W_DARRAY(logz),W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
