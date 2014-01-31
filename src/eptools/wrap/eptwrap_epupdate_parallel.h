/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_PARALLEL
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_EPUPDATE_PARALLEL_H
#define EPTWRAP_EPUPDATE_PARALLEL_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_epupdate_parallel(int ain,int aout,W_IARRAY(potids),
				 W_IARRAY(numpot),W_DARRAY(parvec),
				 W_IARRAY(parshrd),W_DARRAY(cmu),
				 W_DARRAY(crho),W_IARRAY(updind),
				 W_IARRAY(rstat),W_DARRAY(alpha),W_DARRAY(nu),
				 W_DARRAY(logz),W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
