/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMARGINALS
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_FACT_COMPMARGINALS_H
#define EPTWRAP_FACT_COMPMARGINALS_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_fact_compmarginals(int ain,int aout,int n,int m,
				  W_IARRAY(rp_rowind),W_IARRAY(rp_colind),
				  W_DARRAY(rp_bvals),W_DARRAY(rp_pi),
				  W_DARRAY(rp_beta),W_DARRAY(margpi),
				  W_DARRAY(margbeta),W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
