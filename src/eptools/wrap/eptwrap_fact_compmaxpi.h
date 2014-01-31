/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMAXPI
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_FACT_COMPMAXPI_H
#define EPTWRAP_FACT_COMPMAXPI_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_fact_compmaxpi(int ain,int aout,int n,int m,W_IARRAY(rp_rowind),
			      W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			      W_DARRAY(rp_pi),W_DARRAY(rp_beta),int sd_k,
			      W_IARRAY(sd_subind),int sd_subexcl,
			      W_IARRAY(sd_numvalid),W_IARRAY(sd_topind),
			      W_DARRAY(sd_topval),W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
