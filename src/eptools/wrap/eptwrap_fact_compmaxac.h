/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMAXAC
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_FACT_COMPMAXAC_H
#define EPTWRAP_FACT_COMPMAXAC_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_fact_compmaxac(int ain,int aout,int n,int m,W_IARRAY(rp_rowind),
			      W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			      W_DARRAY(rp_pi),W_DARRAY(rp_beta),
			      W_IARRAY(rp_tauind),W_DARRAY(rp_a),W_DARRAY(rp_c),
			      int sda_k,int sdc_k,W_IARRAY(sda_numvalid),
			      W_IARRAY(sda_topind),W_DARRAY(sda_topval),
			      W_IARRAY(sdc_numvalid),W_IARRAY(sdc_topind),
			      W_DARRAY(sdc_topval),W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
