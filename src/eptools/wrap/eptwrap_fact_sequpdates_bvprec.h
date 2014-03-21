/* -------------------------------------------------------------------
 * EPTWRAP_FACT_SEQUPDATES_BVPREC
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_FACT_SEQUPDATES_BVPREC_H
#define EPTWRAP_FACT_SEQUPDATES_BVPREC_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_fact_sequpdates_bvprec(int ain,int aout,int n,int m,
				      W_IARRAY(updjind),W_IARRAY(pm_potids),
				      W_IARRAY(pm_numpot),W_DARRAY(pm_parvec),
				      W_IARRAY(pm_parshrd),
				      W_ARRAY(pm_annobj,void*),
				      W_IARRAY(rp_rowind),W_IARRAY(rp_colind),
				      W_DARRAY(rp_bvals),W_DARRAY(rp_pi),
				      W_DARRAY(rp_beta),W_IARRAY(rp_tauind),
				      W_DARRAY(rp_a),W_DARRAY(rp_c),
				      W_DARRAY(margpi),W_DARRAY(margbeta),
				      W_DARRAY(marga),W_DARRAY(margc),
				      double piminthres,double aminthres,
				      double cminthres,double dampfact,
				      W_IARRAY(sd_numvalid),W_IARRAY(sd_topind),
				      W_DARRAY(sd_topval),
				      W_IARRAY(sda_numvalid),
				      W_IARRAY(sda_topind),W_DARRAY(sda_topval),
				      W_IARRAY(sdc_numvalid),
				      W_IARRAY(sdc_topind),
				      W_DARRAY(sdc_topval),W_IARRAY(sd_subind),
				      int sd_subexcl,W_IARRAY(rstat),
				      W_DARRAY(delta),W_DARRAY(sd_dampfact),
				      int* sd_nupd,int* sd_nrec,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
