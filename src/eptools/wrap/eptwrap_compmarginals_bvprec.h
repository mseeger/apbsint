/* -------------------------------------------------------------------
 * EPTWRAP_COMPMARGINALS_BVPREC
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_COMPMARGINALS_BVPREC_H
#define EPTWRAP_COMPMARGINALS_BVPREC_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_compmarginals_bvprec(int ain,int aout,W_IARRAY(rp_tauind),
				    W_DARRAY(rp_a),W_DARRAY(rp_c),
				    W_DARRAY(marga),W_DARRAY(margc),
				    W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
