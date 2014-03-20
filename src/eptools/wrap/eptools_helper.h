/* -------------------------------------------------------------------
 * Helper functions for EPTOOLS wrapper functions
 * Declarations
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_HELPER_H
#define EPTOOLS_HELPER_H

// Macros for wrapper functions
#include "src/main.h"
#include "src/eptools/wrap/eptools_helper_macros.h"
#include "src/eptools/wrap/eptools_helper_basic.h"

// Helper functions

class PotentialManager;
class FactorizedEPRepresentation;

void createPotentialManager(W_IARRAY(potids),W_IARRAY(numpot),W_DARRAY(parvec),
			    W_IARRAY(parshrd),W_ARRAY(annobj,void*),
			    Handle<PotentialManager>& potMan,W_ERRORARGS);

void createFactEPRepres(int numN,int numM,W_IARRAY(rp_rowind),
			W_IARRAY(rp_colind),W_DARRAY(rp_bvals),W_DARRAY(rp_pi),
			W_DARRAY(rp_beta),
			Handle<FactorizedEPRepresentation>& epRepr,W_ERRORARGS);

void createFactEPRepres_bvprec(int numN,int numM,W_IARRAY(rp_rowind),
			       W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			       W_DARRAY(rp_pi),W_DARRAY(rp_beta),
			       W_IARRAY(rp_tauind),W_DARRAY(rp_a),
			       W_DARRAY(rp_c),
			       Handle<FactorizedEPRepresentation>& epRepr,
			       W_ERRORARGS);

#endif
