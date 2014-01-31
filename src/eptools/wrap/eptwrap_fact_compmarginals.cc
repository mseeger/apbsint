/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMARGINALS
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * EP with factorized Gaussian backbone.
 * Compute marginals on variables from EP (message) parameters, overwrite
 * MARGPI, MARGBETA.
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array]
 * - RP_BETA:     " [double array]
 * - MARGPI:      Marginal pi parameters written here
 * - MARGBETA:    Marginal beta parameters written here
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_fact_compmarginals.h"
#include "src/eptools/FactorizedEPRepresentation.h"

void eptwrap_fact_compmarginals(int ain,int aout,int n,int m,
				W_IARRAY(rp_rowind),W_IARRAY(rp_colind),
				W_DARRAY(rp_bvals),W_DARRAY(rp_pi),
				W_DARRAY(rp_beta),W_DARRAY(margpi),
				W_DARRAY(margbeta),W_ERRORARGS)
{
  Handle<FactorizedEPRepresentation> epRepr;

  try {
    /* Read arguments */
    if (ain!=9)
      W_RETERROR(2,"Need 9 input arguments");
    if (aout!=0)
      W_RETERROR(2,"No return arguments");
    W_CHKSIZE(margpi,n,"MARGPI");
    W_CHKSIZE(margbeta,n,"MARGBETA");
    createFactEPRepres(n,m,W_ARR(rp_rowind),W_ARR(rp_colind),W_ARR(rp_bvals),
		       W_ARR(rp_pi),W_ARR(rp_beta),epRepr,W_ERRARGS);
    /* Compute marginals */
    epRepr->compMarginals(margbeta,margpi);
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
