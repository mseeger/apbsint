/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMARGINALS_BVPREC
 *
 * Extension of EPTWRAP_FACT_COMPMARGINALS for models including
 * bivariate precision potentials.
 * Marginals are computed for x variables (MARGPI, MARGBETA) and tau
 * variables (MARGA, MARGC).
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array]
 * - RP_BETA:     " [double array]
 * - RP_TAUIND    " [int32 array]
 * - RP_A:        " [double array]
 * - RP_C:        " [double array]
 * - MARGPI:      Marginal pi parameters written here
 * - MARGBETA:    " (beta)
 * - MARGA:       " (a)
 * - MARGC:       " (c)
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_fact_compmarginals_bvprec.h"
#include "src/eptools/FactorizedEPRepresentation.h"

void eptwrap_fact_compmarginals_bvprec(int ain,int aout,int n,int m,
				       W_IARRAY(rp_rowind),W_IARRAY(rp_colind),
				       W_DARRAY(rp_bvals),W_DARRAY(rp_pi),
				       W_DARRAY(rp_beta),W_IARRAY(rp_tauind),
				       W_DARRAY(rp_a),W_DARRAY(rp_c),
				       W_DARRAY(margpi),W_DARRAY(margbeta),
				       W_DARRAY(marga),W_DARRAY(margc),
				       W_ERRORARGS)
{
  Handle<FactorizedEPRepresentation> epRepr;

  try {
    /* Read arguments */
    if (ain!=14)
      W_RETERROR(2,"Need 14 input arguments");
    if (aout!=0)
      W_RETERROR(2,"No return arguments");
    W_CHKSIZE(margpi,n,"MARGPI");
    W_CHKSIZE(margbeta,n,"MARGBETA");
    createFactEPRepres_bvprec(n,m,W_ARR(rp_rowind),W_ARR(rp_colind),
			      W_ARR(rp_bvals),W_ARR(rp_pi),W_ARR(rp_beta),
			      W_ARR(rp_tauind),W_ARR(rp_a),W_ARR(rp_c),epRepr,
			      W_ERRARGS);
    int numk=epRepr->numPrecVariables();
    if (numk==0)
      W_RETERROR(1,"Must have bivariate precision potentials");
    /* Compute marginals */
    W_CHKSIZE(marga,numk,"MARGA");
    W_CHKSIZE(margc,numk,"MARGC");
    epRepr->compTauMarginals(marga,margc);
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
