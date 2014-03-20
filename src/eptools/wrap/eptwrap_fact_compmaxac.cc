/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMAXAC
 *
 * Same as EPTWRAP_FACT_COMPMAXPI, but for Gamma a, c parameters
 * corresponding to bivariate precision potentials. Uses
 * 'FactEPMaximumAValues', 'FactEPMaximumCValues' in place of
 * 'FactEPMaximumPiValues'.
 *
 * NOTE: RP_ROWIND, RP_COLIND, RP_BVALS, RP_PI, RP_BETA are not used,
 * but have to be passed in order to create a
 * 'FactorizedEPRepresentation' object.
 *
 * Input:
 * - N:            Number of variables
 * - M:            Number of factors
 * - RP_ROWIND:    Factorized EP representation [int32 array]
 * - RP_COLIND:    " [int32 array]
 * - RP_BVALS:     " [double array]
 * - RP_PI:        " [double array]
 * - RP_BETA:      " [double array]
 * - RP_TAUIND     " [int32 array]
 * - RP_A:         " [double array]
 * - RP_C:         " [double array]
 * - SDA_K:        Value K for max a (must be >1)
 * - SDC_K:        Value K for max c (must be >1)
 *
 * Return:
 * - SDA_NUMVALID: Max a data structure [int32 array]
 * - SDA_TOPIND:   " [int32 array]
 * - SDA_TOPVAL:   " [double array]
 * - SDC_NUMVALID: Max c data structure [int32 array]
 * - SDC_TOPIND:   " [int32 array]
 * - SDC_TOPVAL:   " [double array]
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_fact_compmaxac.h"
#include "src/eptools/FactorizedEPRepresentation.h"
#include "src/eptools/FactEPMaximumAValues.h"
#include "src/eptools/FactEPMaximumCValues.h"

void eptwrap_fact_compmaxac(int ain,int aout,int n,int m,W_IARRAY(rp_rowind),
			    W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			    W_DARRAY(rp_pi),W_DARRAY(rp_beta),
			    W_IARRAY(rp_tauind),W_DARRAY(rp_a),W_DARRAY(rp_c),
			    int sda_k,int sdc_k,W_IARRAY(sda_numvalid),
			    W_IARRAY(sda_topind),W_DARRAY(sda_topval),
			    W_IARRAY(sdc_numvalid),W_IARRAY(sdc_topind),
			    W_DARRAY(sdc_topval),W_ERRORARGS)
{
  int i;
  Handle<FactorizedEPRepresentation> epRepr;
  ArrayHandle<int> sda_numvalidA,sda_topindA;
  ArrayHandle<int> sdc_numvalidA,sdc_topindA;
  ArrayHandle<double> sda_topvalA,sdc_topvalA;
  Handle<FactEPMaximumAValues> epMaxA;
  Handle<FactEPMaximumCValues> epMaxC;

  try {
    /* Read arguments */
    if (ain!=12)
      W_RETERROR(2,"Need 12 input arguments");
    if (aout!=6)
      W_RETERROR(2,"Need 6 return arguments");
    if (sda_k<=1)
      W_RETERROR(1,"SDA_K: Must be >1");
    if (sdc_k<=1)
      W_RETERROR(1,"SDC_K: Must be >1");
    /* Representation */
    createFactEPRepres_bvprec(n,m,W_ARR(rp_rowind),W_ARR(rp_colind),
			      W_ARR(rp_bvals),W_ARR(rp_pi),W_ARR(rp_beta),
			      W_ARR(rp_tauind),W_ARR(rp_a),W_ARR(rp_c),epRepr,
			      W_ERRARGS);
    int numk=epRepr->numPrecVariables();
    if (numk==0)
      W_RETERROR(1,"Must have bivariate precision potentials");
    /* Return arguments */
    W_CHKSIZE(sda_numvalid,numk,"SDA_NUMVALID");
    i=numk*(sda_k+1);
    W_CHKSIZE(sda_topind,i,"SDA_TOPIND");
    W_CHKSIZE(sda_topval,i,"SDA_TOPVAL");
    W_MASKARRAY(sda_numvalid);
    W_MASKARRAY(sda_topind);
    W_MASKARRAY(sda_topval);
    W_CHKSIZE(sdc_numvalid,numk,"SDC_NUMVALID");
    i=numk*(sdc_k+1);
    W_CHKSIZE(sdc_topind,i,"SDC_TOPIND");
    W_CHKSIZE(sdc_topval,i,"SDC_TOPVAL");
    /* Create max data structures */
    for (i=0; i<numk; i++)
      sda_numvalid[i]=sdc_numvalid[i]=1;
    try {
      epMaxA.changeRep(new FactEPMaximumPiValues(epRepr,sda_k,sda_numvalidA,
						 sda_topindA,sda_topvalA));
      epMaxA->recompute(); // Recompute from scratch
    } catch (StandardException ex) {
      W_RETERROR_ARGS(1,"Cannot create FactEPMaximumAValues (selective damping):\n%s",ex.msg());
    }
    try {
      epMaxC.changeRep(new FactEPMaximumPiValues(epRepr,sdc_k,sdc_numvalidA,
						 sdc_topindA,sdc_topvalA));
      epMaxC->recompute(); // Recompute from scratch
    } catch (StandardException ex) {
      W_RETERROR_ARGS(1,"Cannot create FactEPMaximumCValues (selective damping):\n%s",ex.msg());
    }
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
