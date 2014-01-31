/* -------------------------------------------------------------------
 * EPTWRAP_FACT_COMPMAXPI
 *
 * EP with factorized Gaussian backbone.
 * Computes top-K values in 'FactEPMaximumPiValues' data structure from
 * scratch ('FactEPMaximumPiValues::recompute').
 * This data structure is used for selective damping, see
 * EPTOOLS_FACT_SEQUPDATES.
 *
 * If SD_SUBIND is given, it is a subset of 0:(M-1), sorted in
 * ascending order. See 'FactEPMaximumPiValues', fields 'subInd' and
 * 'subExcl'
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array]
 * - RP_BETA:     " [double array]
 * - SD_K:        Value K (must be >1)
 * - SD_SUBIND    See above. Optional [int32 array]
 * - SD_SUBEXCL   ". Def.: 0 [int]
 *
 * Return:
 * - SD_NUMVALID: Max pi data structure [int32 array]
 * - SD_TOPIND:   " [int32 array]
 * - SD_TOPVAL:   " [double array]
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_fact_compmaxpi.h"
#include "src/eptools/FactorizedEPRepresentation.h"
#include "src/eptools/FactEPMaximumPiValues.h"

void eptwrap_fact_compmaxpi(int ain,int aout,int n,int m,W_IARRAY(rp_rowind),
			    W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			    W_DARRAY(rp_pi),W_DARRAY(rp_beta),int sd_k,
			    W_IARRAY(sd_subind),int sd_subexcl,
			    W_IARRAY(sd_numvalid),W_IARRAY(sd_topind),
			    W_DARRAY(sd_topval),W_ERRORARGS)
{
  int i;
  Handle<FactorizedEPRepresentation> epRepr;
  ArrayHandle<int> sd_numvalidA,sd_topindA,sd_subindA;
  ArrayHandle<double> sd_topvalA;
  Handle<FactEPMaximumPiValues> epMaxPi;

  try {
    /* Read arguments */
    if (ain<8 || ain>10)
      W_RETERROR(2,"Wrong number of input arguments");
    if (aout!=3)
      W_RETERROR(2,"Need 3 return arguments");
    if (sd_k<=1)
      W_RETERROR(1,"SD_K: Must be >1");
    if (ain<10)
      sd_subexcl=0;
    if (ain>8) {
      if (nsd_subind==0 || nsd_subind>m)
	W_RETERROR(1,"SD_SUBIND: Wrong size");
    } else
      sd_subind=0;
    /* Representation */
    createFactEPRepres(n,m,W_ARR(rp_rowind),W_ARR(rp_colind),W_ARR(rp_bvals),
		       W_ARR(rp_pi),W_ARR(rp_beta),epRepr,W_ERRARGS);
    /* Return arguments */
    W_CHKSIZE(sd_numvalid,n,"SD_NUMVALID");
    i=n*(sd_k+1);
    W_CHKSIZE(sd_topind,i,"SD_TOPIND");
    W_CHKSIZE(sd_topval,i,"SD_TOPVAL");
    /* Create max pi data structure */
    W_MASKARRAY(sd_numvalid);
    W_MASKARRAY(sd_topind);
    W_MASKARRAY(sd_topval);
    if (sd_subind!=0)
      W_MASKARRAY(sd_subind);
    for (i=0; i<n; i++)
      sd_numvalid[i]=1; // Just to make constructor happy
    try {
      epMaxPi.changeRep(new FactEPMaximumPiValues(epRepr,sd_k,sd_numvalidA,
						  sd_topindA,sd_topvalA,
						  sd_subindA,sd_subexcl!=0));
      epMaxPi->recompute(); // Recompute from scratch
    } catch (StandardException ex) {
      W_RETERROR_ARGS(1,"Cannot create FactEPMaximumPiValues (selective damping):\n%s",ex.msg());
    }
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
