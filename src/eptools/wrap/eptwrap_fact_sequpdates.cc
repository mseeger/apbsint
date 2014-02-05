/* -------------------------------------------------------------------
 * EPTWRAP_FACT_SEQUPDATES
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * EP with factorized Gaussian backbone. Run a number of updates on
 * potentials, one of the after (sequential updating).
 * Operates on ([I]: Input, [I/O]: Input/output, some vectors are
 * overwritten):
 * - Potential manager [I]. PM_POTIDS, PM_NUMPOT, PM_PARVEC, PM_PARSHRD,
 *   PM_ANNOBJ, see EPTOOLS_EPUPDATE_PARALLEL comments.
 *   NOTE: Need full potential manager, even if updates only done on
 *   subset of potentials.
 * - Representation [I/O]. Structure and content of coupling factor
 *   B [I], EP (message) parameters [I/O]
 * - Marginals on variables [I/O]. MARGPI, MARGBETA
 * - Support data structure for selective damping mechanism [I/O].
 *   SD_NUMVALID, SD_TOPIND, SD_TOPVAL
 *
 * There are M potentials (factors), N variables. We run EP updates on
 * potentials in UPDJIND, one after the other. Messages and marginals
 * are factorized Gaussians, given by natural parameters (pi, beta).
 * An update modifies marginals MARGPI, MARGBETA and EP (message)
 * parameters RP_PI, RP_BETA. If DAMPFACT>0, the update is damped. There
 * may also be selective damping (see below).
 * Updates can fail for various reasons. They are either skipped, or
 * selective damping is applied:
 * - Cavity marginal undefined: If pi < PIMINTHRES/2
 * - New marginal undefined: If pi < PIMINTHRES/2
 * RSTAT is return status for each update. Codes defined in
 * 'FactorizedEPDriver':
 * - 0 [updSuccess]: Update successful
 * - 1 [updCavityInvalid]: Cavity marginal undefined. Update skipped
 * - 2 [updNumericalError]: Local EP update raises error. Update skipped
 * - 3 [updMarginalsInvalid]: New marginals undefined. Update skipped
 * - 4 [updCavCondSkipped]: Selective damping requires skipping
 * DELTA is relative change in moments for each non-skipped update, or 0
 * for skipped ones. The entry is the maximum relative difference for
 * means and stddevs (before and after update).
 *
 * Representation:
 * Consists of RP_ROWIND, RP_COLIND, RP_BVALS, RP_PI, RP_BETA. Details in
 * 'FactorizedEPRepresentation' comments. Internal representation
 * automatically compiled by Matlab code, see EPT.BFACT_INTREPRES.
 * RP_PI, RP_BETA (EP parameters) are I/O, their content is overwritten.
 *
 * Selective damping (optional):
 * SD_NUMVALID, SD_TOPIND, SD_TOPVAL, SD_SUBIND, SD_SUBEXCL. Details in
 * technical report, 'FactorizedEPDriver', 'FactEPMaximumPiValues'
 * comments. Idea is to ensure that for all EP parameters and marginals:
 * pi >= PIMINTHRES. This is a precondition (not checked). If the
 * condition is violated after an update, we apply the
 * minimum amount of extra damping (may be on top of DAMPFACT). In the
 * extreme case, the update is skipped (RSTAT: updCavCondSkipped).
 * Effective damping factor used for each update can be returned in
 * SD_DAMPFACT.
 * Requires 'FactEPMaximumPiValues' data structure, represented by SD_XXX
 * variables (I/O, content is overwritten).
 * SD_NUPD, SD_NREC return statistics about this datastructure (number of
 * update calls and block recomputations).
 *
 * Input:
 * - N:           Number of variables
 * - M:           Number of factors
 * - UPDJIND:     Update on these potentials, in order [int32 array]
 * - PM_POTIDS:   Potential manager [int32 array]
 * - PM_NUMPOT:   " [int32 array]
 * - PM_PARVEC:   " [double array]
 * - PM_PARSHRD:  " [int32 array]
 * - PM_ANNOBJ:   " [void* array]
 * - RP_ROWIND:   Factorized EP representation [int32 array]
 * - RP_COLIND:   " [int32 array]
 * - RP_BVALS:    " [double array]
 * - RP_PI:       " [double array; I/O]
 * - RP_BETA:     " [double array; I/O]
 * - MARGPI:      Variable marginals [I/O]
 * - MARGBETA:    " [I/O]
 * - PIMINTHRES:  See above. Positive
 * - DAMPFACT:    Damping factor, in [0,1). Optional, def. is 0
 * - SD_NUMVALID: Selective damping. Optional [int32 array; I/O]
 * - SD_TOPIND:   " [int32 array; I/O]
 * - SD_TOPVAL:   " [double array; I/O]
 * - SD_SUBIND    " [int32 array]
 * - SD_SUBEXCL   ". Def.: false
 *
 * Return:
 * - RSTAT:       Return stati for each update. Optional [int32]
 * - DELTA:       See above. Optional
 * - SD_DAMPFACT: See above. Optional, only if selective damping
 * - SD_NUPD:     " [int32]
 * - SD_NREC:     " [int32]
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_fact_sequpdates.h"
#include "src/eptools/FactorizedEPDriver.h"
#include "src/eptools/FactEPMaximumPiValues.h"

void eptwrap_fact_sequpdates(int ain,int aout,int n,int m,W_IARRAY(updjind),
			     W_IARRAY(pm_potids),W_IARRAY(pm_numpot),
			     W_DARRAY(pm_parvec),W_IARRAY(pm_parshrd),
			     W_ARRAY(pm_annobj,void*),W_IARRAY(rp_rowind),
			     W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			     W_DARRAY(rp_pi),W_DARRAY(rp_beta),W_DARRAY(margpi),
			     W_DARRAY(margbeta),double piminthres,
			     double dampfact,W_IARRAY(sd_numvalid),
			     W_IARRAY(sd_topind),W_DARRAY(sd_topval),
			     W_IARRAY(sd_subind),int sd_subexcl,
			     W_IARRAY(rstat),W_DARRAY(delta),
			     W_DARRAY(sd_dampfact),int* sd_nupd,int* sd_nrec,
			     W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain<16 || ain>22)
      W_RETERROR(2,"Wrong number of input arguments");
    if (aout>5)
      W_RETERROR(2,"Too many return arguments");
    if (n<1) W_RETERROR(1,"N wrong");
    if (m<1) W_RETERROR(1,"M wrong");
    if (nupdjind==0)
      W_RETERROR(1,"UPDJIND must not be empty");
    Interval<int> ivM(0,m-1,IntVal::ivClosed,IntVal::ivClosed);
    if (ivM.check(updjind,nupdjind)!=0)
      W_RETERROR(1,"UPDJIND: Entries of out range");
    //printMsgStdout("Point 1");
    /* Potential manager */
    Handle<PotentialManager> potMan;
    createPotentialManager(W_ARR(pm_potids),W_ARR(pm_numpot),W_ARR(pm_parvec),
			   W_ARR(pm_parshrd),W_ARR(pm_annobj),potMan,W_ERRARGS);
    if (potMan->size()!=m)
      W_RETERROR(1,"PM_*: Potential manager has wrong size");
    /* Representation of B */
    //printMsgStdout("Point 2");
    Handle<FactorizedEPRepresentation> epRepr;
    createFactEPRepres(n,m,W_ARR(rp_rowind),W_ARR(rp_colind),W_ARR(rp_bvals),
		       W_ARR(rp_pi),W_ARR(rp_beta),epRepr,W_ERRARGS);
    /* Variable marginals */
    ArrayHandle<double> margpiA,margbetaA;
    //printMsgStdout("Point 3");
    W_CHKSIZE(margpi,n,"MARGPI");
    W_CHKSIZE(margpi,n,"MARGBETA");
    W_MASKARRAY(margpi);
    W_MASKARRAY(margbeta);
    if (piminthres<=0.0)
      W_RETERROR(1,"PIMINTHRES must be positive");
    int sd_k=0; // K of selective damping (0 if not active)
    ArrayHandle<int> sd_numvalidA,sd_topindA,sd_subindA;
    ArrayHandle<double> sd_topvalA;
    if (ain>16) {
      if (dampfact<0.0 || dampfact>=1.0)
	W_RETERROR(1,"DAMPFACT: Out of range");
      if (ain>17) {
	// Selective damping
	//printMsgStdout("Point 4");
	if (ain<20)
	  W_RETERROR(1,"Need all SD_XXX or none");
	W_CHKSIZE(sd_numvalid,n,"SD_NUMVALID");
	W_MASKARRAY(sd_numvalid);
	sd_k = (nsd_topind/n)-1;
	if (sd_k<=0 || nsd_topind!=n*(sd_k+1))
	  W_RETERROR(1,"SD_TOPIND: Invalid size");
	W_MASKARRAY(sd_topind);
	W_CHKSIZE(sd_topval,nsd_topind,"SD_TOPVAL");
	W_MASKARRAY(sd_topval);
	//printMsgStdout("Point 5");
	if (ain>20) {
	  if (nsd_subind==0 || nsd_subind>m)
	    W_RETERROR(1,"SD_SUBIND: Wrong size");
	  W_MASKARRAY(sd_subind);
	  if (ain==21)
	    sd_subexcl=0;
	}
      }
    } else
      dampfact=0.0;
    /* Return arguments: Default values and check sizes */
    if (aout<5) {
      sd_nrec=0;
      if (aout<4) {
	sd_nupd=0;
	if (aout<3) {
	  sd_dampfact=0;
	  if (aout<2) {
	    delta=0;
	    if (aout==0)
	      rstat=0;
	  }
	}
      }
    }
    if (aout>0) {
      W_CHKSIZE(rstat,nupdjind,"RSTAT");
      if (aout>1) {
	W_CHKSIZE(delta,nupdjind,"DELTA");
	if (aout>2) {
	  if (sd_k==0)
	    W_RETERROR(1,"Cannot return SD_XXX");
	  W_CHKSIZE(sd_dampfact,nupdjind,"SD_DAMPFACT");
	}
      }
    }
    /* Create max_pi data structure (only if selective damping) */
    Handle<FactEPMaximumPiValues> epMaxPi;
    //printMsgStdout("Point 6");
    if (sd_k>0) {
      try {
	//sprintf(errMsg,"MEX: n=%d,K=%d,numvalid=%d,topind=%d,topval=%d",numN,
	//      sd_k,sd_numvalid.size(),sd_topind.size(),sd_topval.size());
	//printMsgStdout(errMsg);
	epMaxPi.changeRep(new FactEPMaximumPiValues(epRepr,sd_k,
						    sd_numvalidA,sd_topindA,
						    sd_topvalA,sd_subindA,
						    sd_subexcl));
      } catch (StandardException ex) {
	W_RETERROR_ARGS(1,"Cannot create FactEPMaximumPiValues (selective damping):\n%s",ex.msg());
      } catch (...) {
	W_RETERROR(1,"Cannot create FactEPMaximumPiValues (selective damping): Unspecified exception");
      }
    }
    /* Create EP driver */
    Handle<FactorizedEPDriver> epDriver;
    //printMsgStdout("Point 7");
    try {
      epDriver.changeRep(new FactorizedEPDriver(potMan,epRepr,margbetaA,
						margpiA,piminthres,epMaxPi));
    } catch (StandardException ex) {
      W_RETERROR_ARGS(1,"Cannot create FactorizedEPDriver:\n%s",ex.msg());
    } catch (...) {
      W_RETERROR(1,"Cannot create FactorizedEPDriver: Unspecified exception");
    }

    /* Main loop over updates */
    for (int i=0; i<nupdjind; i++) {
      int j=updjind[i];
      //sprintf(errstr,"i=%d, j=%d",i,j);
      //printMsgStdout(errstr);
      int irstat=epDriver->sequentialUpdate(j,dampfact,(delta!=0)?(delta+i):0,
					    (sd_dampfact!=0)?(sd_dampfact+i):0);
      if (rstat!=0) rstat[i]=irstat;
      if (irstat!=FactorizedEPDriver::updSuccess && delta!=0)
	delta[i]=0.0;
      if (irstat!=FactorizedEPDriver::updSuccess && sd_dampfact!=0)
	sd_dampfact[i]=1.0;
    }
    //printMsgStdout("Point 9");
    if (sd_nupd!=0) {
      int inrec;
      epMaxPi->getStats(*sd_nupd,inrec);
      if (sd_nrec!=0) *sd_nrec=inrec;
    }
    //printMsgStdout("Point 10");
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s", ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
