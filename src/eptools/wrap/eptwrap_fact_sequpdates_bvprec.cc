/* -------------------------------------------------------------------
 * EPTWRAP_FACT_SEQUPDATES_BVPREC
 *
 * Same as EPTWRAP_FACT_SEQUPDATES, but in the presence of bivariate
 * precision potentials. There can be standard univariate potentials,
 * followed by >=1 precision potentials.
 *
 * There are K precision variables tau_k. The final M_prec of M
 * potentials are bivariate precision, where 1 <= M_prec <= M.
 * Representation extended by RP_TAUIND (mapping j -> k), RP_A, RP_C
 * for [a_jk], [c_jk], see 'FactorizedEPRepresentation'.
 * Marginals [a_k], [c_k] in MARGA, MARGC.
 * {A|C}MINTHRES same role as PIMINTHRES, for a|c instead of pi.
 *
 * Selective damping: Mechanism for pi is extended to a|c. SD_NUMVALID
 * applies for all. SD{A|C}_xxx plays role of SD_xxx.
 * ATTENTION: SD_SUBIND, SD_SUBEXCL are for SD w.r.t. pi only, the a|c
 * mechanism runs over all precision potentials.
 * The return values SD_Nxxx are sums over all SD mechanisms.
 *
 * Input:
 * - N:            Number of variables x_i
 * - M:            Number of factors
 * - UPDJIND:      Update on these potentials, in order [int32 array]
 * - PM_POTIDS:    Potential manager [int32 array]
 * - PM_NUMPOT:    " [int32 array]
 * - PM_PARVEC:    " [double array]
 * - PM_PARSHRD:   " [int32 array]
 * - PM_ANNOBJ:    " [void* array]
 * - RP_ROWIND:    Factorized EP representation [int32 array]
 * - RP_COLIND:    " [int32 array]
 * - RP_BVALS:     " [double array]
 * - RP_PI:        " [double array; I/O]
 * - RP_BETA:      " [double array; I/O]
 * - RP_TAUIND     " [int32 array]
 * - RP_A:         " [double array; I/O]
 * - RP_C:         " [double array; I/O]
 * - MARGPI:       Variable marginals [I/O]
 * - MARGBETA:     " [I/O]
 * - MARGA:        " [I/O]
 * - MARGC:        " [I/O]
 * - PIMINTHRES:   See above. Positive
 * - AMINTHRES:    "
 * - CMINTHRES:    "
 * - DAMPFACT:     Damping factor, in [0,1). Optional, def. is 0
 * - SD_NUMVALID:  Selective damping. Optional [int32 array; I/O]
 * - SD_TOPIND:    " [int32 array; I/O]
 * - SD_TOPVAL:    " [double array; I/O]
 * - SDA_NUMVALID: " [int32 array; I/O]
 * - SDA_TOPIND:   " [int32 array; I/O]
 * - SDA_TOPVAL:   " [double array; I/O]
 * - SDC_NUMVALID: " [int32 array; I/O]
 * - SDC_TOPIND:   " [int32 array; I/O]
 * - SDC_TOPVAL:   " [double array; I/O]
 * - SD_SUBIND:    " [int32 array]
 * - SD_SUBEXCL:   ". Def.: false
 *
 * Return:
 * - RSTAT:        Return stati for each update. Optional [int32]
 * - DELTA:        See above. Optional
 * - SD_DAMPFACT:  See above. Optional, only if selective damping
 * - SD_NUPD:      " [int32]
 * - SD_NREC:      " [int32]
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_fact_sequpdates_bvprec.h"
#include "src/eptools/FactorizedEPDriver.h"
#include "src/eptools/FactEPMaximumPiValues.h"
#include "src/eptools/FactEPMaximumAValues.h"
#include "src/eptools/FactEPMaximumCValues.h"

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
				    W_DARRAY(sd_topval),W_IARRAY(sda_numvalid),
				    W_IARRAY(sda_topind),W_DARRAY(sda_topval),
				    W_IARRAY(sdc_numvalid),W_IARRAY(sdc_topind),
				    W_DARRAY(sdc_topval),W_IARRAY(sd_subind),
				    int sd_subexcl,W_IARRAY(rstat),
				    W_DARRAY(delta),W_DARRAY(sd_dampfact),
				    int* sd_nupd,int* sd_nrec,W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain<23 || ain>35)
      W_RETERROR(2,"Wrong number of input arguments");
    if (aout>5)
      W_RETERROR(2,"Too many return arguments");
    if (n<1) W_RETERROR(1,"N must be positive");
    if (m<1) W_RETERROR(1,"M must be positive");
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
    createFactEPRepres_bvprec(n,m,W_ARR(rp_rowind),W_ARR(rp_colind),
			      W_ARR(rp_bvals),W_ARR(rp_pi),W_ARR(rp_beta),
			      W_ARR(rp_tauind),W_ARR(rp_a),W_ARR(rp_c),epRepr,
			      W_ERRARGS);
    int numk=epRepr->numPrecVariables();
    if (numk==0)
      W_RETERROR(1,"Must have bivariate precision potentials");
    /* Variable marginals */
    ArrayHandle<double> margpiA,margbetaA,margaA,margcA;
    //printMsgStdout("Point 3");
    W_CHKSIZE(margpi,n,"MARGPI");
    W_CHKSIZE(margpi,n,"MARGBETA");
    W_CHKSIZE(marga,numk,"MARGA");
    W_CHKSIZE(margc,numk,"MARGC");
    W_MASKARRAY(margpi);
    W_MASKARRAY(margbeta);
    W_MASKARRAY(marga);
    W_MASKARRAY(margc);
    if (piminthres<=0.0)
      W_RETERROR(1,"PIMINTHRES must be positive");
    if (aminthres<=0.0)
      W_RETERROR(1,"AMINTHRES must be positive");
    if (cminthres<=0.0)
      W_RETERROR(1,"CMINTHRES must be positive");
    int sd_k=0,sda_k=0,sdc_k=0; // K for max data structure (0: not active)
    ArrayHandle<int> sd_numvalidA,sd_topindA,sd_subindA;
    ArrayHandle<int> sda_numvalidA,sda_topindA;
    ArrayHandle<int> sdc_numvalidA,sdc_topindA;
    ArrayHandle<double> sd_topvalA,sda_topvalA,sdc_topvalA;
    if (ain>23) {
      if (dampfact<0.0 || dampfact>=1.0)
	W_RETERROR(1,"DAMPFACT: Out of range");
      if (ain>24) {
	// Selective damping
	//printMsgStdout("Point 4");
	if (ain<27)
	  W_RETERROR(1,"Need all SD_TOPxxx");
	W_CHKSIZE(sd_numvalid,n,"SD_NUMVALID");
	W_MASKARRAY(sd_numvalid);
	sd_k = (nsd_topind/n)-1;
	if (sd_k<=0 || nsd_topind!=n*(sd_k+1))
	  W_RETERROR(1,"SD_TOPIND: Invalid size");
	W_MASKARRAY(sd_topind);
	W_CHKSIZE(sd_topval,nsd_topind,"SD_TOPVAL");
	W_MASKARRAY(sd_topval);
	if (ain>27) {
	  if (ain<30)
	    W_RETERROR(1,"Need all SDA_xxx");
	  W_CHKSIZE(sda_numvalid,numk,"SDA_NUMVALID");
	  W_MASKARRAY(sda_numvalid);
	  sda_k = (nsda_topind/numk)-1;
	  if (sda_k<=0 || nsda_topind!=numk*(sda_k+1))
	    W_RETERROR(1,"SDA_TOPIND: Invalid size");
	  W_CHKSIZE(sda_topval,nsda_topind,"SDA_TOPVAL");
	  W_MASKARRAY(sda_topind);
	  W_MASKARRAY(sda_topval);
	  if (ain>30) {
	    if (ain<33)
	      W_RETERROR(1,"Need all SDC_xxx");
	    W_CHKSIZE(sdc_numvalid,numk,"SDC_NUMVALID");
	    W_MASKARRAY(sdc_numvalid);
	    sdc_k = (nsdc_topind/numk)-1;
	    if (sdc_k<=0 || nsdc_topind!=numk*(sdc_k+1))
	      W_RETERROR(1,"SDC_TOPIND: Invalid size");
	    W_CHKSIZE(sdc_topval,nsdc_topind,"SDC_TOPVAL");
	    W_MASKARRAY(sdc_topind);
	    W_MASKARRAY(sdc_topval);
	    if (ain>33) {
	      if (nsd_subind==0 || nsd_subind>m)
		W_RETERROR(1,"SD_SUBIND: Wrong size");
	      W_MASKARRAY(sd_subind);
	      if (ain==34)
		sd_subexcl=0;
	    }
	  }
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
    /* Selective damping: Create max data structures */
    Handle<FactEPMaximumPiValues> epMaxPi;
    Handle<FactEPMaximumAValues> epMaxA;
    Handle<FactEPMaximumCValues> epMaxC;
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
    if (sda_k>0) {
      try {
	epMaxA.changeRep(new FactEPMaximumAValues(epRepr,sda_k,sda_numvalidA,
						  sda_topindA,sda_topvalA));
      } catch (StandardException ex) {
	W_RETERROR_ARGS(1,"Cannot create FactEPMaximumAValues (selective damping):\n%s",ex.msg());
      } catch (...) {
	W_RETERROR(1,"Cannot create FactEPMaximumAValues (selective damping): Unspecified exception");
      }
    }
    if (sdc_k>0) {
      try {
	epMaxC.changeRep(new FactEPMaximumCValues(epRepr,sdc_k,sdc_numvalidA,
						  sdc_topindA,sdc_topvalA));
      } catch (StandardException ex) {
	W_RETERROR_ARGS(1,"Cannot create FactEPMaximumCValues (selective damping):\n%s",ex.msg());
      } catch (...) {
	W_RETERROR(1,"Cannot create FactEPMaximumCValues (selective damping): Unspecified exception");
      }
    }
    /* Create EP driver */
    Handle<FactorizedEPDriver> epDriver;
    //printMsgStdout("Point 7");
    try {
      epDriver.changeRep(new FactorizedEPDriver(potMan,epRepr,margbetaA,
						margpiA,margaA,margcA,
						piminthres,aminthres,cminthres,
						epMaxPi,epMaxA,epMaxC));
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
      int nupd=0,nrec=0,tu,tr;
      epMaxPi->getStats(tu,tr); nupd+=tu; nrec+=tr;
      if (sda_k>0) {
	epMaxA->getStats(tu,tr); nupd+=tu; nrec+=tr;
      }
      if (sdc_k>0) {
	epMaxC->getStats(tu,tr); nupd+=tu; nrec+=tr;
      }
      *sd_nupd=nupd;
      if (sd_nrec!=0) *sd_nrec=nrec;
    }
    //printMsgStdout("Point 10");
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s", ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
