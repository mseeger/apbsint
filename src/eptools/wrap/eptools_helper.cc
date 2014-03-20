/* -------------------------------------------------------------------
 * Helper functions for EPTOOLS wrapper functions
 * Implementation
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/potentials/PotManagerFactory.h"
#include "src/eptools/FactorizedEPRepresentation.h"

/*
 * Parses arguments POTIDS, NUMPOT, PARVEC, PARSHRD and creates a potential
 * manager.
 */
void createPotentialManager(W_IARRAY(potids),W_IARRAY(numpot),W_DARRAY(parvec),
			    W_IARRAY(parshrd),W_ARRAY(annobj,void*),
			    Handle<PotentialManager>& potMan,W_ERRORARGS)
{
  ArrayHandle<int> potidsA,numpotA,parshrdA;
  ArrayHandle<double> parvecA;
  ArrayHandle<void*> annobjA;

  W_CHKSIZE(numpot,npotids,"NUMPOT");
  W_MASKARRAY(potids);
  W_MASKARRAY(numpot);
  W_MASKARRAY(parvec);
  W_MASKARRAY(parshrd);
  W_MASKARRAY(annobj);
  try {
    potMan.changeRep(PotManagerFactory::create(potidsA,numpotA,parvecA,
					       parshrdA,annobjA));
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Cannot create potential manager:\n%s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Cannot create potential manager: Unspecified exception");
  }
}

/*
 * Creates 'FactorizedEPRepresentation' for a model with standard univariate
 * potentials only (argument group 'atypeUnivariate').
 */
void createFactEPRepres(int numN,int numM,W_IARRAY(rp_rowind),
			W_IARRAY(rp_colind),W_DARRAY(rp_bvals),W_DARRAY(rp_pi),
			W_DARRAY(rp_beta),
			Handle<FactorizedEPRepresentation>& epRepr,W_ERRORARGS)
{
  ArrayHandle<int> rp_rowindA,rp_colindA;
  ArrayHandle<double> rp_bvalsA,rp_piA,rp_betaA;

  W_CHKSIZE(rp_pi,nrp_bvals,"RP_PI");
  W_CHKSIZE(rp_beta,nrp_bvals,"RP_BETA");
  W_MASKARRAY(rp_rowind);
  W_MASKARRAY(rp_colind);
  W_MASKARRAY(rp_bvals);
  W_MASKARRAY(rp_pi);
  W_MASKARRAY(rp_beta);
  try {
    epRepr.changeRep(new FactorizedEPRepresentation(numN,numM,rp_rowindA,
						    rp_colindA,rp_bvalsA,
						    rp_betaA,rp_piA));
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Cannot create B representation:\n%s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Cannot create B representation: Unspecified exception");
  }
}

/*
 * Creates 'FactorizedEPRepresentation' for a model with some bivariate
 * precision potentials (argument group 'atypeBivarPrec').
 */
void createFactEPRepres_bvprec(int numN,int numM,W_IARRAY(rp_rowind),
			       W_IARRAY(rp_colind),W_DARRAY(rp_bvals),
			       W_DARRAY(rp_pi),W_DARRAY(rp_beta),
			       W_IARRAY(rp_tauind),W_DARRAY(rp_a),
			       W_DARRAY(rp_c),
			       Handle<FactorizedEPRepresentation>& epRepr,
			       W_ERRORARGS)
{
  ArrayHandle<int> rp_rowindA,rp_colindA,rp_tauindA;
  ArrayHandle<double> rp_bvalsA,rp_piA,rp_betaA,rp_aA,rp_cA;

  W_CHKSIZE(rp_pi,nrp_bvals,"RP_PI");
  W_CHKSIZE(rp_beta,nrp_bvals,"RP_BETA");
  W_CHKSIZE(rp_c,nrp_a,"RP_C");
  W_MASKARRAY(rp_rowind);
  W_MASKARRAY(rp_colind);
  W_MASKARRAY(rp_bvals);
  W_MASKARRAY(rp_pi);
  W_MASKARRAY(rp_beta);
  W_MASKARRAY(rp_tauind);
  W_MASKARRAY(rp_a);
  W_MASKARRAY(rp_c);
  try {
    epRepr.changeRep(new FactorizedEPRepresentation(numN,numM,rp_rowindA,
						    rp_colindA,rp_bvalsA,
						    rp_betaA,rp_piA,rp_aA,
						    rp_cA,rp_tauindA));
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Cannot create B representation:\n%s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Cannot create B representation: Unspecified exception");
  }
}
