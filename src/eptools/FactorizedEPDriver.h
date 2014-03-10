/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class FactorizedEPDriver
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_FACTORIZEDEPDRIVER_H
#define EPTOOLS_FACTORIZEDEPDRIVER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/PotentialManager.h"
#include "src/eptools/FactEPRepresBivarPrec.h"
#include "src/eptools/FactEPMaximumPiValues.h"
#include "src/eptools/FactEPMaximumAValues.h"
#include "src/eptools/FactEPMaximumCValues.h"

//BEGINNS(eptools)
#define MAXRELDIFF(a,b) (fabs((a)-(b))/std::max(fabs(a),std::max(fabs(b),1e-8)))

  /**
   * Driver for expectation propagation with factorized backbone. Two
   * different cases are supported here (selected at construction):
   * - Univariate potentials, inference over x. All potentials must be
   *   in argument group 'atypeUnivariate'. Representation of type
   *   'FactorizedEPRepresentation'
   * - Bivariate precision potentials, inference over x and tau. All
   *   potentials must be in argument group 'atypeBivarPrec'. Representation
   *   of type 'FactEPRepresBivarPrec'
   * <p>
   * Univariate potentials:
   *   p(x) = Z^-1 prod_j t_j(s_j),  s = B x, x[0:n-1], s[0:m-1]
   * {t_j(s_j)} given by 'epPots' (argument group 'atypeUnivariate'). B,
   * message parameters beta, pi maintained in 'epRepr'
   * ('FactorizedEPRepresentation').
   * Marginal moments represented by [beta_i], [pi_i] ('margBeta',
   * 'margPi').
   * Relationship between marginals and message parameters:
   *   pi_i = sum_j pi_ji,   beta_i = sum_j beta_ji
   * NOTE: A model may consist of several drivers and representations. In
   * particular, if it contains both univariate and bivariate potentials, two
   * drivers are involved. pi_i, beta_i are sums over message parameters of
   * all parts.
   * <p>
   * Bivariate potentials over precision parameters:
   *   p(x,tau) = Z^-1 prod_j t_j(s_j,tau_k(j)),  s = B x,
   *              x[0:n-1], s[0:m-1], tau[0:K-1]
   * {t_j(.)} given by 'epPots' (argument group 'atypeBivarPrec'). B, [k(j)],
   * message parameters beta, pi, a, c maintained in 'epRepr'
   * ('FactEPRepresBivarPrec').
   * Marginal moments represented by [beta_i], [pi_i], [a_k], [c_k]
   * ('margBeta', 'margPi', 'margA', 'margC').
   * Relationship between marginals and message parameters:
   *   a_k = sum_j a_jk,   c_k = sum_j c_jk
   * NOTE: Even if the model consists of several drivers and representations,
   * only one of them can maintain potentials over tau.
   * <p>
   * Sequential EP update:
   * Main service 'sequentialUpdate'. Involves computing cavity marginals,
   * doing local EP update on t_j(.), computing new EP (message) parameters,
   * apply damping (optional), update marginals. The typical reaction to
   * things going wrong is to skip the update (but see "Selective damping"
   * below). Return status:
   * - updSuccess: All fine
   * - updCavityInvalid: Require pi_{-ji} >= eps/2 for all i in V_j, where
   *   eps=='piMinThres' (passed at construction)
   *   Also: a_{-jk} >= 0.5*'aMinThres', c_{-jk} >= 0.5*'cMinThres'
   * - updNumericalError: 'EPScalarPotential::compMoments' returns false
   * - updMarginalsInvalid: Require pi_i >= eps/2 for all i in V_j, for the
   *   new marginals after the (damped) update.
   *   Also: a_k >= 0.5*'aMinThres', c_k >= 0.5*'cMinThres'.
   * <p>
   * Selective damping/skipping:
   * This is done iff the corresponding 'MaximumValuesService' objects are
   * given ('epMaxPi' for pi, 'epMaxA' for a, 'epMaxC' for c). We ensure that
   * after the update:
   *   pi_i - max_k pi_ki >= eps,  pi_i >=eps    for all i, [*]
   * given that this condition holds before the update as well. Here,
   * eps=='piMinThres'. The 'epMaxPi' object maintains max_k pi_ki for all
   * i, it is updated here as well.
   * Details given in the note. We first an update with damping factor
   * 'dampFact'. If this violates [*], we determine the smallest damping
   * factor for which [*] is true. This may be 1, in which case the update
   * is skipped.
   * Same for a's (c's) with 'aMinThres' ('cMinThres') respectively. We
   * use the smallest damping factor s.t. all constraints are fulfilled.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactorizedEPDriver
  {
  public:
    // Constants

    static const int updSuccess         =0;
    static const int updCavityInvalid   =1;
    static const int updNumericalError  =2;
    static const int updMarginalsInvalid=3;
    static const int updCavCondSkipped  =4;

  protected:
    // Members

    Handle<PotentialManager> epPots;       // Potential manager
    Handle<FactorizedEPRepresentation> epRepr;
    FactEPRepresBivarPrec* epReprPrec;     // Only for bivar. prec. pots.
    ArrayHandle<double> margBeta,margPi;
    double piMinThres;
    Handle<FactEPMaximumPiValues> epMaxPi; // Selective damping (optional)
    ArrayHandle<double> margA,margC;       // Only for bivar. prec. pots.
    double aMinThres,cMinThres;
    Handle<FactEPMaximumAValues> epMaxA;
    Handle<FactEPMaximumCValues> epMaxC;
    ArrayHandle<double> buffVec;

  public:
    // Public methods

    /**
     * Constructor (univariate potentials). All potentials in 'pepPots' must
     * be in argument group 'EPScalarPotential::atypeUnivariate'.
     *
     * @param pepPots
     * @param pepRepr
     * @param pmargBeta
     * @param pmargPi
     * @param ppiMinThres
     * @param pepMaxPi    Optional
     */
    FactorizedEPDriver(const Handle<PotentialManager>& pepPots,
		       const Handle<FactorizedEPRepresentation>& pepRepr,
		       const ArrayHandle<double>& pmargBeta,
		       const ArrayHandle<double>& pmargPi,double ppiMinThres,
		       const Handle<FactEPMaximumPiValues>& pepMaxPi=
		       HandleZero<FactEPMaximumPiValues>::get()) :
      epPots(pepPots),epRepr(pepRepr),epReprPrec(0),margBeta(pmargBeta),
      margPi(pmargPi),piMinThres(ppiMinThres),epMaxPi(pepMaxPi) {
      int numN=pepRepr->numVariables();

      if (ppiMinThres<=0.0 || pmargBeta.size()!=numN || pmargPi.size()!=numN)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (pepPots->size()!=
	  pepPots->numArgumentGroup(EPScalarPotential::atypeUnivariate))
	throw InvalidParameterException(EXCEPT_MSG("Potentials must be in group 'atypeUnivariate'"));
    }

    /**
     * Constructor (bivariate precision potentials). All potentials in
     * 'pepPots' must be in argument group
     * 'EPScalarPotential::atypeBivarPrec'.
     *
     * @param pepPots
     * @param pepRepr
     * @param pmargBeta
     * @param pmargPi
     * @param pmargA
     * @param pmargC
     * @param ppiMinThres
     * @param paMinThres
     * @param pcMinThres
     * @param pepMaxPi    Optional
     * @param pepMaxA     "
     * @param pepMaxC     "
     */
    FactorizedEPDriver(const Handle<PotentialManager>& pepPots,
		       const Handle<FactEPRepresBivarPrec>& pepRepr,
		       const ArrayHandle<double>& pmargBeta,
		       const ArrayHandle<double>& pmargPi,
		       const ArrayHandle<double>& pmargA,
		       const ArrayHandle<double>& pmargC,double ppiMinThres,
		       double paMinThres,double pcMinThres,
		       const Handle<FactEPMaximumPiValues>& pepMaxPi=
		       HandleZero<FactEPMaximumPiValues>::get(),
		       const Handle<FactEPMaximumAValues>& pepMaxA=
		       HandleZero<FactEPMaximumAValues>::get(),
		       const Handle<FactEPMaximumCValues>& pepMaxC=
		       HandleZero<FactEPMaximumCValues>::get()) :
      epPots(pepPots),epRepr(pepRepr),epReprPrec(pepRepr.p()),
      margBeta(pmargBeta),margPi(pmargPi),margA(pmargA),margC(pmargC),
      piMinThres(ppiMinThres),aMinThres(paMinThres),cMinThres(pcMinThres),
      epMaxPi(pepMaxPi),epMaxA(pepMaxA),epMaxC(pepMaxC) {
      int numN=pepRepr->numVariables(),numK=pepRepr->numPrecVars();

      if (ppiMinThres<=0.0 || pmargBeta.size()!=numN || pmargPi.size()!=numN ||
	  paMinThres<=0.0 || pcMinThres<=0.0 || pmargA.size()!=numK ||
	  pmargC.size()!=numK)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (pepPots->size()!=
	  pepPots->numArgumentGroup(EPScalarPotential::atypeBivarPrec))
	throw InvalidParameterException(EXCEPT_MSG("Potentials must be in group 'atypeBivarPrec'"));
    }

    virtual ~FactorizedEPDriver() {}

    virtual int numVariables() const {
      return epRepr->numVariables();
    }

    virtual int numPotentials() const {
      return epRepr->numPotentials();
    }

    virtual int numPrecVars() const {
      if (epReprPrec==0)
	throw WrongStatusException(EXCEPT_MSG(""));
      return epReprPrec->numPrecVars();
    }

    virtual const PotentialManager& getEPPotentials() const {
      return *epPots;
    }

    virtual const ArrayHandle<double>& getMarginalsBeta() const {
      return margBeta;
    }

    virtual const ArrayHandle<double>& getMarginalsPi() const {
      return margPi;
    }

    virtual const ArrayHandle<double>& getMarginalsA() const {
      if (epReprPrec==0)
	throw WrongStatusException(EXCEPT_MSG(""));
      return margA;
    }

    virtual const ArrayHandle<double>& getMarginalsC() const {
      if (epReprPrec==0)
	throw WrongStatusException(EXCEPT_MSG(""));
      return margC;
    }

    /**
     * Runs sequential EP update on potential t_j(.). See header comment.
     * If selective damping is active ('epMaxXXX'!=0), the effective damping
     * factor can be returned via 'effDamp'. If equal to 'dampFact', no
     * selective damping was necessary. In the worst case, the update is
     * skipped (ret. status 'updCavCondSkipped').
     * In 'delta', we return the maximum of relative change in mean and
     * stddev. on s_j (not on x), or on s_j and tau_k(j).
     *
     * @param j        Potential index to update on
     * @param dampFact Damping factor in [0,1). Def.: 0 (no damping)
     * @param delta    S.a. Optional
     * @param effDamp  S.a. Optional
     * @return         Return status ('updSuccess' for success)
     */
    virtual int sequentialUpdate(int j,double dampFact=0.0,double* delta=0,
				 double* effDamp=0);
  };

  // Inline methods

  /*
   * Arrays: XX is 'beta', 'pi'
   * - vjInd:  V_j
   * - bP:     b_ji
   * - XXP:    XX_ji, EP parameters (overwritten only at end)
   * - mXXP:   XX_i, marginals (overwritten only at end)
   * - cXXP:   First cavity, then updated EP pars
   * - mprXXP: First updated EP pars (without damping), then
   *           new XX_i, marginals
   * Required, because an update can be skipped until the very end.
   */
  inline int FactorizedEPDriver::sequentialUpdate(int j,double dampFact,
						  double* delta,double* effDamp)
  {
    int i,ii,vjSz,k;
    double temp,temp2,cH,cRho,bval,nu,alpha,cPi,cBeta,pi,beta,tilPi,tilBeta,
      prPi,prBeta,kappa,thres2=0.5*piMinThres,mH,mRho,cA,cC,hatA,hatC,prA,prC;
    const int* vjInd;
    const double* bP;
    double* betaP,*piP,*cBetaP,*cPiP,*mBetaP,*mPiP,*mprBetaP,*mprPiP,*aP,*cP;
    double inp[4],ret[4];
    bool hasBVPrec=(epReprPrec!=0);
    char debMsg[200]; // DEBUG!

    MYASS((!hasBVPrec && epPots->getPot(j).getArgumentGroup()==
	   EPScalarPotential::atypeUnivariate) ||
	  (hasBVPrec && epPots->getPot(j).getArgumentGroup()==
	   EPScalarPotential::atypeBivarPrec));
    if (dampFact<0.0 || dampFact>=1.0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    // Access to data for j. Temporary arrays
    //printMsgStdout("epRepr->accessRow");
    epRepr->accessRow(j,vjSz,vjInd,bP,betaP,piP);
    //printMsgStdout("OK");
    if (hasBVPrec)
      k=epReprPrec->accessTauRow(j,aP,cP); // k=k(j)
    mBetaP=margBeta.p(); mPiP=margPi.p();
    if (buffVec.size()<4*vjSz) {
      //sprintf(debMsg,"buffVec.changeRep: %d",4*vjSz);
      //printMsgStdout(debMsg);
      buffVec.changeRep(4*vjSz);
      //printMsgStdout("OK");
    }
    cBetaP=buffVec.p(); cPiP=cBetaP+vjSz;
    mprBetaP=cPiP+vjSz; mprPiP=mprBetaP+vjSz;
    // Compute cavity marginals
    // The marginal moments on s_j ('mH', 'mRho') are required to compute
    // '*delta' below
    //printMsgStdout("sequentialUpdate: Compute cavities"); // DEBUG!
    cH=cRho=mH=mRho=0.0;
    for (ii=0; ii<vjSz; ii++) {
      i=vjInd[ii];
      //sprintf(debMsg,"ii=%d,i=%d,mPi=%f,mBeta=%f,pi=%f,beta=%f",ii,i,
      //      mPiP[i],mBetaP[i],piP[ii],betaP[ii]); // DEBUG!
      //printMsgStdout(debMsg);
      if ((cPiP[ii]=cPi=mPiP[i]-piP[ii])<thres2)
	return updCavityInvalid; // EP update failed
      cBetaP[ii]=cBeta=mBetaP[i]-betaP[ii];
      bval=bP[ii]; temp=bval/cPi;
      cRho+=bval*temp;
      cH+=temp*cBeta;
      temp=bval/mPiP[i];
      mRho+=bval*temp;
      mH+=temp*mBetaP[i];
    }
    if (hasBVPrec) {
      if ((cA=margA[k]-(*aP))<0.5*aMinThres)
	return updCavityInvalid; // EP update failed
      if ((cC=margC[k]-(*cP))<0.5*cMinThres)
	return updCavityInvalid; // EP update failed
    }
    //sprintf(debMsg,"Cav. marg.: cH=%f,cRho=%f",cH,cRho); // DEBUG!
    //printMsgStdout(debMsg);
    // Local EP update
    inp[0]=cH; inp[1]=cRho;
    if (hasBVPrec) {
      inp[2]=cA; inp[3]=cC;
    }
    if (!epPots->getPot(j).compMoments(inp,ret)) {
      // DEBUG:
      if (!hasBVPrec)
	sprintf(debMsg,"UUPS: j=%d, cH=%f,cRho=%f",j,cH,cRho);
      else
	sprintf(debMsg,"UUPS: j=%d, cH=%f,cRho=%f,cA=%f,cC=%f",j,cH,cRho,cA,
		cC);
      printMsgStdout(debMsg);
      return updNumericalError; // EP update failed
    }
    alpha=ret[0]; nu=ret[1];
    if (hasBVPrec) {
      // New marginal a, c parameters (without damping)
      hatA=ret[2]; hatC=ret[3];
    }
    // Compute new EP parameters without damping (to 'mprXXP'). If
    // selective damping is active, we also determine the effective damping
    // factor (overwrites 'dampFact')
    //printMsgStdout("sequentialUpdate: 2"); // DEBUG!
    for (ii=0; ii<vjSz; ii++) {
      // Undamped EP update
      i=vjInd[ii];
      bval=bP[ii]; // b_{ji}
      pi=piP[ii]; beta=betaP[ii]; // pi_{ji}, beta_{ji}
      cPi=cPiP[ii]; cBeta=cBetaP[ii]; // pi_{-ji}, beta_{-ji}
      //sprintf(debMsg,"ii=%d,i=%d,bval=%f,cPi=%f,cBeta=%f",ii,i,bval,cPi,
      //      cBeta); // DEBUG!
      //printMsgStdout(debMsg);
      // 'tilPi', 'tilBeta': tilde{pi}_{ji}, tilde{beta}_{ji}, EP updates
      // without damping
      if (fabs(bval)>1e-6) {
	// |b_ji| large enough: Simpler equations
	// 'temp2' is pi_{-ji}/b_ji. 'temp' is (pi_ji)'/(pi_{-ji} nu_j)
	temp2=cPi/bval;
	if ((temp=temp2/bval-nu)<1e-10) {
	  // DEBUG
	  sprintf(debMsg,
		  "UUPS: j=%d, cH=%f, cRho=%f, alpha=%f, nu=%f\n"
		  "      b=%f, denom=%f",j,cH,cRho,alpha,nu,bval,temp);
	  printMsgStdout(debMsg);
	  return updNumericalError; // EP update failed
	}
	temp=1.0/temp; // e_ji
	tilPi=temp*cPi*nu;
	tilBeta=temp*(cBeta*nu+temp2*alpha);
      } else {
	// Very small but non-zero |b_ji| (will probably never happen)
	// 'temp' is b_ji / (pi_{-ji} - b_ji^2 nu_j)
	if ((temp=cPi-nu*bval*bval)<1e-10)
	  return updNumericalError; // EP update failed
	temp=bval/temp;
	tilPi=temp*bval*nu*cPi;
	tilBeta=temp*(cBeta*bval*nu+cPi*alpha);
      }
      mprPiP[ii]=tilPi; mprBetaP[ii]=tilBeta; // Intermed. storage
      if (!(epMaxPi==0) && tilPi<pi) {
	// Selective damping to ensure that pi_{-ki} >= eps for all k,i
	kappa=epMaxPi->getMaxValue(i); // kappa_i
	if (kappa<=0.0) {
	  sprintf(debMsg,"ERROR(maxPi,j=%d,i=%d): kappa_i=%f (negative)",j,i,
		  kappa);
	  printMsgStdout(debMsg);
	  return updNumericalError;
	}
	// Value for eta:
	eta=1.0-std::min((mPiP[i]-kappa-piMinThres)/(pi-tilPi),1.0);
	if (eta>=0.98) {
	  // EP update has to be skipped
	  if (effDamp!=0) *effDamp=1.0;
	  return updCavCondSkipped;
	}
	if (kappa==pi) {
	  // This should not happen often. Have to ensure that new kappa_i is
	  // positive. If this is not the case, the update is skipped.
	  // ATTENTION: If this case happens frequently, have to choose better
	  // response, f.ex. increasing 'eta' in small steps.
	  prPi=eta*pi+(1.0-eta)*tilPi; // pi_{ji}' for current 'eta'
	  piP[ii]=prPi;
	  epMaxPi->update(i,j,prPi);
	  kappa=epMaxPi->getMaxValue(i); // kappa_i'
	  piP[ii]=pi; // Back to old state
	  epMaxPi->update(i,j,pi);
	  if (kappa<=0.0) {
	    // Assuming this case almost never happens, we just skip the update
	    printMsgStdout("UUPS(pi selective damping; skipping update due to negative kappa)");
	    if (effDamp!=0) *effDamp=1.0;
	    return updCavCondSkipped;
	  }
	}
	dampFact=std::max(dampFact,eta);
      }
    }
    if (hasBVPrec) {
      // Selective damping
      prA=hatA-cA; prC=hatC-cC;
      if (!(epMaxA==0) && prA<*aP) {
	// Selective damping to ensure that a_{-jk} >= 'aMinThres' for all j,k
	kappa=epMaxA->getMaxValue(k); // kappa_k
	if (kappa<=0.0) {
	  sprintf(debMsg,"ERROR(maxA,j=%d,k=%d): kappa_k=%f (negative)",j,k,
		  kappa);
	  printMsgStdout(debMsg);
	  return updNumericalError;
	}
	// Value for eta:
	eta=1.0-std::min((margA[k]-kappa-aMinThres)/(*aP-prA),1.0);
	if (eta>=0.98) {
	  // EP update has to be skipped
	  if (effDamp!=0) *effDamp=1.0;
	  return updCavCondSkipped;
	}
	if (kappa==*aP) {
	  // Ensure that new kappa_k is positive
	  prA+=eta*(*aP-prA); // a_{jk}' for current 'eta'
	  temp=*aP; *aP=prA;
	  epMaxA->update(k,j,prA);
	  kappa=epMaxA->getMaxValue(k); // kappa_k'
	  *aP=temp; // Back to old state
	  epMaxA->update(k,j,temp);
	  if (kappa<=0.0) {
	    // Assuming this case almost never happens, we just skip the update
	    printMsgStdout("UUPS(a selective damping; skipping update due to negative kappa)");
	    if (effDamp!=0) *effDamp=1.0;
	    return updCavCondSkipped;
	  }
	}
	dampFact=std::max(dampFact,eta);
      }
      if (!(epMaxC==0) && prC<*cP) {
	// Selective damping to ensure that c_{-jk} >= 'cMinThres' for all j,k
	kappa=epMaxC->getMaxValue(k); // kappa_k
	if (kappa<=0.0) {
	  sprintf(debMsg,"ERROR(maxC,j=%d,k=%d): kappa_k=%f (negative)",j,k,
		  kappa);
	  printMsgStdout(debMsg);
	  return updNumericalError;
	}
	// Value for eta:
	eta=1.0-std::min((margC[k]-kappa-cMinThres)/(*cP-prC),1.0);
	if (eta>=0.98) {
	  // EP update has to be skipped
	  if (effDamp!=0) *effDamp=1.0;
	  return updCavCondSkipped;
	}
	if (kappa==*cP) {
	  // Ensure that new kappa_k is positive
	  prC+=eta*(*cP-prC); // c_{jk}' for current 'eta'
	  temp=*cP; *cP=prC;
	  epMaxC->update(k,j,prC);
	  kappa=epMaxC->getMaxValue(k); // kappa_k'
	  *cP=temp; // Back to old state
	  epMaxC->update(k,j,temp);
	  if (kappa<=0.0) {
	    // Assuming this case almost never happens, we just skip the update
	    printMsgStdout("UUPS(c selective damping; skipping update due to negative kappa)");
	    if (effDamp!=0) *effDamp=1.0;
	    return updCavCondSkipped;
	  }
	}
	dampFact=std::max(dampFact,eta);
      }
    }
    if (effDamp!=0) *effDamp=dampFact;
    // Determine new EP parameters with damping (overwrite 'cXXP') and
    // new marginals (to 'mprXXP'). This is done because the update can
    // still fail ('updMarginalsInvalid')
    for (ii=0; ii<vjSz; ii++) {
      pi=piP[ii]; beta=betaP[ii]; // Current parameters
      cPi=cPiP[ii]; cBeta=cBetaP[ii]; // Cavity
      prPi=mprPiP[ii]; prBeta=mprBetaP[ii]; // Undamped update
      //sprintf(debMsg,"ii=%d,prPi=%f,prBeta=%f",ii,prPi,prBeta); // DEBUG!
      //printMsgStdout(debMsg);
      // Damping
      if (dampFact>0.0) {
	prPi+=dampFact*(pi-prPi);
	prBeta+=dampFact*(beta-prBeta);
      }
      // New marginals (overwrite 'mprXXP') and EP parameters (overwrite
      // 'cXXP')
      if ((mprPiP[ii]=cPi+prPi)<thres2)
	return updMarginalsInvalid; // EP update failed
      mprBetaP[ii]=cBeta+prBeta;
      cPiP[ii]=prPi; cBetaP[ii]=prBeta; // New EP parameters
    }
    // New EP parameters and marginals for Gamma parameters: Write back
    if (hasBVPrec) {
      prA=hatA-cA; prC=hatC-cC;
      if (dampFact>0.0) {
	prA+=dampFact*(*aP-prA);
	prC+=dampFact*(*cP-prC);
      }
      if (cA+prA<0.5*aMinThres || cC+prC<0.5*cMinThres)
	return updMarginalsInvalid; // EP update failed
      *aP=prA; *cP=prC;
      margA[k]=cA+prA; margC[k]=cC+prC;
      if (!(epMaxA==0))
	epMaxA->update(k,j,prA);
      if (!(epMaxC==0))
	epMaxC->update(k,j,prC);
      // HIER: '*delta'!!
      // Need mean and stddev. of Gamma!
    }
    // Update succeeded: Write back new EP parameters and marginals
    double mprH=0.0,mprRho=0.0; // For '*delta'
    for (ii=0; ii<vjSz; ii++) {
      i=vjInd[ii];
      //sprintf(debMsg,"ii=%d,i=%d,beta'=%f,pi'=%f,mbeta'=%f,mpi'=%f",ii,i,
      //      cBetaP[ii],cPiP[ii],mprBetaP[ii],mprPiP[ii]); // DEBUG!
      //printMsgStdout(debMsg);
      betaP[ii]=cBetaP[ii]; piP[ii]=cPiP[ii]; // New EP pars
      mBetaP[i]=mprBetaP[ii]; mPiP[i]=mprPiP[ii]; // New marginals
      // For '*delta':
      bval=bP[ii]; temp=bval/mprPiP[ii];
      mprRho+=bval*temp;
      mprH+=temp*mprBetaP[ii];
      if (!(epMaxPi==0))
	epMaxPi->update(i,j,piP[ii]); // Update max-pi object
    }
    if (delta!=0) {
      mRho=sqrt(mRho); mprRho=sqrt(mprRho);
      *delta=std::max(MAXRELDIFF(mH,mprH),MAXRELDIFF(mRho,mprRho));
    }

    return updSuccess;
  }

#undef MAXRELDIFF
//ENDNS

#endif
