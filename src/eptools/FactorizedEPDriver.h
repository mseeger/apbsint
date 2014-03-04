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
#include "src/eptools/FactorizedEPRepresentation.h"
#include "src/eptools/FactEPMaximumPiValues.h"

//BEGINNS(eptools)
#define MAXRELDIFF(a,b) (fabs((a)-(b))/std::max(fabs(a),std::max(fabs(b),1e-8)))

  /**
   * Driver for expectation propagation with factorized backbone. All
   * potentials must be in argument group
   * 'EPScalarPotential::atypeUnivariate'.
   * <p>
   * For details, see technical note.
   *   p(x) = Z^-1 prod_j t_j(s_j),  s = B x, x[0:n-1], s[0:m-1]
   * {t_j(s_j)} given by 'epPots'. B m-by-n sparse matrix,
   * representation by contiguous rows (below). EP parameters
   * (factor-to-variable messages) beta, pi same size as B nonzeros.
   * Sequential EP update on one t_j(s_j), updates all marginals
   * and beta_j, pi_j.
   * <p>
   * Sparse factor B and EP (message) parameters maintained in
   * 'epRepr' ('FactorizedEPRepresentation').
   * Marginal moments represented by [beta_i], [pi_i] ('margBeta',
   * 'margPi'). Must have at all times:
   *   pi_i = sum_j pi_ji,   beta_i = sum_j beta_ji
   * NOTE: 'FactorizedEPRepresentation::compMarginals' computes marginals
   * from EP parameters.
   * <p>
   * Sequential EP update:
   * Main service 'sequentialUpdate'. Involves computing cavity marginals,
   * doing local EP update on t_j(s_j), computing new EP (message) parameters,
   * apply damping (optional), update marginals. The typical reaction to
   * things going wrong is to skip the update (but see "Selective damping"
   * below). Return status:
   * - updSuccess: All fine
   * - updCavityInvalid: Require pi_{-ji} >= eps/2 for all i in V_j, where
   *   eps=='piMinThres' (passed at construction)
   * - updNumericalError: 'EPScalarPotential::compMoments' returns false
   * - updMarginalsInvalid: Require pi_i >= eps/2 for all i in V_j, for the
   *   new marginals after the (damped) update
   * <p>
   * Selective damping/skipping:
   * This is done iff 'epMaxPi' ('FactEPMaximumPiValues') is given. We
   * ensure that after the update:
   *   pi_i - max_k pi_ki >= eps,  pi_i >=eps    for all i, [*]
   * given that this condition holds before the update as well. Here,
   * eps=='piMinThres'. The 'epMaxPi' object maintains max_k pi_ki for all
   * i, it is updated here as well.
   * Details given in the note. We first an update with damping factor
   * 'dampFact'. If this violates [*], we determine the smallest damping
   * factor for which [*] is true. This may be 1, in which case the update
   * is skipped.
   *
   * NOTE: None of this should be required for models with log-concave
   * potentials. For other models, have to see how this works.
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
    ArrayHandle<double> margBeta,margPi;
    double piMinThres;
    Handle<FactEPMaximumPiValues> epMaxPi; // Selective damping (optional)
    ArrayHandle<double> buffVec;

  public:
    // Public methods

    /**
     * Constructor. All potentials in 'pepPots' must be in argument group
     * 'EPScalarPotential::atypeUnivariate'.
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
      epPots(pepPots),epRepr(pepRepr),margBeta(pmargBeta),margPi(pmargPi),
      piMinThres(ppiMinThres),epMaxPi(pepMaxPi) {
      int numN=pepRepr->numVariables();

      if (ppiMinThres<=0.0 || pmargBeta.size()!=numN || pmargPi.size()!=numN)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (pepPots->size()!=
	  pepPots->numArgumentGroup(EPScalarPotential::atypeUnivariate))
	throw InvalidParameterException(EXCEPT_MSG("Potentials must be in group 'atypeUnivariate'"));
    }

    virtual ~FactorizedEPDriver() {}

    virtual int numVariables() const {
      return epRepr->numVariables();
    }

    virtual int numPotentials() const {
      return epRepr->numPotentials();
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

    /**
     * Runs sequential EP update on potential t_j(s_j). See header comment.
     * If selective damping is active ('epMaxPi'!=0), the effective damping
     * factor can be returned via 'effDamp'. If equal to 'dampFact', no
     * selective damping was necessary. In the worst case, the update is
     * skipped (ret. status 'updCavCondSkipped').
     * In 'delta', we return the relative change in mean and stddev. on
     * s_j (not on x).
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
    int i,ii,vjSz;
    double temp,temp2,cPi,cBeta,cH,cRho,bval,nu,alpha,prBeta,prPi,
      thres2=0.5*piMinThres,mH,mRho;
    const int* vjInd;
    const double* bP;
    double* betaP,*piP,*cBetaP,*cPiP,*mBetaP,*mPiP,*mprBetaP,*mprPiP;
    double inp[2],ret[2];
    char debMsg[200]; // DEBUG!

    if (dampFact<0.0 || dampFact>=1.0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    // Access to data for j. Temporary arrays
    //printMsgStdout("epRepr->accessRow");
    epRepr->accessRow(j,vjSz,vjInd,bP,betaP,piP);
    //printMsgStdout("OK");
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
    //sprintf(debMsg,"Cav. marg.: cH=%f,cRho=%f",cH,cRho); // DEBUG!
    //printMsgStdout(debMsg);
    // Local EP update
    inp[0]=cH; inp[1]=cRho;
    if (!epPots->getPot(j).compMoments(inp,ret)) {
      sprintf(debMsg,"UUPS: j=%d, cH=%f,cRho=%f",j,cH,cRho); // DEBUG
      printMsgStdout(debMsg);
      return updNumericalError; // EP update failed
    }
    alpha=ret[0]; nu=ret[1];
    // Compute new EP parameters without damping (to 'mprXXP'). If
    // selective damping is active, we also determine the effective damping
    // factor (overwrites 'dampFact')
    //printMsgStdout("sequentialUpdate: 2"); // DEBUG!
    for (ii=0; ii<vjSz; ii++) {
      // Undamped EP update
      i=vjInd[ii];
      bval=bP[ii];
      cPi=cPiP[ii]; cBeta=cBetaP[ii];
      //sprintf(debMsg,"ii=%d,i=%d,bval=%f,cPi=%f,cBeta=%f",ii,i,bval,cPi,
      //      cBeta); // DEBUG!
      //printMsgStdout(debMsg);
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
	prPi=temp*cPi*nu;
	prBeta=temp*(cBeta*nu+temp2*alpha);
     } else {
	// Very small but non-zero |b_ji| (will probably never happen)
	// 'temp' is b_ji / (pi_{-ji} - b_ji^2 nu_j)
	if ((temp=cPi-nu*bval*bval)<1e-10)
	  return updNumericalError; // EP update failed
	temp=bval/temp;
	prPi=temp*bval*nu*cPi;
	prBeta=temp*(cBeta*bval*nu+cPi*alpha);
      }
      mprPiP[ii]=prPi; mprBetaP[ii]=prBeta; // Intermed. storage
      if (!(epMaxPi==0) && prPi<piP[ii]) {
	//printMsgStdout("sequentialUpdate: HAEH???"); // DEBUG!
	// Selective damping to ensure that pi_{-ki} >= eps for all k,i
	temp2=epMaxPi->getMaxBeta(i); // kappa_i
	temp=std::min((mPiP[i]-std::max(temp2,0.0)-piMinThres)/(piP[ii]-prPi),
		      1.0);
	if (temp<=0.02) {
	  // EP update has to be skipped
	  if (effDamp!=0) *effDamp=1.0;
	  return updCavCondSkipped;
	}
	dampFact=std::max(dampFact,1.0-temp);
      }
    }
    if (effDamp!=0) *effDamp=dampFact;
    // Determine new EP parameters with damping (overwrite 'cXXP') and
    // new marginals (to 'mprXXP'). This is done because the update can
    // still fail ('updMarginalsInvalid')
    for (ii=0; ii<vjSz; ii++) {
      cPi=cPiP[ii]; cBeta=cBetaP[ii]; // Cavity
      prPi=mprPiP[ii]; prBeta=mprBetaP[ii]; // Undamped update
      //sprintf(debMsg,"ii=%d,prPi=%f,prBeta=%f",ii,prPi,prBeta); // DEBUG!
      //printMsgStdout(debMsg);
      // Damping
      if (dampFact>0.0) {
	//printMsgStdout("sequentialUpdate: HAEH2???"); // DEBUG!
	prPi+=dampFact*(piP[ii]-prPi);
	prBeta+=dampFact*(betaP[ii]-prBeta);
      }
      // New marginals (overwrite 'mprXXP') and EP parameters (overwrite
      // 'cXXP')
      if ((mprPiP[ii]=cPi+prPi)<thres2)
	return updMarginalsInvalid; // EP update failed
      mprBetaP[ii]=cBeta+prBeta;
      cPiP[ii]=prPi; cBetaP[ii]=prBeta; // New EP parameters
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
