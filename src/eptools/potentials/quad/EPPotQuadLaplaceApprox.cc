/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class EPPotQuadLaplaceApprox
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/quad/EPPotQuadLaplaceApprox.h"
#include "src/eptools/potentials/SpecfunServices.h"

//BEGINNS(eptools)
  EPPotQuadLaplaceApprox::EPPotQuadLaplaceApprox(const Handle<QuadPotProximal>& qpot,const Handle<QuadratureServices>& qserv) :
    EPPotQuadrature(qpot),quadServ(qserv),qpotProx(qpot.p())
  {
    if (!qpot->hasSecondDerivatives())
      throw InvalidParameterException(EXCEPT_MSG("Need 2nd derivatives"));
    // Check interval and waypoints
    double a,b;
    bool aInf,bInf;
    ArrayHandle<double> wayPts;
    qpot->getInterval(a,aInf,b,bInf,wayPts);
    if (!aInf && !bInf && b<=a)
      throw InvalidParameterException(EXCEPT_MSG("Interval [a,b] must not be empty"));
    int wsz=wayPts.size();
    if (qpot->hasWayPoints() && wsz>0) {
      double temp,prev=wayPts[0];
      for (int i=1; i<wsz; i++,prev=temp)
	if ((temp=wayPts[i])<=prev)
	  throw InvalidParameterException(EXCEPT_MSG("Waypoint list must be increasing"));
      if ((!aInf && a>=wayPts[0]) || (!bInf && b<=wayPts[wsz-1]))
	throw InvalidParameterException(EXCEPT_MSG("Waypoints must lie in (a,b)"));
    }
    // Integrand function
    intFuncPars.qpot=qpot.p();
    intFuncPars.k=0;
    intFunc.function=&EPPotQuadLaplaceApprox_intFunc;
    intFunc.params=(void*) &intFuncPars;
  }

  /*
   * TODO: Verbosity for debugging. And meaningful reaction to errors (right
   * now, we probably fail too often).
   */
  bool EPPotQuadLaplaceApprox::compMoments(double cmu,double crho,
					   double& alpha,double& nu,
					   double* logz,double eta) const
  {
    int i,wsz;
    double a,b,sstar,sigma;
    bool aInf,bInf,isCritical;
    ArrayHandle<double> wayPts;

    if (crho<1e-14 || eta<1e-10 || eta>1.0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    // Determine mode of integrand
    if (!qpotProx->proximal(cmu,eta*crho,sstar))
      return false; // Update fails if mode search not successful
    // Interval [a,b] and waypoints. Can we use 2nd derivative at 'sstar'?
    qpotProx->getInterval(a,aInf,b,bInf,wayPts);
    wsz=qpotProx->hasWayPoints()?wayPts.size():0;
    isCritical=false;
    if (!aInf && fabs(sstar-a)<1e-5)
      isCritical=true;
    else if (!bInf && fabs(sstar-b)<1e-5)
      isCritical=true;
    else if (wsz>0) {
      for (i=0; i<wsz; i++)
	if (fabs(sstar-wayPts[i])<1e-5) {
	  isCritical=true; break;
	}
    }
    // Configure integrand function (except for sigma). Have to do this here,
    // so can use 'getD2H' (does not depend on sigma)
    intFuncPars.h=cmu;
    intFuncPars.rho=crho;
    intFuncPars.eta=eta;
    intFuncPars.sstar=sstar;
    intFuncPars.k=0; // Z (0th order) is first
    if (isCritical)
      sigma=sqrt(crho); // Sitting on critical point: Fallback
    else {
      sigma=intFuncPars.getD2H(sstar); // h''(s_*)
      if (sigma<-(1e-10))
	// UUPS: Not really a minimum point! Fallback
	sigma=sqrt(crho);
      else
	sigma=1.0/sqrt(sigma+(1e-8)); // Beware of small curvature
    }
    // Finalize initialization and transform integration interval
    intFuncPars.sigma=sigma;
    intFuncPars.init(); // Precomputations
    if (!aInf) a=(a-sstar)/sigma;
    if (!bInf) b=(b-sstar)/sigma;
    for (i=0; i<wsz; i++)
      wayPts[i]=(wayPts[i]-sstar)/sigma;
    // Run quadrature calls. We first estimate the normalization constant
    // Z_til after mode normalization, then 1st and 2nd moment
    // TODO: Verbosity! React to errors appropriately.
    double ztil,ex1,ex2;
    if (quadServ->quad(intFunc,a,aInf,b,bInf,ztil,qpotProx->hasWayPoints(),
		       wayPts)!=0)
      return false; // Quadrature failure
    if (ztil<(1e-12))
      return false; // Z_til too small (failure of mode normalization?)
    if (logz!=0)
      *logz = log(ztil)-intFuncPars.hsstar+log(sigma)-
	0.5*(log(crho)+SpecfunServices::m_ln2pi);
    // We fold Z_til into the 1st and 2nd moment computation by subtracting
    // log Z_til from h(sstar) in 'intFuncPars'.
    // NOTE: Do not call 'intFuncPars.init()', would overwrite 'hsstar'!
    intFuncPars.hsstar-=log(ztil);
    intFuncPars.k=1;
    if (quadServ->quad(intFunc,a,aInf,b,bInf,ex1,qpotProx->hasWayPoints(),
		       wayPts)!=0)
      return false; // Quadrature failure
    intFuncPars.k=2;
    if (quadServ->quad(intFunc,a,aInf,b,bInf,ex2,qpotProx->hasWayPoints(),
		       wayPts)!=0)
      return false; // Quadrature failure
    // Can alpha, nu be estimated more directly? Here, we compute them from
    // E[x], E[x^2], expectation w.r.t. p_hat.
    alpha=(sigma*ex1+sstar-cmu)/crho;
    ex2-=ex1*ex1; // Variance
    nu=(1.0-ex2*sigma*sigma/crho)/crho;
  }
//ENDNS
