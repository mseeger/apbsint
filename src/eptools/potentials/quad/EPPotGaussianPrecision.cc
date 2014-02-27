/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class EPPotGaussianPrecision
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/quad/EPPotGaussianPrecision.h"

//BEGINNS(eptools)
  bool EPPotGaussianPrecision::compMoments(double cmu,double crho,double ca,
					   double cc,double& alpha,double& nu,
					   double& hata,double& hatc,
					   double* logz,double eta) const
  {
    int i,wsz,verbose=quadServ->getVerbose();
    double temp,limA,vstar,sigma;

    if (eta!=1.0)
      throw NotImplemException(EXCEPT_MSG(""));
    if (crho<1e-14 || ca<1e-14 || cc<1e-14)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (verbose>0)
      cout << "EPPotGaussianPrecision::compMoments: cmu=" << cmu << ",crho="
	   << crho << ",ca=" << ca << ",cc=" << cc << endl;
    // Prepare integrand function
    intFuncPars.a=ca;
    intFuncPars.cdrho=cc/crho;
    temp=cmu-yscal;
    intFuncPars.xi=temp*temp/crho;
    intFuncPars.l=0;
    intFuncPars.init();
    // Laplace transformation if ca>0.5
    bool doLaplace=(ca>0.5001);
    if (doLaplace) {
      // Determine mode of integrand
      double gamma=2.0*intFuncPars.cdrho;
      double x0,x1,x2;
      // DEBUG:
      if (gamma<1e-12)
	cout << "EPPotGaussianPrecision::compMoments: Small gamma=" << gamma
	     << " (numerical issues!)" << endl;
      if (SpecfunServices::rootsCubicPolynomial(2.0*(gamma-ca+1.0)/gamma,(gamma+intFuncPars.xi-4.0*ca+3.0)/gamma,(1.0-2.0*ca)/gamma,x0,x1,x2)==3) {
	// Either x0 or x2
	if (x2<=0.0) {
	  cout << "ERROR(EPPotGaussianPrecision::compMoments): All cubic roots are negative!" << endl;
	  return false;
	}
	if (x0<=0.0)
	  x0=x2;
	else {
	  // All roots are positive. Pick the one with smaller h(v) value
	  if (intFuncPars.getH(x2)<intFuncPars.getH(x0))
	    x0=x2;
	}
      } else if (x0<=0.0) {
	cout << "ERROR(EPPotGaussianPrecision::compMoments): Cubic root is negative!" << endl;
	return false;
      }
      vstar=x0; // Mode of integrand
      if (verbose>0)
	cout << "  v_* = " << vstar << endl;
      // Set 'sigma' depending on curvature at mode
      sigma=intFuncPars.getD2H(vstar); // h''(v_*)
      if (sigma<-(1e-10))
	// UUPS: Not really a minimum point! Fallback
	sigma=1.0;
      else
	sigma=1.0/sqrt(sigma+(1e-8)); // Beware of small curvature
      if (verbose>0)
	cout << "  sigma = " << sigma << endl;
    } else {
      // No transformation
      vstar=0.0; sigma=1.0;
    }
    // HIER!
    // Finalize initialization and transform integration interval
    intFuncPars.sigma=sigma;
    intFuncPars.init(); // Precomputations
    if (verbose>0)
      cout << "  sigma=" << sigma << endl;
    if (!aInf) a=(a-sstar)/sigma;
    if (!bInf) b=(b-sstar)/sigma;
    for (i=0; i<wsz; i++)
      wayPts[i]=(wayPts[i]-sstar)/sigma;
    // Run quadrature calls. We first estimate the normalization constant
    // Z_til after mode normalization, then 1st and 2nd moment
    // TODO: Verbosity! React to errors appropriately.
    double ztil,ex1,ex2;
    if (quadServ->quad(intFunc,a,aInf,b,bInf,ztil,qpotProx->hasWayPoints(),
		       wayPts)!=0) {
      if (verbose>0)
	cout << "  Quad(k=0) fails" << endl;
      return false; // Quadrature failure
    }
    if (ztil<(1e-12)) {
      if (verbose>0)
	cout << "  Z_til too small (" << ztil << ")" << endl;
      return false; // Z_til too small (failure of mode normalization?)
    }
    if (logz!=0)
      *logz = log(ztil)-intFuncPars.hsstar+log(sigma)-
	0.5*(log(crho)+SpecfunServices::m_ln2pi);
    // We fold Z_til into the 1st and 2nd moment computation by subtracting
    // log Z_til from h(sstar) in 'intFuncPars'.
    // NOTE: Do not call 'intFuncPars.init()', would overwrite 'hsstar'!
    intFuncPars.hsstar-=log(ztil);
    intFuncPars.k=1;
    if (quadServ->quad(intFunc,a,aInf,b,bInf,ex1,qpotProx->hasWayPoints(),
		       wayPts)!=0) {
      if (verbose>0)
	cout << "  Quad(k=1) fails" << endl;
      return false; // Quadrature failure
    }
    intFuncPars.k=2;
    if (quadServ->quad(intFunc,a,aInf,b,bInf,ex2,qpotProx->hasWayPoints(),
		       wayPts)!=0) {
      if (verbose>0)
	cout << "  Quad(k=2) fails" << endl;
      return false; // Quadrature failure
    }
    // Can alpha, nu be estimated more directly? Here, we compute them from
    // E[x], E[x^2], expectation w.r.t. p_hat.
    alpha=(sigma*ex1+sstar-cmu)/crho;
    ex2-=ex1*ex1; // Variance
    nu=(1.0-ex2*sigma*sigma/crho)/crho;

    return true;
    
  }
//ENDNS
