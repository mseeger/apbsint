/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class EPPotGaussianPrecision
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/quad/EPPotGaussianPrecision.h"

//BEGINNS(eptools)
  bool EPPotGaussianPrecision::compMoments(const double* inp,double* ret,
					   double* logz,double eta) const
  {
    int verbose=quadServ->getVerbose();
    double temp,vstar,sigma;
    double cmu=inp[0],crho=inp[1],ca=inp[2],cc=inp[3];

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
    // Quadrature calls: log Z_tilde and kappa moments
    // The integral are over [0,infty) (before transformation), and the
    // integrand is smooth. It has a singularity at 0 iff 'ca'<1/2.
    double lztil,ex1,ex2,hvstar,limA;
    intFuncPars.vstar=vstar;
    intFuncPars.sigma=sigma;
    hvstar=doLaplace?intFuncPars.getH(vstar):0.0;
    intFuncPars.off=hvstar;
    limA=-vstar/sigma;
    if (quadServ->quad(intFunc,limA,false,limA,true,lztil,true)!=0) {
      if (verbose>0)
	cout << "  Quad(lztil, l=0) fails" << endl;
      return false; // Quadrature failure
    }
    if (lztil<(1e-12)) {
      if (verbose>0)
	cout << "  Z_til too small (" << lztil << ")" << endl;
      return false; // Z_til too small (failure of mode normalization?)
    }
    lztil=log(lztil)-hvstar; // log Z_til
    if (logz!=0)
      *logz = lztil-0.5*(log(crho)+SpecfunServices::m_ln2pi);
    intFuncPars.off=-lztil; // New offset is -log Z_til
    intFuncPars.l=1;
    if (quadServ->quad(intFunc,limA,false,limA,true,ex1,true)!=0) {
      if (verbose>0)
	cout << "  Quad(ex1, l=1) fails" << endl;
      return false; // Quadrature failure
    }
    intFuncPars.l=2;
    if (quadServ->quad(intFunc,limA,false,limA,true,ex2,true)!=0) {
      if (verbose>0)
	cout << "  Quad(ex2, l=2) fails" << endl;
      return false; // Quadrature failure
    }
    ret[0]=ex1*(yscal-cmu)/crho; // alpha
    ret[1]=(ex1-intFuncPars.xi*(ex2-ex1*ex1))/crho; // nu
    // tau moments -> a_hat, c_hat
    intFuncPars.l=0; // Reset
    intFuncPars.off=-lztil;
    intFuncPars.a=ca+1.0;
    intFuncPars.init();
    if (quadServ->quad(intFunc,limA,false,limA,true,ex1,true)!=0) {
      if (verbose>0)
	cout << "  Quad(ex_tau1) fails" << endl;
      return false; // Quadrature failure
    }
    ex1*=(ca/cc);
    if (ex1<(1e-12)) {
      if (verbose>0)
	cout << "  E[tau] too small (" << ex1 << ")" << endl;
      return false;
    }
    intFuncPars.off=-log(ex1);
    intFuncPars.a=ca+2.0;
    intFuncPars.init();
    if (quadServ->quad(intFunc,limA,false,limA,true,ex2,true)!=0) {
      if (verbose>0)
	cout << "  Quad(ex_tau2) fails" << endl;
      return false; // Quadrature failure
    }
    ex2*=((ca+1.0)/cc);
    if (ex2-ex1<(1e-12)) {
      if (verbose>0)
	cout << "  x2-x1 too small (" << ex2-ex1 << ")" << endl;
      return false;
    }
    ret[3]=1.0/(ex2-ex1); // hat_c
    ret[2]=ex1*ret[3]; // hat_a

    return true;
  }
//ENDNS
