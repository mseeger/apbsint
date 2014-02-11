/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class QuadPotProximalNewton
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/quad/QuadPotProximalNewton.h"
#include "lhotse/optimize/OneDimSolver.h"

//BEGINNS(eptools)
  QuadPotProximalNewton::QuadPotProximalNewton(double pacc,double pfacc,
					       int pverb) :
    acc(pacc),facc(pfacc),verbose(pverb)
  {
    if (pacc<=0.0 || pfacc<=0.0 || pverb<0)
      throw InvalidParameterException(EXCEPT_MSG(""));
  }

  bool QuadPotProximalNewton::proximal(double h,double rho,double& sstar) const
  {
    if (proxFun==0)
      proxFun.changeRep(new QuadPotProximalNewton_Func1D(this));
    proxFun->setPars(h,rho);
    try {
      // Initial bracket
      double bL,bR;
      initBracket(h,rho,bL,bR);
      int brRight=(bR>bL)?OneDimSolver::brackRightRegular:
	OneDimSolver::brackRightInfinite;
      if (verbose>0) {
	cout << "  QuadPotProximalNewton: Bracket=[" << bL << ",";
	if (bR>bL)
	  cout << bR << "]" << endl;
	else
	  cout << "infty)" << endl;
      }
      // Run Newton solver
      sstar = OneDimSolver::newton(proxFun.p(),bL,bR,acc,facc,brRight,0.0,
				   "QuadPotProximalNewton");
    } catch (...) {
      return false; // Exception thrown in 'OneDimSolver::newton'
    }

    return true;
  }
//ENDNS
