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
  QuadPotProximalNewton::QuadPotProximalNewton(double pacc,double pfacc) :
    acc(pacc),facc(pfacc)
  {
    if (pacc<=0.0 || pfacc<=0.0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    proxFun.changeRep(new QuadPotProximalNewton_Func1D(this));
  }

  bool QuadPotProximalNewton::proximal(double h,double rho,double& sstar)
  {
    proxFun->setPars(h,rho);
    try {
      // Initial bracket
      double bL,bR;
      initBracket(h,rho,bL,bR);
      int brRight=(bR>bL)?OneDimSolver::brackRightRegular:
	OneDimSolver::brackRightInfinite;
      // Run Newton solver
      sstar = OneDimSolver::newton(proxFun.p(),bL,bR,acc,facc,brRight,0.0,
				   "QuadPotProximalNewton");
    } catch (...) {
      return false; // Exception thrown in 'OneDimSolver::newton'
    }

    return true;
  }
//ENDNS
