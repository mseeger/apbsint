/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class EPPotPoissonExpRate
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/quad/EPPotPoissonExpRate.h"
#include "lhotse/optimize/OneDimSolver.h"

//BEGINNS(eptools)
  /*
   * Bracket [L,R]
   *   f(s) = exp(s) + s - a
   * (1) If a<=1: L = a-exp(a), R = a
   * (2) If a>1: R = log(a) gives f(R) = log(a) > 0. We choose
   * L = (1-u) log(a), which gives
   * f(L) = a (a^{-u} - 1) + (1-u) log(a). Setting
   * a (a^{-u} - 1) + log(a) equal to 0 gives
   * u = -log(1 - (log a)/a)/log(a) and
   * f(L) = -u log(a) < 0. Here, u is roughly 1/a as a gets large, so the
   * bracket size behaves as log(a)/a.
   */
  bool EPPotPoissonExpRate::proximal(double h,double rho,double& sstar) const
  {
    double ascal,bL,bR;

    if (rho<(1e-16))
      throw InvalidParameterException(EXCEPT_MSG(""));
    proxFun->setA(ascal=h+yscal*rho+log(rho));
    if (ascal<=1.001) {
      bL=ascal-exp(ascal); bR=ascal;
    } else {
      bR=log(ascal);
      bL=bR-log1p(-bR/ascal);
    }
    try {
      // Run Newton solver
      sstar = OneDimSolver::newton(proxFun.p(),bL,bR,acc,facc,
				   OneDimSolver::brackRightRegular,0.0,
				   "EPPotPoissonExpRate") - log(rho);
    } catch (...) {
      return false; // Exception thrown in 'OneDimSolver::newton'
    }

    return true;
  }
//ENDNS
