/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotPoissonExpRate
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTPOISSONEXPRATE_H
#define EPTOOLS_EPPOTPOISSONEXPRATE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadPotProximal.h"
#include "src/eptools/potentials/quad/EPPotPoissonCommon.h"
#include "lhotse/optimize/FuncOneDim.h"

//BEGINNS(eptools)
  /**
   * Represents
   *   f(s) = e^s + s - a
   * Used to drive 'OneDimSolver'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotPoissonExpRate_Func1D : public FuncOneDim
  {
  protected:
    // Members

    double ascal;

  public:
    // Public methods

    EPPotPoissonExpRate_Func1D() : ascal(1.0) {}

    void setA(double pa) {
      ascal=pa;
    }

    bool hasDerivative() const {
      return true;
    }

    void eval(double x,double* f,double* df) {
      double temp=exp(x);
      *f=temp+x-ascal;
      *df=temp+1.0;
    }
  };

  /**
   * Poisson potential with exponential rate function:
   *   t(s) = (y!)^-1 lam(s)^y exp(-lam(s)),  y in N,
   *   lam(s) = exp(s)
   * Parameters: y (nonneg. int.)
   * <p>
   * We need numerical quadrature for this potential. 'proximal' requires a
   * very simple 1D Newton optimization, which is implemented here (we use
   * 'OneDimSolver', but time could be saved by running few hardcoded Newton
   * steps).
   * <p>
   * ATTENTION: If 'SpecfunServices::logGamma' is not implemented, we drop
   * the constant (y!)^-1 in front.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotPoissonExpRate : public QuadPotProximal,public EPPotPoissonCommon
  {
  protected:
    // Members

    double acc,facc;
    mutable Handle<EPPotPoissonExpRate_Func1D> proxFun;

  public:
    // Public methods

    /**
     * @param py    Value for y
     * @param pacc  See 'OneDimSolver::newton'. >0
     * @param pfacc "
     */
    EPPotPoissonExpRate(double py,double pacc,double pfacc) :
      EPPotPoissonCommon(py),acc(pacc),facc(pfacc) {
      if (pacc<=0.0 || pfacc<=0.0)
	throw InvalidParameterException(EXCEPT_MSG(""));
      proxFun.changeRep(new EPPotPoissonExpRate_Func1D());
    }

    bool isLogConcave() const {
      return true;
    }

    bool hasFirstDerivatives() const {
      return true;
    }

    bool hasSecondDerivatives() const {
      return true;
    }

    double eval(double s,double* dl=0,double* ddl=0) const {
      double temp=exp(s);

      if (dl!=0) (*dl)=temp-yscal;
      if (ddl!=0) (*ddl)=temp;
      return temp-s*yscal+logYFact;
    }

    bool proximal(double h,double rho,double& sstar) const;
  };
//ENDNS

#endif
