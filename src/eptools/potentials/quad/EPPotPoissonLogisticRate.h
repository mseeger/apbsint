/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotPoissonLogisticRate
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTPOISSONLOGISTICRATE_H
#define EPTOOLS_EPPOTPOISSONLOGISTICRATE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadPotProximalNewton.h"
#include "src/eptools/potentials/quad/EPPotPoissonCommon.h"

//BEGINNS(eptools)
  /**
   * Poisson potential with logistic rate function:
   *   t(s) = (y!)^-1 lam(s)^y exp(-lam(s)),  y in N,
   *   lam(s) = log(1 + exp(s))
   * Parameters: y (nonneg. int.)
   * <p>
   * We need numerical quadrature for this potential. We use the generic
   * Newton implementation.
   * <p>
   * ATTENTION: If 'SpecfunServices::logGamma' is not implemented, we drop
   * the constant (y!)^-1 in front.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotPoissonLogisticRate : public QuadPotProximalNewton,
				   public EPPotPoissonCommon
  {
  public:
    // Public methods

    /**
     * @param py    Value for y
     * @param pacc  See 'OneDimSolver::newton'. >0
     * @param pfacc "
     * @param pverb Passed to 'QuadPotProximalNewton'. Def.: 0
     */
    EPPotPoissonLogisticRate(double py,double pacc,double pfacc,int pverb=0) :
      QuadPotProximalNewton(pacc,pfacc,pverb),EPPotPoissonCommon(py) {}

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
      double temp,sig,lam,sgdlm;

      if (s>=0.0) {
	temp=exp(-s);
	sig=1.0/(1.0+temp);
	lam=s+log1p(temp);
	sgdlm=sig/lam;
      } else {
	temp=exp(s);
	sig=temp/(1.0+temp);
	lam=log1p(temp);
	sgdlm=(s>-10.0)?(sig/lam):(1.0/(1.0+temp));
      }
      if (dl!=0)
	(*dl) = sig-yscal*sgdlm;
      if (ddl!=0) {
	temp=1.0-sig;
	(*ddl) = sig*temp+yscal*sgdlm*(sgdlm-temp);
      }

      return lam-yscal*log(lam)+logYFact;
    }

    void initBracket(double h,double rho,double& l,double& r) const;
  };

  /*
   * Details in technical report
   */
  inline void EPPotPoissonLogisticRate::initBracket(double h,double rho,
						    double& l,double& r) const
  {
    int i;
    double a,sga;
    double acand[] = {2.20, 1.39, 0.85, 0.41, 0.0};

    if (rho<(1e-16))
      throw InvalidParameterException(EXCEPT_MSG(""));
    l=h-rho;
    for (i=0; i<5; i++) {
      a=acand[i]; sga=1.0/(1.0+exp(-a));
      r=h-sga*rho;
      if (yscal>0.0)
	r=0.5*(r+sqrt(r*r+4.0*yscal*rho));
      if (r>a)
	break;
    }
  }
//ENDNS

#endif
