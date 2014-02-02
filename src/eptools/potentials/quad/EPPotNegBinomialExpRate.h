/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotNegBinomialExpRate
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTNEGBINOMIALEXPRATE_H
#define EPTOOLS_EPPOTNEGBINOMIALEXPRATE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadPotProximalNewton.h"
#include "src/eptools/potentials/quad/EPPotNegBinomialCommon.h"

//BEGINNS(eptools)
  /**
   * Negative binomial potential with exponential rate function:
   *   t(s)   = C (1-p(s))^r p(s)^y,   y in N, r>0,
   *   p(s)   = lam(s)/(r + lam(s)),
   *   C      = Gamma(r+y)/(Gamma(y+1) Gamma(r))
   *   lam(s) = exp(s)
   * Parameters: y (nonneg. int.), r (positive)
   * <p>
   * We need numerical quadrature for this potential. We use the generic
   * Newton implementation.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotNegBinomialExpRate : public QuadPotProximalNewton,
				  public EPPotNegBinomialCommon
  {
  public:
    // Public methods

    /**
     * @param py    Value for y
     * @param pr    Value for r
     * @param pacc  See 'OneDimSolver::newton'. >0
     * @param pfacc "
     */
    EPPotNegBinomialExpRate(double py,double pr,double pacc,double pfacc) :
      QuadPotProximalNewton(pacc,pfacc),EPPotPoissonCommon(py,pr) {}

    bool isLogConcave() const {
      return true;
    }

    bool hasFirstDerivatives() const {
      return true;
    }

    bool hasSecondDerivatives() const {
      return true;
    }

    // HIER!!
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

      return lam-yscal*lam;
    }

    void initBracket(double h,double rho,double& l,double& r) const;
  };

  /*
   * Details in technical report
   */
  inline void EPPotNegBinomialExpRate::initBracket(double h,double rho,
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
