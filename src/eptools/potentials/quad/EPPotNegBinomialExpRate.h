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
   * lam(s) = exp(s).
   * Parameters: y (nonneg. int.), r (positive)
   * <p>
   * We need numerical quadrature for this potential. We use the generic
   * Newton implementation.
   * <p>
   * ATTENTION: If 'SpecfunServices::logGamma' is not implemented, we
   * use C=1 as constant in front.
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
     * @param pverb Passed to 'QuadPotProximalNewton'. Def.: 0
     */
    EPPotNegBinomialExpRate(double py,double pr,double pacc,double pfacc,
			    int pverb=0) :
      QuadPotProximalNewton(pacc,pfacc,pverb),EPPotNegBinomialCommon(py,pr) {}

    bool isLogConcave() const {
      return true;
    }

    bool hasFirstDerivatives() const {
      return true;
    }

    bool hasSecondDerivatives() const {
      return true;
    }

    double eval(double s,double* dl=0,double* ddl=0) const;

    void initBracket(double h,double rho,double& l,double& r) const {
      if (rho<(1e-16))
	throw InvalidParameterException(EXCEPT_MSG(""));
      l=h-rscal*rho; r=h+yscal*rho;
    }
  };

  inline double EPPotNegBinomialExpRate::eval(double s,double* dl,
					      double* ddl) const
  {
    double lgr=log(rscal),sig,temp,ret;

    if (s>=lgr) {
      temp=exp(lgr-s);
      sig=1.0/(1.0+temp);
      ret=rscal*s+(rscal+yscal)*log1p(temp);
    } else {
      temp=exp(s-lgr);
      sig=temp/(1.0+temp);
      ret=-yscal*s+(rscal+yscal)*(lgr+log1p(temp));
    }
    if (dl!=0)
      (*dl) = (yscal+rscal)*sig;
    if (ddl!=0)
      (*ddl) = (yscal+rscal)*sig*(1.0-sig);

    return ret;
  }
//ENDNS

#endif
