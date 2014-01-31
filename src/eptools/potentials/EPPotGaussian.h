/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotGaussian
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTGAUSSIAN_H
#define EPTOOLS_EPPOTGAUSSIAN_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"
#include "src/eptools/potentials/SpecfunServices.h"

//BEGINNS(eptools)
  /**
   * Gaussian potential:
   *   t(s) = N(s | y,ssq) = N(y | s,ssq)
   * Parameters: y (mean), ssq (variance)
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotGaussian : public EPScalarPotential
  {
  protected:
    double yscal,ssq;

  public:
    EPPotGaussian(double py=0.0,double pssq=1.0) : yscal(py),ssq(pssq) {
      if (pssq<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    int numPars() const {
      return 2;
    }

    void getPars(double* pv) const {
      pv[0]=yscal; pv[1]=ssq;
    }

    void setPars(const double* pv) {
      if (!isValidPars(pv))
	throw InvalidParameterException(EXCEPT_MSG(""));
      yscal=pv[0]; ssq=pv[1];
    }

    bool isValidPars(const double* pv) const {
      return (pv[1]>=1e-13);
    }

    bool suppFractional() const {
      return true;
    }

    bool isLogConcave() const {
      return true;
    }

    /*
     * About fractional:
     *   t(s)^eta = N(s | y,ssq/eta) eta^{-1/2}
     *   log Z(ssq,eta) = log Z(ssq/eta,1) - 0.5 log(eta)
     */
    bool compMoments(double cmu,double crho,double& alpha,double& nu,
		     double* logz=0,double eta=1.0) const {
      double temp;

      if (crho<=0.0 || eta>1.0 || eta<=0.0)
	return false;
      nu=1.0/(crho+ssq/eta);
      alpha=nu*(temp=yscal-cmu);
      if (logz!=0)
	*logz=-0.5*(nu*temp*temp-log(nu)+SpecfunServices::M_LN2PI+log(eta));

      return true;
    }
  };
//ENDNS

#endif
