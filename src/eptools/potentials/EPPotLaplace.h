/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotLaplace
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTLAPLACE_H
#define EPTOOLS_EPPOTLAPLACE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPPotQuantileRegress.h"

//BEGINNS(eptools)
  /**
   * Laplace (double exponential) potential:
   *   t(s)  = (tau/2) exp( -tau |y-s| )
   * Parameters: y, tau>0.
   * <p>
   * NOTE: Special case of 'EPPotQuantileRegress', just calls static method
   * there.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotLaplace : public EPScalarPotential
  {
  protected:
    // Members

    double yscal,tau;

  public:
    // Public methods

    EPPotLaplace(double py=0.0,double ptau=1.0) {
      setY(py); setTau(ptau);
    }

    virtual double getTau() const {
      return tau;
    }

    virtual void setTau(double ptau) {
      if (ptau<1e-12)
	throw InvalidParameterException(EXCEPT_MSG(""));
      tau=ptau;
    }

    virtual double getY() const {
      return yscal;
    }

    virtual void setY(double py) {
      yscal=py;
    }

    int numPars() const {
      return 2;
    }

    void getPars(double* pv) const {
      pv[0]=getY(); pv[1]=getTau();
    }

    void setPars(const double* pv) {
      setY(pv[0]); setTau(pv[1]);
    }

    bool isValidPars(const double* pv) const {
      return (pv[1]>=1e-12);
    }

    bool suppFractional() const {
      return true;
    }

    bool isLogConcave() const {
      return true;
    }

    /*
     * t(s)^eta corresponds to C times 'EPPotQuantileRegress' with
     *   kappa = 1/2, xi = 2 eta tau, C = (tau/2)^eta.
     */
    bool compMoments(double cmu,double crho,double& alpha,double& nu,
		     double* logz=0,double eta=1.0) const {
      if (crho<1e-14 || eta<1e-10 || eta>1.0)
	throw InvalidParameterException(EXCEPT_MSG(""));

      bool ret =
	EPPotQuantileRegress::compMomentsInt(cmu,crho,2.0*eta*tau,yscal,0.5,
					     alpha,nu,logz);
      if (ret && logz!=0)
	(*logz) += eta*log(0.5*tau);

      return ret;
    }
  };
//ENDNS

#endif
