/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotDebugQuadLaplace
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTDEBUGQUADLAPLACE_H
#define EPTOOLS_EPPOTDEBUGQUADLAPLACE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadPotProximal.h"

//BEGINNS(eptools)
  /**
   * DEBUG code. Implements Laplace potential ('EPPotLaplace') as
   * 'QuadPotProximal'. This is to test the quadrature implementation in the
   * presence of a waypoint (which is y).
   *   t(s)  = (tau/2) exp( -tau |y-s| )
   * Parameters: y, tau>0.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotDebugQuadLaplace : public QuadPotProximal
  {
  protected:
    // Members

    double yscal,tau;

  public:
    // Public methods

    EPPotDebugQuadLaplace(double py=0.0,double ptau=1.0) {
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

    bool isLogConcave() const {
      return true;
    }

    bool hasFirstDerivatives() const {
      return true;
    }

    bool hasSecondDerivatives() const {
      return true;
    }

    bool hasWayPoints() const {
      return true;
    }

    /**
     * The integration interval is all of R. y is a waypoint, l(s) is not
     * differentiable there.
     */
    void getInterval(double& a,bool& aInf,double& b,bool& bInf,
		     ArrayHandle<double>& wayPts) const {
      aInf=bInf=true;
      if (wayPts.size()!=1)
	wayPts.changeRep(1);
      wayPts[0]=yscal;
    }

    double eval(double s,double* dl=0,double* ddl=0) const {
      double ret;
      // NOTE: Does not complain if s==yscal, treats it like s>yscal
      if (s>=yscal) {
	ret=tau*(s-yscal)-log(0.5*tau);
	if (dl!=0) *dl=tau;
      } else {
	ret=tau*(yscal-s)-log(0.5*tau);
	if (dl!=0) *dl=-tau;
      }
      if (ddl!=0) *ddl=0.0;

      return ret;
    }

    /*
     * This is the usual l_1 proximal map:
     * If x = s-y: argmin_x kappa |x| + 0.5 (x-mu)^2,
     * kappa = rho*tau, mu = h-y.
     * The solution x_* is soft shrinkage of mu by kappa.
     * NOTE: This maps s_* = y for all h close to y, so we sit on the
     * waypoint then (where l(s) is not differentiable).
     */
    bool proximal(double h,double rho,double& sstar) const {
      double mu=h-yscal,kap=rho*tau;

      sstar=yscal+((mu>kap)?(mu-kap):((mu<-kap)?(mu+kap):0.0));
      return true;
    }
  };
//ENDNS

#endif
