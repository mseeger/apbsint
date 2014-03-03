/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPPotQuadLaplaceApprox
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTQUADLAPLACEAPPROX_H
#define EPTOOLS_EPPOTQUADLAPLACEAPPROX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/EPPotQuadrature.h"
#include "src/eptools/potentials/quad/QuadPotProximal.h"
#include "src/eptools/potentials/quad/QuadratureServices.h"

//BEGINNS(eptools)
  /**
   * Information passed to 'EPPotQuadLaplaceApprox_intFunc'.
   * The integrand is
   *   g(x) = x^k exp( h(s_*) - h(s_* + sigma*x) ),
   *   h(s) = eta*l(s) + (s-h)^2/(2*rho)
   * Here, l(s) = -log t(s) is represented by the 'QuadraturePotential'
   * object.
   * <p>
   * NOTE: h(s) here lacks the additive constant 0.5*log(2*pi*rho) in
   * the TR.
   */
  class EPPotQuadLaplaceApprox_intFuncParams
  {
  public:
    QuadraturePotential* qpot;
    double h,rho,eta;
    double sstar,sigma;
    double hsstar;
    int k;

    /**
     * Must be called whenever parameters have been changed. Throws exception
     * if parameters invalid, does precomputations.
     */
    void init() {
      if (qpot==0 || rho<(1e-16) || eta<=0.0 || eta>1.0 || sigma<(1e-16) ||
	  k<0 || k>2)
	throw InvalidParameterException(EXCEPT_MSG(""));
      hsstar=getH(sstar);
    }

    double getH(double s) const {
      double temp=s-h;
      return eta*qpot->eval(s)+0.5*temp*temp/rho;
    }

    double getG(double x) const {
      double ret=exp(hsstar-getH(sstar+sigma*x));

      switch (k) {
      case 1:
	ret*=x;
	break;
      case 2:
	ret*=(x*x);
	break;
      }

      return ret;
    }

    /**
     * Second derivative h''(s). This does not depend on 'sstar', 'sigma'.
     */
    double getD2H(double s) const {
      double temp;
      qpot->eval(s,0,&temp);
      return eta*temp+1.0/rho;
    }
  };

  /**
   * Integrand function g(x) passed to quadrature routine. See
   * 'EPPotQuadLaplaceApprox_intFuncParams' comments.
   */
  inline double EPPotQuadLaplaceApprox_intFunc(double x,void* params)
  {
    return ((const EPPotQuadLaplaceApprox_intFuncParams*) params)->getG(x);
  }

  /**
   * Implementation of 'EPPotQuadrature': EP update service via numerical
   * quadrature. The integration variable is transformed by way of a
   * Laplace approximation.
   * <p>
   * Details are given in the technical report. The quadrature potential
   * object 'quadPot' must be of type 'QuadPotProximal', and 2nd derivatives
   * must be given. If there are points in (a,b) where l(s) = -log t(s) is
   * not twice continuously differentiable, they should be passed as waypoints
   * (even though we do not require 'quadPot->hasWayPoints' to return true: we
   * assume no waypoints then).
   *
   * We first determine the mode of the integrand for Z by 'proximal'. The
   * value of the integrand there is pulled outside, which hopefully counters
   * underflow. We also transform the integration variable using the 2nd
   * derivative of l(s) there. This is not done (and we standardize using the
   * variance rho) if the minimum point mu of l(s) is equal or very close to a
   * waypoint or one of a, b. The normalized and transformed integrand g(x)
   * is represented by 'EPPotQuadLaplaceApprox_intFuncParams'. It is passed
   * to the quadrature code in 'quadServ' for computing 0th, 1st, 2nd moments.
   * <p>
   * NOTE: Using the Laplace transformation together with sophisticated
   * adaptive quadrature code is probably overkill. But it can be combined
   * with cheap non-adaptive quadrature (e.g., Gauss-Hermite). If the
   * transformed integral is
   *   int_a^b g(x) d x,
   * the idea is that g(x) is "close to" N(x|0,1), so Gauss-Hermite could
   * be applied to
   *   int_a^b [g(x)/N(x|0,1)] N(x|0,1) d x,
   * where N(x|0,1) is the standardized weight function, and g(x)/N(x|0,1)
   * is (hopefully) well-approximated by a low-order polynomial.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotQuadLaplaceApprox : public EPPotQuadrature
  {
  protected:
    // Members

    QuadPotProximal* qpotProx;            // 'quadPot' with correct type
    Handle<QuadratureServices> quadServ;  // Quadrature services
    quad_function intFunc;                // Represents integrand g(x)
    mutable EPPotQuadLaplaceApprox_intFuncParams intFuncPars;

  public:
    // Public methods

    /**
     * Constructor
     *
     * @param qpot  Quadrature potential (must be 'QuadPotProximal')
     * @param qserv Quadrature services
     */
    EPPotQuadLaplaceApprox(const Handle<QuadPotProximal>& qpot,
			   const Handle<QuadratureServices>& qserv);

    bool suppFractional() const {
      return true; // Fractional EP generally supported
    }

    /**
     * Right now, we return with failure if the proximal map computation
     * fails. This could be replaced by a fallback, say evaluate the
     * integrand at cmu and another dedicated place, normalizing by the
     * maximum for these, and transforming by crho.
     * <p>
     * We also return with failure if any of the quadrature service calls
     * returns with status !=0. Again, this may be too stringent.
     */
    bool compMoments(const double* inp,double* ret,double* logz=0,
		     double eta=1.0) const;
  };
//ENDNS

#endif
