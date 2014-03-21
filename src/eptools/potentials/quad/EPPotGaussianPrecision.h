/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotGaussianPrecision
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTGAUSSIANPRECISION_H
#define EPTOOLS_EPPOTGAUSSIANPRECISION_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"
#include "src/eptools/potentials/SpecfunServices.h"
#include "src/eptools/potentials/quad/QuadratureServices.h"

//BEGINNS(eptools)
  /**
   * Represents integrand function g(x) and its parameters:
   *   g(x) = exp( off - h_l(v_* + sigma*x) ),
   * where h_l(v) depends on a, c/rho, xi and l:
   *   h_l(v) = -log f_l(kappa,xi) - log G(v|a,c/rho),
   *   kappa = v/(1+v)
   */
  class EPPotGaussianPrecision_intFuncParams
  {
  public:
    double a,cdrho,xi;
    int l;
    double vstar,sigma;
    double off,cnst;

    /**
     * Must be called whenever parameters a, c/rho or xi are changed.
     */
    void init() {
      if (a<(1e-16) || cdrho<(1e-16) || xi<0.0 || l<0 || l>2)
	throw InvalidParameterException(EXCEPT_MSG(""));
      cnst=0.5*xi-a*log(cdrho)+SpecfunServices::logGamma(a);
    }

    double getH(double v) const {
      int dl=(double) l;

      return (dl+0.5)*log1p(v)-(a+dl-0.5)*log(v)-0.5*xi/(1.0+v)+cdrho*v+cnst;
    }

    double getG(double x) const {
      return exp(off-getH(vstar+sigma*x));
    }

    /**
     * Second derivative h_0''(v). This does not depend on 'cdrho', 'l'.
     */
    double getD2H(double v) const {
      double temp=v+1.0;

      return (a-0.5)/v/v-(0.5+xi/temp)/temp/temp;
    }
  };

  /**
   * Integrand function g(x) passed to quadrature routine. See
   * 'EPPotGaussianPrecision_intFuncParams' comments.
   */
  inline double EPPotGaussianPrecision_intFunc(double x,void* params)
  {
    return ((const EPPotGaussianPrecision_intFuncParams*) params)->getG(x);
  }

  /**
   * Gaussian potential of input s and precision tau:
   *   t(s,tau) = N(s | y, tau^-1)
   * Parameters: y.
   * Argument group: 'atypeBivarPrec'.
   * <p>
   * We need quadrature in order to integrate w.r.t. tau, done by 'quadServ'.
   * If bounded above, the integrand is transformed by a Laplace
   * approximation (see also 'EPPotQuadLaplaceApprox'), the corr. mode can
   * be solved for analytically (requires roots of cubic equation).
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotGaussianPrecision : public EPScalarPotential
  {
  protected:
    // Members

    double yscal;
    Handle<QuadratureServices> quadServ; // Quadrature services
    quad_function intFunc;               // Represents integrand g(x)
    mutable EPPotGaussianPrecision_intFuncParams intFuncPars;

  public:
    // Public methods

    /**
     * Constructor
     *
     * @param qserv Quadrature services
     * @param py    Init. value for y. Def.: 0
     */
    EPPotGaussianPrecision(const Handle<QuadratureServices>& qserv,
			   double py=0.0) : quadServ(qserv) {
      setY(py);
      // Integrand function
      intFuncPars.l=0;
      intFunc.function=&EPPotGaussianPrecision_intFunc;
      intFunc.params=(void*) &intFuncPars;
    }

    virtual double getY() const {
      return yscal;
    }

    virtual void setY(double py) {
      yscal=py;
    }

    int numPars() const {
      return 1;
    }

    void getPars(double* pv) const {
      pv[0]=getY();
    }

    void setPars(const double* pv) {
      setY(pv[0]);
    }

    bool isValidPars(const double* pv) const {
      return true;
    }

    bool suppFractional() const {
      return false;
    }

    bool isLogConcave() const {
      return false;
    }

    static int getArgumentGroup_static() {
      return atypeBivarPrec;
    }

    int getArgumentGroup() const {
      return EPPotGaussianPrecision::getArgumentGroup_static();
    }

    /**
     * - inp: [cmu,crho,ca,cc]
     * - ret: [alpha,nu,ahat,chat]
     * Integration over tau requires numerical quadrature over [0,infty).
     * If 'ca'>1/2, we apply a Laplace approximation to transform and
     * normalize the integrand (see 'EPPotQuadLaplaceApprox'). Mode finding
     * is analytically tractable (requires roots of cubic equation).
     * If 'ca'<=1/2, the integrand's mode is at 0 (left boundary), a
     * singularity if 'ca'<1/2. In this case, we do not apply a
     * transformation.
     */
    bool compMoments(const double* inp,double* ret,double* logz=0,
		     double eta=1.0) const;
  };
//ENDNS

#endif
