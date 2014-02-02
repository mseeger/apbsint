/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class QuadratureServices
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_QUADRATURESERVICES_H
#define EPTOOLS_QUADRATURESERVICES_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/default.h"

//BEGINNS(eptools)
  /**
   * Type for integrand functions. Additional parameters can be passed via
   * 'params'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  typedef struct {
    double (*function) (double x,void* params);
    void* params;
  } quad_function;

  /**
   * Base class for numerical quadrature services required by subclasses of
   * 'EPPotQuadrature'.
   * <p>
   * An instance of this class can be shared by several potential objects,
   * as long as the corresponding quadrature calls do not happen in parallel.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class QuadratureServices
  {
  public:
    // Public methods

    /**
     * @return Does 'quad' return an estimate of the final absolute error?
     */
    virtual bool hasAbsErrorEstimate() const = 0;

    /**
     * Main quadrature service:
     *   I = int_a^b f(x) d x.
     * a is -infty if 'aInf'==true ('a' ignored then), same with b, +infty,
     * and 'bInf'.
     * 'wayPts' (optional) is an array of "critical" points, where f(x) is
     * singular, discontinuous, not differentiable, etc. If given, the
     * entries must be increasing, and a < 'wayPts[0]', 'wayPts[end-1]' < b
     * (the array does not include a, b).
     * <p>
     * An absolute error estimate is written to 'abserr' if given, and if
     * 'hasAbsErrorEstimate' returns true, otherwise it is ignored.
     * The return status is 0 if the computation succeeded (the semantics
     * depend on the implementation: typically, estimates of absolute and/or
     * relative error fall below thresholds). If 'errmsg' is given and the
     * ret. status is != 0, an error message is written there (semantics
     * depend on implementation).
     *
     * @param fun    Integrand function f
     * @param a      Left interval boundary a (ignored if 'aInf'==true)
     * @param aInf   Is a = -infty?
     * @param b      Right interval boundary b (ignored if 'bInf'==true)
     * @param bInf   Is b = +infty?
     * @param ival   Integral value I returned here
     * @param wayPts S.a. Can be 0
     * @param abserr S.a. Def.: 0
     * @param errmsg S.a. Def.: 0
     * @return       Return status (0: Success)
     */
    virtual int quad(quad_function& fun,double a,bool aInf,double b,bool bInf,
		     double& ival,const ArrayHandle<double>& wayPts=
		     ArrayHandleZero<double>::get(),double* abserr=0,
		     string* errmsg=0);
  };
//ENDNS

#endif
