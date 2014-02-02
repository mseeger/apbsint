/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class QuadraturePotential
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_QUADRATUREPOTENTIAL_H
#define EPTOOLS_QUADRATUREPOTENTIAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalPotentialBase.h"

//BEGINNS(eptools)
  /**
   * Base class for services required by numerical quadrature implementations
   * of the 'EPScalarPotential' interface.
   * <p>
   * The 'eval' service returns the value of l(s) = -log t(s), as well as
   * (optional) 1st and 2nd derivative.
   *
   * 'getInterval' specifies the integration interval [a,b]. Here, a can be
   * -infty and/or b can be +infty. Optionally, the method returns an
   * increasing list of waypoints s_i (which can be empty):
   * - a < s_1 < ... < s_K < b
   * - l(s) is smooth in any open subinterval, but may be discontinuous,
   *   nondifferentiable, even singular at any waypoint or a, b
   * This information is returned if 'hasWayPoints' returns true.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class QuadraturePotential : public virtual EPScalPotentialBase
  {
  public:
    // Public methods

    /**
     * @return Can 1st derivatives be requested in 'eval'?
     */
    virtual bool hasFirstDerivatives() const = 0;

    /**
     * If this returns true, so must be 'hasFirstDerivatives'.
     *
     * @return Can 2nd derivatives be requested in 'eval'?
     */
    virtual bool hasSecondDerivatives() const = 0;

    /**
     * See header comment and 'getInterval'.
     *
     * @return Are waypoints returned in 'getInterval'?
     */
    virtual bool hasWayPoints() const = 0;

    /**
     * If l(s) = -log t(s), we return l(s) and (optionally) its 1st and 2nd
     * derivative. The optional return arguments are ignored if respective
     * 'hasXXXDerivatives' returns false. Passing 0 means they are not
     * requested.
     *
     * @param s   Argument
     * @param dl  First derivative l'(s) ret. here. Optional
     * @param ddl Second derivative l''(s) ret. here. Optional
     * @return    Function value l(s)
     */
    virtual double eval(double s,double* dl=0,double* ddl=0) const = 0;

    /**
     * Returns endpoints of integration interval [a,b] in 'a', 'b'.
     * If a==-infty, 'aInf'=true and 'a' is ignored. Likewise with b,
     * +infty, 'bInf'.
     * If 'hasWayPoints' returns true, an increasing list of waypoints
     * s_i is returned in 'wayPts' (see header comment). Note that a, b
     * are excluded. The list can be empty.
     * If 'hasWayPoints' returns false, 'wayPts' is ignored.
     *
     * @param a      Left interval boundary a (ignored if 'aInf'=true)
     * @param aInf   Is a = -infty?
     * @param b      Right interval boundary b (ignored if 'bInf'=true)
     * @param bInf   Is b = +infty?
     * @param wayPts S.a. Only if 'hasWayPoints' returns true
     */
    virtual void getInterval(double& a,bool& aInf,double& b,bool& bInf,
			     ArrayHandle<double>& wayPts) const = 0;
  };
//ENDNS

#endif
