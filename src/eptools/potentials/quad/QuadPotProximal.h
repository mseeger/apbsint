/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class QuadPotProximal
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_QUADPOTPROXIMAL_H
#define EPTOOLS_QUADPOTPROXIMAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadraturePotential.h"

//BEGINNS(eptools)
  /**
   * Extends 'QuadraturePotential' by the proximal map service. This is
   * required by certain quadrature implementations which transform variables
   * so to avoid underflow for the normalization constant.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class QuadPotProximal : public virtual QuadraturePotential
  {
  public:
    // Public methods

    /**
     * If l(s) = -log t(s), the proximal map (in our context) is:
     *   s_* = argmin_s rho l(s) + (1/2) (s - h)^2
     * If l(s) is convex ('isLogConcave' returns true), this is uniquely
     * solvable by 1D convex minimization.
     * NOTE: Typically, h is the cavity mean h{-}, while rho is eta*rho{-},
     * rho{-} the cavity variance, eta the fractional parameter.
     *
     * @param h     Parameter
     * @param rho   Parameter (positive)
     * @param sstar Result s_* ret. here
     * @return      Successful?
     */
    virtual bool proximal(double h,double rho,double& sstar) const = 0;
  };
//ENDNS

#endif
