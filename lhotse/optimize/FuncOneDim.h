//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: optimize
 * Desc.:  Header abstract class FuncOneDim
 * ------------------------------------------------------------------- */

#ifndef FUNCONEDIM_H
#define FUNCONEDIM_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/optimize/default.h"

//BEGINNS(optimize)
  /**
   * Simple abstract class to represent a function from a subset of \R into
   * \R, together with its first derivative. This is used by the 1-D solvers
   * in class 'OneDimSolver'.
   * All subclasses must implement the virtual methods def. here. Most of them
   * will have additional methods, e.g. to chance other parameters than the
   * eval. point the function depends upon.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FuncOneDim
  {
  public:
    // Destructor
    virtual ~FuncOneDim() {}

    // Public methods

    /**
     * If this method returns true, then 'eval' will return valid deriv.
     * values via 'df', otherwise the return of 'df' is undefined.
     *
     * @return See above
     */
    virtual bool hasDerivative() const = 0;

    /**
     * Returns function value and derivative at a point 'x' in the domain. If
     * 'x' is not in the domain, a 'InvalidParameterException' should be
     * thrown.
     * NOTE: The deriv. is returned only if 'hasDerivative' returns true,
     * otherwise 'df' is not touched.
     *
     * @param x   Eval. point
     * @param f   Function value ret. here
     * @param df  Derivative ret. here (but see above!)
     */
    virtual void eval(double x,double* f,double* df) = 0;

    /**
     * In additional to what 'eval' does,
     *   f(x) / f'(x)
     * is returned in '*rat'. This is done only if 'hasDerivative' returns
     * true.
     * NOTE: Default implementation just calls 'eval'. Overwrite this if
     * stability issues are likely!
     *
     * @param x   Eval. point
     * @param f   Function value ret. here
     * @param df  Derivative ret. here (but see above!)
     * @param rat S.a.
     */
    virtual void evalStable(double x,double* f,double* df,double* rat) {
      eval(x,f,df);
      if (hasDerivative()) *rat=(*f)/(*df);
    }
  };
//ENDNS

#endif
