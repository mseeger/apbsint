//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: optimize
 * Desc.:  Header class OneDimSolver
 * ------------------------------------------------------------------- */

#ifndef ONEDIMSOLVER_H
#define ONEDIMSOLVER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * TODO:
 * - Solver which does not need derivatives
 */

#include "lhotse/optimize/default.h"
#include "lhotse/optimize/FuncOneDim.h"

//BEGINNS(optimize)
  /**
   * Collects static methods which numerically solve f(x)=0 for x from a
   * subset of \R. f(x) is encoded as a instance of a subclass of 'FuncOneDim'.
   * The bracket-based methods will usually guarantee that a solution is found
   * within a given initial bracket, up to a given accuracy, given that f(x)
   * is continuously differentiable, and f(l)*f(r)<0 for the initial bracket
   * [l,r].
   * <p>
   * Typically, subclasses of 'FuncOneDim' require the derivative to be
   * computed alongside with the function value. Methods which do not need this
   * information, will simply ignore the derivative value.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class OneDimSolver
  {
  public:
    // Constants

    static const int brackRightRegular =0;
    static const int brackRightBound   =1;
    static const int brackRightInfinite=2;

    // Public static methods

    /**
     * One-dimensional Newton solver. The function object f=='func' has to
     * return derivatives.
     * <p>
     * The method maintains a bracket [l,r], containing a zero, l=='l',
     * r=='r', l<r, and f(l)*f(r)<0. 'r' need not be passed (see below).
     * Once the bracket size shrinks below 'acc', one of the bracket ends is
     * returned as solution (the one where f has been evaluated last
     * recently). Also, if |f(x)| is below 'facc', x is returned.
     * The algorithm is started from l. At each iteration, a Newton step is
     * attempted from the last recent point (l or r). If this falls out of
     * [l,r], a bisection step is done. If a Newton step is taken, but the
     * bracket shrinks by less than 15%, the next step is a bisection.
     * <p>
     * Not passing the right bracket end r:
     * If 'brRight'!='brackRightRegular', no right bracket end has to be
     * given. In this case, if 'r' > 'l', the first step is taken to 'r',
     * otherwise to 'l'+'acc' (very conservative). We first try to
     * to find a right bracket end r, then continue as above.
     * NOTE: This search for r may fail for some continuously differentiable
     * functions. You are safe if the sign of f'(x) is constant for x>l and the
     * opposite of the sign of f(l), for example if f(x), x>l, is increasing,
     * and f(l)<0.
     *
     * If R>l is known s.t. f(x) becomes unbounded for x -> R from the left
     * (with a sign opposite of f(l)), use 'brRight'=='brackRightBound' pass R
     * via 'boundR'. We assume that f(x) can be eval. for x <= R-'acc'. Still,
     * we do not set r=R-'acc', trying to avoid evaluating f there. If f does
     * not change sign in [l,R-'acc'], the algorithm fails.
     * If no such R is known, pass 'brRight'=='brackRightInfinite'. 'boundR' is
     * ignored in this case.
     *
     * The first step is to x='r', or to x='l'+'acc' (see above). Say, f(l)<0.
     * We have l<x with f(l)<0, we are done if f(x)>0, or |f(x)|<'facc'.
     * Suppose that f(x)<0. We fit a quadratic q with q'(x) = f'(x),
     * q(x) = f(x), q(l) = f(l). If this is convex, we solve q(x) = 0 for the
     * next x, l becomes the old x. If q is concave, we use the line l with
     * l(x) = f(x), l'(x) = f'(x), and solve l(x) = 0 (a Newton step). Note
     * that the quadratic step (if done), is shorter than the Newton step.
     * In general, if R is given, the step to x is capped from the right, so
     * that R-'acc'-x does not shrink by a factor >9/10.
     * <p>
     * NOTE: If this method is used to minimize a convex function (passing its
     * derivative for f), it may be that only a right bracket end r is known.
     * In that case, pass f(-x) for 'func', and l=-r.
     * <p>
     * NOTE: Code similar to NR routine 'rtsafe', but corrected. Can be
     * improved upon (f.ex., by Brent's method).
     * <p>
     * f(l), f'(l), and f(l)/f'(l) may be passed via 'fl', 'df', 'rat' (all
     * or none). This saves the first eval. of f at l.
     *
     * @param func      Function f
     * @param l         Initial left bracket end
     * @param r         Initial right bracket end; or first step to find such
     *                  an end
     * @param acc       Accuracy in argument
     * @param facc      Accuracy in function value
     * @param brRight   See above. Def: 'brackRightRegular'
     * @param boundR    Required iff 'brRight'=='brackRightBound'. See above
     * @param fl        Value f(l). Optional
     * @param df        Value f'(l). Optional
     * @param rat       Value f(l)/f'(l). Optional
     * @param debName   Name for debug messages. Def.: 0 (no debug messages)
     *                  (only if HAVE_DEBUG defined)
     * @return          Solution
     */
    static double newton(FuncOneDim* func,double l,double r,double acc,
			 double facc,int brRight,double boundR,double fl,
			 double df,double rat,const char* debName=0);

    static double newton(FuncOneDim* func,double l,double r,double acc,
			 double facc,int brRight=brackRightRegular,
			 double boundR=0.0,const char* debName=0) {
      double fl,df,rat;

      if (!func->hasDerivative())
	throw InvalidParameterException("'func' must return derivatives!");
      // rat = fl/df
      func->evalStable(l,&fl,&df,&rat);

      return newton(func,l,r,acc,facc,brRight,boundR,fl,df,rat,debName);
    }
  };
//ENDNS

#endif
