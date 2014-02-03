//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: optimize
 * Desc.:  Definition of class OneDimSolver
 * ------------------------------------------------------------------- */

#include "lhotse/optimize/OneDimSolver.h"

//BEGINNS(optimize)
  const int OneDimSolver::brackRightRegular ;
  const int OneDimSolver::brackRightBound   ;
  const int OneDimSolver::brackRightInfinite;

  // Public static methods

  /*
   * Losely related to NR routine 'rtsafe', which is incorrect. We correct the
   * mistakes and also allow for a missing right bracket end.
   */
#define MAXIT 100 // for 'newton'
#define MY_SIGN(x) (((x)>=0.0)?1:-1)
  double OneDimSolver::newton(FuncOneDim* func,double l,double r,double acc,
			      double facc,int brRight,double boundR,double fl,
			      double df,double rat,const char* debName)
  {
    int j,lsgn;
    double f,df2,rat2,dx,temp,rts,olds,alpha;
    bool nextBisect,didNewton,numerErr;

    if (!func->hasDerivative())
      throw InvalidParameterException("'func' must return derivatives!");
    lsgn=MY_SIGN(fl);
    if (fabs(fl)<facc) return l;
    if (brRight==brackRightRegular) {
      if (l>=r) throw InvalidParameterException(EXCEPT_MSG(""));
      func->eval(r,&f,&df2);
      if (fabs(f)<facc) return r;
      if (lsgn==MY_SIGN(f)) {
	cout << "l=" << l << ": fl=" << fl << endl
	     << "r=" << r << ": fr=" << f << endl;
	throw InvalidParameterException("Root must be bracketed in [l,r]");
      }
    } else {
      // Find right bracket end r (see header comment)
      // At the end, r > l, f(r) and f(l) have opposite signs. l may move to
      // the right as well, but does not change sign. 'fl', 'df', 'rat'
      // contain f(l) and f'(l).
      bool isBound=(brRight==brackRightBound);
      j=0;
#ifdef HAVE_DEBUG
      if (debName!=0)
	cout << debName << ": Finding bracket: l=" << l << ",fl=" << fl
	     << ",dl=" << df << endl;
#endif
      if ((dx=r-l)<=0.0) dx=acc; // First step
      if (isBound && l+dx>boundR-acc)
	throw InvalidParameterException("Initial step violates 'boundR'");
      for (;;) {
	if (j++>MAXIT)
	  throw NumericalException("OneDimSolver::newton failed: Maximum number of iterations exceeded");
	do {
	  rts=l+dx;
	  numerErr=false;
	  try {
	    func->evalStable(rts,&f,&df2,&rat2);
	  } catch (...) {
	    if (dx==acc)
	      throw NumericalException("OneDimSolver::newton failed: Cannot find right bracket end!");
	    dx=acc; numerErr=true; // try again with step size 'acc'
	  }
	} while (numerErr);
#ifdef HAVE_DEBUG
	if (debName!=0)
	  cout << debName << ":   rts=" << rts << ",f=" << f << ",df=" << df2
	       << endl;
#endif
	if (fabs(f)<facc) return rts;
	else if (MY_SIGN(f)!=lsgn) break; // OK, found right bracket end
	// Try quadratic or linear (Newton) step. If neither applies, just use
	// 'dx' again
	alpha=(fl-f)/(l-rts)-df2;
	if (MY_SIGN(alpha)==lsgn && fabs(alpha)>10.0*facc*(rts-l)) {
	  // OK, use quadratic approx.
	  alpha/=(l-rts);
	  dx=(lsgn==-1)?(0.5*(sqrt(df2*df2-4.0*alpha*f)-df2)/alpha):
	    (-0.5*(sqrt(df2*df2-4.0*alpha*f)+df2)/alpha);
	} else if (rat2<0.0)
	  // Use linear approx. (Newton step)
	  dx=-rat2;
	if (isBound && dx>(temp=0.9*(boundR-acc-rts)))
	  dx=temp; // Cap the step
	if (dx<acc) {
	  dx=acc; // Minimum step size
	  if (isBound && rts+dx>boundR-acc)
	    throw NumericalException("OneDimSolver::newton failed: Cannot find right bracket end!");
	}
	l=rts; // new left bracket end
	fl=f; df=df2; rat=rat2;
#ifdef HAVE_DEBUG
	if (debName!=0)
	  cout << debName << ":   l=" << l << ",fl=" << fl << ",dl=" << df
	       << endl;
#endif
      }
      r=rts; // right bracket end detected
    }

    // From here on, we have a bracket [l,r].
    // Start with rts==l (value, derivs. in fl, df, rat), Newton step
    // NOTE: The function values at l, r are not maintained (not needed).
    if ((olds=r-l)<acc) return l;
    if (MY_SIGN(fl)!=lsgn) throw InternalException(EXCEPT_MSG("")); // Sanity
    rts=l; f=fl;
    nextBisect=false; // Start with Newton step
    for (j=0; j<=MAXIT; j++) {
      // We are at rts, which may be equal to l or r (bracket ends).
      // Function value and deriv. are in f, df, rat.
      // We choose a bisection step if either 'nextBisect' is true or the
      // Newton step falls out of the bracket.
      // If we take the Newton step and find out that by doing so, the
      // bracket shrunk by a fraction less than 0.15, we set 'nextBisect'
      // to true, which leads to bisection in the next turn.
#ifdef HAVE_DEBUG
      if (debName!=0)
	cout << debName << ": [l=" << l << ",r=" << r << "]" << endl;
#endif
      temp=rts-rat;
      if (nextBisect || temp<=l || temp>=r) {
	// Bisection step
	rts=0.5*(l+r);
	didNewton=false;
#ifdef HAVE_DEBUG
	if (debName!=0)
	  cout << debName << ":   Bisect: rts=" << rts << endl;
#endif
      } else {
	// Newton step
	rts=temp;
	didNewton=true;
#ifdef HAVE_DEBUG
	if (debName!=0)
	  cout << debName << ":   Newton: rts=" << rts << endl;
#endif
      }
      func->evalStable(rts,&f,&df,&rat);
      if (fabs(f)<facc) return rts;
      else if (MY_SIGN(f)==lsgn) l=rts;
      else r=rts;
      if ((temp=r-l)<acc) return rts;
#ifdef HAVE_DEBUG
      if (debName!=0)
	cout << debName << ":   f(rts)=" << f << ",df=" << df << endl;
#endif
      nextBisect=didNewton && (temp>0.85*olds); // Bisection step next?
      olds=temp;
    }
    throw NumericalException("OneDimSolver::newton failed: Maximum number of iterations exceeded");
  }
#undef MY_SIGN
#undef MAXIT
//ENDNS
