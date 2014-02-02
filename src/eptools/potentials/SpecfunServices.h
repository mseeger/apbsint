/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class SpecfunServices
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_SPECFUNSERVICES_H
#define EPTOOLS_SPECFUNSERVICES_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/default.h"

//BEGINNS(eptools)
#define ERF_CODY_LIMIT1 0.6629
#define ERF_CODY_LIMIT2 5.6569
  /**
   * Collects static methods for computing certain special functions.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SpecfunServices
  {
  public:
    // Constants

    static const double m_ln2pi  = 1.83787706640934533908193770913;
    static const double m_ln2    = 0.69314718055994530941723212146;
    static const double m_sqrtpi = 1.77245385090551602729816748334;
    static const double m_sqrt2  = 1.41421356237309504880168872421;

    // Static methods

    /**
     * @param z Argument
     * @return  log N(z|0,1)
     */
    static double logPdfNormal(double z);

    /**
     * @param z Argument
     * @return  Phi(z), c.d.f. of N(0,1)
     */
    static double cdfNormal(double z);

    /**
     * If Phi(z) denotes the c.d.f. of N(0,1), this method computes
     * log Phi(z).
     * NOTE: The technical report defines
     *   F(x) = log(1 - Phi(x)).
     * This method computes F(-z).
     *
     * @param z Argument
     * @return  log Phi(z)
     */
    static double logCdfNormal(double z);

    /**
     * If Phi(z) denotes the c.d.f. of N(0,1), this method computes
     *   f(z) = (d/dz) log Phi(z) = N(z)/Phi(z).
     * NOTE: The technical report defines the hazard function
     *   h(x) = N(x)/(1 - Phi(x)).
     * This method computes h(-z).
     *
     * @param z Argument
     * @return  (d/dz) log Phi(z)
     */
    static double derivLogCdfNormal(double z);

    /**
     * Computes natural log of Gamma(z) for z>0. Note that if z is a
     * natural number, then z! = Gamma(z+1).
     *
     * @param z Argument (positive)
     * @return  log Gamma(z)
     */
    static double logGamma(double z);

  protected:
    // Internal static methods

    /**
     * For x >= ERF_CODY_LIMIT1, define Q(x) by
     *   1 - Phi(x) approx N(x) x^{-1} Q(x).
     * We compute Q(x) according to
     *   Cody
     *   Rational Chebyshev approximation to the error function
     * This is done differently for x >= ERF_CODY_LIMIT2 and
     * ERF_CODY_LIMIT1 <= x < ERF_CODY_LIMIT2.
     * NOTE: Q(x) -> 1 for x->infty.
     */
    static double erfRationalHelper(double x);

    /**
     * Implements rational function R_3(y),  y = x^2/2,
     * which is used if 0 <= x < ERF_CODY_LIMIT1. In this range:
     *   Phi(x) approx (1 + (x/sqrt(2)) R_3(x^2/2))/2
     * See
     *   Cody
     *   Rational Chebyshev approximation to the error function
     */
    static double erfRationalHelperR3(double y);
  };

  // Inline methods

  inline double SpecfunServices::logPdfNormal(double z)
  {
    return -0.5*(m_ln2pi+z*z);
  }

  inline double SpecfunServices::cdfNormal(double x)
  {
    double res;

    if (fabs(z)<ERF_CODY_LIMIT1) {
      // Part 3 approximation:
      // Phi(z) approx (1 + y R_3(y^2))/2, y = z/sqrt(2)
      res=0.5*(1.0 + (z/M_SQRT2)*erfRationalHelperR3(0.5*z*z));
    } else {
      // Part 1 or 2 approximation:
      // Phi(z) approx N(z) Q(-z)/(-z), z<0
      // NOTE: The case z >= ERF_CODY_LIMIT1 is uncritical, we could even use
      // a cheaper approximation then
      if (z<0.0)
	res=pdfNormal(z)*erfRationalHelper(-z)/(-z);
      else
	res=1.0-pdfNormal(z)*erfRationalHelper(z)/z;
    }

    return res;
  }

  inline double SpecfunServices::logCdfNormal(double z)
  {
    double res;

    if (fabs(z)<ERF_CODY_LIMIT1) {
      // Part 3 approximation:
      // Phi(z) approx (1 + y R_3(y^2))/2, y = z/sqrt(2)
      res=log1p((z/m_sqrt2)*erfRationalHelperR3(0.5*z*z))-m_ln2;
    } else {
      // Part 1 or 2 approximation:
      // Phi(z) approx N(z) Q(-z)/(-z), z<0
      // NOTE: The case z >= ERF_CODY_LIMIT1 is uncritical, we could even use
      // a cheaper approximation then
      if (z<0.0)
	res=logPdfNormal(z)-log(-z)+log(erfRationalHelper(-z));
      else
	res=log1p(-exp(logPdfNormal(z))*erfRationalHelper(z)/z);
    }

    return res;
  }

  inline double SpecfunServices::derivLogCdfNormal(double z)
  {
    double res;

    if (fabs(z)<ERF_CODY_LIMIT1) {
      // Part 3 approximation:
      // Phi(z) approx (1 + y R_3(y^2))/2, y = z/sqrt(2)
      res=2.0*exp(logPdfNormal(z))/
	(1.0 + (z/m_sqrt2)*erfRationalHelperR3(0.5*z*z));
    } else {
      // Part 1 or 2:
      // Phi(z) approx N(z) Q(-z)/(-z), z<0
      // NOTE: The case z >= ERF_CODY_LIMIT1 is uncritical, we could even use
      // a cheaper approximation then
      if (z<0.0)
	res=-z/erfRationalHelper(-z);
      else {
	double temp=exp(logPdfNormal(z));
	res=temp/(1.0-temp*erfRationalHelper(z)/z);
      }
    }

    return res;
  }

  inline double SpecfunServices::erfRationalHelper(double x)
  {
    int i;
    double res,den,y;

    MYASS(x>0.0);
    if (x>=ERF_CODY_LIMIT2) {
      // x/sqrt(2) >= 4
      // Q(x)   = 1 + sqrt(pi) y R_1(y),
      // R_1(y) = poly(p_j,y) / poly(q_j,y),   y = 2/x^2
      // Ordering of arrays: 4,3,2,1,0,5 (only for numerator p_j; q_5=1)
      // ATTENTION: The p_j are negative of the entries here
      double p[]={3.05326634961232344e-1, 3.60344899949804439e-1,
		  1.25781726111229246e-1, 1.60837851487422766e-2,
		  6.58749161529837803e-4, 1.63153871373020978e-2};
      double q[]={2.56852019228982242,    1.87295284992346047,
		  5.27905102951428412e-1, 6.05183413124413191e-2,
		  2.33520497626869185e-3};
      y=2.0/x/x;
      res=y*p[5]; den=y;
      for (i=0; i<4; i++) {
	res=(res+p[i])*y; den=(den+q[i])*y;
      }
      // Minus, because p[j] values have to be negated
      res=1.0-m_sqrtpi*y*(res+p[4])/(den+q[4]);
    } else {
      // x/sqrt(2) < 4, x/sqrt(2) >= 0.469
      // Q(x)   = sqrt(pi) y R_2(y),
      // R_2(y) = poly(p_j,y) / poly(q_j,y),   y = x/sqrt(2)
      // Ordering of arrays: 7,6,5,4,3,2,1,0,8 (only p_8; q_8=1)
      double p[]={5.64188496988670089e-1, 8.88314979438837594,
		  6.61191906371416295e+1, 2.98635138197400131e+2,
		  8.81952221241769090e+2, 1.71204761263407058e+3,
		  2.05107837782607147e+3, 1.23033935479799725e+3,
		  2.15311535474403846e-8};
      double q[]={1.57449261107098347e+1, 1.17693950891312499e+2,
		  5.37181101862009858e+2, 1.62138957456669019e+3,
		  3.29079923573345963e+3, 4.36261909014324716e+3,
		  3.43936767414372164e+3, 1.23033935480374942e+3};
      y=x/m_sqrt2;
      res=y*p[8]; den=y;
      for (i=0; i<7; i++) {
	res=(res+p[i])*y; den=(den+q[i])*y;
      }
      res=m_sqrtpi*y*(res+p[7])/(den+q[7]);
    }

    return res;
  }

  inline double SpecfunServices::erfRationalHelperR3(double y)
  {
    int i;
    double nom,den;

    MYASS(y>=0.0);
    // R_3(y) = poly(p_j,y) / poly(q_j,y)
    // Ordering of arrays: 3,2,1,0,4 (only for p_5; q_5=1)
    double p[] = {3.16112374387056560,    1.13864154151050156e+2,
		  3.77485237685302021e+2, 3.20937758913846947e+3,
		  1.85777706184603153e-1};
    double q[] = {2.36012909523441209e+1, 2.44024637934444173e+2,
		  1.28261652607737228e+3, 2.84423683343917062e+3};
    nom=y*p[4]; den=y;
    for (i=0; i<3; i++) {
      nom=(nom+p[i])*y; den=(den+q[i])*y;
    }

    return (nom+p[3])/(den+q[3]);
  }

#ifndef HAVE_LIBGSL
#include "src/eptools/potentials/SpecfunServices_basic.h"
#else
#include "src/eptools/potentials/SpecfunServices_workaround.h"
#endif

#undef ERF_CODY_LIMIT2
#undef ERF_CODY_LIMIT1
//ENDNS

#endif
