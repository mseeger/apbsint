/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPPotPoissonCommon
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTPOISSONCOMMON_H
#define EPTOOLS_EPPOTPOISSONCOMMON_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadraturePotential.h"
#include "src/eptools/potentials/SpecfunServices.h"

//BEGINNS(eptools)
  /**
   * Some common code for Poisson potential classes with different rate
   * functions lam(s):
   *   t(s) = (y!)^-1 lam(s)^y exp(-lam(s)),  y in N,
   * Parameters: y (nonneg. int.)
   * <p>
   * ATTENTION: If 'SpecfunServices::logGamma' is not implemented, we drop
   * the constant (y!)^-1 in front.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotPoissonCommon : public virtual QuadraturePotential
  {
  protected:
    // Members

    double yscal;
    double logYFact; // log(y!)

  public:
    // Public methods

    /**
     * @param py    Value for y
     */
    EPPotPoissonCommon(double py) {
      setY(py);
    }

    virtual double getY() const {
      return yscal;
    }

    virtual void setY(double py) {
      if (!isValidPars(&py))
	throw InvalidParameterException(EXCEPT_MSG(""));
      yscal=py;
      setLogYFact();
    }

    int numPars() const {
      return 1;
    }

    void getPars(double* pv) const {
      pv[0]=yscal;
    }

    void setPars(const double* pv) {
      setY(pv[0]);
    }

    bool isValidPars(const double* pv) const {
      int i=(int) ceil(pv[0]);

      return (i>=0 && ((double) i)==pv[0]);
    }

    bool hasWayPoints() const {
      return true;
    }

    /**
     * The integration interval is all of R, and l(s) is smooth everywhere.
     * NOTE: Assumes that the log rate function is smooth.
     */
    void getInterval(double& a,bool& aInf,double& b,bool& bInf,
		     ArrayHandle<double>& wayPts) const {
      aInf=bInf=true;
      wayPts.changeRep(0); // no way points
    }

  protected:
    /**
     * Computes log(y!) = log Gamma(y+1), stores in 'logYFact'.
     * ATTENTION: If 'SpecfunServices::logGamma' is not implemented, we
     * set 'logYFact' to 0.
     */
    void setLogYFact() {
      try {
	logYFact=SpecfunServices::logGamma(yscal+1.0);
      } catch (NotImplemException ex) {
	logYFact=0.0;
      }
    }
  };
//ENDNS

#endif
