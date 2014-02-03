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

//BEGINNS(eptools)
  /**
   * Some common code for Poisson potential classes with different rate
   * functions lam(s):
   *   t(s) = lam(s)^y exp(-lam(s)),  y in N,
   * Parameters: y (nonneg. int.)
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotPoissonCommon : public virtual QuadraturePotential
  {
  protected:
    // Members

    double yscal;

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

      return (i>=0 && ((double) i)==py);
    }

    bool hasWayPoints() const {
      return true;
    }

    /**
     * The integration interval is all of R, and l(s) is smooth everywhere.
     */
    void getInterval(double& a,bool& aInf,double& b,bool& bInf,
		     ArrayHandle<double>& wayPts) const {
      aInf=bInf=true;
      wayPts.changeRep(0); // no way points
    }
  };
//ENDNS

#endif
