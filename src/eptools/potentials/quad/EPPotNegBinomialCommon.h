/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPPotNegBinomialCommon
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTNEGBINOMIALCOMMON_H
#define EPTOOLS_EPPOTNEGBINOMIALCOMMON_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadraturePotential.h"
#include "src/eptools/potentials/SpecfunServices.h"

//BEGINNS(eptools)
  /**
   * Some common code for negative binomial potential classes with
   * different rate functions lam(s):
   *   t(s) = C (1-p(s))^r p(s)^y,   y in N, r>0,
   *   p(s) = lam(s)/(r + lam(s)),
   *   C    = Gamma(r+y)/(Gamma(y+1) Gamma(r))
   * Parameters: y (nonneg. int.), r (positive)
   * <p>
   * ATTENTION: If 'SpecfunServices::logGamma' is not implemented, we
   * use C=1 as constant in front.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotNegBinomialCommon : public virtual QuadraturePotential
  {
  protected:
    // Members

    double yscal,rscal;
    double logConst;    // log C(y,r) + r log r

  public:
    // Public methods

    /**
     * @param py Value for y
     * @param pr Value for r
     */
    EPPotNegBinomialCommon(double py,double pr) : yscal(0.0),rscal(1.0) {
      setY(py);
      setR(pr);
    }

    virtual double getY() const {
      return yscal;
    }

    virtual double getR() const {
      return rscal;
    }

    virtual void setY(double py) {
      double pv[2];

      pv[0]=py; pv[1]=rscal;
      setPars(pv);
    }

    virtual void setR(double pr) {
      double pv[2];

      pv[0]=yscal; pv[1]=pr;
      setPars(pv);
    }

    int numPars() const {
      return 2;
    }

    void getPars(double* pv) const {
      pv[0]=yscal; pv[1]=rscal;
    }

    void setPars(const double* pv) {
      if (!isValidPars(pv))
	throw InvalidParameterException(EXCEPT_MSG(""));
      yscal=pv[0]; rscal=pv[1];
      update();
    }

    bool isValidPars(const double* pv) const {
      int i=(int) ceil(pv[0]);

      return (i>=0 && ((double) i)==pv[0] && pv[1]>1e-12);
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
    void update() {
      // Recompute log C
      try {
	logConst=SpecfunServices::logGamma(rscal+yscal)-
	  SpecfunServices::logGamma(yscal+1.0)-
	  SpecfunServices::logGamma(rscal);
      } catch (NotImplemException ex) {
	logConst=0.0;
      }
      logConst+=rscal*log(rscal);
    }
  };
//ENDNS

#endif
