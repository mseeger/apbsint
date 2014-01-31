/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotQuantileRegress
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTQUANTILEREGRESS_H
#define EPTOOLS_EPPOTQUANTILEREGRESS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"
#include "src/eptools/potentials/SpecfunServices.h"

//BEGINNS(eptools)
  /**
   * Quantile regression potential:
   *   t(s)  = tt(xi (y - s)),
   *   tt(r) = exp(-kappa [r]_+ - (1-kappa) [-r]_+),
   *   [r]_+ = x I{x>0}
   * Parameters: y, xi>0, kappa in (0,1).
   * <p>
   * NOTE: The Laplace potential (see 'EPPotLaplace') is a special
   * case, which is implemented calling a static method here.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotQuantileRegress : public EPScalarPotential
  {
  protected:
    // Members

    double yscal,xi,kappa;

  public:
    // Public methods

    EPPotQuantileRegress(double py=0.0,double pxi=1.0,double pkappa=0.5) :
      kappa(pkappa) {
      if (kappa<=0.0 || kappa>=1.0)
	throw InvalidParameterException(EXCEPT_MSG(""));
       setY(py); setXi(pxi);
    }

    virtual double getXi() const {
      return xi;
    }

    virtual void setXi(double pxi) {
      if (pxi<1e-12)
	throw InvalidParameterException(EXCEPT_MSG(""));
      xi=pxi;
    }

    virtual double getY() const {
      return yscal;
    }

    virtual void setY(double py) {
      yscal=py;
    }

    virtual double getKappa() const {
      return kappa;
    }

    int numPars() const {
      return 3;
    }

    void getPars(double* pv) const {
      pv[0]=getY(); pv[1]=getXi(); pv[2]=kappa;
    }

    void setPars(const double* pv) {
      if (pv[2]<=0.0 || pv[2]>=1.0)
	throw InvalidParameterException(EXCEPT_MSG(""));
      setY(pv[0]); setXi(pv[1]); kappa=pv[2];
    }

    bool isValidPars(const double* pv) const {
      return (pv[1]>=1e-12 && pv[2]>0.0 && pv[2]<1.0);
    }

    bool suppFractional() const {
      return true;
    }

    bool isLogConcave() const {
      return true;
    }

    bool compMoments(double cmu,double crho,double& alpha,double& nu,
		     double* logz=0,double eta=1.0) const {
      if (crho<1e-14 || eta<1e-10 || eta>1.0)
	throw InvalidParameterException(EXCEPT_MSG(""));

      return compMomentsInt(cmu,crho,xi*eta,yscal,kappa,alpha,nu,logz);
    }

    /**
     * Implements 'compMoments'. Called by 'EPPotLaplace' as well.
     * NOTE: No fractional parameter 'eta'. Multiply this parameter into xi.
     */
    static bool compMomentsInt(double cmu,double crho,double xi,double yscal,
			       double kappa,double& alpha,double& nu,
			       double* logz);
  };

  inline bool
  EPPotQuantileRegress::compMomentsInt(double cmu,double crho,double xi,
				       double yscal,double kappa,double& alpha,
				       double& nu,double* logz)
  {
    double temp,kapc,hh,hr,rhor,sqrhor,argf,li01,li02,logi0,q;

    kapc = 1.0-kappa;
    hh = yscal-cmu;
    hr = xi*hh; rhor=xi*xi*crho; // Moments h_r, rho_r
    sqrhor = xi*sqrt(crho);
    argf = kappa*sqrhor-hr/sqrhor;
    // 'SpecfunServices:logCdfNormal' implements F(-x)
    li01 = 0.5*kappa*(kappa*rhor-2*hr) + SpecfunServices::logCdfNormal(-argf);
    li02 = 0.5*kapc*(kapc*rhor+2*hr) +
      SpecfunServices::logCdfNormal(argf-sqrhor);
    if (li01>=li02) {
      temp = exp(li02-li01);
      logi0 = li01+log1p(temp);
      q = temp/(1.0+temp);
    } else {
      temp = exp(li01-li02);
      logi0 = li02+log1p(temp);
      q = 1.0/(1.0+temp);
    }
    if (logz!=0) *logz=logi0;
    alpha = xi*(kappa-q);
    nu = xi*xi*(exp(-0.5*(hh*hh/crho+SpecfunServices::m_ln2pi)-logi0)/sqrhor
		- q*(1.0-q));

    return true;
  }
//ENDNS

#endif
