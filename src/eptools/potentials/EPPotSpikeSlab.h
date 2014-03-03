/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotSpikeSlab
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTSPIKESLAB_H
#define EPTOOLS_EPPOTSPIKESLAB_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"
#include "src/eptools/potentials/SpecfunServices.h"

//BEGINNS(eptools)
  /**
   * Basic spike and slab potential (Gaussian slab):
   *   t(s) = (1-p) delta_0(s) + p N(s|0,v),   v>0, c=log(p/(1-p)).
   * Parameters: c, v>0.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotSpikeSlab : public EPScalarPotential
  {
  protected:
    // Members

    double lpscal,vscal; // c = log(p/(1-p)), v>0

  public:
    // Public methods

    /**
     * Constructor
     *
     * @param plp Logit log(p/(1-p)). Def.: 0
     * @param pv  Slab variance. Def.: 1
     */
    EPPotSpikeSlab(double plp=0.0,double pv=1.0) {
      setVariance(pv);
      setLogitP(plp);
    }

    int numPars() const {
      return 2;
    }

    virtual double getVariance() const {
      return vscal;
    }

    virtual void setVariance(double pv) {
      if (pv<(1e-12))
	throw InvalidParameterException(EXCEPT_MSG(""));
      vscal=pv;
    }

    virtual double getLogitP() const {
      return lpscal;
    }

    virtual void setLogitP(double plp) {
      lpscal=plp;
    }

    void getPars(double* pv) const {
      pv[0]=lpscal; pv[1]=vscal;
    }

    void setPars(const double* pv) {
      if (!isValidPars(pv))
	throw InvalidParameterException(EXCEPT_MSG(""));
      setVariance(pv[1]);
      setLogitP(pv[0]);
    }

    bool isValidPars(const double* pv) const {
      return (pv[1]>=(1e-12));
    }

    bool suppFractional() const {
      return false;
    }

    bool isLogConcave() const {
      return false;
    }

    bool compMoments(const double* inp,double* ret,double* logz=0,
		     double eta=1.0) const;

  protected:
    // Internal methods

    /**
     * Does job of 'compMoments', but based on the unnormalized cavity
     * marginal exp( beta{-} s - 0.5 pi{-} s^2) instead of N(s|mu{-},rho{-}).
     * This means that the value log Z_hat is based on
     *   Z_hat = int t(s) exp( beta{-} s - 0.5 pi{-} s^2) d s,
     * it has to be corrected in 'compMoments'.
     * The method formally works as long as
     *   1 + pi{-} v >= 1e-16.
     * In contrast, 'compMoments' requires pi{-} to be bounded away from 0.
     *
     * @param cbeta Cavity parameter beta{-}
     * @param cpi   Cavity parameter pi{-}
     * @param alpha Value alpha ret. here
     * @param nu    Value nu ret. here
     * @param logzh Value log Z_hat ret. here. Optional
     * @return      Success?
     */
    bool compMomentsInt(double cbeta,double cpi,double& alpha,double& nu,
			double* logzh) const;
  };

  inline bool
  EPPotSpikeSlab::compMoments(const double* inp,double* ret,double* logz,
			      double eta) const
  {
    double cpi,cbeta,cmu=inp[0],crho=inp[1];
    bool rstat;

    if (eta!=1.0)
      throw NotImplemException("Fractional updates not implemented");
    if (crho<(1e-16))
      throw NumericalException(EXCEPT_MSG(""));
    cpi=1.0/crho; cbeta=cmu/crho;
    if (rstat=compMomentsInt(cbeta,cpi,ret[0],ret[1],logz)) {
      if (logz!=0)
	*logz-=0.5*(cbeta*cmu+log(crho)+SpecfunServices::m_ln2pi);
    }

    return rstat;
  }

  /*
   * Adapted from 'EPPotSpikeSlab' in module 'epscal/potentials', but using
   * the simplification implemented in 'EPPotGaussMixture'.
   */
  inline bool
  EPPotSpikeSlab::compMomentsInt(double cbeta,double cpi,double& alpha,
				 double& nu,double* logzh) const
  {
    double temp,temp2,bmsq,rho2,r2,z2m1;

    // Natural parameters of "cavity" (may be undefined)
    if (1.0+cpi*vscal<(1e-16))
      return false;
    bmsq=cbeta*cbeta;
    rho2=vscal/(1.0+cpi*vscal);
    temp=lpscal+0.5*(rho2*bmsq-log1p(cpi*vscal)); // log(Z_2/(1-p))
    // r_2 = Z_2/Z (note that Z_1 = 1-p):
    r2=1.0/(1.0+(temp2=exp(-temp)));
    z2m1=-rho2*cpi; // z_2 - 1
    if (logzh!=0) {
      *logzh=log1p(temp2)+temp-log1p(exp(lpscal)); // log Z_hat
      //cout << "temp=" << temp << ",logz=" << *logzh << endl;
    }
    temp=1.0+r2*z2m1; // A_til
    alpha=-cbeta*temp;
    nu=temp*cpi-bmsq*r2*(1.0-r2)*z2m1*z2m1;

    return true;
  }
//ENDNS

#endif
