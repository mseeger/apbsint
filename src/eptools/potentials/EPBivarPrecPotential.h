/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPBivarPrecPotential
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPBIVARPRECPOTENTIAL_H
#define EPTOOLS_EPBIVARPRECPOTENTIAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalPotentialBase.h"

//BEGINNS(eptools)
  /**
   * Expectation propagation: Interface for potential t(s,tau), where
   * s, tau are scalars, and tau>0 is a precision variable (inverse
   * variance).
   * <p>
   * Fractional EP:
   * If 'suppFractional' returns true, fractional EP updates are
   * supported. In this case, t(s,tau) is replaced by t(s,tau)^eta,
   * where eta in (0,1] is a parameter passed to 'compMoments' or other
   * methods. Cavity moments must be computed accordingly.
   * <p>
   * Default constructor, annotations: This is similar to 'EPScalarPotential',
   * see comments there.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPBivarPrecPotential : public virtual EPScalPotentialBase
  {
  public:
    /**
     * See header comment.
     *
     * @return Do we support fractional EP?
     */
    virtual bool suppFractional() const {
      return false; // not supp. by def.
    }

    /**
     * Given the cavity marginals
     *   q{-}(s,tau) = N(s|mu{-},rho{-}), Gamma(tau|a{-},c{-})
     * where
     *   Gamma(tau|a,c) propto tau^{a-1} e^{-c tau}
     * is the Gamma familty, the tilted distribution is
     *
     *   P_hat(s,tau) = Z^{-1} t(s,tau)^eta q{-}(s,tau)
     *
     * Here, eta==1 by default, it can be in (0,1) if 'suppFractional'
     * returns true.
     * If hmu, hrho denote mean, variance of P_hat(s), this method
     * returns log Z, alpha, nu, s.t.
     *
     *   hmu  = mu{-} + alpha rho{-},
     *   hrho = rho{-} (1 - nu rho{-})
     *
     * It also returns a_hat, c_hat s.t. Gamma(tau|a_hat,c_hat) has the
     * same mean, variance as P_hat(tau).
     * The computation can fail for numerical reasons, in which case false
     * is returned (also if eta<1 and 'suppFractional' returns false).
     * New values for EP parameters beta, pi (corr. to s) can be computed
     * as:
     *
     *   pi'   = nu/(1 - nu rho{-}) + (1-eta) pi,
     *   beta' = (nu mu{-} + alpha)/(1 - nu rho{-}) + (1-eta) beta
     *
     * @param cmu   Mean of cavity marginal q{-}(s)
     * @param crho  Variance of cavity marginal q{-}(s)
     * @param ca    Shape parameter of cavity marginal q{-}(tau)
     * @param cc    Scale parameter of cavity marginal q{-}(tau)
     * @param alpha Result for s posterior
     * @param nu    "
     * @param hata  Result for tau posterior
     * @param hatc  "
     * @param logz  log Z ret. here. Optional
     * @param eta   S.a. Def.: 1
     * @return      Success?
     */
    virtual bool compMoments(double cmu,double crho,double ca,double cc,
			     double& alpha,double& nu,double& hata,double& hatc,
			     double* logz=0,double eta=1.0) const = 0;
  };
//ENDNS

#endif
