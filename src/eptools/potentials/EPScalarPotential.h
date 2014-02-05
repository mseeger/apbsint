/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPScalarPotential
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPSCALARPOTENTIAL_H
#define EPTOOLS_EPSCALARPOTENTIAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalPotentialBase.h"

//BEGINNS(eptools)
  /**
   * Expectation propagation: Interface for potential t(s), s a scalar
   * variable. This is the base class for EP update services
   * <p>
   * Fractional EP:
   * If 'suppFractional' returns true, fractional EP updates are
   * supported. In this case, t(s) is replaced by t(s)^eta, where
   * eta in (0,1] is a parameter passed to 'compMoments' or other
   * methods. Cavity moments must be computed accordingly.
   * <p>
   * NOTE: More general, we could compute EP updates without determining
   * the cavity marginal. For certain potentials, this can be done even
   * if the cavity marginal is not defined, and it could be numerically
   * more robust.
   * <p>
   * Default constructor:
   * A subclass must also provide a default constructor, which uses some
   * valid set of parameters. This is called by
   * 'EPPotentialFactory::createDefault'. This default constructor must only
   * depend on construction parameters (see below): typical subclasses do not
   * have such, and their default constructor is parameter-free.
   * ATTENTION: This call has to be done explicitly in 'EPPotentialFactory',
   * since a constructor cannot be virtual.
   * This is required to break a chicken-and-egg situation: We need an
   * existing object in order to check a parameter setting via 'isValidPars'.
   * <p>
   * Annotations:
   * A 'EPScalarPotential' type can be annotated (see 'EPPotentialFactory',
   * 'PotManagerFactory'), meaning that it contains a member which is not
   * created on-demand (and which may be shared).
   * For an annotated type, all constructors must accept a void*, which is
   * cast to the annotation type (depends on 'EPScalarPotential' subclass).
   * In particular, if this argument is 0, an exception with a meaningful
   * text must be raised.
   * ATTENTION: This mechanism is unsafe. A non-zero value cannot be properly
   * checked for validity.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPScalarPotential : public virtual EPScalPotentialBase
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
     * Given the cavity marginal N(s|mu{-},rho{-}), the tilted
     * distribution is
     *
     *   P_hat(s) = Z^{-1} t(s)^eta N(s|mu{-},rho{-})
     *
     * Here, eta==1 by default, it can be in (0,1) if 'suppFractional'
     * returns true.
     * If hmu, hrho denote mean, variance of P_hat(s), this method
     * returns log Z, alpha and nu, s.t.
     *
     *   hmu  = mu{-} + alpha rho{-},
     *   hrho = rho{-} (1 - nu rho{-})
     *
     * The computation can fail for numerical reasons, in which case false
     * is returned (also if eta<1 and 'suppFractional' returns false).
     * New values for EP parameters beta, pi can be computed as:
     *
     *   pi'   = nu/(1 - nu rho{-}) + (1-eta) pi,
     *   beta' = (nu mu{-} + alpha)/(1 - nu rho{-}) + (1-eta) beta
     *
     * @param cmu    Cavity marginal mean
     * @param crho   Cavity marginal variance
     * @param alpha  Value alpha ret. here
     * @param nu     Value nu ret. here
     * @param logz   log Z ret. here. Optional
     * @param eta    S.a. Def.: 1
     * @return       Success?
     */
    virtual bool compMoments(double cmu,double crho,double& alpha,double& nu,
			     double* logz=0,double eta=1.0) const = 0;
  };
//ENDNS

#endif
