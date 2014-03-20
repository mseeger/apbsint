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
   * Expectation propagation: Interface for potential t(.).
   * By default, t(s) depends on a scalar variable s, but see "argument
   * groups" below. Base class for EP update services
   * <p>
   * Fractional EP:
   * If 'suppFractional' returns true, fractional EP updates are
   * supported. In this case, t(.) is replaced by t(.)^eta, where
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
   * <p>
   * Argument groups:
   * Potential types (subclasses of this class) are grouped according to the
   * input and return arguments to 'compMoments' (major service). In general,
   * 'compMoments' receives an input vector 'inp' and returns a return vector
   * 'ret'. Structure and semantics of these vectors depends on the argument
   * group.
   *
   * Currently supported argument groups:
   * - atypeUnivariate: Standard univariate t(s), s scalar.
   * - atypeBivarPrec: Bivariate t(s,tau), s scalar, tau>0 precision parameter
   * See 'compMoments' comments for input/return vectors.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPScalarPotential : public virtual EPScalPotentialBase
  {
  public:
    // Potential argument groups
    static const int atypeUnivariate=0;
    static const int atypeBivarPrec =1;

    /**
     * See header comment.
     *
     * @return Do we support fractional EP?
     */
    virtual bool suppFractional() const {
      return false; // not supp. by def.
    }

    static int getArgumentGroup_static() {
      return atypeUnivariate; // Default
    }

    /**
     * @return Potential argument group
     */
    virtual int getArgumentGroup() const {
      return getArgumentGroup_static();
    }

    /**
     * Local EP update. Maps input vector 'inp' (cavity moments) to return
     * vector. Details depend on argument group.
     * <p>
     * Argument group 'atypeUnivariate' (default):
     * Given the cavity marginal N(s|mu{-},rho{-}), the tilted
     * distribution is
     *
     *   P_hat(s) = Z^{-1} t(s)^eta N(s|mu{-},rho{-})
     *
     * The input vector 'inp' is [mu{-},rho{-}].
     * Here, eta==1 by default, it can be in (0,1) if 'suppFractional'
     * returns true.
     * If hmu, hrho denote mean, variance of P_hat(s), this method returns
     * [alpha,nu] in 'ret', s.t.
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
     * <p>
     * Argument group 'atypeBivarPrec':
     * The potential is t(s,tau), tau>0 a precision parameter. Given cavity
     * marginals
     *   q{-}(s,tau) = N(s|mu{-},rho{-}), Gamma(tau|a{-},c{-})
     * where
     *   Gamma(tau|a,c) propto tau^{a-1} e^{-c tau}
     * is the Gamma familty, the tilted distribution is
     *
     *   P_hat(s,tau) = Z^{-1} t(s,tau)^eta q{-}(s,tau)
     *
     * The input vector 'inp' is [mu{-},rho{-},a{-},c{-}]. The return vector
     * 'ret' is [alpha,nu,a_hat,c_hat], where alpha, nu as above, and
     * Gamma(a_hat,c_hat) has same mean, variance as P_hat(tau).
     *
     * @param inp    Input vector (cavity marginal)
     * @param ret    Return vector
     * @param logz   log Z ret. here. Optional
     * @param eta    S.a. Def.: 1
     * @return       Success?
     */
    virtual bool compMoments(const double* inp,double* ret,double* logz=0,
			     double eta=1.0) const = 0;
  };
//ENDNS

#endif
