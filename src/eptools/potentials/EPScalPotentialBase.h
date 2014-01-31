/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPScalPotentialBase
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPSCALPOTENTIALBASE_H
#define EPTOOLS_EPSCALPOTENTIALBASE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/default.h"

//BEGINNS(eptools)
  /**
   * Expectation propagation: Interface for potential t(s), s a scalar
   * variable.
   * <p>
   * This base class defines basic services to be implemented for every
   * potential. Higher-order service, in particular pertaining to EP
   * updates, are delegated to 'EPScalarPotential'.
   * <p>
   * A potential is configured by 'numPars' parameters (type double),
   * this can be zero.
   * <p>
   * Construction parameters:
   * Some subclasses require parameters in the default constructor, and also
   * to determine the number of parameters in 'numberPars'. These parameters
   * are construction parameters. They must form the prefix of the parameter
   * vector. 'numberConstPars' returns the number of construction parameters
   * (def. implementation: 0).
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPScalPotentialBase
  {
  public:
    // Public methods

    virtual ~EPScalPotentialBase() {}

    /**
     * @return Number of parameters (can be zero)
     */
    virtual int numPars() const = 0;

    /**
     * See header comment.
     *
     * @return Number of construction parameters
     */
    virtual int numConstPars() const {
      return 0; // Default implementation: No construction parameters
    }

    /**
     * Parameter vector (of size 'numPars') written to 'pv'.
     *
     * @param pv Parameter vector written here
     */
    virtual void getPars(double* pv) const = 0;

    /**
     * Restrictions may apply on different values. A
     * 'InvalidParameterException' is thrown if 'pv' violates these.
     *
     * @param pv New parameter vector
     */
    virtual void setPars(const double* pv) = 0;

    /**
     * @param pv Parameter vector
     * @return   Is configuration 'pv' valid?
     */
    virtual bool isValidPars(const double* pv) const = 0;

    /**
     * Log-concavity simplifies EP algoritms to some extent.
     *
     * @return Is log t(s) (generalized) concave?
     */
    virtual bool isLogConcave() const = 0;
  };
//ENDNS

#endif
