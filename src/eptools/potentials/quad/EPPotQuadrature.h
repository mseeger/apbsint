/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class EPPotQuadrature
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTQUADRATURE_H
#define EPTOOLS_EPPOTQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/quad/QuadraturePotential.h"
#include "src/eptools/potentials/EPScalarPotential.h"

//BEGINNS(eptools)
  /**
   * Base class of implementations of 'EPScalarPotential' which use
   * numerical quadrature routines in order to provide EP update services.
   * To this end, a 'QuadraturePotential' object has to be maintained here
   * in 'quadPot'.
   * All 'EPScalPotentialBase' services are mapped back to this object as
   * well.
   * <p>
   * NOTE: Any further services useful for (almost) all quadrature
   * implementations should be added here.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotQuadrature : public EPScalarPotential
  {
  protected:
    // Members

    Handle<QuadraturePotential> quadPot;

  public:
    // Public methods

    /**
     * Constructor
     *
     * @param qpot Quadrature potential services
     */
    EPPotQuadrature(const Handle<QuadraturePotential>& qpot) : quadPot(qpot) {}

    int numPars() const {
      return quadPot->numPars();
    }

    int numConstPars() const {
      return quadPot->numConstPars();
    }

    void getPars(double* pv) const {
      quadPot->getPars(pv);
    }

    void setPars(const double* pv) {
      quadPot->setPars(pv);
    }

    bool isValidPars(const double* pv) const {
      return quadPot->isValidPars(pv);
    }

    bool isLogConcave() const {
      return quadPot->isLogConcave();
    }
  };
//ENDNS

#endif
