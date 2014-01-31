/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header abstract class PotentialManager
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_POTENTIALMANAGER_H
#define EPTOOLS_POTENTIALMANAGER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"

//BEGINNS(eptools)
  /**
   * Manager for a set of scalar potentials t_j(s_j), of type
   * 'EPScalarPotential'.
   * <p>
   * The service 'getPot' returns the potential object for an index j.
   * A typical implementation has to serve 'getPot' being called by
   * a sequential loop over all or a subset of potentials.
   * <p>
   * ATTENTION: At present, implementations are typically not thread-safe.
   * 'getPot' returns a reference to 'EPScalarPotential'. This object must
   * not be used once 'getPot' is called again.
   * TODO: Thread-safe variant.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class PotentialManager
  {
  public:
    // Public methods

    virtual ~PotentialManager() {}

    /**
     * @return Number of potentials
     */
    virtual int size() const = 0;

    /**
     * NOTE: The returned object should be read-accessed only. In particular,
     * 'setPars' must not be used: the object returned is typically a temp.
     * copy anyway.
     *
     * @param j Potential index
     * @return  Potential object t_j(s_j)
     */
    virtual const EPScalarPotential& getPot(int j) const = 0;
  };
//ENDNS

#endif
