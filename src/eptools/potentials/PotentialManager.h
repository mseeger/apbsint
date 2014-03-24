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
   * Manager for a set of potentials t_j(.), of type 'EPScalarPotential'.
   * <p>
   * The service 'getPot' returns the potential object for an index j.
   * A typical implementation has to serve 'getPot' being called by
   * a sequential loop over all or a subset of potentials.
   * <p>
   * A PM may contain potentials of different argument groups (see
   * 'EPScalarPotential'). If it contains bivariate precision potentials
   * ('EPScalarPotential::atypeBivarPrec'), these must form a contiguous
   * suffix. The index of the first such potential is determined as
   *   'size() - numArgumentGroup(EPScalarPotential::atypeBivarPrec)'.
   * ATTENTION: Implementations have to ensure that this constaint holds
   * true.
   * <p>
   * ATTENTION: Typical implementations are not thread-safe, in that calls
   * to 'getPot' re-use the same 'EPScalarPotential' object, instead of
   * creating a new one.
   * If a thread-safe variant is required, the return type must be changed
   * from 'const EPScalarPotential&' to 'Handle<EPScalarPotential>'.
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
     * Each potential belongs to an argument group (see 'EPScalarPotential').
     *
     * @param atype Argument group
     * @return      Number of potential in argument group 'atype'
     */
    virtual int numArgumentGroup(int atype) const = 0;

    /**
     * NOTE: The returned object should be read-accessed only. In particular,
     * 'setPars' must not be used: the object returned is typically a temp.
     * copy anyway.
     * ATTENTION: Not thread-safe.
     *
     * @param j Potential index
     * @return  Potential object t_j(.)
     */
    virtual const EPScalarPotential& getPot(int j) const = 0;
  };
//ENDNS

#endif
