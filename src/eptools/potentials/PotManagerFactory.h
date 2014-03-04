/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class PotManagerFactory
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_POTMANAGERFACTORY_H
#define EPTOOLS_POTMANAGERFACTORY_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/PotentialManager.h"

//BEGINNS(eptools)
  /**
   * Factory method to create 'PotentialManager' object from a compressed
   * description.
   * <p>
   * The object is either a 'DefaultPotManager' or a 'ContainerPotManager'
   * of 'DefaultPotManager' children. The simple description in terms of
   * flat vectors is for the Matlab MEX and Python interfaces.
   * <p>
   * Checking internal representation for errors:
   * The 'create' service does not check the validity of potential
   * parameter vectors, and its exceptions are not in general meaningful
   * error messages.
   * The 'checkRepres' service runs an exhaustive number of checks on the
   * representation normally passed to 'create'. Its error messages can be
   * returned to the user.
   * NOTE: This is a compromise to implementing full-blown exception handling.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class PotManagerFactory
  {
  public:
    // Public static methods

    /**
     * Creates potential manager, type 'ContainerPotManager' with K
     * 'DefaultPotManager' children (or a single 'DefaultPotManager', if
     * K==1). 'potIDs' (size K) contains potential IDs
     * ('EPPotentialFactory'), 'numPot' (size K) contains numbers N_k of
     * potentials per 'DefaultPotManager'. 'parVec' (double array) and
     * 'parShrd' (int array) are concatenations of the corr. arrays passed
     * to the 'DefaultPotManager' constructor.
     * 'annObj' (size K) contains void* to annotation objects. Entries are
     * ignored for types without annotations (should be 0 to be safe), but
     * are mandatory for annotated types.
     * ATTENTION: Unsafe! A non-zero pointer cannot be checked for validity.
     * This could crash the program or worse.
     * <p>
     * NOTE: The K 'EPScalarPotential' objects created here are
     * default-constructed (see 'EPPotentialFactory::createDefault'). If the
     * potential type has construction parameters, these must be in the prefix
     * of 'parVec', and the corr. 'parShrd' entries must be 1.
     * <p>
     * NOTE: We do not check the representation for validity, unless this
     * hinders the creating of the potential manager. In particular, the
     * content of 'parVec' can be invalid (except for construction parameters,
     * s.a.).
     * Call 'checkRepres' for a representation before using it with 'create'.
     * If 'checkPure'==true, all potentials must be in the same argument
     * group (see 'EPScalarPotential').
     * <p>
     * ATTENTION: Implementation unsafe, does not copy content or handles
     * of 'parVec', 'parShrd', just uses pointers. If content is changed
     * or deallocated, the PM becomes invalid.
     * ==> Use only to create temporary potential managers, not to be
     *     kept around.
     *
     * @param potIDs
     * @param numPot
     * @param parVec
     * @param parShrd
     * @param annObj
     * @param checkPure S.a. Def.: false
     * @return          New potential manager
     */
    static PotentialManager* create(const ArrayHandle<int>& potIDs,
				    const ArrayHandle<int>& numPot,
				    const ArrayHandle<double>& parVec,
				    const ArrayHandle<int>& parShrd,
				    const ArrayHandle<void*>& annObj,
				    bool checkPure=false);

    /**
     * Check representation (as passed to 'create') for validity. If an
     * error is detected, a 'InvalidParameterException' with a meaningful
     * message is thrown.
     * The K parts corr. to children are referred to as "blocks". Potentials
     * are numbered relative to a block. If K==1, blocks are not mentioned.
     * The parameter vector constellation for each potential is checked.
     * Block or potential positions are 0-floor. Pass 'posoff'=1 to make
     * them 1-floor in the error message (e.g., when using this in a MEX
     * function).
     * <p>
     * Construction parameters: If the potential type for a block has
     * construction parameters, they must form the prefix of the corr.
     * 'parVec' part, and the corr. 'parShrd' entries must all be 1.
     * <p>
     * If 'checkPure'==true, all potentials must be in the same argument
     * group (see 'EPScalarPotential') for the representation to be valid.
     *
     * @param potIDs
     * @param numPot
     * @param parVec
     * @param parShrd
     * @param annObj
     * @param posoff    S.a. Def.: 0
     * @param checkPure S.a. Def.: false
     */
    static void checkRepres(const ArrayHandle<int>& potIDs,
			    const ArrayHandle<int>& numPot,
			    const ArrayHandle<double>& parVec,
			    const ArrayHandle<int>& parShrd,
			    const ArrayHandle<void*>& annObj,int posoff=0,
			    bool checkPure=false);
  };
//ENDNS

#endif
