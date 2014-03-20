/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotentialFactory
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTENTIALFACTORY_H
#define EPTOOLS_EPPOTENTIALFACTORY_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"

//BEGINNS(eptools)
  /**
   * Provides factory method for supported 'EPScalarPotential'
   * subclasses. All these have to be registered here with a unique
   * ID.
   * <p>
   * Registration is static and compile-time, nothing fancy. The
   * subclass 'EPPotentialNamedFactory' also maintains the
   * association
   *   Name (string) -> ID (nonneg. int)
   * Here, Name is exposed towards interface and must not change,
   * while IDs may change internally.
   * <p>
   * NOTE: Can this be done better? The implementation has to modified
   * for every new 'EPScalarPotential' subclass, and each static method
   * here has to be modified.
   * Would be better if subclasses registered themselves here, passing
   * information and function pointers (would not work for constructors).
   * ==> ???
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotentialFactory
  {
  protected:
    // Internal constants (potential IDs)

    static const int potGaussian    =0;
    static const int potLaplace     =1;
    static const int potProbit      =2;
    static const int potHeaviside   =3;
    static const int potExponential =4;
    static const int potQuantRegress=5;
    static const int potGaussMixture=6;
    static const int potSpikeSlab   =7;
    static const int potLast        =7;
#ifdef HAVE_WORKAROUND
#include "src/eptools/potentials/EPPotentialFactory_workaround.h"
#endif

  public:
    // Public static methods

    static bool isValidID(int pid) {
#ifndef HAVE_WORKAROUND
      return (pid>=0 && pid<=potLast);
#else
      return (pid>=0 && pid<=potLast) || (pid>=potFirst_workaround &&
					  pid<=potLast_workaround);
#endif
    }

    /**
     * See 'EPScalarPotential' comments.
     *
     * @param pid Potential ID
     * @return    Argument group ID ('EPScalarPotential::atypeXXX')
     */
    static int getArgumentGroup(int pid);

    /**
     * Creates 'EPScalarPotential' object of correct type, given ID.
     * In 'pv', a valid initial parameter vector has to be passed. Use
     * 'createDefault' for default construction.
     * <p>
     * In 'annot', a void* to an annotation can be passed. This is ignored
     * by types which do not have annotations, but is mandatory for types
     * which do.
     * NOTE: The validity of a non-zero 'annot' is not checked. Passing
     * an invalid can lead to a crash.
     *
     * @param pid   Potential ID
     * @param pv    Initial parameter vector
     * @param annot S.a. Def.: 0
     * @return      New object
     */
    static EPScalarPotential* create(int pid,const double* pv,void* annot=0);

    /**
     * Creates 'EPScalarPotential' object of correct type, given ID. The
     * default constructor of the type is called. The type may require
     * construction parameters (see 'EPScalarPotential' header comments),
     * in which case 'pv' must point to these. Otherwise, 'pv' is ignored.
     * <p>
     * In 'annot', a void* to an annotation can be passed. This is ignored
     * by types which do not have annotations, but is mandatory for types
     * which do.
     * NOTE: The validity of a non-zero 'annot' is not checked. Passing
     * an invalid can lead to a crash.
     *
     * @param pid   Potential ID
     * @param pv    Construction parameters. Def.: 0
     * @param annot S.a. Def.: 0
     * @return      New object
     */
    static EPScalarPotential* createDefault(int pid,const double* pv=0,
					    void* annot=0);

#ifdef HAVE_WORKAROUND
  protected:
    static int getArgumentGroup_workaround(int pid);

    static EPScalarPotential* create_workaround(int pid,const double* pv,
						void* annot);

    static EPScalarPotential* createDefault_workaround(int pid,const double* pv,
						       void* annot);
  public:
#endif
  };
//ENDNS

#endif
