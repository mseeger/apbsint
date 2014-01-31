/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class DefaultPotManager
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_DEFAULTPOTMANAGER_H
#define EPTOOLS_DEFAULTPOTMANAGER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/PotentialManager.h"

//BEGINNS(eptools)
  /**
   * Generic implementation of EP potential manager based on single
   * 'EPScalarPotential' object.
   * <p>
   * The object 'epPot' must be passed upon construction. All N=='num'
   * potentials t_j(s_j) are represented by 'epPot', but can have
   * different parameters. Namely, each of the parameters can be
   * shared by all t_j(s_j), or be different for each potential.
   *
   * Parameter values are given by the flat vector 'parVec' and the
   * index vector 'parOff'. The first value for the k-th parameter is
   * at 'parVec[parOff[k]]'. The difference 'parOff[k+1]-parOff[k]' is
   * either 1 (shared; 'parShrd[k]' true) or N (individual;
   * 'parShrd[k]' false).
   * <p>
   * ATTENTION: This implementation is not thread-safe. 'epPot' is
   * used by all 'getPot' calls, and the object is reconfigured
   * accordingly.
   * <p>
   * TODO: Currently, parameter values are fixed upon construction.
   * Should allow them to be modified later on.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DefaultPotManager : public PotentialManager
  {
  protected:
    // Members

    mutable Handle<EPScalarPotential> epPot;  // Potential object
    int num;                                  // Number of potentials N
    ArrayHandle<double> parVec;               // See header comment
    ArrayHandle<int> parOff;                  // "
    ArrayHandle<int> parShrd;                 // Is parameter shared?
    mutable ArrayHandle<double> tmpVec;

  public:
    // Public methods

    /**
     * Constructor. See header comment for 'parVec', 'parOff'. The parameter
     * values are checked for validness.
     *
     * @param peppot     Potential object
     * @param pnum       Number of potentials N
     * @param ppvec      Value for 'parVec'
     * @param ppshd      Value for 'parShrd'
     * @param checkValid Check whether parameters are valid? Def.: true
     */
    DefaultPotManager(const Handle<EPScalarPotential>& peppot,int pnum,
		      const ArrayHandle<double>& ppvec,
		      const ArrayHandle<int>& ppshd,bool checkValid=true);

    int size() const {
      return num;
    }

    const EPScalarPotential& getPot(int j) const {
      if (j<0 || j>=size()) throw OutOfRangeException(EXCEPT_MSG(""));
      if (parOff.size()>0) {
	getPotPars(j,tmpVec.p());
	epPot->setPars(tmpVec.p());
      }

      return *epPot;
    }

  protected:
    // Internal methods

    /**
     * @param j   Potential index
     * @param arr Parameters written here
     */
    void getPotPars(int j,double* arr) const {
      for (int i=0; i<parOff.size(); i++)
	arr[i]=parVec[parOff[i]+(parShrd[i]?0:j)];
    }
  };
//ENDNS

#endif
