/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class ContainerPotManager
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_CONTAINERPOTMANAGER_H
#define EPTOOLS_CONTAINERPOTMANAGER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/PotentialManager.h"

//BEGINNS(eptools)
  /**
   * Container class for 'PotentialManager'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class ContainerPotManager : public PotentialManager
  {
  protected:
    // Members

    ArrayHandle<Handle<PotentialManager> > pmArr;
    ArrayHandle<int> startPos;

  public:
    // Public methods

    /**
     * Constructor. The array 'parr' is copied here.
     *
     * @param parr Array of child objects
     */
    ContainerPotManager(const ArrayHandle<Handle<PotentialManager> >& parr) {
      int num=parr.size(),i,off,sz,nprec;
      bool havePrec=false;

      if (num==0) throw InvalidParameterException(EXCEPT_MSG(""));
      startPos.changeRep(num);
      // Also have to check that potentials of group 'atypeBivarPrec' come
      // last (if at all). 'havePrec' is false until such potentials are
      // detected.
      for (i=off=0; i<num; i++) {
	startPos[i]=off; sz=parr[i]->size();
	off+=sz;
	nprec=pmArr[i]->numArgumentGroup(EPScalarPotential::atypeBivarPrec);
	if (havePrec && nprec<sz)
	  throw InvalidParameterException(EXCEPT_MSG("'atypeBivarPrec' potentials must form suffix"));
	havePrec=(nprec>0);
      }
      pmArr.copy(parr);
    }

    int size() const {
      int i=pmArr.size()-1;

      return startPos[i]+pmArr[i]->size();
    }

    int numArgumentGroup(int atype) const {
      int ret=0;

      for (int i=0; i<pmArr.size(); i++)
	ret+=pmArr[i]->numArgumentGroup(atype);

      return ret;
    }

    const EPScalarPotential& getPot(int j) const {
      int i,ic;

      if (j<0 || j>=size()) throw OutOfRangeException(EXCEPT_MSG(""));
      i=getRelPos(j,ic);

      return pmArr[ic]->getPot(i);
    }

  protected:
    // Internal methods

    /**
     * @param j  Potential index
     * @param ic Child index ret. here
     * @return   Position within 'pmArr[ic]'
     */
    int getRelPos(int j,int& ic) const {
      int num=pmArr.size();

      MYASS(j>=0 && j<size());
      for (ic=0; ic<num; ic++)
	if (j<startPos[ic]) break;

      return j-startPos[--ic];
    }
  };
//ENDNS

#endif
