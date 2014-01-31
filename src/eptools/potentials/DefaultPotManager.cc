/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class DefaultPotManager
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/DefaultPotManager.h"

//BEGINNS(eptools)
  DefaultPotManager::DefaultPotManager(const Handle<EPScalarPotential>& peppot,
				       int pnum,
				       const ArrayHandle<double>& ppvec,
				       const ArrayHandle<int>& ppshd,
				       bool checkValid) :
    epPot(peppot),num(pnum),parVec(ppvec),parShrd(ppshd)
  {
    int i,np=ppshd.size(),off;

    if (pnum<=0 || np!=peppot->numPars())
      throw InvalidParameterException(EXCEPT_MSG(""));
    parOff.changeRep(np);
    for (i=0,off=0; i<np; i++) {
      parOff[i]=off;
      off+=(ppshd[i]?1:pnum);
    }
    if (ppvec.size()!=off) throw InvalidParameterException(EXCEPT_MSG(""));
    tmpVec.changeRep(np);
    // Check whether all parameters are valid
    if (checkValid && np>0) {
      for (i=0; i<pnum; i++) {
	getPotPars(i,tmpVec.p());
	if (!peppot->isValidPars(tmpVec.p()))
	  throw InvalidParameterException(EXCEPT_MSG(""));
      }
    }
  }
//ENDNS
