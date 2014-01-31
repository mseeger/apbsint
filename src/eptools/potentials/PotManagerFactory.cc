/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class PotManagerFactory
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/PotManagerFactory.h"
#include "src/eptools/potentials/DefaultPotManager.h"
#include "src/eptools/potentials/ContainerPotManager.h"
#include "src/eptools/potentials/EPPotentialFactory.h"

//BEGINNS(eptools)
  PotentialManager* PotManagerFactory::create(const ArrayHandle<int>& potIDs,
					      const ArrayHandle<int>& numPot,
					      const ArrayHandle<double>& parVec,
					      const ArrayHandle<int>& parShrd)
  {
    int i,j,k,numk=potIDs.size();
    ArrayHandle<Handle<PotentialManager> > parr;
    ArrayHandle<double> pvecMsk;
    ArrayHandle<int> shrdMsk;
    double* pvecP=parVec.p();
    int* shrdP=parShrd.p();

    if (numPot.size()!=numk)
      throw InvalidParameterException(EXCEPT_MSG(""));
    for (k=0; k<numk; k++)
      if (!EPPotentialFactory::isValidID(potIDs[k]) || numPot[k]<=0)
	throw InvalidParameterException(EXCEPT_MSG(""));
    if (numk>1)
      parr.changeRep(numk);
    for (k=0; k<numk; k++) {
      // Create 'DefaultPotManager' object
      //sprintf(debStr,"k=%d",k); printMsgStdout(debStr); // DEBUG!
      int pid=potIDs[k],npot=numPot[k];
      // Have to create instance in order to request number of parameters
      // NOTE: We can pass 'pvecP' for construction parameters (if any),
      // without having to prepare a parameter vector or even knowing the
      // size. These parameters must form the prefix
      Handle<EPScalarPotential> epPot;
      try {
	epPot.changeRep(EPPotentialFactory::createDefault(pid,pvecP));
      } catch (StandardException ex) {
	throw InvalidParameterException(EXCEPT_MSG("Cannot create potential object"));
      }
      int numConstPars=epPot->numConstPars();
      int npar=epPot->numPars(); // Could be 0
      if (numConstPars>0) {
	// Checks for construction parameters
	if (npar<numConstPars) throw InvalidParameterException(EXCEPT_MSG(""));
	for (i=0; i<numConstPars; i++)
	  if (!shrdP[i]) throw InvalidParameterException(EXCEPT_MSG(""));
      }
      if (npar>0) {
	shrdMsk.changeRep(shrdP,npar,false);
	shrdP+=npar;
	for (i=j=0; i<npar; i++)
	  j+=(shrdMsk[i]?1:npot); // Length of 'parVec' part
	pvecMsk.changeRep(pvecP,j,false);
	pvecP+=j;
      } else {
	// Potential has no parameters
	shrdMsk.changeRep(0); pvecMsk.changeRep(0);
      }
      // ATTENTION: 'shrdMsk', 'pvecMsk' do not own their buffers, and
      // they do not copy 'parShrd', 'parVec' content!
      PotentialManager* pmanP=
	new DefaultPotManager(epPot,npot,pvecMsk,shrdMsk,false);
      if (numk==1)
	return pmanP; // Single 'DefaultPotManager'
      else
	parr[k].changeRep(pmanP);
    }
    MYASS(numk>1);

    return new ContainerPotManager(parr);
  }

  void PotManagerFactory::checkRepres(const ArrayHandle<int>& potIDs,
				      const ArrayHandle<int>& numPot,
				      const ArrayHandle<double>& parVec,
				      const ArrayHandle<int>& parShrd,
				      int posoff)
  {
    int k,numk=potIDs.size();
    ArrayHandle<double> pvecMsk,tmpVec;
    ArrayHandle<int> shrdMsk,parOff;
    double* pvecP=parVec.p();
    int* shrdP=parShrd.p();
    char errStr[100];

    if (numPot.size()!=numk)
      throw InvalidParameterException("NUMPOT has wrong size");
    for (k=0; k<numk; k++) {
      if (!EPPotentialFactory::isValidID(potIDs[k])) {
	sprintf(errStr,"Block %d: POTIDS entry invalid",k+posoff);
	throw InvalidParameterException(errStr);
      }
      if (numPot[k]<=0) {
	sprintf(errStr,"Block %d: NUMPOT entry must be positive",k+posoff);
	throw InvalidParameterException(errStr);
      }
    }
    for (k=0; k<numk; k++) {
      int i,j,pid=potIDs[k],npot=numPot[k];
      // Have to create instance in order to request number of parameters
      Handle<EPScalarPotential> epPot;
      try {
	epPot.changeRep(EPPotentialFactory::createDefault(pid,pvecP));
      } catch (StandardException ex) {
	sprintf(errStr,"Block %d: Cannot create potential object (%s)",
		k+posoff,ex.msg());
	throw InvalidParameterException(errStr);
      }
      int numConstPars=epPot->numConstPars();
      int npar=epPot->numPars(); // Could be 0
      if (numConstPars>0) {
	// Checks for construction parameters
	if (npar<numConstPars) {
	  sprintf(errStr,"Block %d: Need %d construction parameters",k+posoff,
		  numConstPars);
	  throw InvalidParameterException(errStr);
	}
	for (i=0; i<numConstPars; i++)
	  if (!shrdP[i]) {
	    sprintf(errStr,
		    "Block %d: PARSHRD invalid for construction parameters",
		    k+posoff);
	    throw InvalidParameterException(errStr);
	  }
      }
      if (npar>0) {
	if (parShrd.p()+parShrd.size()<shrdP+npar)
	  throw InvalidParameterException("PARSHRD too short");
	shrdMsk.changeRep(shrdP,npar,false);
	shrdP+=npar;
	if (parOff.size()<npar) {
	  // Helper arrays for parameter validity checks
	  parOff.changeRep(npar);
	  tmpVec.changeRep(npar);
	}
	for (i=j=0; i<npar; i++) {
	  parOff[i]=j;
	  j+=(shrdMsk[i]?1:npot);
	}
	if (parVec.p()+parVec.size()<pvecP+j)
	  throw InvalidParameterException("PARVEC too short");
	pvecMsk.changeRep(pvecP,j,false);
	pvecP+=j;
	// Loop over potentials: Check validity
	for (i=0; i<npot; i++) {
	  // Assemble parameter vector (see 'DefaultPotManager::getPotPars')
	  for (j=0; j<npar; j++)
	    tmpVec[j]=pvecMsk[parOff[j]+(shrdMsk[j]?0:i)];
	  if (!epPot->isValidPars(tmpVec.p())) {
	    if (numk>1)
	      sprintf(errStr,"Potential %d in block %d: Invalid parameters",
		      i+posoff,k+posoff);
	    else
	      sprintf(errStr,"Potential %d: Invalid parameters",i+posoff);
	    throw InvalidParameterException(errStr);
	  }
	}
      }
    }
    if (parShrd.p()+parShrd.size()>shrdP)
      throw InvalidParameterException("PARSHRD too long");
    if (parVec.p()+parVec.size()>pvecP)
      throw InvalidParameterException("PARVEC too long");
  }
//ENDNS
