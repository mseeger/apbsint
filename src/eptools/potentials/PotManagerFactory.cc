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
					      const ArrayHandle<int>& parShrd,
					      const ArrayHandle<void*>& annObj)
  {
    int i,j,k,numk=potIDs.size(),atype;
    ArrayHandle<Handle<PotentialManager> > parr;
    ArrayHandle<double> pvecMsk;
    ArrayHandle<int> shrdMsk;
    double* pvecP=parVec.p();
    int* shrdP=parShrd.p();
    bool hasBVPrec=false;

    if (numPot.size()!=numk || annObj.size()!=numk)
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
      //cout << "ManFactory: k=" << k << endl; // DEBUG!
      Handle<EPScalarPotential> epPot;
      try {
	epPot.changeRep(EPPotentialFactory::createDefault(pid,pvecP,annObj[k]));
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
      atype=epPot->getArgumentGroup();
      if (hasBVPrec && atype!=EPScalarPotential::atypeBivarPrec)
	throw InvalidParameterException(EXCEPT_MSG(""));
      hasBVPrec=(atype==EPScalarPotential::atypeBivarPrec);
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
				      const ArrayHandle<void*>& annObj,
				      int posoff,
				      const ArrayHandle<int>& tauInd)
  {
    int k,numk=potIDs.size(),atype,numBVPrec=0;
    ArrayHandle<double> pvecMsk,tmpVec;
    ArrayHandle<int> shrdMsk,parOff;
    double* pvecP=parVec.p();
    int* shrdP=parShrd.p();
    ArrayHandle<char> errStr;

    if (numPot.size()!=numk || annObj.size()!=numk)
      throw InvalidParameterException("NUMPOT or ANNOBJ have wrong size");
    for (k=0; k<numk; k++) {
      if (!EPPotentialFactory::isValidID(potIDs[k])) {
	errStr.changeRep(50);
	sprintf(errStr.p(),"Block %d: POTIDS entry invalid",k+posoff);
	throw InvalidParameterException(errStr.p());
      }
      if (numPot[k]<=0) {
	errStr.changeRep(59);
	sprintf(errStr.p(),"Block %d: NUMPOT entry must be positive",k+posoff);
	throw InvalidParameterException(errStr.p());
      }
    }
    for (k=0; k<numk; k++) {
      int i,j,pid=potIDs[k],npot=numPot[k];
      // Have to create instance in order to request number of parameters
      Handle<EPScalarPotential> epPot;
      try {
	epPot.changeRep(EPPotentialFactory::createDefault(pid,pvecP,
							  annObj[k]));
      } catch (StandardException ex) {
	errStr.changeRep(65+strlen(ex.msg()));
	sprintf(errStr.p(),"Block %d: Cannot create potential object (%s)",
		k+posoff,ex.msg());
	throw InvalidParameterException(errStr.p());
      }
      int numConstPars=epPot->numConstPars();
      int npar=epPot->numPars(); // Could be 0
      if (numConstPars>0) {
	// Checks for construction parameters
	if (npar<numConstPars) {
	  errStr.changeRep(71);
	  sprintf(errStr.p(),"Block %d: Need %d construction parameters",
		  k+posoff,numConstPars);
	  throw InvalidParameterException(errStr.p());
	}
	for (i=0; i<numConstPars; i++)
	  if (!shrdP[i]) {
	    errStr.changeRep(73);
	    sprintf(errStr.p(),
		    "Block %d: PARSHRD invalid for construction parameters",
		    k+posoff);
	    throw InvalidParameterException(errStr.p());
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
	    errStr.changeRep(74);
	    if (numk>1)
	      sprintf(errStr.p(),"Potential %d in block %d: Invalid parameters",
		      i+posoff,k+posoff);
	    else
	      sprintf(errStr.p(),"Potential %d: Invalid parameters",i+posoff);
	    throw InvalidParameterException(errStr.p());
	  }
	}
      }
      atype=epPot->getArgumentGroup();
      if (atype==EPScalarPotential::atypeBivarPrec)
	numBVPrec+=npot;
      else if (numBVPrec>0)
	throw InvalidParameterException("Potentials of group 'atypeBivarPrec' must come last");
    }
    if (parShrd.p()+parShrd.size()>shrdP)
      throw InvalidParameterException("PARSHRD too long");
    if (parVec.p()+parVec.size()>pvecP)
      throw InvalidParameterException("PARVEC too long");
    // Check 'tauInd' (if given)
    if (numBVPrec==0) {
      if (!(tauInd==0))
	throw InvalidParameterException("TAUIND only together with 'atypeBivarPrec' potentials");
    } else {
      checkBVPrecTauInd(tauInd,numBVPrec);
    }
  }

  void PotManagerFactory::checkBVPrecTauInd(const ArrayHandle<int>& tauInd,
					    int numBVPrec)
  {
    if (numBVPrec<=0)
      throw InvalidParameterException("No bivariate precision potentials");
    if (tauInd==0)
      throw InvalidParameterException("TAUIND must be given");
    int dimK=tauInd[numBVPrec],i,j,sz;
    if (dimK<=0 || tauInd.size()!=2*numBVPrec+dimK+2)
      throw InvalidParameterException("TAUIND wrong size");
    ArrayHandle<bool> karr(dimK);
    for (k=0; k<dimK; k++) karr[k]=false;
    for (j=sz=0; j<numBVPrec; j++) {
      k=tauInd[j];
      if (k<0 || k>=dimK)
	throw InvalidParameterException("TAUIND wrong");
      if (!karr[k]) {
	karr[k]=true; sz++;
      }
    }
    if (sz<dimK)
      throw InvalidParameterException("TAUIND: Every k value must occur at least once");
    for (k=0; k<dimK; k++) {
      j=tauInd[k+numBVPrec+1];
      sz=tauInd[k+numBVPrec+2]-j;
      if (sz<1 || j<numBVPrec+dimK+2 || j+sz>tauInd.size())
	throw InvalidParameterException("TAUIND wrong");
      if (!Range::isIncreasing(tauInd.p()+j,sz) ||
	  tauInd[j]<0 || tauInd[j+sz-1]>=numBVPrec)
	throw InvalidParameterException("TAUIND wrong");
      for (i=0; i<sz; i++)
	if (tauInd[tauInd[i+j]]!=k)
	  throw InvalidParameterException("TAUIND wrong: Forward and inverse different");
    }
  }
//ENDNS
