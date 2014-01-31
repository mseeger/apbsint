/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class FactorizedEPRepresentation
 * ------------------------------------------------------------------- */

#include "src/eptools/FactorizedEPRepresentation.h"

//BEGINNS(eptools)
  FactorizedEPRepresentation::FactorizedEPRepresentation(int pnumN,int pnumM,const ArrayHandle<int>& prowInd,const ArrayHandle<int>& pcolInd,const ArrayHandle<double>& pbmatVals,const ArrayHandle<double>& pbetaVals,const ArrayHandle<double>& ppiVals) :
    numN(pnumN),numM(pnumM),rowInd(prowInd),colInd(pcolInd),
    bmatVals(pbmatVals),betaVals(pbetaVals),piVals(ppiVals)
  {
    int j,sz,off,nnz=pbmatVals.size();

    if (numN==0 || numM==0 || pbetaVals.size()!=nnz || ppiVals.size()!=nnz ||
	prowInd.size()<=numM+1 || pcolInd.size()<=numN+1)
      throw InvalidParameterException(EXCEPT_MSG(""));
    // Run some basic checks
    if (prowInd[numM]!=nnz || prowInd[0]!=0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    for (j=0; j<numM; j++) {
      off=prowInd[j]; sz=prowInd[j+1]-off;
      // NOTE: Zero rows are not allowed!
      if (sz<=0 || sz>numN)
	throw InvalidParameterException(EXCEPT_MSG(""));
    }
    if (pcolInd[numN]!=2*nnz+numN+1 || pcolInd[0]!=numN+1)
    for (j=0; j<numN; j++) {
      off=pcolInd[j]; sz=pcolInd[j+1]-off;
      if (sz%2==1)
	throw InvalidParameterException(EXCEPT_MSG(""));
      sz/=2;
      // NOTE: Zero columns are allowed
      if (sz<0 || sz>numM)
	throw InvalidParameterException(EXCEPT_MSG(""));
    }
  }
//ENDNS
