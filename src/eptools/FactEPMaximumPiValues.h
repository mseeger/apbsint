/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class FactEPMaximumPiValues
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_FACTEPMAXIMUMPIVALUES_H
#define EPTOOLS_FACTEPMAXIMUMPIVALUES_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/MaximumValuesService.h"
#include "src/eptools/FactorizedEPRepresentation.h"

//BEGINNS(eptools)
  /**
   * Specialization of 'MaximumValuesService' to max_k pi_ki, where the
   * factor group (coupling factor B) and the pi values are maintained
   * by a 'FactorizedEPRepresentation' object.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactEPMaximumPiValues : public MaximumValuesService
  {
  protected:
    // Additional members

    Handle<FactorizedEPRepresentation> epRepr;

  public:
    // Public methods

    /**
     * Constructor. Consistency of 'ptopVal' with 'pepRepr' is not checked.
     *
     * @param pepRepr   EP representation
     * @param pmaxSize  K
     * @param pnumValid Entries must be in 1:pmaxSize
     * @param ptopInd
     * @param ptopVal
     * @param psubInd   Optional
     * @param psubExcl  Def.: false
     */
    FactEPMaximumPiValues(const Handle<FactorizedEPRepresentation>& pepRepr,
			  int pmaxSize,const ArrayHandle<int>& pnumValid,
			  const ArrayHandle<int>& ptopInd,
			  const ArrayHandle<double>& ptopVal,
			  const ArrayHandle<int>& psubInd=
			  ArrayHandleZero<int>::get(),bool psubExcl=false) :
      MaximumValuesService(pepRepr->numVariables(),pepRepr->numPotentials(),
			   pmaxSize,pnumValid,ptopInd,ptopVal,psubInd,psubExcl),
      epRepr(pepRepr) {}

    int numVariables() const {
      return epRepr->numVariables();
    }

    int numFactors() const {
      return epRepr->numPotentials();
    }

    int getFactorValues(int i,const int*& vind,const int*& jind,
			const double*& xarr) const {
      // Maps directly to 'FactorizedEPRepresentation::accessCol'
      double* bP,*betaP;

      return epRepr->accessCol(i,vind,jind,bP,betaP,xarr)
    }
  };
//ENDNS

#endif
