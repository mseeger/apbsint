/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class FactEPMaximumCValues
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_FACTEPMAXIMUMCVALUES_H
#define EPTOOLS_FACTEPMAXIMUMCVALUES_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/MaximumValuesService.h"
#include "src/eptools/FactorizedEPRepresentation.h"

//BEGINNS(eptools)
  /**
   * Specialization of 'MaximumValuesService' to max_j c_jk, where the
   * structure j -> k and the c values (Gamma parameters) are maintained
   * by a 'FactorizedEPRepresentation' object.
   * <p>
   * NOTE: The index j over potentials is 0-based, it runs over bivariate
   * precision potentials only.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactEPMaximumCValues : public MaximumValuesService
  {
  protected:
    // Additional members

    Handle<FactorizedEPRepresentation> epRepr;

  public:
    // Public methods

    /**
     * Constructor. Consistency of 'ptopVal' with 'pepRepr' is not checked.
     * 'pepRepr' must contain precision potentials.
     *
     * @param pepRepr   EP representation (bivariate precision potentials)
     * @param pmaxSize  K
     * @param pnumValid Entries must be in 1:pmaxSize
     * @param ptopInd
     * @param ptopVal
     * @param psubInd   Optional
     * @param psubExcl  Def.: false
     */
    FactEPMaximumCValues(const Handle<FactorizedEPRepresentation>& pepRepr,
			 int pmaxSize,const ArrayHandle<int>& pnumValid,
			 const ArrayHandle<int>& ptopInd,
			 const ArrayHandle<double>& ptopVal,
			 const ArrayHandle<int>& psubInd=
			 ArrayHandleZero<int>::get(),bool psubExcl=false) :
      MaximumValuesService(pepRepr->numPrecVariables(),
			   pepRepr->numBVPrecPotentials(),pmaxSize,pnumValid,
			   ptopInd,ptopVal,psubInd,psubExcl),epRepr(pepRepr) {
      if (pepRepr->numPrecVariables()==0)
	throw WrongStatusException(EXCEPT_MSG(""));
    }

    int numVariables() const {
      return epRepr->numPrecVariables();
    }

    int numFactors() const {
      return epRepr->numBVPrecPotentials();
    }

    int getFactorValues(int i,const int*& vind,const int*& jind,
			const double*& xarr) const {
      // Maps directly to 'FactorizedEPRepresentation::accessTauCol'
      const double* aP;
      int sz=epRepr->accessTauCol(i,vind,aP,xarr);
      jind=vind;

      return sz;
    }
  };
//ENDNS

#endif
