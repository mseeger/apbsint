/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class FactEPMaximumAValues
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_FACTEPMAXIMUMAVALUES_H
#define EPTOOLS_FACTEPMAXIMUMAVALUES_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/MaximumValuesService.h"
#include "src/eptools/FactEPRepresBivarPrec.h"

//BEGINNS(eptools)
  /**
   * Specialization of 'MaximumValuesService' to max_j a_jk, where the
   * structure j -> k and the a values (Gamma parameters) are maintained
   * by a 'FactEPRepresBivarPrec' object.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactEPMaximumAValues : public MaximumValuesService
  {
  protected:
    // Additional members

    Handle<FactEPRepresBivarPrec> epRepr;

  public:
    // Public methods

    /**
     * Constructor. Consistency of 'ptopVal' with 'pepRepr' is not checked.
     *
     * @param pepRepr   EP representation (bivariate precision potentials)
     * @param pmaxSize  K
     * @param pnumValid Entries must be in 1:pmaxSize
     * @param ptopInd
     * @param ptopVal
     * @param psubInd   Optional
     * @param psubExcl  Def.: false
     */
    FactEPMaximumAValues(const Handle<FactEPRepresBivarPrec>& pepRepr,
			 int pmaxSize,const ArrayHandle<int>& pnumValid,
			 const ArrayHandle<int>& ptopInd,
			 const ArrayHandle<double>& ptopVal,
			 const ArrayHandle<int>& psubInd=
			 ArrayHandleZero<int>::get(),bool psubExcl=false) :
      MaximumValuesService(pepRepr->numPrecVars(),pepRepr->numPotentials(),
			   pmaxSize,pnumValid,ptopInd,ptopVal,psubInd,psubExcl),
      epRepr(pepRepr) {}

    int numVariables() const {
      return epRepr->numPrecVars();
    }

    int numFactors() const {
      return epRepr->numPotentials();
    }

    int getFactorValues(int i,const int*& vind,const int*& jind,
			const double*& xarr) const {
      // Maps directly to 'FactEPRepresBivarPrec::accessTauCol'
      const double* cP;
      int sz=epRepr->accessTauCol(i,vind,xarr,cP);
      jind=vind;

      return sz;
    }
  };
//ENDNS

#endif
