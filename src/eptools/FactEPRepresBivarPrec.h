/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class FactEPRepresBivarPrec
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_FACTEPREPRESBIVARPREC_H
#define EPTOOLS_FACTEPREPRESBIVARPREC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/FactorizedEPRepresentation.h"

//BEGINNS(eptools)
  /**
   * Representation (coupling factor, message parameters) for part of
   * model consisting of bivariate precision parameter potentials.
   * Extends 'FactorizedEPRepresentation'.  
   * <p>
   * The model part is
   *   prod_j t_j(s_j,tau_k(j)),   s = B x,  j=0:(m-1),
   * the potentials belong to argument group
   * 'EPScalarPotential::atypeBivarPrec'. The part B, beta, pi is maintained
   * by the superclass.
   *
   * [tau_k] are precision variables, their exponential family is Gamma (see
   * 'EPScalarPotential::compMoments'; parameters a>0, c>0). j -> k is a
   * function, maintained by the flat index 'tauInd':
   * - Index k(j) [m]
   * - For each k=0:(K-1): Start offset of J_k = {j | k(j)==k} [K]
   * - Dummy entry (start offset of J_K if it existed) [1]
   * - J_k, k=0:(K-1), each ascending order [m]
   * Moreover, message (Gamma) parameters are kept in 'aVals', 'cVals' (flat,
   * size m). 'accessTauRow', 'accessTauCol' are access methods indexed by
   * j, k resp.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactEPRepresBivarPrec : public FactorizedEPRepresentation
  {
  protected:
    // Additional members

    int numK;                        // Number of tau_k variables
    ArrayHandle<int> tauInd;         // j <-> k (header comment)
    ArrayHandle<double> aVals,cVals; // Parameters of tau messages

  public:
    // Public methods

    /**
     * Constructor.
     * NOTE: Arrays are not copied, but referred to.
     *
     * @param pnumN     Number x_i variables n
     * @param pnumM     Number potentials m
     * @param prowInd
     * @param pcolInd
     * @param pbmatVals
     * @param pbetaVals
     * @param ppiVals
     * @param pnumK     Number tau_k variables K
     * @param ptauInd
     * @param paVals
     * @param pcVals
     */
    FactEPRepresBivarPrec(int pnumN,int pnumM,const ArrayHandle<int>& prowInd,
			  const ArrayHandle<int>& pcolInd,
			  const ArrayHandle<double>& pbmatVals,
			  const ArrayHandle<double>& pbetaVals,
			  const ArrayHandle<double>& ppiVals,int pnumK,
			  const ArrayHandle<int>& ptauInd,
			  const ArrayHandle<double>& paVals,
			  const ArrayHandle<double>& pcVals) :
      FactorizedEPRepresentation(pnumN,pnumM,prowInd,pcolInd,pbmatVals,
				 pbetaVals,ppiVals),
      numK(pnumK),tauInd(ptauInd),aVals(paVals),cVals(pcVals) {
      int j,k,sz;

      if (pnumK==0 || ptauInd.size()!=2*pnumM+pnumK+1 ||
	  paVals.size()!=pnumM || pcVals.size()!=pnumM)
	throw InvalidParameterException(EXCEPT_MSG(""));
      // Run some basic checks
      for (j=0; j<pnumM; j++)
	if (ptauInd[j]<0 || ptauInd[j]>=pnumK)
	  throw InvalidParameterException(EXCEPT_MSG(""));
      for (k=0; k<pnumK; k++) {
	j=ptauInd[k+pnumM];
	sz=ptauInd[k+pnumM+1]-j;
	if (sz<1 || j<pnumM+pnumK || j+sz>ptauInd.size())
	  throw InvalidParameterException(EXCEPT_MSG(""));
	if (!Range::isIncreasing(ptauInd.p()+j,sz) ||
	    ptauInd[j]<0 || ptauInd[j+sz-1]>=pnumM)
	  throw InvalidParameterException(EXCEPT_MSG(""));
      }
    }

    virtual int numPrecVars() const {
      return numK;
    }

    /**
     * Access to precision parameter data for potential j.
     * NOTE: Use this for write access to EP parameters.
     *
     * @param j     Potential index
     * @param aP    Address for a_jk
     * @param cP    Address for c_jk
     * @return      k==k(j)
     */
    virtual int accessTauRow(int j,double*& aP,double*& cP) {
      if (j<0 || j>=numM) throw InvalidParameterException(EXCEPT_MSG(""));
      int k=tauInd[j];
      aP=aVals.p()+k; cP=cVals.p()+k;

      return k;
    }

    /**
     * Access to precision parameter data for variable tau_k. 'jInd'
     * is J_k, which is also the index into 'aP', 'cP'.
     *
     * @param k    Variable index
     * @param jInd Support index J_k (ascending)
     * @param aP   a message parameters (flat array)
     * @param cP   c message parameters (flat array)
     * @return     Size |J_k|
     */
    virtual int accessTauCol(int k,const int*& jInd,const double*& aP,
			     const double*& cP) {
      if (k<0 || k>=numK) throw InvalidParameterException(EXCEPT_MSG(""));
      int off=tauInd[numM+k];
      int sz=tauInd[numM+k+1]-off;
      jInd=tauInd.p()+off;
      aP=aVals.p(); cP=cVals.p();

      return sz;
    }

    /**
     * Compute parameters of Gamma marginals on [tau_k] from message
     * parameters 'aVals', 'cVals'.
     * If 'increm'==true, the marginals are added to 'margA', 'margC'.
     *
     * @param margA  Marginal pars. a ret. here
     * @param margC  Marginal pars. c ret. here
     * @param increm Incremental? Def.: false
     */
    virtual void compTauMarginals(double* margA,double* margC,
				  bool increm=false);
  };

  // Inline methods

  inline void
  FactEPRepresBivarPrec::compTauMarginals(double* margA,double* margC,
					  bool increm)
  {
    int k,j,jj,sz;
    double mA,mC;
    const double* aP,*cP;
    const int* jInd;

    for (k=0; k<numK; k++) {
      sz=accessTauCol(k,viInd,jInd,aP,cP);
      for (j=0,mA=mC=0.0; j<sz; j++) {
	jj=jInd[j];
	mA+=aP[jj]; mC+=cP[jj];
      }
      if (!increm) {
	margA[i]=mA; margC[i]=mC;
      } else {
	margA[i]+=mA; margC[i]+=mC;
      }
    }
  }
//ENDNS

#endif
