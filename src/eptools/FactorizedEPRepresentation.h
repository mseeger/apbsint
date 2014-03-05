/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class FactorizedEPRepresentation
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_FACTORIZEDEPREPRESENTATION_H
#define EPTOOLS_FACTORIZEDEPREPRESENTATION_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/default.h"

//BEGINNS(eptools)
  /**
   * Represents coupling factor B (sparsity pattern and content) for
   * expectation propagation with factorized backbone. The EP
   * (message) parameters (same size and structure as B) are also
   * maintained here.
   * <p>
   * j=0:(m-1) indexes potentials, i=0:(n-1) indexes variables.
   * 'accessRow' provides access to data for potential j:
   * - Support index V_j (ascending): Nonzeros in B(j,:) are at V_j
   * - Corr. nonzero entries of B(j,:) and EP pars. pi(j,:), beta(j,:),
   *   as flat arrays
   * NOTE: We do not allow B(j,:) to have no nonzeros, for any j. This
   * would lead to failed updates on the j-th potential.
   *
   * 'accessCol' provides access to data for variable i:
   * - Support index V_i (ascending): Nonzeros in B(:,i) are at V_i
   * - Index J_i and array pointers bP, piP, betaP, so that
   *   nonzeros of B(:,i) are bP(J_i), ...
   * NOTE: B(:,i) may be all zero, even though this does not make sense
   * for inference (variable i should be eliminated).
   * <p>
   * Internal representation:
   * This class only maintains and operates on arrays, which are passed
   * from outside upon construction. This "unsafe" interface is to
   * support Matlab MEX, where these arrays are passed from and back to
   * Matlab for each call. Flat arrays:
   * - 'bmatVals': Nonzeros of B
   * - 'betaVals', 'piVals': EP (message) parameters beta_ji, pi_ji
   *
   * Row index:
   * Compatible to basic sparse matrix format. 'rowInd' has two parts,
   * 0:m and (m+1):end. For j=0:(m-1), 'rowInd[j]' is offset into
   * 'bmatVals' (values for row j), 'rowInd[j]+(m+1)' is offset into
   * (2nd part of) 'rowInd' (support index V_j). The size |V_j| is given
   * by 'rowInd[j+1]-rowInd[j]', and 'rowInd[m]' is the number of
   * nonzeros.
   *
   * Column index:
   * 'colInd' has two parts, 0:n and (n+1):end. For i=0:(n-1),
   * 'colInd[i]' is offset into (2nd part of) 'colInd'. The block for
   * column i contains V_i (ascending), followed by J_i (index into
   * 'bmatVals'). Both have size sz_i, and 2*sz_i is given by
   * 'colInd[i+1]-colInd[i]'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactorizedEPRepresentation
  {
  protected:
    // Members

    int numN,numM;                        // Number variables, potentials
    ArrayHandle<int> rowInd;              // "
    ArrayHandle<int> colInd;              // "
    ArrayHandle<double> bmatVals;
    ArrayHandle<double> betaVals,piVals;

  public:
    // Public methods

    /**
     * Constructor. Some simple checks are performed on the data structures,
     * but in particular the consistency between 'rowInd' and 'colInd' is not
     * checked here.
     * NOTE: Arrays are not copied, but referred to. The content of EP parameter
     * vectors is overwritten.
     *
     * @param pnumN     Number variables n
     * @param pnumM     Number potentials m
     * @param prowInd
     * @param pcolInd
     * @param pbmatVals
     * @param pbetaVals
     * @param ppiVals
     */
    FactorizedEPRepresentation(int pnumN,int pnumM,
			       const ArrayHandle<int>& prowInd,
			       const ArrayHandle<int>& pcolInd,
			       const ArrayHandle<double>& pbmatVals,
			       const ArrayHandle<double>& pbetaVals,
			       const ArrayHandle<double>& ppiVals) :
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

    virtual ~FactorizedEPRepresentation() {}

    virtual int numVariables() const {
      return numN;
    }

    virtual int numPotentials() const {
      return numM;
    }

    /**
     * Access to data for potential j.
     * NOTE: Use this for write access to EP parameters.
     * The nonzeros of B(j,:) form a contiguous part in a flat array, same
     * for beta, pi. We return the start offset into this flat array for
     * j.
     *
     * @param j     Potential index
     * @param vjSz  Size |V_j| ret. here
     * @param vjInd Support index V_j
     * @param bP    Nonzeros of B(j,:)
     * @param betaP Nonzeros of beta(j,:)
     * @param piP   Nonzeros of pi(j,:)
     * @return      Offset, s.a.
     */
    virtual int accessRow(int j,int& vjSz,const int*& vjInd,const double*& bP,
			  double*& betaP,double*& piP);

    /**
     * Access to data for variable i.
     *
     * @param i     Variable index
     * @param viInd Support index V_i
     * @param jiInd Index J_i into flat arrays
     * @param bP    Nonzeros of B (flat array)
     * @param betaP Nonzeros of beta (flat array)
     * @param piP   Nonzeros of pi (flat array)
     * @return      Size |V_i|, number nonzeros B(:,i)
     */
    virtual int accessCol(int i,const int*& viInd,const int*& jiInd,
			  const double*& bP,const double*& betaP,
			  const double*& piP);

    /**
     * Compute Gaussian marginals on variables from 'betaVals', 'piVals'.
     * If 'increm'==true, the marginals are added to 'margBeta', 'margPi'.
     *
     * @param margBeta Marginal pars. beta ret. here
     * @param margPi   Marginal pars. pi ret. here
     * @param increm   Incremental? Def.: false
     */
    virtual void compMarginals(double* margBeta,double* margPi,
			       bool increm=false);
  };

  // Inline methods

  inline int
  FactorizedEPRepresentation::accessRow(int j,int& vjSz,const int*& vjInd,
					const double*& bP,double*& betaP,
					double*& piP)
  {
    int jOff;

    if (j<0 || j>=numM) throw InvalidParameterException(EXCEPT_MSG(""));
    jOff=rowInd[j];
    vjSz=rowInd[j+1]-jOff;
    bP=bmatVals.p()+jOff;
    betaP=betaVals.p()+jOff; piP=piVals.p()+jOff;
    vjInd=rowInd.p()+(jOff+numM+1);

    return jOff;
  }

  inline int
  FactorizedEPRepresentation::accessCol(int i,const int*& viInd,
					const int*& jiInd,const double*& bP,
					const double*& betaP,const double*& piP)
  {
    int iOff,viSz;

    if (i<0 || i>=numN) throw InvalidParameterException(EXCEPT_MSG(""));
    iOff=colInd[i];
    viSz=(colInd[i+1]-iOff)>>1;
    viInd=colInd.p()+iOff;
    jiInd=colInd.p()+(iOff+viSz);
    bP=bmatVals.p();
    betaP=betaVals.p(); piP=piVals.p();

    return viSz;
  }

  inline void
  FactorizedEPRepresentation::compMarginals(double* margBeta,double* margPi,
					    bool increm)
  {
    int i,j,jj,viSz;
    double mBeta,mPi;
    const double* bP,*betaP,*piP;
    const int* viInd,*jiInd;

    for (i=0; i<numVariables(); i++) {
      viSz=accessCol(i,viInd,jiInd,bP,betaP,piP);
      for (j=0,mBeta=mPi=0.0; j<viSz; j++) {
	jj=jiInd[j];
	mPi+=piP[jj]; mBeta+=betaP[jj];
      }
      if (!increm) {
	margPi[i]=mPi; margBeta[i]=mBeta;
      } else {
	margPi[i]+=mPi; margBeta[i]+=mBeta;
      }
    }
  }
//ENDNS

#endif
