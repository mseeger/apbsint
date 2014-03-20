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
#include "src/eptools/potentials/PotManagerFactory.h"

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
   * <p>
   * Bivariate precision potentials:
   * If the model contains bivariate precision potentials (see
   * 'EPScalarPotential', group 'atypeBivarPrec'), the representation is
   * extended. The representation here is independent of the potential
   * manager (see 'PotManagerFactory'). The message parameters are 'aVals',
   * 'cVals' (Gamma parameters, see 'EPScalarPotential'). Each such potential
   * j is associated with one precision parameter k=k(j). The flat index
   * 'tauInd' contains the assignment j -> k and its inverse, where both j
   * and k are 0-based.
   * NOTE: These potentials always come last, so j is converted to general
   * potential position by adding m - m_prec, where m_prec is the number of
   * precision potentials. 'tauInd':
   * - Index k(j) [m_prec]
   * - Number K of tau_k entries [1]
   * - For each k=0:(K-1): Start offset of J_k = {j | k(j)==k} [K]
   * - Dummy entry (start offset of J_K if it existed) [1]
   * - J_k, k=0:(K-1), each ascending order [m_prec]
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
    int numK;                             // Number of precision variables
    ArrayHandle<double> aVals,cVals;      // Only if precision potentials
    ArrayHandle<int> tauInd;              // "

  public:
    // Public methods

    /**
     * Constructor (default: no bivariate precision potentials). Some simple
     * checks are performed on the data structures, but in particular the
     * consistency between 'rowInd' and 'colInd' is not checked here.
     * NOTE: Arrays are not copied, but referred to. The content of EP
     * parameter vectors is overwritten.
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
      bmatVals(pbmatVals),betaVals(pbetaVals),piVals(ppiVals),numK(0)
    {
      checkInternalRepres(pnumN,pnumM,prowInd,pcolInd,pbmatVals,pbetaVals,
			  ppiVals);
    }

    /**
     * Constructor (bivariate precision potentials). See header comment
     * for additional parts. 'ptauInd' is checked by
     * 'PotManagerFactory::checkBVPrecTauInd'.
     *
     * @param pnumN     Number variables n
     * @param pnumM     Number potentials m
     * @param prowInd
     * @param pcolInd
     * @param pbmatVals
     * @param pbetaVals
     * @param ppiVals
     * @param paVals
     * @param pcVals
     * @param ptauInd
     */
    FactorizedEPRepresentation(int pnumN,int pnumM,
			       const ArrayHandle<int>& prowInd,
			       const ArrayHandle<int>& pcolInd,
			       const ArrayHandle<double>& pbmatVals,
			       const ArrayHandle<double>& pbetaVals,
			       const ArrayHandle<double>& ppiVals,
			       const ArrayHandle<double>& paVals,
			       const ArrayHandle<double>& pcVals,
			       const ArrayHandle<int>& ptauInd) :
      numN(pnumN),numM(pnumM),rowInd(prowInd),colInd(pcolInd),
      bmatVals(pbmatVals),betaVals(pbetaVals),piVals(ppiVals),aVals(paVals),
      cVals(pcVals),tauInd(ptauInd)
    {
      checkInternalRepres(pnumN,pnumM,prowInd,pcolInd,pbmatVals,pbetaVals,
			  ppiVals);
      int numBVPrec=paVals.size();
      if (numBVPrec>numM)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (numBVPrec==0 || pcVals.size()!=numBVPrec)
	throw InvalidParameterException(EXCEPT_MSG(""));
      PotManagerFactory::checkBVPrecTauInd(ptauInd,numBVPrec);
      numK=ptauInd[numBVPrec];
    }

  private:
    void checkInternalRepres(int pnumN,int pnumM,
			     const ArrayHandle<int>& prowInd,
			     const ArrayHandle<int>& pcolInd,
			     const ArrayHandle<double>& pbmatVals,
			     const ArrayHandle<double>& pbetaVals,
			     const ArrayHandle<double>& ppiVals);
  public:

    virtual ~FactorizedEPRepresentation() {}

    virtual int numVariables() const {
      return numN;
    }

    virtual int numPotentials() const {
      return numM;
    }

    virtual int numBVPrecPotentials() const {
      return aVals.size();
    }

    virtual int numPrecVariables() const {
      return numK;
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

    /**
     * Only if bivar. prec. potentials.
     * Access to precision parameter data for potential j.
     * NOTE: j is absolute potential index, not index into (0-based)
     * 'aVals', ... See header comment.
     * NOTE: Use this for write access to EP parameters.
     *
     * @param j  Potential index
     * @param aP Address for a_jk
     * @param cP Address for c_jk
     * @return   k==k(j)
     */
    virtual int accessTauRow(int j,double*& aP,double*& cP);

    /**
     * Access to precision parameter data for variable tau_k. 'jInd'
     * is J_k, which is also the index into 'aP', 'cP'.
     * NOTE: 'jInd' entries are 0-based, so do not corr. to absolute
     * potential positions (see header comment).
     *
     * @param k    Variable index
     * @param jInd Support index J_k (ascending)
     * @param aP   a message parameters (flat array)
     * @param cP   c message parameters (flat array)
     * @return     Size |J_k|
     */
    virtual int accessTauCol(int k,const int*& jInd,const double*& aP,
			     const double*& cP);

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

  inline  void
  FactorizedEPRepresentation::checkInternalRepres(int pnumN,int pnumM,const ArrayHandle<int>& prowInd,const ArrayHandle<int>& pcolInd,const ArrayHandle<double>& pbmatVals,const ArrayHandle<double>& pbetaVals,const ArrayHandle<double>& ppiVals)
  {
    int j,sz,off,nnz=pbmatVals.size();

    if (pnumN==0 || pnumM==0 || pbetaVals.size()!=nnz ||
	ppiVals.size()!=nnz || prowInd.size()<=pnumM+1 ||
	pcolInd.size()<=pnumN+1)
      throw InvalidParameterException(EXCEPT_MSG(""));
    // Run some basic checks
    if (prowInd[pnumM]!=nnz || prowInd[0]!=0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    for (j=0; j<pnumM; j++) {
      off=prowInd[j]; sz=prowInd[j+1]-off;
      // NOTE: Zero rows are not allowed!
      if (sz<=0 || sz>pnumN)
	throw InvalidParameterException(EXCEPT_MSG(""));
    }
    if (pcolInd[pnumN]!=2*nnz+pnumN+1 || pcolInd[0]!=pnumN+1)
      for (j=0; j<pnumN; j++) {
	off=pcolInd[j]; sz=pcolInd[j+1]-off;
	if (sz%2==1)
	  throw InvalidParameterException(EXCEPT_MSG(""));
	sz/=2;
	// NOTE: Zero columns are allowed
	if (sz<0 || sz>pnumM)
	  throw InvalidParameterException(EXCEPT_MSG(""));
      }
  }

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

  inline int
  FactorizedEPRepresentation::accessTauRow(int j,double*& aP,double*& cP)
  {
    if (numK==0) throw WrongStatusException(EXCEPT_MSG(""));
    int startPos=numM-aVals.size();
    if (j<startPos || j>=numM)
      throw InvalidParameterException(EXCEPT_MSG(""));
    j-=startPos;
    aP=aVals.p()+j; cP=cVals.p()+j;

    return tauInd[j];
  }

  inline int
  FactorizedEPRepresentation::accessTauCol(int k,const int*& jInd,
					   const double*& aP,const double*& cP)
  {
    if (numK==0) throw WrongStatusException(EXCEPT_MSG(""));
    if (k<0 || k>=numK) throw InvalidParameterException(EXCEPT_MSG(""));
    int numBVPrec=aVals.size(); // m_prec
    int off=tauInd[k+numBVPrec+1];
    int sz=tauInd[k+numBVPrec+2]-off;
    jInd=tauInd.p()+off;
    aP=aVals.p(); cP=cVals.p();

    return sz;
  }

  inline void
  FactorizedEPRepresentation::compTauMarginals(double* margA,double* margC,
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
	margA[k]=mA; margC[k]=mC;
      } else {
	margA[k]+=mA; margC[k]+=mC;
      }
    }
  }
//ENDNS

#endif
