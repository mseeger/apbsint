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

#include <algorithm>
#include "src/eptools/FactorizedEPRepresentation.h"

//BEGINNS(eptools)
  /**
   * Data structure for selective damping mechanism in
   * 'FactorizedEPDriver'.
   * <p>
   * Goal is to track max_k pi_ki for every variable i. The pi_ki change
   * all the time. Should be simple, and max should have to be recomputed
   * seldomly (needs run over V_i, which can be large).
   *
   * For K=='maxSize', we keep up to K entries (pi_ki,k) for each i. The
   * valid entries (between 1 and K) correspond to the maximum ones, sorted
   * in descending order. When a new (pi_ji,j) comes in, this list is
   * updated, whereby it can shrink by one entry. If it becomes empty, it
   * is recomputed (for i only).
   *
   * The pi_ji and the structure of the B coupling factor are maintained
   * in 'epRepr' ('FactorizedEPRepresentation').
   * <p>
   * Top-K lists maintained in arrays 'topVal', 'topInd' (size n*(K+1))
   * each, entries for i start at i*(K+1), first 'numValid[i]' are valid.
   * The last entry is a dummy entry.
   * NOTE: 'numValid[i]' must not be 0 for any i. This means that every
   * variable i must be touched by at least one potential j not excluded
   * by 'subInd'.
   *
   * NOTE: We use simple flat vectors here, instead of structs. This is to
   * simplify transfer with Matlab without copying.
   * <p>
   * Maximum over subset:
   * If 'subInd' is given, max_k pi_ki does not run over all k. If
   * 'subExcl'==false, max_k runs over 'subInd'. If 'subExcl'==true, max_k
   * runs over the complement of 'subInd'. 'subInd' must be sorted in
   * ascending order.
   * NOTE: Involves binary search over 'subInd' for every 'update' and
   * 'recompute' call, so choose 'subExcl' to keep 'subInd' small.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FactEPMaximumPiValues
  {
  protected:
    // Members

    Handle<FactorizedEPRepresentation> epRepr;
    int maxSize;
    ArrayHandle<int> numValid;
    ArrayHandle<int> topInd;
    ArrayHandle<double> topVal;
    ArrayHandle<int> subInd;
    bool subExcl;
    int statNUpd,statNRec;

  public:
    // Public methods

    /**
     * Constructor. Consistency of 'ptopVal' with 'pepRepr' is not checked.
     * To compute top-K lists from scratch, pass arbitrary arrays, then
     * call 'recompute'.
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
      epRepr(pepRepr),maxSize(pmaxSize),numValid(pnumValid),topInd(ptopInd),
      topVal(ptopVal),subInd(psubInd),subExcl(psubExcl) {
      int numN=pepRepr->numVariables();

      if (pmaxSize<1 || pnumValid.size()!=numN ||
	  ptopInd.size()!=numN*(pmaxSize+1) || ptopVal.size()!=ptopInd.size())
	throw InvalidParameterException(EXCEPT_MSG(""));
      Interval<int> ivK(1,pmaxSize,IntVal::ivClosed,IntVal::ivClosed);
      if (ivK.check(pnumValid.p(),numN)!=0)
	throw InvalidParameterException(EXCEPT_MSG("pnumValid: Entries out of range"));
      if (!(psubInd==0)) {
	int sz=psubInd.size(),numM=pepRepr->numPotentials();
	if (!Range::isIncreasing(psubInd.p(),sz))
	  throw InvalidParameterException(EXCEPT_MSG("psubInd must be sorted in ascending order"));
	if (psubInd[0]<0 || psubInd[sz-1]>=numM)
	  throw InvalidParameterException(EXCEPT_MSG("psubInd: Out of range"));
	if ((!psubExcl && sz<pmaxSize) || (psubExcl && numM-sz<pmaxSize))
	  throw InvalidParameterException(EXCEPT_MSG("psubInd: Too small"));
      }
      resetStats();
    }

    virtual ~FactEPMaximumPiValues() {}

    /**
     * Recompute top-K list for variable i. If i is not used, all top-K
     * lists are recomputed.
     * NOTE: 'numValid[i]' must be >0 afterwards (exception is thrown
     * otherwise).
     *
     * @param i Variable index. Optional
     */
    virtual void recompute(int i);

    virtual void recompute() {
      for (int i=0; i<epRepr->numVariables(); i++)
	recompute(i);
    }

    /**
     * @param i Variable index
     * @return  max_k pi_ki
     */
    virtual double getMaxBeta(int i) const {
      return topVal[i*(maxSize+1)];
    }

    /**
     * Notification that pi_ji has a new value. This must already have been
     * written back ('epRepr').
     * NOTE: Call this method directly after writing back single pi_ji
     * values.
     * j must not be excluded by 'subInd'. This is not checked!
     *
     * @param i   Variable index
     * @param j   Potential index
     * @param val New value pi_ji (passed for convenience)
     */
    virtual void update(int i,int j,double val);

    /**
     * Returns statistics collected so far.
     *
     * @param nupd Number of calls to 'update'
     * @param nrec Number of 'recompute' calls done by 'update'
     */
    virtual void getStats(int& nupd,int& nrec) const {
      nupd=statNUpd; nrec=statNRec;
    }

    virtual void resetStats() {
      statNUpd=statNRec=0;
    }

  protected:
    // Helper methods

    /**
     * Insert entry (val,j) into top-K list for i. Assumes that j is not
     * in 'topInd' for i, and that j is not excluded by 'subInd'.
     */
    void insertEntry(int i,int j,double val);

    /**
     * Check whether j is in 'topInd' for i. If so, remove the corr. entry.
     * This may leave the top-K list for i empty. Assumes that j is not
     * excluded by 'subInd'.
     *
     * @return Was j found (and therefore removed)?
     */
    bool removeEntry(int i,int j);
  };

  // Inline methods

  inline void FactEPMaximumPiValues::recompute(int i)
  {
    int j,jj,k,viSz;
    const double* bP,*betaP,*piP;
    const int* viInd,*jiInd;

    viSz=epRepr->accessCol(i,viInd,jiInd,bP,betaP,piP);
    numValid[i]=0;
    for (k=0; k<viSz; k++) {
      jj=jiInd[k]; j=viInd[k];
      // Skip j if excluded by 'subInd'
      if (!(subInd==0) &&
	  std::binary_search(subInd.p(),subInd.p()+subInd.size(),j)==subExcl)
	continue;
      insertEntry(i,j,piP[jj]);
    }
    if (numValid[i]==0)
      throw WrongStatusException(EXCEPT_MSG("Cannot have numValid[i]==0. Representation invalid now!"));
  }

  inline void FactEPMaximumPiValues::update(int i,int j,double val)
  {
    if (i<0 || j<0 || i>=epRepr->numVariables() ||
	j>=epRepr->numPotentials())
      throw InvalidParameterException(EXCEPT_MSG(""));
    //if (!(subInd==0) &&
    //std::binary_search(subInd.p(),subInd.p()+subInd.size(),j)==subExcl)
    //  throw InvalidParameterException(EXCEPT_MSG("Cannot update for excluded j"));
    if (val<=topVal[i*(maxSize+1)+numValid[i]-1]) {
      // New pi_ji smaller than other list entries
      if (removeEntry(i,j)) {
	// If top-K list is empty: Have to recompute
	if (numValid[i]==0) {
	  recompute(i);
	  statNRec++;
	}
      }
    } else {
      // New entry has to be inserted into top-K list
      removeEntry(i,j); // Remove entry for j, if present
      insertEntry(i,j,val);
    }
    statNUpd++;
  }

  inline void FactEPMaximumPiValues::insertEntry(int i,int j,double val)
  {
    int k,num=numValid[i],cpj;
    double cpv;
    int* tiP;
    double* tvP;

    k=i*(maxSize+1);
    tiP=topInd.p()+k; tvP=topVal.p()+k;
    if (num==maxSize && val<=tvP[maxSize-1])
      return; // 'val' smaller than all others
    for (k=0; k<num && val<=tvP[k]; k++);
    // Needs dummy entry at end:
    for (; k<=num; k++) {
      cpv=tvP[k]; cpj=tiP[k];
      tvP[k]=val; tiP[k]=j;
      val=cpv; j=cpj;
    }
    if (num<maxSize)
      numValid[i]++; // Increase list size
  }

  inline bool FactEPMaximumPiValues::removeEntry(int i,int j)
  {
    int k,num=numValid[i];
    int* tiP;
    double* tvP;

    MYASS(num>0);
    k=i*(maxSize+1);
    tiP=topInd.p()+k; tvP=topVal.p()+k;
    for (k=0; k<num && tiP[k]!=j; k++);
    if (k==num)
      return false; // j not in list
    // Needs dummy entry at end:
    for (; k<num-1; k++) {
      tiP[k]=tiP[k+1]; tvP[k]=tvP[k+1];
    }
    numValid[i]--;

    return true;
  }
//ENDNS

#endif
