/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header virtual class MaximumValuesService
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_MAXIMUMVALUESSERVICE_H
#define EPTOOLS_MAXIMUMVALUESSERVICE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <algorithm>

//BEGINNS(eptools)
  /**
   * Data structure for maintaining maximum values in the context of a
   * bipartite graph. Used to drive the selective damping mechanism in
   * factorized EP.
   * <p>
   * Let n variables be indexed by i, m factors indexed by k. The
   * structure is a bipartite factor graph (variable and factor nodes).
   * Each factor k is linked to >=1 variable, each variable is linked
   * to >=1 factor (factor nodes can be restricted, see 'subInd').
   * For variable i, denote by V_i the set of factors connected to i.
   * Each link (k,i) is associated a variable x_ki, these values change
   * all the time (the graph structure is fixed).
   * The object tracks max_k x_ki for each variable i.
   *
   * For K=='maxSize', we keep up to K entries (x_ki,k) for each i. The
   * valid entries (between 1 and K) correspond to the maximum ones, sorted
   * in descending order. When a new (x_ki,k) comes in, this list is
   * updated, whereby it can shrink by one entry. If it becomes empty, it
   * is recomputed (for i only).
   *
   * Read access to graph structure and values (V_i, {x_ki}) is via pure
   * virtual methods 'numVariables', 'numFactors', 'getFactorValues'. These
   * must be implemented by subclasses (f.ex., 'FactEPMaximumPiValues').
   * <p>
   * Top-K lists maintained in arrays 'topVal', 'topInd' (size n*(K+1))
   * each, entries for i start at i*(K+1), first 'numValid[i]' are valid.
   * The last entry is a dummy entry.
   * NOTE: 'numValid[i]' must not be 0 for any i. This means that every
   * variable i must be touched by at least one factor k not excluded
   * by 'subInd'.
   * <p>
   * Maximum over subset of factors:
   * If 'subInd' is given, max_k x_ki does not run over all k. If
   * 'subExcl'==false, max_k runs over 'subInd'. If 'subExcl'==true, max_k
   * runs over the complement of 'subInd'. 'subInd' must be sorted in
   * ascending order.
   * NOTE: Involves binary search over 'subInd' for every 'update' and
   * 'recompute' call, so choose 'subExcl' to keep 'subInd' small.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MaximumValuesService
  {
  protected:
    // Members

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
     * Constructor. Arrays are not copied, no consistency checks.
     * To compute top-K lists from scratch, pass arbitrary arrays, then
     * call 'recompute'.
     *
     * @param pn        Number of variables n
     * @param pm        Number of factors m
     * @param pmaxSize  K
     * @param pnumValid Entries must be in 1:pmaxSize
     * @param ptopInd
     * @param ptopVal
     * @param psubInd   Optional
     * @param psubExcl  Def.: false
     */
    MaximumValuesService(int pn,int pm,int pmaxSize,
			 const ArrayHandle<int>& pnumValid,
			 const ArrayHandle<int>& ptopInd,
			 const ArrayHandle<double>& ptopVal,
			 const ArrayHandle<int>& psubInd=
			 ArrayHandleZero<int>::get(),bool psubExcl=false) :
      maxSize(pmaxSize),numValid(pnumValid),topInd(ptopInd),topVal(ptopVal),
      subInd(psubInd),subExcl(psubExcl) {
      if (pmaxSize<1 || pnumValid.size()!=pn ||
	  ptopInd.size()!=pn*(pmaxSize+1) || ptopVal.size()!=ptopInd.size())
	throw InvalidParameterException(EXCEPT_MSG(""));
      Interval<int> ivK(1,pmaxSize,IntVal::ivClosed,IntVal::ivClosed);
      if (ivK.check(pnumValid.p(),pn)!=0)
	throw InvalidParameterException(EXCEPT_MSG("pnumValid: Entries out of range"));
      if (!(psubInd==0)) {
	int sz=psubInd.size();
	if (!Range::isIncreasing(psubInd.p(),sz))
	  throw InvalidParameterException(EXCEPT_MSG("psubInd must be sorted in ascending order"));
	if (psubInd[0]<0 || psubInd[sz-1]>=pm)
	  throw InvalidParameterException(EXCEPT_MSG("psubInd: Out of range"));
	if ((!psubExcl && sz<pmaxSize) || (psubExcl && pm-sz<pmaxSize))
	  throw InvalidParameterException(EXCEPT_MSG("psubInd: Too small"));
      }
      resetStats();
    }

    virtual ~MaximumValuesService() {}

    /**
     * Has to be implemented by subclasses
     *
     * @return Number n of variables
     */
    virtual int numVariables() const = 0;

    /**
     * Has to be implemented by subclasses
     *
     * @return Number m of factors
     */
    virtual int numFactors() const = 0;

    /**
     * Has to be implemented by subclasses.
     * For variable index i, V_i contains factors k which connect to i, V_i
     * is sorted in increasing order. J_i is an index into the flat vector
     * 'xarr': if k==V_i[l], then x_ki is in 'xarr[J_i[l]]'. The size of
     * V_i, J_i is returned.
     *
     * @param i    Variable index
     * @param vind Index V_i
     * @param jind Index J_i
     * @param xarr Flat array for x values
     * @return     Length of V_i, J_i
     */
    virtual int getFactorValues(int i,const int*& vind,const int*& jind,
				const double*& xarr) const = 0;

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
      for (int i=0; i<numVariables(); i++)
	recompute(i);
    }

    /**
     * @param i Variable index
     * @return  max_k x_ki
     */
    virtual double getMaxValue(int i) const {
      return topVal[i*(maxSize+1)];
    }

    /**
     * Notification that x_ji has a new value. This must already have been
     * written back (so that 'getFactorValues' is up-2-date).
     * NOTE: Call this method directly after writing back single x_ji
     * values.
     * j must not be excluded by 'subInd'. This is not checked!
     *
     * @param i   Variable index
     * @param j   Factor index
     * @param val New value x_ji (passed for convenience)
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

  inline void MaximumValuesService::recompute(int i)
  {
    int j,jj,k,viSz;
    const double* xP;
    const int* viInd,*jiInd;

    viSz=getFactorValues(i,viInd,jiInd,xP);
    numValid[i]=0;
    for (k=0; k<viSz; k++) {
      jj=jiInd[k]; j=viInd[k];
      // Skip j if excluded by 'subInd'
      if (!(subInd==0) &&
	  std::binary_search(subInd.p(),subInd.p()+subInd.size(),j)==subExcl)
	continue;
      insertEntry(i,j,xP[jj]);
    }
    if (numValid[i]==0)
      throw WrongStatusException(EXCEPT_MSG("Cannot have numValid[i]==0. Representation invalid now!"));
  }

  inline void MaximumValuesService::update(int i,int j,double val)
  {
    if (i<0 || j<0 || i>=numVariables() || j>=numFactors())
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (val<=topVal[i*(maxSize+1)+numValid[i]-1]) {
      // New x_ji smaller than other list entries
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

  inline void MaximumValuesService::insertEntry(int i,int j,double val)
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

  inline bool MaximumValuesService::removeEntry(int i,int j)
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
