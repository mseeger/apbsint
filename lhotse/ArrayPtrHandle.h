//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class ArrayPtrHandle
 * ------------------------------------------------------------------- */

#ifndef ARRAYPTRHANDLE_H
#define ARRAYPTRHANDLE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/ArrayHandle.h"

/**
 * Extends 'ArrayHandle<T*>' by the functionality of owning all or certain
 * array elements (in the sense that the destructor here destroys these
 * elements as well). By default, the objects owns all array elements.
 * If the array should act as mask only, i.e. NOT owning any elements,
 * 'ArrayHandle<T*>' should be used.
 * <p>
 * Deallocation of the representation:
 * By def., this means that the destructors for all array elements are
 * called. However, it is possible upon construction or via 'changeRep'
 * to pass a boolean array of the same length of the array. In this case,
 * an element in the array is considered "our own" iff the corr. entry
 * in the boolean array is true. Now, deallocation means that the constructors
 * of only those elements are called that are "our own". We draw a copy of
 * the boolean array and store it in 'ourOwn'. Like with 'org', handles
 * ref. to the same repr. share this array. If 'ourOwn' is 0, all elements
 * are "our own" (default).
 * <p>
 * NOTE: 'ArrayPtrHandle' objects (or 'ArrayHandle<T*>') are often used to
 * pass arrays masking out elements from other structures. Pass
 * 'ArrayHandle<T*>' in this case, because the "owner" status does not matter.
 * <p>
 * ATTENTION:
 * Deassoc. from a buffer has to be done here explicitely by calling
 * 'deassoc'. If a superclass method is called, it will call its own
 * 'deassoc' and 'deleteRep' which does NOT destroy elements of 'rep'!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class ArrayPtrHandle : public ArrayHandle<T*>
{
protected:
  // Additional members

  ArrayHandle<bool> ourOwn; // see header comment (0 -> own them all)

public:
  // Constructors

  /**
   * Default constructor. Represents the zero handle
   */
  ArrayPtrHandle() : ArrayHandle<T*>() {}

  /**
   * Constructor. 'pp' is an array of T* elements. 'pp' itself is copied here
   * and does not have to be retained. Entries of 'pp' belong to us dep. on
   * the corr. 'oOwn' entries. If 'oOwn'==0, then all 'pp' entries belong to
   * us. A typical creation is:
   *   ArrayHandle<T*> pp(l);
   *   for (i=0; i<l; i++) {
   *     pp[i]=new SUB_OF_T(...,parDepenOn[i],...);
   *   }
   *   ArrayPtrHandle hand(pp);
   * In this example, 'oOwn'==0 and all elements belong to us.
   *
   * @param pp   Array
   * @param oOwn Optional. Array for 'ourOwn' (see header comment)
   */
  explicit ArrayPtrHandle(const ArrayHandle<T*>& pp,const bool* oOwn=0) :
  ArrayHandle<T*>(pp.size()) {
    int l=pp.size();
    if (l>0) {
      memmove(this->rep,pp,l*sizeof(T*));
      if (oOwn!=0) {
	ourOwn.changeRep(l);
	memmove(ourOwn.p(),oOwn,l*sizeof(bool));
      }
    }
  }

  /**
   * Copy constructor.
   *
   * @param r Handle<T> to copy
   */
  ArrayPtrHandle(const ArrayPtrHandle<T>& r) :
  ArrayHandle<T*>(r),ourOwn(r.ourOwn) {}

  /**
   * Copy constructor.
   * Works if T is strict subclass of T2. Used for implicit conversion
   *
   * @param r Handle<T2> to copy
   */
  template<class T2> ArrayPtrHandle(const ArrayPtrHandle<T2>& r) :
  ArrayHandle<T*>(r),ourOwn(r.getOurOwnInternal()) {}

  /**
   * Have to call 'deassoc' here in order to have the repres. be
   * properly destroyed (superclass. destr. would not do this).
   */
  ~ArrayPtrHandle() {
    this->deassoc();
  }

  /**
   * Assignment operator. Note that the representation is NOT duplicated.
   * Works if T is strict subclass of T2.
   *
   * @param r Handle rvalue
   * @return  *this (supports chaining)
   */
  template<class T2> ArrayPtrHandle<T>& operator=(const ArrayPtrHandle<T2>& r)
  {
    if (this!=&r) {
      this->deassoc();
      ArrayHandle<T*>::operator=(r);
      ourOwn=r.getOurOwnInternal();
    }

    return *this;
  }

  /**
   * Assignment operator. Note that the representation is NOT duplicated.
   *
   * @param r Handle rvalue
   * @return  *this (supports chaining)
   */
  ArrayPtrHandle<T>& operator=(const ArrayPtrHandle<T>& r)
  {
    if (this!=&r) {
      this->deassoc();
      ArrayHandle<T*>::operator=(r);
      ourOwn=r.ourOwn;
    }

    return *this;
  }

  /**
   * Substitute for 'dynamic_cast'.
   * If T2 can be dyn-casted to T, we return a conversion of 'r'
   * into ArrayPtrHandle<T>, otherwise return zero.
   *
   * @param r ArrayPtrHandle<T2>
   * @return  ArrayPtrHandle<T>
   */
  template<class T2> static ArrayPtrHandle<T> cast(const ArrayPtrHandle<T2>& r)
  {
    T* pp=DYNCAST(T,r.p());
    if (pp!=0)
      return ArrayPtrHandle<T>(pp,r.size(),r.getOrgInternal(),
			       r.getOurOwnInternal());
    else
      return ArrayPtrHandle<T>(); // zero
  }

  /**
   * Changes representation to 'pp'. The same comments as given in the
   * constructor above, apply to 'pp' here: the objects 'pp' points to
   * must have been newly created, and must not be referred to from somewhere
   * else. The array 'pp' itself is copied here.
   *
   * @param pp   New representation
   * @param oOwn Optional. Array for 'ourOwn' (see header comment)
   */
  void changeRep(const ArrayHandle<T*>& pp,const bool* oOwn=0) {
    int l=pp.size();
    this->deassoc(); ourOwn.changeRep(0);
    if (l>0) {
      ArrayHandle<T*>::changeRep(l);
      memmove(this->rep,pp,l*sizeof(T*));
      if (oOwn!=0) {
	ourOwn.changeRep(l);
	memmove(ourOwn.p(),oOwn,l*sizeof(bool));
      }
    }
  }

  void changeRep(int l) {
    if (l!=0)
      throw InternalException("'ArrayPtrHandle' does not allow default initialisation!");
    this->deassoc(); ourOwn.changeRep(0); // zero handle
  }

  /**
   * FOR INTERNAL USE ONLY!
   * NOTE: Has to be public, we are not friend of 'ArrayPtrHandle<T2*>'.
   *
   * @return 'ourOwn' member
   */
  const ArrayHandle<bool>& getOurOwnInternal() const {
    return ourOwn;
  }

protected:
  // Internal methods

  /**
   * Required by 'cast' method.
   * The newly created array will use watch 'newOrg' (counter is increm. here)
   * instead of creating a new one.
   *
   * @param repP   Repres. pointer
   * @param lenP   Array length
   * @param newOrg Watcher to be used (counter increm. here)
   * @param ownP   'ourOwn' array
   */
  ArrayPtrHandle(T* repP,int lenP,MemWatchBase* newOrg,
		 const ArrayHandle<bool>& ownP) :
  ArrayHandle<T*>(repP,lenP,newOrg),ourOwn(ownP) {}

  /**
   * Also have to destroy array elements, dep. on 'ourOwn'. 'ourOwn' itself
   * is also de-alloc.
   */
  void deleteRep() {
    int i;
    if (ourOwn.isZero()) {
      for (i=0; i<this->len; i++) delete this->rep[i];
    } else {
      for (i=0; i<this->len; i++)
	if (ourOwn[i]) delete this->rep[i];
      ourOwn.changeRep(0);
    }
  }
};

/**
 * Manages a zero handle as static member and returns a reference upon call
 * of 'get'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class ArrayPtrHandleZero
{
protected:
  // Static members

  static ArrayPtrHandle<T> zeroH;

public:
  /**
   * @return Reference to zero handle
   */
  static ArrayPtrHandle<T>& get() {
    return zeroH;
  }
};

template<class T> ArrayPtrHandle<T> ArrayPtrHandleZero<T>::zeroH;

#endif
