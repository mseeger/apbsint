//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class MemWatcher, ArrayHandle
 * ------------------------------------------------------------------- */

#ifndef ARRAYHANDLE_H
#define ARRAYHANDLE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#ifndef HAVE_NO_BLAS
#include "lhotse/matrix/predecl.h" // for friend decl. in 'ArrayHandle'
#include "lhotse/matrix/ArrayUtilsBasic.h"
#endif

/**
 * Every dynamically allocated memory region must be represented by an
 * object of this class. The object maintains a ref. counter which is
 * increased for every reference which is in use to this memory region.
 * The memory region must be deallocated only once this counter drops
 * to zero. This is done in the destructor here.
 * <p>
 * NOTE: Usually, the memory region should be alloc. in the constructor
 * here. This allows the use of non-standard alloc./de-alloc. mechanisms
 * by just replacing code here
 * ==> For compat. reasons, a pointer to a newly alloc. region can be
 *     passed to the constructor. In this case, de-alloc. for this region
 *     MUST work as in the destructor here!
 * <p>
 * Enforcing all this is the responsibility of the participating classes.
 * The major players are 'ArrayHandle' and its subclasses for controlled
 * arrays, and 'BaseVector', 'BaseMatrix' for vector/matrix objects.
 * The main usage is securing mask objects or wrappers which share use
 * of a buffer with other objects.
 * - associate with a buffer: store a pointer to the 'MemWatcher' object
 *   and increase ref. counter: 'incr'
 * - de-assoc. with a buffer: decrease ref. counter: 'decr'. If it drops
 *   to 0, delete the 'MemWatcher' object (which de-alloc. the buffer in its
 *   destructor). Remove pointer to 'MemWatcher' object
 * NOTE: 'BaseVector', 'BaseMatrix' use services of 'MatTimeStamp' for
 * this.
 * <p>
 * NOTE: A user may only control part of the memory region. It typically
 * maintains a separate pointer into the region. Whether it accesses part
 * of the region or not, cannot be checked here.
 * NOTE: Also, the length of the memory region is NOT maintained here.
 * <p>
 * Buffers of 'BaseVector', 'BaseMatrix' objects:
 * A second problem is to secure that mask objects point into valid regions
 * of other objects. This is done using a timestamp mech. ('MatTimeStamp').
 * In order to protect against the ref. object being destroyed, the mech.
 * also use the 'isAssoc' flag here.
 * The flag is typ. false all the time. It is only set to true by the
 * object which owns the buffer (matrix/vector) once it de-assoc. from the
 * buffer. Mask objects will query this flag to decide whether they are
 * still valid. See 'MatTimeStamp'.
 * <p>
 * See doc/memManage.txt for more details.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class MemWatchBase
{
protected:
  int pcount;   // ref. counter
public:
  bool isAssoc; // see above
  //bool debug; // DEBUG!!

public:
  MemWatchBase() : pcount(0),isAssoc(false) {
    //debug=false; // DEBUG!!
  }

  virtual ~MemWatchBase() {
    if (pcount!=0)
      printMsgStdout("ERROR (MemWatchBase::delete): Ref. counter not zero!");
  }

#ifdef DEBUG_TRACKHANDLES
  virtual void incr(int debugCause=0) {
#else
  void incr() {
#endif
    pcount++;
  }

  /**
   * Decrements counter by 1, then returns true iff it is 0
   */
#ifdef DEBUG_TRACKHANDLES
  virtual bool decr(int debugCause=0) {
#else
  bool decr() {
#endif
    --pcount;
    return (pcount==0);
  }

  int getRefCount() const {
    return pcount;
  }
};

template<class T> class MemWatcher : public MemWatchBase
{
protected:
  T* buff;   // pointer to mem. region

public:
  /**
   * Def. constructor. A region of size 'len' is allocated.
   * NOTE: 'len' is not stored in this object!
   * If 'pp' is given, it points to a memory region of size 'len' which is
   * to be owned by the watcher. In this case, 'buff' is init.
   * with 'pp', no memory is alloc.
   *
   * @param len S.a.
   * @param pp  S.a. Def.: 0
   */
  explicit MemWatcher(int len,T* pp=0,uchar debStat=0) : MemWatchBase() {
    if (len<=0)
      throw InvalidParameterException("MemWatcher: 'len' must be positive");
    if (pp!=0) buff=pp;
    else {
      try {
	buff=new T[len];
      } catch (std::bad_alloc ex) {
	string msg("MemWatcher: Cannot allocate block of memory\nByte size: ");
	char sbuff[30];
	sprintf(sbuff,"%ld",((long) len)*((long) sizeof(T)));
	msg+=sbuff;
	throw MemAllocException(msg.c_str());
      }
#ifdef DEBUG_TRACKHANDLES
      debugType temp;
      temp.tag=debStat; temp.sz=len*sizeof(T);
      temp.ptr=this;
      debugMem.insert(pair<int,debugType>((int) buff,temp));
#endif
    }
    incr(); // Incr. ref. counter (set to 1)
  }

  ~MemWatcher() {
#ifdef DEBUG_TRACKHANDLES
    MAP_ITER(int,debugType) it=debugMem.find((int) buff);
    if (it!=debugMem.end())
      debugMem.erase(it);
    else {
      char msg[100];
      sprintf(msg,"*** MemWatcher-destr.: Cannot find %d",(int) buff);
      printMsgStdout(msg);
    }
#endif
    //if (debug) cout << "~MemWatcher" << endl; // DEBUG!!
    delete[] buff;
  }

  T* getBuff() const {
    return buff;
  }

#ifdef DEBUG_TRACKHANDLES
  void incr(int debugCause=0) {
    MemWatchBase::incr(debugCause);
    if (debugCause<0 || debugCause>6)
      throw InternalException("MemWatcher::incr(1)");
    // Increase tag counter (only for tag==8!)
    MAP_ITER(int,debugType) it=debugMem.find((int) buff);
    if (it!=debugMem.end()) {
      if ((*it).second.tag==8)
  	(*it).second.cnt[debugCause]++;
    } else
      throw InternalException("MemWatcher::incr(2)");
  }

  bool decr(int debugCause=0) {
    if (MemWatchBase::decr(debugCause))
      return true;
    else {
      if (debugCause<0 || debugCause>6)
  	throw InternalException("MemWatcher::decr(1)");
      // Decrease tag counter (only for tag==8!)
      MAP_ITER(int,debugType) it=debugMem.find((int) buff);
      if (it!=debugMem.end()) {
  	if ((*it).second.tag==8)
  	  (*it).second.cnt[debugCause]--;
      } else
  	throw InternalException("MemWatcher::decr(2)");
      return false;
    }
  }
#endif
};

/**
 * Provides similar functionality to 'Handle', but for arrays of objects
 * from class T. Not a 'Handle' subclass. Wrapping an array into a 'Handle'
 * does not work!
 * Apart from the advantages of handles, this gives you [] operators with
 * index control, and the 'size' method.
 * NOTE: ALL array allocations must be run through 'ArrayHandle', to assure
 * that non-standard mem. managers are supported!
 * <p>
 * NOTE: The repr. pointer is maintained in 'rep' which seems
 * redundant (stored in mem.watcher 'org' as well)
 * ==> 'rep' can be different from 'org', bec. the handle may use part of a
 *     buffer region only. Ex.: 'BaseVector::getFlatBuff'.
 *     NOTE: For handles created here, 'rep' and 'org' are identical and the
 *     handle uses the whole buffer.
 * <p>
 * NOTE: Subtle issue: 'org' is maintained as 'MemWatchBase*' instead of
 * 'MemWatcher<T>*', although it refers to a 'MemWatcher<T>' object.
 * This is because of the polymorphic features: if T is subclass of T2,
 * an ArrayHandle<T> can be assigned an ArrayHandle<T2>. This would not
 * work if 'org' were 'MemWatcher<T2>*'.
 * <p>
 * NOTE: Once created, the size of an 'ArrayHandle' cannot be changed.
 * This is because many handles may share a given buffer region, even
 * different parts of it ('rep' need not point to the first buffer elem.,
 * and 'len' can be smaller than its length!). If we changed the size,
 * other handles could not be notified.
 * <p>
 * Handles not owning their buffers:
 * If 'org'==0 although 'rep' and 'len' are given, the handle does not own
 * its buffer and will not de-alloc. it upon destruction. There is no ref.
 * counting in this case.
 * ==> Not recommended, only to interface foreign code!!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class ArrayHandle
{
#ifndef HAVE_NO_BLAS
  friend class BaseVector<T>; // for protected constructor
  friend class BaseMatrix<T>; // "
#endif

protected:
  // Members

  T* rep;            // repres. pointer
  int len;           // length of array
  MemWatchBase* org; // ref. counter and mem. region

public:
  // Constructors

  /**
   * Default constructor. Represents the zero handle
   */
  ArrayHandle() : rep(0),len(0),org(0) {}

  /**
   * Constructor, for newly created mem. region, starting at 'pp' and of size
   * 'l'.
   * ATTENTION: Use other constructor 'ArrayHandle(int)' if possible!
   * 'pp' must be s.t. it is de-alloc. by 'MemWatcher' destructor!
   * If 'myOwn'==false, the handle does not own the buffer (see header comm.).
   * NOTE: Use
   *   ArrayHandle<T> arr(l);
   * (below) instead of
   *   ArrayHandle<T> arr(new T[l],l);
   *
   * @param pp    Pointer to representation (newly created)
   * @param l     Array length
   * @param myOwn S.a. Def.: true
   */
  explicit ArrayHandle(T* pp,int l,bool myOwn=true) : rep(pp),len(0),org(0) {
    if (pp!=0) {
      if ((len=l)<=0)
	throw InvalidParameterException("ArrayHandle: Array length must be positive");
      // 'org' points to 'MemWatcher<T>' object, although its type is
      // 'MemWatchBase*'
      if (myOwn) org=new MemWatcher<T>(l,pp);
    }
  }

  /**
   * Constructor. The array (of length 'l') is allocated here. The elements
   * are not init. other than what the T def. constructor does.
   * NOTE: Alloc. is done by 'MemWatcher' constructor.
   *
   * @param l Array length
   */
  explicit ArrayHandle(int l) : rep(0),len(0),org(0) {
    if (l<0)
      throw InvalidParameterException("ArrayHandle: Negative array length");
    else if (l>0) {
      len=l;
      // 'org' points to 'MemWatcher<T>' object, although its type is
      // 'MemWatchBase*'
      MemWatcher<T>* watch=new MemWatcher<T>(l);
      rep=watch->getBuff(); org=watch;
    }
  }

  /**
   * Copy constructor
   *
   * @param r Source object
   */
  ArrayHandle(const ArrayHandle<T>& r) : rep(r.rep),len(r.len),org(r.org) {
#ifndef DEBUG_TRACKHANDLES
    if (org!=0) org->incr();
#else
    if (org!=0) org->incr(1);
#endif
  }

  /**
   * Copy constructor
   * Works if T is strict subclass of T2. Used for implicit conversion.
   * ==> If T no subclass of T2, the assignment of 'rep' breaks (at compile
   *     time)
   *
   * @param r ArrayHandle<T2> to copy
   */
  template<class T2> ArrayHandle(const ArrayHandle<T2>& r) :
  rep(r.p()),len(r.size()),org(r.getMemWatch()) {
#ifndef DEBUG_TRACKHANDLES
    if (org!=0) org->incr();
#else
    if (org!=0) org->incr(1);
#endif
  }

  virtual ~ArrayHandle() {
    deassoc();
  }

  /**
   * Element access
   *
   * @param pos Position in array
   * @return    Ref. to array element
   */
  const T& operator[](int pos) const {
    if (pos<0 || pos>=len) throw OutOfRangeException("ArrayHandle: pos");
    return rep[pos];
  }

  T& operator[](int pos) {
    if (pos<0 || pos>=len) throw OutOfRangeException("ArrayHandle: pos");
    return rep[pos];
  }

  /**
   * Assignment operator. Note that the representation is NOT duplicated
   *
   * @param r Handle rvalue
   * @return  *this (supports chaining)
   */
  ArrayHandle<T>& operator=(const ArrayHandle<T>& r) {
    if (this!=&r) {
      deassoc();
      rep=r.rep; len=r.len; org=r.org;
#ifndef DEBUG_TRACKHANDLES
    if (org!=0) org->incr();
#else
    if (org!=0) org->incr(1);
#endif
    }

    return *this;
  }

  /**
   * Assignment operator. Note that the representation is NOT duplicated
   *
   * @param r Handle rvalue
   * @return  *this (supports chaining)
   */
  template<class T2> ArrayHandle<T>& operator=(const ArrayHandle<T2>& r) {
    if (this!=&r) {
      deassoc();
      rep=r.p(); len=r.size(); org=r.getMemWatch();
#ifndef DEBUG_TRACKHANDLES
      if (org!=0) org->incr();
#else
      if (org!=0) org->incr(1);
#endif
    }

    return *this;
  }

  /**
   * More general variant of 'operator='. This handle is assigned to the
   * part 'off':('off'+'sz'-1) of 'r', or to the part 'off':end (whatever
   * shorter.
   * NOTE: This handle and 'r' must be different objects.
   *
   * @param r   Handle rvalue
   * @param off Start position
   * @param sz  Length of part
   */
  template<class T2> void assign(const ArrayHandle<T2>& r,int off,int sz) {
    if (off>r.size()) throw OutOfRangeException(EXCEPT_MSG(""));
    if (this==&r) throw InvalidParameterException(EXCEPT_MSG(""));
    deassoc();
    if (off<r.size()) {
      rep=r.p()+off; len=std::min(sz,r.size()-off); org=r.getMemWatch();
#ifndef DEBUG_TRACKHANDLES
      if (org!=0) org->incr();
#else
      if (org!=0) org->incr(1);
#endif
    }
  }

  /**
   * Copies the content of 'src' into this array. If 'src' and this array have
   * the same size, this buffer is overwritten. Otherwise, this handle is
   * de-assoc., then assoc. with a new buffer of the correct size (even if
   * this buffer is > 'src' buffer).
   * If 'sz' is given, only the prefix of size 'sz' of 'src' is copied.
   *
   * @param src Source array
   * @param sz  S.a. Optional
   */
  void copy(const ArrayHandle<T>& src,int sz=-1) {
    if (sz==-1) sz=src.size();
    else if (sz<1 || sz>src.size())
      throw OutOfRangeException(EXCEPT_MSG(""));
    if (len!=sz)
      changeRep(sz);
    // ATTENTION: Cannot use 'memmove' here, because T may be complex type
    // with overloaded '='!
    for (int i=0; i<len; i++) rep[i]=src.rep[i];
  }

  /**
   * @return Is this the zero handle?
   */
  bool isZero() const {
    return (rep==0);
  }

  /**
   * Comparison operator. Returns true iff this handle and 'a' ref. to the
   * same object (i.e. have same 'rep' and 'len'), or if they are both the
   * zero handle.
   * NOTE: 'a' can be 0 (int) to repr. the zero handle. 'hand==0' is equiv.
   * to 'hand.isZero()'.
   *
   * @param a S.a.
   */
  virtual bool operator==(const ArrayHandle<T>& a) const {
    return (rep==a.rep && len==a.len);
  }

  template<class T2> bool operator==(const ArrayHandle<T2>& a) const {
    return (rep==a.p() && len==a.size());
  }

  virtual bool operator==(int a) const {
    if (a!=0) throw InvalidParameterException("ArrayHandle: Invalid use of param. 'a'!");
    return (rep==0);
  }

  virtual bool operator!=(const ArrayHandle<T>& a) const {
    return (rep!=a.rep || len!=a.len);
  }

  template<class T2> bool operator!=(const ArrayHandle<T2>& a) const {
    return (rep!=a.p() || len!=a.size());
  }

  virtual bool operator!=(int a) const {
    if (a!=0) throw InvalidParameterException("ArrayHandle: Invalid use of param. 'a'!");
    return (rep!=0);
  }

  /**
   * @return Reference count (0 if this is the 0 handle or the handle does
   *         not own its buffer)
   */
  int getRefCount() const {
    return (org!=0)?org->getRefCount():0;
  }

  /**
   * @return Array length, or 0 if zero handle
   */
  int size() const {
    return len;
  }

  /**
   * Substitute for 'dynamic_cast' (DYNCAST).
   * If T2* can be dyn-casted to T*, we return a conversion of 'r'
   * into ArrayHandle<T>, otherwise return zero.
   *
   * @param r ArrayHandle<T2>
   * @return  ArrayHandle<T>
   */
  template<class T2> static ArrayHandle<T> cast(const ArrayHandle<T2>& r) {
    T* pp=DYNCAST(T,r.p());
    if (pp!=0)
      return ArrayHandle<T>(pp,r.size(),r.getMemWatch());
    else
      return ArrayHandle<T>(); // zero
  }

  /**
   * Conversion operator from handle to (T*). Allows for implicit conversion,
   * a handle to be passed as T* argument.
   * Does not work together with 'dynamic_cast' (DYNCAST), use 'p' instead!
   * Use 'cast' to obtain a dyn-cast for a handle.
   */
  operator T*() const {
    return rep;
  }

  /**
   * @return Pointer to repres., or 0 if zero handle
   */
  T* p() const {
    return rep;
  }

  /**
   * Changes representation. The new representation must not been wrapped
   * by any handles! Never use this method with other repr. than those
   * newly created! To change the handle to become the zero handle, pass
   * the 0 pointer.
   * <p>
   * NOTE: Use 'changeRep(int)' whenever possible! 'pp' must be s.t. it is
   * de-alloc. by destructor code of 'MemWatcher'!
   *
   * @param pp    New representation
   * @param l     Array length of new rep.
   * @param myOwn See constructor. Def.: true
   */
  void changeRep(T* pp,int l,bool myOwn=true) {
    deassoc();
    rep=pp;
    if (pp!=0) {
      if ((len=l)<=0) throw InvalidParameterException("ArrayHandle: Invalid array length");
      // 'org' points to 'MemWatcher<T>' object, although its type is
      // 'MemWatchBase*'
      if (myOwn) org=new MemWatcher<T>(l,pp);
    }
  }

  /**
   * Changes representation. Here, the repr. is allocated inside this
   * method.
   *
   * @param l  Array length of new rep.
   */
  void changeRep(int l) {
    if (l<0)
      throw InvalidParameterException("ArrayHandle: Negative array length");
    deassoc();
    len=l;
    if (l>0) {
      // 'org' points to 'MemWatcher<T>' object, although its type is
      // 'MemWatchBase*'
      MemWatcher<T>* watch=new MemWatcher<T>(l);
      rep=watch->getBuff(); org=watch;
    }
  }

  /**
   * Applies unary function object 'f' to 'n' elements of 'a' (starting from
   * 'aoff'), writing the result into this representation, starting from
   * 'off'. If the result does not fit, an 'OutOfRangeException' is thrown.
   *
   * @param a    Source array
   * @param f    Function object (unary)
   * @param n    S.a. Def.: -1 (all of 'a' from 'aoff')
   * @param off  S.a. Def.: 0
   * @param aoff S.a. Def.: 0
   */
  template<class UnOp,class T2> void
  apply1(const ArrayHandle<T2>& a,const UnOp& f,int n=-1,int off=0,
	 int aoff=0) {
#ifndef HAVE_NO_BLAS
    int na=a.size();

    if (n==-1) n=na-aoff;
    else if (n<0) throw OutOfRangeException(EXCEPT_MSG(""));
    if (off<0 || off+n>len || aoff<0 || aoff+n>na)
      throw OutOfRangeException(EXCEPT_MSG(""));
    ArrayUtilsBasic<T>::applyFunc(rep+off,a.p()+aoff,n,f);
#else
    throw NotImplemException("ArrayHandle::apply1: HAVE_NO_BLAS must not be set");
#endif
  }

  /**
   * Applies binary function object 'f' to 'n' elements of 'a', 'b' (starting
   * from 'aoff', 'boff' resp.), writing the result into this representation,
   * starting from 'off'. If the result does not fit, an
   * 'OutOfRangeException' is thrown.
   *
   * @param a    Source array
   * @param b    Source array
   * @param f    Function object (binary)
   * @param n    S.a. Def.: -1 (all of 'a', 'b', from 'aoff', 'boff')
   * @param off  S.a. Def.: 0
   * @param aoff S.a. Def.: 0
   * @param boff S.a. Def.: 0
   */
  template<class BinOp,class T2,class T3> void
  apply2(const ArrayHandle<T2>& a,const ArrayHandle<T3>& b,const BinOp& f,
	 int n=-1,int off=0,int aoff=0,int boff=0) {
#ifndef HAVE_NO_BLAS
    int na=a.size(),nb=b.size();

    if (n==-1) {
      n=na-aoff;
      if (n!=nb-boff) throw OutOfRangeException(EXCEPT_MSG(""));
    } else if (n<0) throw OutOfRangeException(EXCEPT_MSG(""));
    if (off<0 || off+n>len || aoff<0 || aoff+n>na || boff<0 || boff+n>nb)
      throw OutOfRangeException(EXCEPT_MSG(""));
    ArrayUtilsBasic<T>::applyBinFunc(rep+off,a.p()+aoff,b.p()+boff,n,f);
#else
    throw NotImplemException("ArrayHandle::apply1: HAVE_NO_BLAS must not be set");
#endif
  }

  /**
   * FOR INTERNAL USE ONLY!
   * NOTE: Has to be public, because we are not friend of 'ArrayHandle<T2>'.
   *
   * @return Mem. watcher 'org' used here (or 0 if there is none)
   */
  MemWatchBase* getMemWatch() const {
    return org;
  }

protected:
  // Internal methods

  /**
   * Required by 'cast' method, also by 'BaseVector' (which is friend).
   * 'BaseVector' uses 'MemWatcher' objects as well to guard mem. regions.
   * The newly created array will use watch 'newOrg' (counter is increm. here)
   * instead of creating a new one.
   * NOTE: If 'newOrg'==0, the new handle does not own its buffer (see
   * header comm.)
   *
   * @param repP   Repres. pointer
   * @param lenP   Array length
   * @param newOrg Watcher to be used (counter increm. here)
   */
  ArrayHandle(T* repP,int lenP,MemWatchBase* newOrg) : rep(repP),len(lenP),
  org(newOrg) {
#ifndef DEBUG_TRACKHANDLES
    if (newOrg!=0) newOrg->incr();
#else
    if (newOrg!=0) newOrg->incr(1);
#endif
  }

  /**
   * Destroy representation. Note that 'org' is deleted in 'deassoc'
   * which leads to de-alloc. of the underlying mem. region. Additional
   * clean-ups have to be done here ('org' and the mem. region still exist
   * here!), see 'ArrayPtrHandle' for example.
   */
  virtual void deleteRep() {}

  /**
   * If 'org'!=0, we decrease the ref. counter by 1. If it drops to 0, we
   * destroy 'org' which de-alloc. the mem. region.
   */
  void deassoc() {
#ifdef DEBUG_TRACKHANDLES
    if (org!=0 && org->decr(1)) {
#else
    if (org!=0 && org->decr()) {
#endif
      deleteRep();
      // Although 'org' is 'MemWatchBase*', it really refers to an
      // 'MemWatcher<T2>' object, where T2 is the original element type
      // (T is subclass of T2). The virtual destructor does the right
      // thing.
      delete org; // de-alloc. mem. region
    }
    rep=0; len=0; org=0;
  }
};

/**
 * Manages a zero handle as static member and returns a reference upon call
 * of 'get'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class ArrayHandleZero
{
protected:
  // Static members

  static ArrayHandle<T> zeroH;

public:
  /**
   * @return Reference to zero handle
   */
  static ArrayHandle<T>& get() {
    return zeroH;
  }
};

template<class T> ArrayHandle<T> ArrayHandleZero<T>::zeroH;

#endif
