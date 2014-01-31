//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class Handle
 * ------------------------------------------------------------------- */

#ifndef HANDLE_H
#define HANDLE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Helper class for 'Handle'. Manages ref. counter and 'ourOwn' flag without
 * being a template class.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class HandleHelper
{
public:
  int pcount;  // ref. counter
  bool ourOwn; // does rep. belong to us (do we have to dealloc.)?

  HandleHelper(bool oOwn=true) : pcount(1),ourOwn(oOwn) {}
};

/**
 * Handles (smart pointers) behave much like pointers, but objects they
 * refer to are destroyed once no more handle refers to them.
 * The object a handle ref. to is called representation. The zero handle
 * is the only handle without a representation.
 * <p>
 * Each representation is associated with a ref. counter ('HandleHelper')
 * storing the number of handles assoc. with the object. The counter is
 * changed together with changes to its handles. Once it drops to 0, the
 * representation is destroyed.
 * <p>
 * Efficiency:
 * Handles are small objects to which the fixed memory manager
 * 'FixedMemManager' applies. Except in very special cases involving
 * structures with very many pointers, handles are as efficient as
 * pointers.
 * <p>
 * Rules for using handles (make your life much easier!):
 * - Use handles to newly alloc. objects by default.
 *   LHOTSE code returns handles to newly created objects. In very special
 *   cases, the pointer can be freed from the handle using 'deassoc' (DO NOT
 *   DO THIS).
 *     Handle<T> hobj(new T(...)); hobj->meth(...);
 *   instead of
 *     T* pobj=new T(...); pobj->meth(...);
 * - Use dynamic allocation and handles for all objects which are later on
 *   referred to by others.
 *     Handle<T> hobj(new T(...));
 *     // ...
 *     i_ref_to_my_arg(hobj);
 *     // ...
 *     // 'hobj' goes out of scope, but if 'i_ref_to_my_arg' still has a
 *     // handle ref., the object is not dealloc.
 *   instead of
 *     T obj(...);
 *     // ...
 *     i_ref_to_my_arg(obj);
 *     // ...
 *     // 'obj' goes out of scope and is dealloc., even though
 *     // 'i_ref_to_my_arg' still has a ref. BOOM!
 * - NEVER create an object non-dynamically, then wrap it into a handle!
 *     T obj(...);
 *     Handle<T> hobj(&obj); // NO!!!!
 * - NEVER wrap on object into two different handles! In general, do not
 *   store a pointer to a representation apart from the handle. To have
 *   handles ha, hb ref. to a, wrap a into ha, then: 'hb=ha' (like copying
 *   a pointer).
 *     Handle<T> hobj(new T(...));
 *     T* pobj=hobj.p(); // DON'T DO THIS
 *     Handle<T> hobj2(pobj); // NO!!!
 * - Never pass pointers to methods, pass handles instead. The argument type
 *   should be 'const Handle<T>&'. This allows polymorphism via the conversion
 *   operator defined here: suppose B is subclass of A and a Handle<B> obj. hb
 *   is passed for a 'const Handle<A>&'. The conversion operator will create a
 *   temp. Handle<A> object with the repr. of hb. This works iff B* can be ass.
 *   to A*.
 *   NOTE: The const here does not mean that a->f(...) fails for the argument
 *   if 'f' is not const, just that the Handle object 'a' itself must not be
 *   changed (e.g., 'a.changeRep(...)' would fail).
 * - Use handles as safe references inside objects.
 *     Handle<T> hobj;
 *   Even if the object 'hobj' refers to goes out of its scope, it will not
 *   be dealloc. before 'hobj' itself is.
 * - A handle is used like a pointer. The operator '->' is implemented, as is
 *   operator*. 'p' returns the pointer to the repres. In most cases,
 *   'func(hobj)' is equivalent to 'func(hobj.p())', if 'func' requires a
 *   pointer arg., but 'DYNCAST(T2,hobj)' does not work: use
 *   'DYNCAST(T2,hobj.p())'.
 * - Do not use 'Handle' for arrays! Use 'ArrayHandle'
 * <p>
 * The zero handle:
 * The equivalent of a zero pointer is the zero handle: it justs represents
 * nothing. The default constructor constructs a zero handle. An existing
 * non-zero handle can be made a zero handle using 'changeRep' and passing 0.
 * A handle can be tested for zero using 'isZero', or the expression
 * 'hobj==0'.
 * <p>
 * Handles not owning their representation:
 * Sometimes, objects have to be wrapped in handles even if they have not
 * been allocated dynamically, because the method to be called requires a
 * handle. This class supports this by the 'myOwn' flag.
 * If this is false, we do not own the repres. and do not touch it even if
 * the ref. count drops to 0.
 * NOTE: Such handles can jeopardise the whole idea behind 'Handle'. Imagine
 * an object being ref. to by owning and non-owning handles, and the handles
 * being destroyed in an ordering which makes a non-owning the last. The
 * repres. would not be destroyed then: a memory leak.
 * ==> Use non-owning handles only temporarily within a fixed score. Do NOT
 *     store them as members of a class, for example. AVOID if possible.
 * <p>
 * Deassociating a handle:
 * If a is ref. to by the single handle ha (ref. count is 1), the assoc. can
 * be removed by
 *   aptr=ha.p();
 *   ha.deassoc();
 * ha becomes the zero handle, but importantly a is not destroyed.
 * NOTE: Use deassociation only in very special situations, namely if object
 * a is inserted into a structure which guarantees the cleanup, and where
 * storing handles for the cleanup would be inefficient.
 * <p>
 * Handles are like pointers:
 * Let B be subclass of A, ha Handle<A>, hb Handle<B>:
 * - 'ha=hb' works because 'operator=' is a template member
 * - if f(const Handle<A>& arg) is a function/method, then 'f(hb)'
 *   works because a template member conversion operator creates a temp.
 *   Handle<A> copy of hb
 * - 'cast' is the analogue of DYNCAST for pointers:
 *   'Handle<B>::cast(ha)' returns a Handle<B> copy of ha iff ha's repres.
 *   can be dyn-casted to B*, otherwise a zero handle
 *   NOTE: 'cast' can have the effect that handles of different types T
 *   refer to the same repres. (using the same counter object) if the type
 *   of the repres. is a common superclass. OK, because the destructor is
 *   polymorphic
 * - '->' behaves polymorphic as with pointers
 * NOTE: These methods req. protected members access between Handle<T>
 * objects of different types. Since we don't know how to declare such
 * a friend relation, we have to make the req. access methods public:
 * 'getXXXInternal'. DO NOT USE THEM!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class Handle
{
protected:
  // Variables

  T* rep;            // pointer to representation of Handle
  HandleHelper* org; // ref. counter and 'ourOwn' flag

public:
  // Constructors

  /**
   * Default constructor. Represents the zero pointer
   */
  Handle() : rep(0),org(0) {}

  /**
   * Constructor.
   * NOTE: This constructor is declared 'explicit'. A normal declaration
   * allows implicit use, which leads to really bad bugs! For example, if
   * a function 'func(const Handle<T>& h,...)' is called passing a T* for
   * 'h' instead of a handle, this constructor would be invoked to wrap
   * the T* into a handle OWNING IT. Upon return of 'func', the handle is
   * destroyed, destroying the T* object as well. This is usually not
   * intended.
   *
   * @param pp   Pointer to representation (newly created by new)
   * @param oOwn Do we own the repres., i.e. are we to destroy it in the end?
   *             Def.: true
   */
  explicit Handle(T* pp,bool oOwn=true) : rep(pp),org(0) {
    if (pp!=0) org=new HandleHelper(oOwn); // ref. count init. with 1
  }

  /**
   * Copy constructor
   *
   * @param r Handle<T> to copy
   */
  Handle(const Handle<T>& r) : rep(r.rep),org(r.org) {
    if (org!=0) org->pcount++;
  }

  /**
   * Copy constructor
   * Works if T2 is strict subclass of T, but NOT for T == T2 (who knows
   * why??). Used for implicit conversion
   *
   * @param r Handle<T2> to copy
   */
  template<class T2> Handle(const Handle<T2>& r) :
    rep(r.p()),org(r.getOrgInternal()) {
    if (org!=0) org->pcount++;
  }

  /**
   * The representation is deallocated iff the ref. counter is zero
   */
  virtual ~Handle() {
    decrRefCount();
  }

  /**
   * This operators renders the handle the same behaviour as a pointer
   *
   * @return Pointer to underlying reference
   */
  T* operator->() const {
    if (rep==0)
      throw WrongStatusException(EXCEPT_MSG("This is the zero handle"));
    return rep;
  }

  /**
   * Assignment operator. Note that the representation is NOT duplicated.
   *
   * @param r Handle rvalue
   * @return  *this (supports chaining)
   */
  Handle<T>& operator=(const Handle<T>& r) {
    if (this!=&r) {
      decrRefCount();
      rep=r.rep; org=r.org;
      if (org!=0) org->pcount++;
    }

    return *this;
  }

  /**
   * Assignment operator. Note that the representation is NOT duplicated.
   * Works if T2 is strict subclass of T, but NOT for T == T2 (who knows
   * why??)
   *
   * @param r Handle rvalue
   * @return  *this (supports chaining)
   */
  template<class T2> Handle<T>& operator=(const Handle<T2>& r) {
    decrRefCount();
    rep=r.p(); org=r.getOrgInternal();
    if (org!=0) org->pcount++;

    return *this;
  }

  /**
   * Just to allow resetting a handle to zero by
   *   hand=0;
   * which is equivalent to
   *  hand.changeRep(0);
   * If the int argument is !=0, an exception is thrown.
   *
   * @return  *this (supports chaining)
   */
  Handle<T>& operator=(int dummy) {
    if (dummy!=0) throw InvalidParameterException(EXCEPT_MSG(""));
    decrRefCount();
    org=0; rep=0;
  }

  /**
   * @return Is this a zero pointer?
   */
  bool isZero() const {
    return (rep==0);
  }

  /**
   * Comparison operator. Returns true iff this handle and 'a' ref. to the
   * same object, or if they are both the zero handle.
   * NOTE: 'a' can be 0 (int) to repr. the zero handle. 'hand==0' is equiv.
   * to 'hand.isZero()'.
   *
   * @param a S.a.
   */
  bool operator==(const Handle<T>& a) const {
    return (rep==a.rep);
  }

  template<class T2> bool operator==(const Handle<T2>& a) const {
    return (rep==a.rep);
  }

  bool operator==(int a) const {
    if (a!=0) throw InvalidParameterException(EXCEPT_MSG(""));
    return (rep==0);
  }

  bool operator!=(const Handle<T>& a) const {
    return (rep!=a.rep);
  }

  template<class T2> bool operator!=(const Handle<T2>& a) const {
    return (rep!=a.rep);
  }

  bool operator!=(int a) const {
    if (a!=0) throw InvalidParameterException(EXCEPT_MSG(""));
    return (rep!=0);
  }

  /**
   * @return Reference count (0 if this is the 0 handle)
   */
  int getRefCount() const {
    return (rep!=0)?org->pcount:0;
  }

  /**
   * Conversion operator from handle to (T*). Allows for implicit conversion,
   * a handle to be passed as T* argument.
   * ==> Does not work together with DYNCAST, use 'p' instead.
   * Why don't we define a member template converting to T2* using a dynamic
   * cast? This would restrict usage of T to polymorphic types!
   */
  operator T*() const {
    return rep;
  }

  /**
   * Substitute for 'dynamic_cast' (DYNCAST).
   * If T2 can be dyn-casted to T, we return a conversion of 'r'
   * into Handle<T>, otherwise return zero.
   *
   * @param r Handle<T2>
   * @return  Handle<T>
   */
  template<class T2> static Handle<T> cast(const Handle<T2>& r) {
    T* pp=DYNCAST(T,r.p());
    if (pp!=0)
      return Handle<T>(pp,r.getOrgInternal());
    else
      return Handle<T>(); // zero
  }

  /**
   * Returns reference to representation.
   * NOTE: Since references, just as pointers, work together with
   * polymorphism, no problems arise here. For example, if the true
   * type of (*rep) is a subclass of T, the operator will still return
   * a reference to this object, and will NOT make a copy to convert to
   * T.
   *
   * @return Reference to representation
   */
  T& operator*() const {
    if (rep==0)
      throw WrongStatusException(EXCEPT_MSG("This is the zero handle"));
    return *rep;
  }

  /**
   * Returns pointer to representation (or 0 if the handle is the zero
   * handle). Note that for a handle h, h->meth(...) is a shorthand
   * for h.p()->meth(...), the former is clearer and should be preferred.
   * This method should be used if the pointer itself is required.
   * <p>
   * Use 'p' to make clear that the rep. pointer is required, or in
   * situations where the implicit conversion to T* does not work (e.g.
   * in a DYNCAST).
   *
   * @return Pointer to representation
   */
  T* p() const {
    return rep;
  }

  /**
   * Changes representation. The new representation must not have been
   * wrapped by other handles! Never use this method with other repr. than
   * those newly created! To change the handle to become the zero handle,
   * pass 0.
   * <p>
   * NOTE: The best rule to avoid problems is to use 'changeRep' ONLY if
   * 'pp' points to an object just created (pref., 'pp' should be a temp.
   * pointer, s.t. its value cannot be stored anywhere else as in this
   * handle) or 'pp'==0.
   * NOTE: To make this handle ref. to the repr. of another one, use
   * 'operator='.
   * <p>
   * NOTE: If 'oOwn'==false, 'pp' can be owned by other handles, but make
   * sure that this handle is deassoc. from 'pp' before 'pp' is destroyed!
   *
   * @param pp   New representation
   * @param oOwn Do we own the repres., i.e. do we destroy it if the counter
   *             drops to 0? Def.: true
   */
  virtual void changeRep(T* pp,bool oOwn=true) {
    if (pp==rep) return;
    decrRefCount();
    rep=pp;
    if (pp!=0) org=new HandleHelper(oOwn);
    else org=0;
  }

  /**
   * Works only if the ref. counter is 1. Deassociates this handle from its
   * representation (which is not destroyed).
   * NOTE: Use in special situations only.
   */
  virtual void deassoc() {
    if (org==0 || org->pcount!=1)
      throw WrongStatusException(EXCEPT_MSG("Cannot be deassociated"));
    delete org;
    org=0; rep=0; // zero
  }

  /**
   * FOR INTERNAL USE ONLY!
   * Need this, because 'Handle<T>' for different T cannot be made friends.
   *
   * @return Watcher member 'org'
   */
  HandleHelper* getOrgInternal() const {
    return org;
  }

protected:
  // Internal methods

  /**
   * Creates handle with given repres. 'repP' and watcher 'orgP'.
   * The ref. counter in the latter is increased.
   * Used by 'cast' method.
   *
   * @param repP S.a.
   * @param orgP S.a.
   */
  Handle(T* repP,HandleHelper* orgP) : rep(repP),org(orgP) {
    if (orgP!=0) orgP->pcount++;
  }

  /**
   * This method simply destroys the representation behind 'rep', without
   * consulting the ref. counter. It is a helper for 'decrRefCount'.
   * The def. implementation just calls 'delete' of T.
   * ==> Overwrite in subclasses for which this is not appropriate
   */
  virtual void deleteRep() {
    if (rep==0)
      printMsgStdout("INTERNAL ERROR (Handle::deleteRep): Called with zero pointer!");
    else {
      delete rep; rep=0; // def. implementation
    }
  }

  /**
   * If 'org'!=0, we decrease the ref. counter by 1. If it drops to 0, we
   * destroy 'org'. In this case, if 'ourOwn'==true, we also destroy the
   * repres. (calling 'deleteRep').
   */
  void decrRefCount() {
    if (org!=0 && (--(org->pcount))==0) {
      if (org->ourOwn) deleteRep();
      delete org;
      rep=0; org=0;
    }
  }
};

/**
 * Manages a zero handle as static member and returns a reference upon call
 * of 'get'.
 * NOTE: For some reason, the compiler does not tolerate a static
 * 'Handle<T>' pointer within 'Handle', so we have to use this separate class.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class HandleZero
{
protected:
  // Static members

  static Handle<T> zeroH;

public:
  /**
   * @return Reference to zero handle
   */
  static Handle<T>& get() {
    return zeroH;
  }
};

template<class T> Handle<T> HandleZero<T>::zeroH;

#endif
