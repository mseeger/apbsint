//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class NullaryFunc, BinderOnly
 * ------------------------------------------------------------------- */

#ifndef NULLARYFUNC_H
#define NULLARYFUNC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <functional>

/**
 * Similar to STL 'unary_function', but for nullary function (no arguments).
 * Declares operator() (without argument).
 * NOTE: Makes sense if there is an internal state which makes the function
 * non-constant (example: PRN generator, see 'Generator').
 * <p>
 * 'mybindarg' extends STL 'bind1st' to map 'unary_function' ->
 * 'NullaryFunc'. 'ptr_0fun' extends STL 'ptr_fun' to 'NullaryFunc'.
 * 'mybind1st' is special variants of STL 'bind1st', in that a
 * 'binary_function' and a 'NullaryFunc' are mapped to a 'unary_function'.
 * 'mybind2nd' is the same, but the 'NullaryFunc' determines the 2nd
 * argument.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class Res> class NullaryFunc
{
public:
  typedef Res result_type;

  /**
   * Function eval. operator
   */
  virtual result_type operator()() const = 0;
};

/**
 * Same as STL's 'bind1st', but maps 'unary_function' to 'NullaryFunc'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class UnOp> class BinderOnly :
  public NullaryFunc<typename UnOp::result_type>
{
protected:
  UnOp op;
  typename UnOp::argument_type arg;

public:
  BinderOnly(const UnOp& x,const typename UnOp::argument_type& v) :
    op(x),arg(v) {}

  typename UnOp::result_type operator()() const {
    return op(arg);
  }
};

template<class UnOp,class T>
inline BinderOnly<UnOp> mybindarg(const UnOp& op,const T& v)
{
  typedef typename UnOp::argument_type ArgType;
  return BinderOnly<UnOp>(op,ArgType(v));
}

/**
 * Like STL's 'bind1st', but argument is bound to nullary function, not
 * to constant.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp,class NullOp> class MyBinder1st :
  public std::unary_function<typename BinOp::second_argument_type,
			     typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1Type;
  typedef typename BinOp::second_argument_type T2Type;
  typedef typename BinOp::result_type T3Type;
  BinOp f;
  NullOp op0;

public:
  MyBinder1st(const BinOp& x,const NullOp& y) : f(x),op0(y) {}

  T3Type operator()(const T2Type& arg) const {
    return f((T1Type) op0(),arg);
  }
};

template<class BinOp,class NullOp>
inline MyBinder1st<BinOp,NullOp> mybind1st(const BinOp& f,const NullOp& op)
{
  return MyBinder1st<BinOp,NullOp>(f,op);
}

/**
 * Like STL's 'bind2nd', but argument is bound to nullary function, not
 * to constant.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp,class NullOp> class MyBinder2nd :
  public std::unary_function<typename BinOp::first_argument_type,
			     typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1Type;
  typedef typename BinOp::second_argument_type T2Type;
  typedef typename BinOp::result_type T3Type;
  BinOp f;
  NullOp op0;

public:
  MyBinder2nd(const BinOp& x,const NullOp& y) : f(x),op0(y) {}

  T3Type operator()(const T1Type& arg) const {
    return f(arg,(T2Type) op0());
  }
};

template<class BinOp,class NullOp>
inline MyBinder2nd<BinOp,NullOp> mybind2nd(const BinOp& f,const NullOp& op)
{
  return MyBinder2nd<BinOp,NullOp>(f,op);
}

/**
 * Like 'mybindarg', but argument is bound to nullary function, not
 * to constant.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class UnOp,class NullOp> class MyBinderNull :
  public NullaryFunc<typename UnOp::result_type>
{
protected:
  typedef typename UnOp::argument_type ArgType;
  typedef typename UnOp::result_type ResType;
  UnOp f;
  NullOp op0;

public:
  MyBinderNull(const UnOp& x,const NullOp& y) : f(x),op0(y) {}

  ResType operator()() const {
    return f((ArgType) op0());
  }
};

template<class UnOp,class NullOp>
inline MyBinderNull<UnOp,NullOp> mybindnull(const UnOp& f,const NullOp& op)
{
  return MyBinderNull<UnOp,NullOp>(f,op);
}

/**
 * Wraps pointer to nullary function
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class Res> class Ptr2NullFunc : public NullaryFunc<Res>
{
protected:
  typedef Res (*FuncType)();
  FuncType func;

public:
  Ptr2NullFunc(FuncType ptr) : func(ptr) {}

  Res operator()() const {
    return (*func)();
  }
};

// Was 'ptr2nullfunc' before
template<class Res>
inline Ptr2NullFunc<Res> ptr_0fun(Res (*func)())
{
  return Ptr2NullFunc<Res>(func);
}

#endif
