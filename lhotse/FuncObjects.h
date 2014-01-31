//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Global function objects
 * ------------------------------------------------------------------- */

#ifndef FUNCOBJECTS_H
#define FUNCOBJECTS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <functional>

/*
 * Collects definition of global function objects extending the STL
 * range of available function objects, predicates, and adapters.
 *
 * Some facts about function objects and adapters (STL):
 * - f(x,y) binary: 'bind1st(f,x_0)' -> f(x_0,.) as unary
 * - f(x,y) binary: 'bind2nd(f,y_0)' -> f(.,y_0) as unary
 * - 'ptr_fun(p)': unary_function or binary_function from pointer to func. 'p'
 * - negators 'not1' (unary_function), 'not2' (binary_function) bind
 *   predicates
 *
 * Our extensions here:
 * - Casting: 'cast1st', 'cast2nd' casts 1st/2nd argument
 *   from T. 'cast1', 'cast2' cast result argument from unary/binary func.
 *   to T (binders)
 *   Example: 'BaseVector::count'
 *   NOTE: Many function object adapters are written in a way which do
 *   required casts automatically.
 * - Plug in same argument: 'equal_args' converts bin. func. f into unary
 *   func. x -> f(x,x).
 *   Example: x -> x^2 is given by 'equal_args(std::multiplies<T>())'
 * - Extending normal binary func. to binary func. on pair<...>:
 *   'pair1st', 'pair2nd'
 *   Example: 'ArrayUtils::sortInd'
 * - Function composition: 'compose11', 'compose12', 'compose21',
 *   'compose22'
 *   See comments on function composition below.
 * - Reverse arguments of bin. func.: 'revargs'
 *   Realizes (x,y) -> f(y,x) from f.
 *   Useful together with 'compose2<x>' variant where g2 is not given. See
 *   comments on function composition below.
 * - STL map as function: 'map_fun'
 * - Printing to stdout as unary operator: 'print_op'
 * - Binary min/max: 'BinFuncMax', 'BinFuncMin'
 *   NOTE: 'BinFuncMax' is equivalent to 'ptr_fun(std::max<T>())',
 *   'BinFuncMin' is equivalent to 'ptr_fun(std::min<T>())'
 *
 * ATTENTION: Old code used 'UnaryFunc', 'BinaryFunc' together with
 * binders 'mybind1st', 'mybind2nd', 'ptr2unafunc', 'ptr2binfunc'.
 * ==> Has to be replaced (does not compile anymore)!
 * In old code, 'CumulBinder' and 'cumul_fun' were defined. Use
 * 'AccumulBinder' and 'accum_fun' instead.
 */

/*
 * On the composition of elementary function objects
 *
 * The idea of function objects is to avoid error-prone loops over
 * containers, using the STL algorithms instead. These are configured by
 * unary function objects f : x -> f(x), binary function objects
 * g : (x,y) -> g(x,y).
 * The major LHOTSE classes do support function objects, for example
 * 'apply1', 'apply2' in 'BaseVector', 'BaseMatrix'. See also the concept
 * of accumulators in 'AccumulFunc.h'. If you come from MATLAB, the concept
 * of applying functions to vectors/matrices is familiar to you.
 *
 * The STL provides elementary arithmetic function objects and predicates.
 * An adapter is used to create a function object from arguments, which can
 * be function objects as well. The STL defines some simple adapters, such
 * as 'bind2nd': 'bind2nd(std::plus<double>(),1.0)' implements the unary
 * function x -> x+1. A pointer to a global function can be wrapped with
 * 'ptr_fun': 'ptr_fun(exp)' implements x -> exp(x).
 *
 * We extend the concept here, by allowing for composition of complex
 * function objects from these elementary ones. We define the following
 * adapters:
 * - compose11(f,g):     x     -> f(g(x))
 * - compose12(f,g):     (x,y) -> f(g(x,y))
 * - compose21(f,g1,g2): x     -> f(g1(x),g2(x))
 * - compose21(f,g1):    x     -> f(g1(x),x) [g2 is identity]
 * - compose22(f,g1,g2): (x,y) -> f(g1(x),g2(y))
 * - compose22(f,g1):    (x,y) -> f(g1(x),y) [g2 is identity]
 * Naming: In 'compose<g><x>', <g> is the number of inner g functions,
 * <x> is the number of args of the resulting func. object.
 * In terms of types of the argument function objects, the types for the
 * resulting func. object are deduced from these. For 'compose21', the
 * argument type of g1 det. the res. argument type. All internal
 * result-to-argument conversions are done using static casts.
 *
 * In 'compose2<x>', if g2 is not given, it is taken to be the identity.
 * The adapter 'revargs' can be used if in this case, g1 should be the
 * identity:
 *   compose21(revargs(f),g2)
 * Or:
 *   revargs(compose22(revargs(f),g2))
 *
 * Examples:
 * x -> 3*x*exp(x) is implemented as:
 *   compose11(bind2nd(std::multiplies<double>(),3.0),
 *             compose21(std::multiplies<double>(),ptr_fun(exp)))
 * x -> x / exp(x) is implemented as:
 *   compose21(revargs(std::divides<double>()),ptr_fun(exp))
 * (x,y) -> x / exp(y) is implemented as:
 *   revargs(compose22(revargs(std::divides<double>()),ptr_fun(exp)))
 *
 * If p1, p2 are unary predicates, the unary predicate p1 && p2 is:
 *   compose21(std::logical_and(),p1,p2)
 * For p1, p2 being binary predicates, use:
 *   compose22(std::logical_and(),p1,p2)
 * Note that the STL negaters std::not1, std::not2 can be obtained as:
 *   std::not1(p) <==> compose11(std::logical_not(),p)
 *   std::not2(p) <==> compose12(std::logical_not(),p)
 * (although the direct use of the STL ones is probably more efficient).
 *
 * Possible drawbacks of using composite function objects:
 * If a CFO is used in a performance-critical computation, it should always
 * be checked whether the compiler actually does a good job in optimizing
 * the code, by comparing it against looping over the container. If the
 * computation is not critical, we would recommend CFOs in order to avoid
 * bugs with loops.
 * Another problem is that compiler errors based on incorrect use of CFOs
 * can be hard to parse. This may be a problem with programmers not used
 * to STL concepts.
 */

/**
 * Cast 1st argument (adapter -> binary function). T is source type.
 * <p>
 * Configured by a binary function f: (T1,T2) -> T3 and type T, represents
 * binary function (T,T2) -> T3 defined as
 *   (x,y) -> f((T1) x,y).
 * Required because most STL binary function objects come with T1 == T2
 * and do not include casts.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T,class BinOp> class Cast1stAdapter :
  public std::binary_function<T,typename BinOp::second_argument_type,
						typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  typedef typename BinOp::result_type T3;
  BinOp f;

public:
  Cast1stAdapter(const BinOp& x) : f(x) {}

  T3 operator()(const T& arg1,const T2& arg2) const {
    return f((T1) arg1,arg2);
  }
};

// NOTE: Template arg. T has to be given explicitly
template<class T,class BinOp> inline Cast1stAdapter<T,BinOp>
cast1st(const BinOp& op)
{
  return Cast1stAdapter<T,BinOp>(op);
}

/**
 * Cast 2nd argument (adapter -> binary function). T is source type.
 * <p>
 * Configured by binary function f: (T1,T2) -> T3 and type T, represents
 * binary function (T1,T) -> T3, def. as:
 *   (x,y) -> f(x,(T2) y).
 * Required because most STL binary function objects come with T1 == T2
 * and do not include casts.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T,class BinOp> class Cast2ndAdapter :
  public std::binary_function<typename BinOp::first_argument_type,T,
			      typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  typedef typename BinOp::result_type T3;
  BinOp f;

public:
  Cast2ndAdapter(const BinOp& x) : f(x) {}

  T3 operator()(const T1& arg1,const T& arg2) const {
    return f(arg1,(T2) arg2);
  }
};

// NOTE: Template arg. T has to be given explicitly
template<class T,class BinOp> inline Cast2ndAdapter<T,BinOp>
cast2nd(const BinOp& op)
{
  return Cast2ndAdapter<T,BinOp>(op);
}

/**
 * Cast result argument (adapter -> unary function). T is result type.
 * <p>
 * Configured by unary function f: T1 -> T2 and type T, represents
 * unary function T1 -> T:
 *   x -> (T) f(x)
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T,class UnOp> class CastResUnAdapter :
  public std::unary_function<typename UnOp::argument_type,T>
{
protected:
  typedef typename UnOp::argument_type T1;
  UnOp f;

public:
  CastResUnAdapter(const UnOp& x) : f(x) {}

  T operator()(const T1& arg) const {
    return (T) f(arg);
  }
};

// NOTE: Template arg. T has to be given explicitly
template<class T,class UnOp> inline CastResUnAdapter<T,UnOp>
cast1(const UnOp& op)
{
  return CastResUnAdapter<T,UnOp>(op);
}

/**
 * Cast result argument (adapter -> binary function). T is result type.
 * <p>
 * Configured by binary function f: (T1,T2) -> T3 and type T, represents:
 *   (x,y) -> (T) f(x,y).
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T,class BinOp> class CastResBinAdapter :
  public std::binary_function<typename BinOp::first_argument_type,
			      typename BinOp::second_argument_type,T>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  BinOp f;

public:
  CastResBinAdapter(const BinOp& x) : f(x) {}

  T operator()(const T1& arg1,const T2& arg2) const {
    return (T) f(arg1,arg2);
  }
};

// NOTE: Template arg. T has to be given explicitly
template<class T,class BinOp> inline CastResBinAdapter<T,BinOp>
cast2(const BinOp& op)
{
  return CastResBinAdapter<T,BinOp>(op);
}

/**
 * Plug in same argument (adapter -> unary function)
 * <p>
 * Configured by binary function f: (T1,T2) -> T3, where cast T1 -> T2 must
 * be possible. Represents unary function T1 -> T3: x -> f(x,x), or better:
 * f(x,(T2) x).
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp> class EqualArgAdapter :
  public std::unary_function<typename BinOp::first_argument_type,
			     typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  typedef typename BinOp::result_type T3;
  BinOp f;

public:
  EqualArgAdapter(const BinOp& x) : f(x) {}

  T3 operator()(const T1& arg) const {
    return f(arg,(T2) arg);
  }
};

template<class BinOp> inline EqualArgAdapter<BinOp>
equal_args(const BinOp& op)
{
  return EqualArgAdapter<BinOp>(op);
}

/**
 * Extending binary function to pair<...> arguments, by applying it to
 * first parts (adapter -> binary function)
 * <p>
 * Configured by a binary function f: (T1,T2) -> T3 and types T4, T5.
 * Represents binary function (pair<T1,T4>,pair<T2,T5>) -> T3, which
 * realizes f on the first parts of the arguments.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T4,class T5,class BinOp> class Pair1stAdapter :
  public std::binary_function<pair<typename BinOp::first_argument_type,T4>,
  pair<typename BinOp::second_argument_type,T5>,typename BinOp::result_type>
{
protected:
  typedef pair<typename BinOp::first_argument_type,T4> TP1;
  typedef pair<typename BinOp::second_argument_type,T5> TP2;
  typedef typename BinOp::result_type T3;
  BinOp f;

public:
  Pair1stAdapter(const BinOp& x) : f(x) {}

  T3 operator()(const TP1& arg1,const TP2& arg2) const {
    return f(arg1.first,arg2.first);
  }
};

// NOTE: Template args. T4, T5 to be given explicitly
template<class T4,class T5,class BinOp> inline Pair1stAdapter<T4,T5,BinOp>
pair1st(const BinOp& op)
{
  return Pair1stAdapter<T4,T5,BinOp>(op);
}

template<class T,class BinOp> inline Pair1stAdapter<T,T,BinOp>
pair1st(const BinOp& op)
{
  return Pair1stAdapter<T,T,BinOp>(op);
}

/**
 * Extending binary function to pair<...> arguments, by applying it to
 * second parts (adapter -> binary function)
 * <p>
 * Configured by a binary function f: (T1,T2) -> T3 and types T4, T5.
 * Represents binary function (pair<T4,T1>,pair<T5,T2>) -> T3, which
 * realizes f on the second parts of the arguments.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T4,class T5,class BinOp> class Pair2ndAdapter :
  public std::binary_function<pair<T4,typename BinOp::first_argument_type>,
  pair<T5,typename BinOp::second_argument_type>,typename BinOp::result_type>
{
protected:
  typedef pair<T4,typename BinOp::first_argument_type> TP1;
  typedef pair<T5,typename BinOp::second_argument_type> TP2;
  typedef typename BinOp::result_type T3;
  BinOp f;

public:
  Pair2ndAdapter(const BinOp& x) : f(x) {}

  T3 operator()(const TP1& arg1,const TP2& arg2) const {
    return f(arg1.second,arg2.second);
  }
};

// NOTE: Template args. T4, T5 to be given explicitly
template<class T4,class T5,class BinOp> inline Pair2ndAdapter<T4,T5,BinOp>
pair2nd(const BinOp& op)
{
  return Pair2ndAdapter<T4,T5,BinOp>(op);
}

template<class T,class BinOp> inline Pair2ndAdapter<T,T,BinOp>
pair2nd(const BinOp& op)
{
  return Pair2ndAdapter<T,T,BinOp>(op);
}

/**
 * Function composition for unary function (adapter -> unary function).
 * Used for 'compose11'.
 * <p>
 * Given f: T2 -> T3, g: T1 -> T2, the composite function x -> f(g(x))
 * [ T1 -> T3 ] is represented.
 * NOTE: The result type of g must be castable to T2.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class UnOpF,class UnOpG> class Compose11Adapter :
  public std::unary_function<typename UnOpG::argument_type,
			     typename UnOpF::result_type>
{
protected:
  typedef typename UnOpG::argument_type T1;
  typedef typename UnOpF::argument_type T2;
  typedef typename UnOpF::result_type T3;
  UnOpF f;
  UnOpG g;

public:
  Compose11Adapter(const UnOpF& of,const UnOpG& og) : f(of),g(og) {}

  T3 operator()(const T1& arg) const {
    return f((T2) g(arg));
  }
};

template<class UnOpF,class UnOpG> inline
Compose11Adapter<UnOpF,UnOpG> compose11(const UnOpF& of,const UnOpG& og)
{
  return Compose11Adapter<UnOpF,UnOpG>(of,og);
}

/**
 * Function composition for binary function (adapter -> binary function).
 * Used for 'compose12'.
 * <p>
 * Given f: T3 -> T4, g: (T1,T2) -> T3, the composite function
 * (x,y) -> f(g(x,y)) [ (T1,T2) -> T4 ] is represented.
 * NOTE: The result type of g must be castable to T3.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class UnOpF,class BinOpG> class Compose12Adapter :
  public std::binary_function<typename BinOpG::first_argument_type,
			      typename BinOpG::second_argument_type,
			      typename UnOpF::result_type>
{
protected:
  typedef typename BinOpG::first_argument_type T1;
  typedef typename BinOpG::second_argument_type T2;
  typedef typename UnOpF::argument_type T3;
  typedef typename UnOpF::result_type T4;
  UnOpF f;
  BinOpG g;

public:
  Compose12Adapter(const UnOpF& of,const BinOpG& og) : f(of),g(og) {}

  T4 operator()(const T1& arg1,const T2& arg2) const {
    return f((T3) g(arg1,arg2));
  }
};

template<class UnOpF,class BinOpG> inline
Compose12Adapter<UnOpF,BinOpG> compose12(const UnOpF& of,const BinOpG& og)
{
  return Compose12Adapter<UnOpF,BinOpG>(of,og);
}

/**
 * Function composition for unary function (adapter -> unary function).
 * Used for 'compose21'.
 * <p>
 * Given f: (T1,T2) -> T3, g1: T -> T1, g2: T -> T2, the composite function
 * x -> f(g1(x),g2(x)) [ T -> T3 ] is represented.
 * NOTE: Type T is taken to be the argument type of g1. The implementation
 * still works if g2: T4 -> T2, so that the cast T -> T4 is OK.
 * Also, the result type of g1 (g2) just has to be castable to T1 (T2).
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp,class UnOp1,class UnOp2> class Compose21Adapter :
  public std::unary_function<typename UnOp1::argument_type,
			     typename BinOp::result_type>
{
protected:
  typedef typename UnOp1::argument_type T;
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  typedef typename BinOp::result_type T3;
  typedef typename UnOp2::argument_type T4;
  BinOp f;
  UnOp1 g1;
  UnOp2 g2;

public:
  Compose21Adapter(const BinOp& of,const UnOp1& og1,const UnOp2& og2) :
    f(of),g1(og1),g2(og2) {}

  T3 operator()(const T& arg) const {
    return f((T1) g1(arg),(T2) g2((T4) arg));
  }
};

/**
 * Same as 'Compose21Adapter', but g2 is the identity.
 * Used for 'compose21' (with 2 arguments).
 * <p>
 * NOTE: Here, T must be castable to T2.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp,class UnOp> class Compose21DefAdapter :
  public std::unary_function<typename UnOp::argument_type,
			     typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  typedef typename BinOp::result_type T3;
  typedef typename UnOp::argument_type T;
  BinOp f;
  UnOp g1;

public:
  Compose21DefAdapter(const BinOp& of,const UnOp& og1) : f(of),g1(og1) {}

  T3 operator()(const T& arg) const {
    return f((T1) g1(arg),(T2) arg);
  }
};

template<class BinOp,class UnOp1,class UnOp2> inline
Compose21Adapter<BinOp,UnOp1,UnOp2> compose21(const BinOp& of,const UnOp1& og1,
					      const UnOp2& og2)
{
  return Compose21Adapter<BinOp,UnOp1,UnOp2>(of,og1,og2);
}

template<class BinOp,class UnOp> inline
Compose21DefAdapter<BinOp,UnOp> compose21(const BinOp& of,const UnOp& og1)
{
  return Compose21DefAdapter<BinOp,UnOp>(of,og1);
}

/**
 * Function composition for binary function (adapter -> binary function)
 * Used for 'compose22'.
 * <p>
 * Given f: (T1,T2) -> T3, g1: T4 -> T1, g2: T5 -> T2, the composite function
 * (x,y) -> f(g1(x),g2(y)) [ (T4,T5) -> T3 ] is represented.
 * NOTE: The implementation still works if the result type of g1 (g2) is
 * castable to T1 (T2).
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp,class UnOp1,class UnOp2> class Compose22Adapter :
  public std::binary_function<typename UnOp1::argument_type,
			      typename UnOp2::argument_type,
			      typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::second_argument_type T2;
  typedef typename BinOp::result_type T3;
  typedef typename UnOp1::argument_type T4;
  typedef typename UnOp2::argument_type T5;
  BinOp f;
  UnOp1 g1;
  UnOp2 g2;

public:
  Compose22Adapter(const BinOp& of,const UnOp1& og1,const UnOp2& og2) :
    f(of),g1(og1),g2(og2) {}

  T3 operator()(const T4& arg1,const T5& arg2) const {
    return f((T1) g1(arg1),(T2) g2(arg2));
  }
};

/**
 * Same as 'Compose22Adapter', but g2 is the identity.
 * Used for 'compose22' (with 2 arguments).
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp,class UnOp> class Compose22DefAdapter :
  public std::binary_function<typename UnOp::argument_type,
			      typename BinOp::second_argument_type,
			      typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type T1;
  typedef typename BinOp::result_type T3;
  typedef typename UnOp::argument_type T4;
  typedef typename BinOp::second_argument_type T5;
  BinOp f;
  UnOp g1;

public:
  Compose22DefAdapter(const BinOp& of,const UnOp& og1) : f(of),g1(og1) {}

  T3 operator()(const T4& arg1,const T5& arg2) const {
    return f((T1) g1(arg1),arg2);
  }
};

template<class BinOp,class UnOp1,class UnOp2> inline
Compose22Adapter<BinOp,UnOp1,UnOp2> compose22(const BinOp& of,const UnOp1& og1,
					      const UnOp2& og2)
{
  return Compose22Adapter<BinOp,UnOp1,UnOp2>(of,og1,og2);
}

template<class BinOp,class UnOp> inline
Compose22DefAdapter<BinOp,UnOp> compose22(const BinOp& of,const UnOp& og1)
{
  return Compose22DefAdapter<BinOp,UnOp>(of,og1);
}

/**
 * Reverses order of arguments to binary function (adapter -> binary
 * function)
 * <p>
 * Given bin. func. f, (x,y) -> f(y,x) is represented.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp> class RevArgsAdapter :
  public std::binary_function<typename BinOp::second_argument_type,
			      typename BinOp::first_argument_type,
			      typename BinOp::result_type>
{
protected:
  typedef typename BinOp::second_argument_type T1;
  typedef typename BinOp::first_argument_type T2;
  typedef typename BinOp::result_type T3;
  BinOp f;

public:
  RevArgsAdapter(const BinOp& of) : f(of) {}

  T3 operator()(const T1& x,const T2& y) const {
    return f(y,x);
  }
};

template<class BinOp> inline
RevArgsAdapter<BinOp> revargs(const BinOp& of)
{
  return RevArgsAdapter<BinOp>(of);
}

/**
 * Map (unary function)
 * <p>
 * Configured by STL map T2 -> T and a dummy element. y = f(x): if x is
 * a key in the map, y is the corr. data element. Otherwise, y is the
 * dummy element.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T2,class T> class UnFuncMap : public std::unary_function<T2,T>
{
protected:
  typedef MAP_TYPE(T2,T) MapType;
  typedef typename MAP_CONSTITER(T2,T) MapIter;
  const MapType* mp;
  T dummy;

public:
  UnFuncMap(const MapType* amp,const T& adummy) :
    mp(amp),dummy(adummy){}

  T operator()(const T2& arg) const {
    MapIter it=mp->find(arg);
    return (it!=mp->end())?(it->second):dummy;
  }
};

template<class T2,class T> inline UnFuncMap<T2,T>
map_fun(const MAP_TYPE(T2,T)* amp,const T& adummy)
{
  return UnFuncMap<T2,T>(amp,adummy);
}

/**
 * Print element (unary operator)
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class UnOperPrint : public std::unary_function<T,void>
{
protected:
  my_string separ;
  ostream& os;

public:
  UnOperPrint(ostream& ostr,const my_string& sep) : os(ostr),separ(sep) {}

  void operator()(const T& arg) {
    os << arg << separ;
  }
};

template<class T> inline UnOperPrint<T>
print_op(ostream& ostr,const my_string& sep=" ")
{
  return UnOperPrint<T>(ostr,sep);
}

/**
 * Maximum (binary function)
 * NOTE: Equivalent to 'ptr_fun(std::max<T>())'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class BinFuncMax : public std::binary_function<T,T,T>
{
public:
  T operator()(const T& arg,const T& arg2) const {
    return (arg<arg2)?arg2:arg;
  }
};

/**
 * Minimum (binary function)
 * NOTE: Equivalent to 'ptr_fun(std::max<T>())'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class BinFuncMin : public std::binary_function<T,T,T>
{
public:
  T operator()(const T& arg,const T& arg2) const {
    return (arg<arg2)?arg:arg2;
  }
};

#endif
