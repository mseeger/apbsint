//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class AccumulFunc
 * ------------------------------------------------------------------- */

#ifndef ACCUMULFUNC_H
#define ACCUMULFUNC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Abstract accumulator superclass.
 * An abstract accumulator defines 'operator()' and 'get'. 'operator()' is
 * fed with an Arg argument, 'get' returns an Res value. The accumulator
 * has an internal state which is updated with each 'operator()', and a
 * repres. of which is returned by 'get'. Example: sum of elements passed
 * to 'operator()'.
 * <p>
 * NOTE: Inner state member(s) have to be decl. mutable, because 'operator()'
 * and 'reset' are const methods (required s.t. const refs of this class
 * can be used).
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class Arg,class Res> class AccumulFunc
{
public:
  typedef Arg argument_type; // like in STL unary_function
  typedef Res result_type;

  /**
   * Accumulation operator
   */
  virtual void operator()(const Arg& arg) const = 0;

  /**
   * Resets state
   */
  virtual void reset() const = 0;

  /**
   * Return state repres.
   */
  virtual Res get() const = 0;
};

/**
 * Binder (adapter -> AccumulFunc<T,T2>)
 * <p>
 * Configured by a binary function f: (T,T2) -> T2 and a starting value s0.
 * Maintains inner state s (type T2). Every time 'operator()' is invoked with
 * arg. a, we execute
 *   s = f(a,s).
 * 'get' returns current value of s, and 'reset' sets s=s0.
 * Example: Pass 'std::plus' for 'f' to obtain vectorized sum.
 * <p>
 * NOTE: If f is (T2,T2) -> T2, but args are of type T, use 'Cast2ndAdapter'
 * (function 'cast1st') to obtain s = f((T2) a,s).
 * Example: Summation of float entries using double summation variable:
 *   dsum=fvec.accumulate(accum_fun(cast1st<float>(std::plus<double>())));
 * <p>
 * NOTE: More complicated accumulators can be realized with composite func.
 * objects (see 'FuncObjects'). For example, an accumulator with
 * f(a,s) = g(h(a),s) can be realized by passing
 *   compose22(g,h)
 * for 'f'.
 * <p>
 * NOTE: Why is the inner state 's' mutable? The methods must be const,
 * while modifying 's'. Otherwise, temp. created objects of this type cannot
 * be passed to other methods as const refs.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class BinOp> class AccumulBinder :
  public AccumulFunc<typename BinOp::first_argument_type,
		     typename BinOp::result_type>
{
protected:
  typedef typename BinOp::first_argument_type Arg;
  typedef typename BinOp::result_type Res;
  typedef typename BinOp::second_argument_type Arg2;
  BinOp f;
  mutable Res s; // inner state
  Res s0;        // init. state

public:
  AccumulBinder(const BinOp& x,const Res& s0val) : f(x),s0(s0val),s(s0val) {}

  void operator()(const Arg& arg) const {
    s=f(arg,(Arg2) s);
  }

  void reset() const {
    s=s0;
  }

  Res get() const {
    return s;
  }
};

template<class BinOp> inline AccumulBinder<BinOp>
accum_fun(const BinOp& op,const typename BinOp::result_type& s0)
{
  return AccumulBinder<BinOp>(op,s0);
}

/*
 * Default Accumulators
 * NOTE: Some default accumulator can be obtained as bindings:
 * - Sum of elements of vector:
 *   s=vec.accumulate(accum_fun(std::plus<T>(),s0));
 * - Cast T -> T2 and sum (f.ex.: float -> double):
 *   s=vec.accumulate(accum_fun(cast1st<T>(std::plus<T2>()),s0));
 * - Count how many entries of vector fulfill a (unary) predicate
 *   s=vec.accumulate(accum_fun(cast2nd<bool>(std::plus<int>()),0),pred());
 *   Here, cast2nd makes 'plus' a function (int,bool) -> int
 * - Sum of elements of f(v) for vector v:
 *   s=vec.accumulate(accum_fun(compose22(std::plus<T>(),f),s0));
 */

/**
 * Max of elements accumulator.
 * Requires > to be defined for T.
 * We also keep a counter 'cnt' which is reset to 0 upon construction and
 * call of 'reset' and increm. for each 'operator()'. If 'operator()' updates
 * 'maxval', 'pos' is set to the actual 'cnt'.
 * The return type is pair<T,int>, ret. the max. value and the value of
 * 'pos'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class AccumMax : public AccumulFunc<T,pair<T,int> >
{
protected:
  mutable T maxval;
  mutable int cnt,pos;
  mutable bool active;

public:
  /**
   * Constructor.
   */
  AccumMax() : cnt(0),active(false) {}

  void reset() const {
    cnt=0; active=false;
  }

  void operator()(const T& arg) const {
    if (!active || arg>maxval) {
      maxval=arg; pos=cnt; active=true;
    }
    cnt++;
  }

  pair<T,int> get() const {
    if (!active) throw WrongStatusException("Accumulator not initialised");
    return pair<T,int>(maxval,pos);
  }
};

/**
 * Same as 'AccumMax', but returns max. value only.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class AccumMaxVal : public AccumulFunc<T,T>
{
protected:
  mutable T maxval;
  mutable bool active;

public:
  /**
   * Constructor.
   */
  AccumMaxVal() : active(false) {}

  void reset() const {
    active=false;
  }

  void operator()(const T& arg) const {
    if (!active || arg>maxval) {
      maxval=arg; active=true;
    }
  }

  T get() const {
    if (!active) throw WrongStatusException("Accumulator not initialised");
    return maxval;
  }
};

/**
 * Same as 'AccumMax', but returns position of max. elem. only.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class AccumMaxPos : public AccumulFunc<T,int>
{
protected:
  mutable T maxval;
  mutable int cnt,pos;
  mutable bool active;

public:
  /**
   * Constructor.
   */
  AccumMaxPos() : cnt(0),active(false) {}

  void reset() const {
    cnt=0; active=false;
  }

  void operator()(const T& arg) const {
    if (!active || arg>maxval) {
      maxval=arg; pos=cnt; active=true;
    }
    cnt++;
  }

  int get() const {
    if (!active) throw WrongStatusException("Accumulator not initialised");
    return pos;
  }
};

/**
 * Min of elements accumulator.
 * Requires > to be defined for T.
 * We also keep a counter 'cnt' which is reset to 0 upon construction and
 * call of 'reset' and increm. for each 'operator()'. If 'operator()' updates
 * 'minval', 'pos' is set to the actual 'cnt'.
 * The return type is pair<T,int>, ret. the min. value and the value of
 * 'pos'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class AccumMin : public AccumulFunc<T,pair<T,int> >
{
protected:
  mutable T minval;
  mutable int cnt,pos;
  mutable bool active;

public:
  /**
   * Constructor.
   */
  AccumMin() : cnt(0),active(false) {}

  void reset() const {
    cnt=0; active=false;
  }

  void operator()(const T& arg) const {
    if (!active || arg<minval) {
      minval=arg; pos=cnt; active=true;
    }
    cnt++;
  }

  pair<T,int> get() const {
    if (!active) throw WrongStatusException("Accumulator not initialised");
    return pair<T,int>(minval,pos);
  }
};

/**
 * Same as 'AccumMin', but returns min. value only.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class AccumMinVal : public AccumulFunc<T,T>
{
protected:
  mutable T minval;
  mutable bool active;

public:
  /**
   * Constructor.
   */
  AccumMinVal() : active(false) {}

  void reset() const {
    active=false;
  }

  void operator()(const T& arg) const {
    if (!active || arg<minval) {
      minval=arg; active=true;
    }
  }

  T get() const {
    if (!active) throw WrongStatusException("Accumulator not initialised");
    return minval;
  }
};

/**
 * Same as 'AccumMin', but returns position of min. elem. only.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class AccumMinPos : public AccumulFunc<T,int>
{
protected:
  mutable T minval;
  mutable int cnt,pos;
  mutable bool active;

public:
  /**
   * Constructor.
   */
  AccumMinPos() : cnt(0),active(false) {}

  void reset() const {
    cnt=0; active=false;
  }

  void operator()(const T& arg) const {
    if (!active || arg<minval) {
      minval=arg; pos=cnt; active=true;
    }
    cnt++;
  }

  int get() const {
    if (!active) throw WrongStatusException("Accumulator not initialised");
    return pos;
  }
};

#endif
