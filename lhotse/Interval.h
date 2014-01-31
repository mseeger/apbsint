//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header classes Interval, DefIVal
 * ------------------------------------------------------------------- */

#ifndef INTERVAL_H
#define INTERVAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/IntVal.h"

/**
 * Represents an interval on the real line. Either of the left/right
 * boundary can be open, closed or infinite. Here, "infinite" means that
 * a boundary check is not done. Intervals are used for range
 * checks (see 'BaseVector::checkBounds').
 * The operators '<' and '==' have to be defined for the type T.
 * <p>
 * Formally, an interval based on left/right boundaries a,b is valid iff:
 * - a or b is infinite, or
 * - a,b not infinite, a < b, or
 * - a,b not infinite, a <= b, a closed or b closed
 * Note that there are valid empty intervals (although they are not useful).
 * The type of the boundaries is coded in 'boundFlag', the 4 lsb's for
 * the lower, the 4 msb's for the upper boundary. For type infinite, the
 * boundary value is ignored.
 * <p>
 * NOTE: The assoc. class 'DefIVal' contains static members for a number
 * of default intervals. It works only for the standard numerical types!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class Interval :
  public IntVal,public std::unary_function<T,bool>
{
protected:
  // Members
  uchar boundFlag; // boundary types (see header and 'ivXXX')
  T lowBnd,uppBnd; // boundaries

public:
  // Constructors

  /**
   * Default constructor. The boundary value is ignored if the corr. type
   * is 'ivInf'.
   *
   * @param lowB    Lower boundary value.
   * @param uppB    Upper boundary value
   * @param lowType Lower b. type. Def.: infinite
   * @param uppType Upper b. type. Def.: infinite
   */
  Interval(const T& lowB,const T& uppB,int lowType=ivInf,int uppType=ivInf) :
  lowBnd(lowB),uppBnd(uppB) {
    if (lowType<0 || lowType>ivLast)
      throw InvalidParameterException("lowType");
    if (uppType<0 || uppType>ivLast)
      throw InvalidParameterException("uppType");
    if (lowType!=ivInf && uppType!=ivInf && !(lowB<uppB)) {
      if ((lowType!=ivClosed && uppType!=ivClosed) || !(lowB==uppB))
	throw InvalidParameterException("Invalid interval range");
    }
    boundFlag=((uchar) uppType) << 4;
    boundFlag|=((uchar) lowType);
  }

  Interval(const Interval& arg) : boundFlag(arg.boundFlag),lowBnd(arg.lowBnd),
      uppBnd(arg.uppBnd) {}

  /**
   * Checks whether 'val' falls into the given interval. If so, 0 is ret.
   * Otherwise, 1 is ret. if 'val' is too small, 2 if 'val' is too large.
   *
   * @param val Value to check
   * @return    See above
   */
  int check(const T& val) const {
    int lowType=((int) boundFlag)&0x0F,uppType=((int) boundFlag)>>4;
    if (lowType!=ivInf && !(lowBnd<val)) {
      if (lowType==ivOpen || !(lowBnd==val)) return 1;
    }
    if (uppType!=ivInf && !(val<uppBnd)) {
      if (uppType==ivOpen || !(uppBnd==val)) return 2;
    }
    return 0;
  }

  /**
   * Makes this class a unary predicate
   *
   * @param elem Element
   * @return     Is element contained in the interval?
   */
  bool operator()(const T& elem) const {
    return (check(elem)==0);
  }

  /**
   * Checks whether all elements of the vector starting at 'vec', size 'sz',
   * falls into interval. If so, 0 is ret. Otherwise, let a be the value of
   * the first violating comp. 1 is ret. if a is too small, 2 if a is too
   * large. The pos. of this element can be ret. in 'pos'
   *
   * @param vec See above
   * @param sz  "
   * @param pos ". Optional
   * @return    "
   */
  int check(const T* vec,int sz,int* pos=0) const {
    int i,ret;
    for (i=0; i<sz; i++,vec++)
      if ((ret=check(*vec))!=0) {
	if (pos!=0) *pos=i;
	return ret;
      }
    return 0;
  }

  /**
   * Same as above, but 'vec' is 'ArrayHandle'.
   *
   * @param vec See above
   * @param pos ". Optional
   * @return    "
   */
  int check(const ArrayHandle<T>& vec,int* pos=0) const {
    return check(vec,vec.size(),pos);
  }
};

/**
 * Contains static members for freq. used standard intervals: positive,
 * non-negative, negative, non-positive.
 * <p>
 * ATTENTION! The implementation works only for the standard numerical types.
 * The corr. intervals are kept as static members, repr. via handles.
 * The static member 'zeroVal' represents the zero value for the type.
 * The intervals are defined relative to 'zeroVal'. 'initInt' does the
 * initialisation (if not already done). All methods call 'init' before they
 * do anything.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class DefIVal
{
protected:
  // Static members
  static bool isInit;
  static T zeroVal;
  static Handle<Interval<T> > ivPos;
  static Handle<Interval<T> > ivNeg;
  static Handle<Interval<T> > ivNonpos;
  static Handle<Interval<T> > ivNonneg;

public:
  // Static methods

  /**
   * For all supported types: Sets 'zeroVal' and calls 'initInt'.
   * Default: throw exception
   */
  static void init();

  /**
   * @return Interval of all positive numbers
   */
  static const Interval<T>& posit() {
    init();
    return *ivPos;
  }

  /**
   * @return Interval of all negative numbers
   */
  static const Interval<T>& negat() {
    init();
    return *ivNeg;
  }

  /**
   * @return Interval of all nonnegative numbers
   */
  static const Interval<T>& nonneg() {
    init();
    return *ivNonneg;
  }

  /**
   * @return Interval of all nonpositive numbers
   */
  static const Interval<T>& nonpos() {
    init();
    return *ivNonpos;
  }
};

template<class T> Handle<Interval<T> > DefIVal<T>::ivPos;
template<class T> Handle<Interval<T> > DefIVal<T>::ivNeg;
template<class T> Handle<Interval<T> > DefIVal<T>::ivNonneg;
template<class T> Handle<Interval<T> > DefIVal<T>::ivNonpos;

template<class T> void DefIVal<T>::init()
{
  if (!isInit) {
    ivPos.changeRep(new Interval<T>(zeroVal,zeroVal,
				    IntVal::ivOpen,
				    IntVal::ivInf));
    ivNeg.changeRep(new Interval<T>(zeroVal,zeroVal,
				    IntVal::ivInf,
				    IntVal::ivOpen));
    ivNonneg.changeRep(new Interval<T>(zeroVal,zeroVal,
				       IntVal::ivClosed,
				       IntVal::ivInf));
    ivNonpos.changeRep(new Interval<T>(zeroVal,zeroVal,
				       IntVal::ivInf,
				       IntVal::ivClosed));
    isInit=true;
  }
}

#endif
