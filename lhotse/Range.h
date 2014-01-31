//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class Range
 * ------------------------------------------------------------------- */

#ifndef RANGE_H
#define RANGE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#ifndef HAVE_NO_BLAS
#include "lhotse/matrix/predecl.h"
#endif

/**
 * Range objects are used to create/re-assign mask objects (mask vectors,
 * mask matrices) or to extract parts of matrices/vectors. They are very slim
 * and usually created temporarily only.
 * <p>
 * A range has size 'sz'. Empty ranges are not allowed, 'sz'>0 (but see below
 * for 'sz'==-1). All entries of a range must be non-neg.
 * There are different kinds of ranges:
 * - flat: from 'start', increment 1.
 * - linear: from 'start', increment 'step'. 'step' can be negative,
 *   but must be != 0.
 * - indexed: using ref. to an index vector (ArrayHandle<int>)
 *   NOTE: The index may contain duplicate entries. Methods using an indexed
 *   range traverse the index from start to end. This may lead to problems if
 *   l-value and r-value in an expression refer to overlapping memory
 *   regions (do not do this!!).
 * A flat or linear range is called literal, and 'end' is
 *   start + (sz-1)*step,
 * where 'sz' is the size of the range.
 * <p>
 * Open range:
 * For a flat (but not for a linear!) range, we can have 'sz'==-1. Then,
 * the range runs until the last element in the context it is applied to
 * (e.g. if applied to create a mask vector from another vector,
 * 'end' is substituted by the position of the last vector element).
 * Such a range is called open.
 * <p>
 * Full range:
 * A range is equal to the full range relative to some size n iff it is flat,
 * 'start'==0 and 'sz'==-1 or 'start'+'sz'==n.
 * A linear (non-flat) range cannot be the full one, even if 'step'==-1.
 * NOTE: The full range with 'sz'==-1 is the only one which can be applied
 * to an empty buffer (see 'checkRange'). The result is empty again ('size'
 * returns 0).
 * <p>
 * Indexed ranges:
 * They can be created from an 'ArrayHandle' or a 'BaseVector<int>'. In both
 * cases, we try to copy a handle only instead of copying the vector content.
 * This fails only if the 'BaseVector<int>' argument is not flat (in which
 * case a flat copy is drawn here).
 * NOTE: Use 'Range' objects only temporarily, create them new when needed.
 * If a indexed range has to persist, create an explicit copy of the
 * underlying index object.
 *
 * A range is a unary predicate which is true iff the argument lies within
 * the range. For an open range, any argument >='start' lies within.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class Range : public std::unary_function<int,bool>
{
public:
  // Constants
  static const int statFlat  =0;
  static const int statLinear=1;
  static const int statIndex =2;

protected:
  // Members
  int status;   // See constants 'statXXX'
  int sz;       // Size of range. -1 => open (flat) range
  int start;    // Valid iff 'status'!='statIndex'
  int step;     // Valid iff 'status'!='statIndex'
  ArrayHandle<int> index; // Used iff 'status'=='statIndex'

public:
  // Constructors

  /**
   * NOTE: The default range is 0,...,'n'-1 if applied to a buffer of length
   * 'n'. If 'pstep' is given and != 1, the range is linear. In this case,
   * 'pend' must not be -1.
   *
   * @param pstart Value for 'start'. Def.: 0
   * @param pend   Value for 'end'. Def.: -1 (open range)
   * @param pstep  Value for 'step'. Def.: 1 (flat range)
   */
  explicit Range(int pstart=0,int pend=-1,int pstep=1);

  /**
   * Constructor for indexed range. The handle is just copied here.
   *
   * @param pindex Index repres. the range
   */
  Range(const ArrayHandle<int>& pindex);

  /**
   * Constructor for indexed range. If 'pindex' has a flat buffer, we obtain
   * an 'ArrayHandle' to this one. Otherwise, a flat copy is drawn. In both
   * cases, the mem. watcher mechanism makes sure that the buffer is not
   * dealloc. prematurely (e.g., it is OK for 'pindex' to be a temp. object).
   * <p>
   * ATTENTION: If no copy of 'pindex' is drawn, changing 'pindex' will change
   * this object!
   *
   * @param pindex Index repres. the range (no copy is drawn!)
   */
#ifndef HAVE_NO_BLAS
  Range(const BaseVector<int>& pindex);
#endif

  /**
   * Copy constructor
   *
   * @param arg Source
   */
  Range(const Range& arg) : status(arg.status),sz(arg.sz),start(arg.start),
  step(arg.step),index(arg.index) {}

  // Public methods

  /**
   * @return status (see constants 'statXXX')
   */
  int getStatus() const {
    return status;
  }

  /**
   * @return Is this a literal range? Otherwise, it's an index
   */
  bool isLiteralRange() const {
    return (status!=statIndex);
  }

  /**
   * @return Is this a flat range (step size 1)?
   */
  bool isFlatRange() const {
    return (status==statFlat);
  }

  /**
   * NOTE: Only a flat range can be the full one.
   *
   * @param len Length of full range
   * @return    Is this the full range?
   */
  bool isFullRange(int len) const {
    return (status==statFlat && start==0 && (sz==-1 || len==sz));
  }

  /**
   * A range is open iff it is flat and 'sz'==-1.
   *
   * @return Is this range open?
   */
  bool isOpen() const {
    return (status==statFlat && sz==-1);
  }

  /**
   * @return Value of 'start'
   */
  int getStart() const {
    if (status==statIndex) throw WrongStatusException(EXCEPT_MSG(""));
    return start;
  }

  /**
   * Returns value of 'end'. If 'sz'==-1, 'n'-1 is returned.
   * NOTE: 'n' is used only if 'sz'==-1. Cannot be used for indexed range.
   *
   * @param n See above
   * @return  "
   */
  int getEnd(int n=0) const {
    if (status==statIndex) throw WrongStatusException(EXCEPT_MSG(""));
    return (sz!=-1)?(start+step*(sz-1)):(n-1);
  }

  /**
   * Cannot be used for indexed range.
   *
   * @return Value of 'step'
   */
  int getStep() const {
    if (status==statIndex) throw WrongStatusException(EXCEPT_MSG(""));
    return step;
  }

  /**
   * @return Ref. to 'index' (or 0 if not indexed range)
   */
  const ArrayHandle<int>& getIndex() const {
    return (status==statIndex)?index:ArrayHandleZero<int>::get();
  }

  /**
   * Returns size of the range. In 'n', the length of the buffer to which
   * the range is to be applied, has to be passed. If 'sz'==-1, we just
   * substitute 'n'-'start'. If 'sz'!=-1, 'n' is not used.
   * NOTE: No checks here for range violations. Use 'checkRange'!
   *
   * @param n See above
   * @return  "
   */
  int size(int n) const {
    return (sz!=-1)?sz:(n-start);
  }

  /**
   * @param pos Position
   * @return    Range index value
   */
  const int operator[](int pos) const {
    if (pos<0) throw OutOfRangeException(EXCEPT_MSG("pos"));
    if (sz==-1) return start+pos;
    else if (pos>=sz) throw OutOfRangeException(EXCEPT_MSG("pos"));
    return (status==statIndex)?index[pos]:(start+step*pos);
  }

  /**
   * Returns largest position of the range, given that it is applied to
   * a buffer with range 0,...,'n'-1.
   * NOTE: 'n' is used only if the range is open.
   *
   * @param n See above
   * @return  "
   */
  int getMaxPos(int n) const;

  /**
   * Checks whether this range violates the buffer range 0,...,'n'-1.
   * NOTE: By definition, the full range does not violate any buffer, even
   * an empty one.
   *
   * @param n See above
   * @return  Does it violate?
   */
  bool checkRange(int n) const;

  /**
   * Checks whether this range is a unique map. This is always true for
   * a literal range, but may be false for an indexed one.
   *
   * @return Is the range a unique map?
   */
  bool isUniqueMap() const;

  /**
   * Extracts this range of index vector 'src', whose size is (at least) 'n',
   * and copies it to 'trg': for a literal range, 'trg[i]=src[start+i*step]'.
   * For an index range, 'trg[i]=src[index[i]]'. If 'doCheck' is true,
   * we first check for violations of the buffer range 0,...,'n'-1. If there
   * is such a violation, nothing is copied and a 'OutOfRangeException'
   * is thrown.
   * NOTE: Use safe version below if possible.
   *
   * @param src     See above
   * @param n       "
   * @param trg     "
   * @param doCheck ". Def.: true
   */
  void mapIndex(const int* src,int n,int* trg,bool doCheck=true) const;

  /**
   * Safe version. 'src' is checked for range violation. If 'trg' is too
   * small, it is realloc. with the correct size.
   *
   * @param src S.a.
   * @param trg S.a.
   */
  void mapIndex(const ArrayHandle<int>& src,ArrayHandle<int>& trg) const;

  /**
   * Extracts subrange 'rng' from this range and returns it as a new object.
   * The new range is flat iff both this and 'rng' are flat, and linear if
   * both this and 'rng' are literal, but at least one is not flat. In any
   * other case, the new range is indexed. The new range is open iff both
   * this and 'rng' are.
   *
   * @param rng See above
   * @return    "
   */
  Range subrange(const Range& rng) const;

  /**
   * Predicate
   *
   * @param arg Value
   * @return    Is 'arg' in the range?
   */
  bool operator()(int arg) const {
    int i;

    switch(status) {
    case statFlat:
      return (arg>=start && (sz==-1 || arg<start+sz));
    case statLinear:
      if ((i=arg-start)<0 || i%step!=0) return false;
      return (i<sz*step);
    default:
      if (arg<0) return false;
      for (i=0; i<sz; i++)
	if (arg==index[i]) return true;
      return false;
    }
  }

  /**
   * @param arg Value
   * @return    Pos. of 'arg' in the range, or -1 (not in range)
   */
  int getPos(int arg) const {
    int i;

    if (arg<0) return -1;
    switch(status) {
    case statFlat:
      return (arg>=start && (sz==-1 || arg<start+sz))?(arg-start):(-1);
    case statLinear:
      return ((i=arg-start)<0 || i%step!=0 || i>=sz*step)?(-1):(i/sz);
    default:
      for (i=0; i<sz; i++)
	if (arg==index[i]) return i;
      return -1;
    }
  }

  /**
   * Returns range which is obtained from this one by adding 'off' to each
   * element. For a flat or linear range, the new one is flat/linear as well,
   * for an indexed one a copy of 'index' is drawn (if 'off'!=0).
   *
   * @param off Offset by which to translate each elem.
   */
  Range translate(int off) const;

  /**
   * Changes range represented by object. Can only be used for flat or
   * linear ranges.
   *
   * @param pstart Value for 'start'. Nonnegative
   * @param pend   Value for 'end'. Def.: -1 (open range)
   * @param pstep  Value for 'step'. Def.: 1 (flat range)
   */
  void reset(int pstart,int pend=-1,int pstep=1);

  // Static methods

  /**
   * Checks whether the index given by 'ind','n' is strictly monotonically
   * increasing.
   * If 'ind' is 'ArrayHandle<int>', 'n' need not be specified.
   *
   * @param ind See above
   * @param n   "
   * @return    Strictly increasing?
   */
  static bool isIncreasing(const int* ind,int n) {
    if (n<0) throw InvalidParameterException(EXCEPT_MSG(""));
    int i,act,prev;

    if (n<2) return true;
    prev=*(ind++);
    if (prev<0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (prev+n-1>ind[n-2]) return false; // quick test
    for (i=1; i<n; i++) {
      act=*(ind++);
      if (act<=prev) return false;
      prev=act;
    }
    return true;
  }

  static bool isIncreasing(const ArrayHandle<int>& ind) {
    return isIncreasing(ind.p(),ind.size());
  }
};

/**
 * Manages a default range (which represents the full range 0,...) as static
 * member and returns a ref. upon 'get'. Useful to pass as def. argument.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class RangeFull
{
protected:
  static Range defR;

public:
  /**
   * @return Reference to def. full range
   */
  static const Range& get() {
    return defR;
  }
};

/**
 * Global function, simply returning 'RangeFull::get()'. Allows to write
 * 'full()' instead of 'RangeFull::get()' to obtain a const Range& to the
 * full range.
 * Use the old form, if 'full' has a different local meaning!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
extern const Range& full();

#endif
