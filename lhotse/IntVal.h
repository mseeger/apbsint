//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header abstract class IntVal
 * ------------------------------------------------------------------- */

#ifndef INTVAL_H
#define INTVAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Abstract superclass of 'Interval'. Required by polymorphic code, in
 * which the type T of 'Interval' is not known.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class IntVal
{
public:
  /* Constants
     NOTE: Must lie within 0,...,15! */
  static const int ivOpen  =0;
  static const int ivClosed=1;
  static const int ivInf   =2;
  static const int ivLast  =2;

public:
  // Public methods

  virtual ~IntVal() {}
};

#endif
