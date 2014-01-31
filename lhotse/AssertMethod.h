//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class AssertMethod
 * ------------------------------------------------------------------- */

#ifndef ASSERTMETHOD_H
#define ASSERTMETHOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

// Macros

/**
 * MYASS(cond)
 * Only if CHECKASSERTS defined.
 * Checks condition 'cond', prints message containing filename 'name' and
 * line number 'line'. Then, if CHECKASSERTS == 1, 'abort' is called to
 * stop the program.
 */
#ifdef CHECKASSERTS
#define MYASS(cond) do {					   \
  if (!(cond)) {						   \
    char* msg = (char*) malloc(sizeof(char)*(strlen(__FILE__)+45));\
    sprintf(msg,"ASSERTION FAILED in '%s' (line %d):",		   \
	    __FILE__,__LINE__);					   \
    printMsgStdout(msg);					   \
    free((void*) msg);						   \
    if ((CHECKASSERTS)==1) abort();				   \
  }								   \
  } while (0)
#else
#define MYASS(cond) do {} while (0)
#endif

/**
 * ASS
 * Should be called by part. classes at entry and leave of every method,
 * and at leave of constructor. Calls 'assertm', passing filename and
 * line number.
 * Empty if CHECKASSERTS not defined.
 */
#ifdef CHECKASSERTS
#define ASS assertm(__FILE__,__LINE__)
#else
#define ASS do {} while (0)
#endif

/**
 * The assertion method 'assertm' can be overwritten in any subclass.
 * The idea is that 'assertm' checks invariances of the object which are
 * pre- and postconditions for every method of the class.
 * 'assertm' can check a subset of plausible conditions only, or none at
 * all.
 * Can also be used for debugging, to track down an inconsistent condition.
 * 'assert' should be called by every method in a participating class, at
 * entry and upon leaving. Use the macro ASS for this.
 * ==> The more methods comply to this, the tighter can a error condition
 *     be bracketed.
 * NOTE: ASS can be called in the constructor as well (end of it).
 * <p>
 * Macros are defined here to facil. coding of 'assertm'. MYASS checks a
 * condition and prints an error message containing line number and file
 * name from where 'assertm' has been called.
 * <p>
 * CHECKASSERTS macro:
 * ASS does not call 'assertm' if the macro CHECKASSERTS is not defined.
 * If it is defined, it contains a status flag:
 * - 0: MYASS just prints an error message
 * - 1: MYASS prints an error message, then triggers a program abortion
 *      by calling 'abort'
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class AssertMethod
{
public:
  /**
   * Assertion method. Check consistency conditions. If one is violated,
   * print descriptive message. Use MYASS to write conditions.
   */
#ifdef CHECKASSERTS
  virtual void assertm(const char* file,int line) const = 0;
#endif
};

#endif
