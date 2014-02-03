//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: optimize
 * Desc.:  Standard exceptions
 * ------------------------------------------------------------------- */

#ifndef EX_OPTIMIZE_H
#define EX_OPTIMIZE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

//BEGINNS(optimize)
  /**
   * Tried to call lastXXX method when YYYNoXXX has been used before (see
   * CritFunc)
   */
  class NoLastAvailException : public StandardException {
  public:
    NoLastAvailException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("NoLastAvailException",mess,file,line) {}
  };

  /**
   * Used to trigger restart of optimizer by the CritFunc. Has to be passed
   * by line search methods and caught by optimizers.
   */
  class RestartException : public StandardException {
  public:
    RestartException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("RestartException",mess,file,line) {}
  };

  /**
   * Used to trigger termination of optimizer by the CritFunc. Has to be passed
   * by line search methods and caught by optimizers. An optimizer should
   * terminate, however signal back to the user that it has been terminated in
   * this way.
   */
  class TerminateException : public StandardException {
  public:
    TerminateException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("TerminateException",mess,file,line) {}
  };
//ENDNS

#endif
