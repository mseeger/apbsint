//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Standard exceptions
 * ------------------------------------------------------------------- */

#ifndef EX_GLOBAL_H
#define EX_GLOBAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef MATLAB_MEX
  /**
   * Error in direct context of calling a method through the Matlab
   * interface. The message is passed back to Matlab.
   */
  class MatIFException : public StandardException {
  public:
    MatIFException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("MatIFException",mess,file,line) {}
  };
#endif

/**
 * Generic exception for invalid parameter given to a method
 */
class InvalidParameterException : public StandardException {
 public:
  InvalidParameterException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("InvalidParameterException",mess,file,line) {}
};

/**
 * Generic exception: A method has been called which cannot be executed given
 * the current status of the object
 */
class WrongStatusException : public StandardException {
 public:
  WrongStatusException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("WrongStatusException",mess,file,line) {}
};

/**
 * Error within FileUtils method, or general file related
 */
class FileUtilsException : public StandardException {
 public:
  FileUtilsException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("FileUtilsException",mess,file,line) {}
};

/**
 * Attempted to load file that has wrong format
 */
class FileFormatException : public StandardException {
public:
  FileFormatException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("FileFormatException",mess,file,line) {}
};

/**
 * Error when parsing the command file, or general parsing related
 */
class ParseException : public StandardException {
public:
  ParseException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("ParseException",mess,file,line) {}
};

/**
 * Unknown key requested from CommandParser
 */
class KeyNotFoundException : public StandardException {
public:
  KeyNotFoundException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("KeyNotFoundException",mess,file,line) {}

  KeyNotFoundException(const KeyNotFoundException& ex) :
    StandardException("KeyNotFoundException",ex.msg()) {}
};

/**
 * Error in 'ArgBlock' (argument block)
 */
class ArgBlockException : public StandardException {
public:
  ArgBlockException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("ArgBlockException",mess,file,line) {}
};

/**
 * The program is run on an architecture which is not supported
 */
class UnsupportedArchitectureException : public StandardException {
 public:
  UnsupportedArchitectureException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("UnsupportedArchitectureException",mess,file,line) {}
};

/**
 * Error due to instability in a numerical routine
 */
class NumericalException : public StandardException {
public:
  NumericalException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("NumericalException",mess,file,line) {}
};

/**
 * Error in memory manager (see 'FixedMemManager')
 */
class MemManagerException : public StandardException {
public:
  MemManagerException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("MemManagerException",mess,file,line) {}
};

/**
 * Type not supported by templatised object
 */
class TypeNotSuppException : public StandardException {
public:
  TypeNotSuppException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("TypeNotSuppException",mess,file,line) {}
};

/**
 * Internal critical error
 */
class InternalException : public StandardException {
public:
  InternalException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("InternalException",mess,file,line) {}
};

/**
 * Access out of range (use instead of STL 'out_of_range'!)
 */
class OutOfRangeException : public StandardException {
public:
  OutOfRangeException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("OutOfRangeException",mess,file,line) {}
};

/**
 * Method/function not implemented
 */
class NotImplemException : public StandardException {
public:
  NotImplemException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("NotImplemException",mess,file,line) {}
};

/**
 * Memory allocation failed
 */
class MemAllocException : public StandardException {
public:
  MemAllocException(const char* mess=0,const char* file=0,int line=0) :
    StandardException("MemAllocException",mess,file,line) {}
};

class NumRecNotAvailableException : public StandardException {
public:
  NumRecNotAvailableException(const char* mess="SORRY: Cannot execute this method.\nYou came across a method which requires (in the moment) code from the\nNumerical Recipes distribution. The authors do not permit the NR source code to be shipped as part of free software.\nIf you do not like this policy, you might consider contacting the authors and tell them about it.\nBTW: The scientific quality of most of the NR routines is very poor, and the few we still use for convenience\nwill be replaced by proper routines in the very near future.",const char* file=0,int line=0) :
    StandardException("NumRecNotAvailableException",mess,file,line) {}
};

#endif
