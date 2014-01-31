//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class DefaultLogs
 * ------------------------------------------------------------------- */

#ifndef DEFAULTLOGS_H
#define DEFAULTLOGS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header
#include "lhotse/LogFile.h"

// Macros

#define ADDGLOBALLOG(str)       ADDLOG(DefaultLogs::global(),str)
#define PRINTGLOBALLOG(str)     PRINTLOG(DefaultLogs::global(),str)
#define AMATGLOBALLOG(txt,mat)  ADDMATLOG(DefaultLogs::global(),txt,mat)
#define PMATGLOBALLOG(txt,mat)  PRINTMATLOG(DefaultLogs::global(),txt,mat)
#define ASCALGLOBALLOG(txt,val) ADDSCALLOG(DefaultLogs::global(),txt,val)
#define PSCALGLOBALLOG(txt,val) PRINTSCALLOG(DefaultLogs::global(),txt,val)
#define AINTGLOBALLOG(txt,val)  ADDINTLOG(DefaultLogs::global(),txt,val)
#define PINTGLOBALLOG(txt,val)  PRINTINTLOG(DefaultLogs::global(),txt,val)
#define GLOBALLOGBUFF           ((DefaultLogs::global()).getBuff())

#define ADDERRORLOG(str)        ADDLOG(DefaultLogs::error(),str)
#define PRINTERRORLOG(str)      PRINTLOG(DefaultLogs::error(),str)
#define AMATERRORLOG(txt,mat)   ADDMATLOG(DefaultLogs::error(),txt,mat)
#define PMATERRORLOG(txt,mat)   PRINTMATLOG(DefaultLogs::error(),txt,mat)
#define ASCALERRORLOG(txt,val)  ADDSCALLOG(DefaultLogs::error(),txt,val)
#define PSCALERRORLOG(txt,val)  PRINTSCALLOG(DefaultLogs::error(),txt,val)
#define AINTERRORLOG(txt,val)   ADDINTLOG(DefaultLogs::error(),txt,val)
#define PINTERRORLOG(txt,val)   PRINTINTLOG(DefaultLogs::error(),txt,val)
#define ERRORLOGBUFF            ((DefaultLogs::error()).getBuff())

#ifdef HAVE_DEBUG
#define ADDDEBUGLOG(str)        ADDLOG(DefaultLogs::debug(),str)
#define PRINTDEBUGLOG(str)      PRINTLOG(DefaultLogs::debug(),str)
#define AMATDEBUGLOG(txt,mat)   ADDMATLOG(DefaultLogs::debug(),txt,mat)
#define PMATDEBUGLOG(txt,mat)   PRINTMATLOG(DefaultLogs::debug(),txt,mat)
#define ASCALDEBUGLOG(txt,val)  ADDSCALLOG(DefaultLogs::debug(),txt,val)
#define PSCALDEBUGLOG(txt,val)  PRINTSCALLOG(DefaultLogs::debug(),txt,val)
#define AINTDEBUGLOG(txt,val)   ADDINTLOG(DefaultLogs::debug(),txt,val)
#define PINTDEBUGLOG(txt,val)   PRINTINTLOG(DefaultLogs::debug(),txt,val)
#define DEBUGLOGBUFF            ((DefaultLogs::debug()).getBuff())
#endif

/**
 * This static class maintains the default logfiles of the system. Currently,
 * these are:
 * - global logfile:  messages of all kinds
 * - error logfile:   warnings and errors
 * - debug logfile:   debug outputs (see also MatlabDebug)
 * - monitor logfile: monitors the optimizers by prot. different criteria
 * Use the macros above to access these logs.
 * <p>
 * ATTENTION: Before accessing any of the default logs, call 'init' and supply
 * the base filename of the task!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class DefaultLogs
{
protected:
  // Static members

  static Handle<LogFile> globalLog;
  static Handle<LogFile> errorLog;
#ifdef HAVE_DEBUG
  static Handle<LogFile> debugLog;
#endif

public:
  // Public static methods

  /**
   * This init. method has to be called before any of the logs is accessed!
   *
   * @param baseName Base filename of the task
   */
  static void init(const my_string& baseName) {
    char buff[200];

    if (globalLog.isZero()) {
      sprintf(buff,"%slhotse-global.log",baseName.c_str());
      globalLog.changeRep(new LogFile(buff));
      sprintf(buff,"%slhotse-error.log",baseName.c_str());
      errorLog.changeRep(new LogFile(buff));
#ifdef HAVE_DEBUG
      sprintf(buff,"%slhotse-debug.log",baseName.c_str());
      debugLog.changeRep(new LogFile(buff));
#endif
    }
  }

  // Public methods

  static LogFile& global() {
    return *globalLog;
  }

  static LogFile& error() {
    return *errorLog;
  }

#ifdef HAVE_DEBUG
  static LogFile& debug() {
    return *debugLog;
  }
#endif
};

#endif
