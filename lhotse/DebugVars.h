//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class DebugVars
 * ------------------------------------------------------------------- */

#ifndef DEBUGVARS_H
#define DEBUGVARS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

class CommandParser; // predecl.

/**
 * Manages static members (as controlled global variables) related to different
 * debugging purposes. Typically, all debugging options are switched off.
 * In order to setup the debugging options, pass a 'CommandParser' object
 * with the command file containing the debug control file parameters (see
 * LHOTSE documentation) to the 'init' method of this class.
 * NOTE: 'init' should be called from the default 'LHOTSE_init' function iff
 * the constant 'HAVE_DEBUG' is defined (see global.h).
 * NOTE: The old 'MATLAB_DEBUG' constant becomes obsolete, since the Matlab
 * debugging facility is now controlled by variables maintained here.
 * Old 'MATLAB_DEBUG' related code can still be commented out by wrapping into
 * 'MATLAB_DEBUG_OLD'.
 * <p>
 * Debug facilities:
 * - Matlab debug (see 'MatlabDebug')
 *   Here, we maintain the base filename to be passed to the 'activate' method.
 *   'init' also controls whether Matlab debug should be activated upon start
 *   of the program or not (default). 'matlabDebugActivate' and
 *   'matlabDebugDeactivate' are wrappers for the corr. 'MatlabDebug' methods.
 *   Control file pars:
 *   - debug-matlab-base-fname
 *   - debug-matlab-activate
 * - Messing up LHOTSE exceptions
 *   If your debugger cannot catch exceptions and track back from where they
 *   came, activate this facility. It leads to some code within the constructor
 *   of 'StandardException' being executed which triggers a segmentation fault.
 *   This CAN be tracked back by the debugger.
 *   NOTE: Does not work for non-LHOTSE exceptions! CF pars:
 *   - debug-messup-exc-activate:
 *     - 0 (def): Not active
 *     - 1: Active, breaks for all 'StandardException' exceptions, except:
 *          KeyNotFoundException
 *   Exceptions can also be messed up selectively. If 'debug-messup-exc-name'
 *   is also given, it contains the name of the exception which should be
 *   messed up (all others are not).
 * - Print exception message at throw time
 *   The message of an exception is typ. not printed at all, or only at
 *   caught time. If this fac. is active, the mess. is printed at throw time,
 *   when the exception object is created. CF pars:
 *   - debug-print-exc-early:
 *     - 0 (def): Not active
 *     - 1: Active, print mess. at throw time
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class DebugVars
{
protected:
  // Members
  static string matDebBaseFName; // Base filename for Matlab debug
  static bool doMessUpExc;       // Do we mess up exceptions?
  static bool doPrintExc;        // Do we print mess. early?
  static string messUpName;

public:
  // Public static methods

  /**
   * Reads initial values for the members maintained here from the
   * 'CommandParser' object 'contFile'.
   *
   * @param contFile Control file object
   */
  static void init(CommandParser& args);

  /**
   * Activates Matlab debug facility, passing 'matDebBaseFName' as base filename.
   */
  static void matlabDebugActivate();

  /**
   * Deactivates Matlab debug facility.
   */
  static void matlabDebugDeactivate();

  /**
   * If 'name' is not empty, exceptions of that name are messed up. If
   * empty, no exceptions are messed up.
   *
   * @param name S.a.
   */
  static void switchMessUpExceptions(const string& name) {
    if (doMessUpExc=(name.length()>0)) {
      messUpName=name;
    }
  }

  /**
   * Returns true iff messing up exceptions feature active.
   * If 'name' is given and 'messupName' given, ret. true only if they
   * are the same.
   * Returns false if 'name' is any of the following:
   * - 'KeyNotFoundException'
   *
   * @param name S.a.
   * @return     S.a.
   */
  static bool doWeMessUpExceptions(const string& name="") {
    if (!doMessUpExc || name=="KeyNotFoundException")
      return false;
    else if (messUpName.length()==0 || name.length()==0)
      return true;
    else
      return (name==messUpName);
  }

  /**
   * @param messUp Shall LHOTSE exceptions print message at throw time?
   */
  static void switchPrintExcEarly(bool var) {
    doPrintExc=var;
  }

  static bool doWePrintExcEarly() {
    return doPrintExc;
  }
};

#endif
