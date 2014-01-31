//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class LogFile
 * ------------------------------------------------------------------- */

#ifndef LOGFILE_H
#define LOGFILE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header

// Macros (for convenience)

#define ADDLOG(logF,str) (logF).add(str)
#define PRINTLOG(logF,str) (logF).print(__FILE__,__LINE__,str)
#define ADDMATLOG(logF,txt,mat) do { \
  (logF).add(txt); (mat).print((logF).getBuff()); (logF).add(""); \
  } while (0)
#define PRINTMATLOG(logF,txt,mat) do { \
  (logF).add(txt); (mat).print((logF).getBuff()); PRINTLOG(logF,""); \
  } while (0)
#define ADDSCALLOG(logF,txt,val) do { \
  sprintf((logF).tempBuff,"%s %10e",txt,(double) val); \
  (logF).add((logF).tempBuff); \
  } while (0)
#define PRINTSCALLOG(logF,txt,val) do { \
  sprintf((logF).tempBuff,"%s %10e",txt,(double) val); \
  PRINTLOG(logF,(logF).tempBuff); \
  } while (0)
#define ADDINTLOG(logF,txt,val) do { \
  sprintf((logF).tempBuff,"%s %d",txt,(int) val); \
  (logF).add((logF).tempBuff); \
  } while (0)
#define PRINTINTLOG(logF,txt,val) do { \
  sprintf((logF).tempBuff,"%s %d",txt,(int) val); \
  PRINTLOG(logF,(logF).tempBuff); \
  } while (0)

/**
 * Maintains log file. Messages to be printed into a logfile can be composed
 * in multiple steps. The instance maintains a buffer to which strings given
 * to the add method are simply appended. The print method appends a final
 * string and prints the whole buffer together with filename and line number
 * of the PRINTLOG macro into the logfile. The buffer variable can also be
 * accessed directly by the getBuff method.
 * <p>
 * The feature of printing filename and line number can be suppressed. In this
 * case, only the strings itself are written into the file, every line term.
 * by newline. The default is to print name and number.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class LogFile
{
protected:
  // Members

  ofstream os;     // Output file stream
  my_string buff;  // String buffer
  bool nameNum;    // Print filename and line number?

public:

  char tempBuff[257]; // Buffer used by macros

  // Constructors

  /**
   * @param fname Filename for logfile
   */
  LogFile(const char* fname) : os(fname),nameNum(true) {
    if (!os) throw FileUtilsException("Cannot create logfile");
    os.precision(20); os.setf(std::ios::scientific);
  }

  // Public methods

  /**
   * @param flag Shall we print filename and line number together with
   *             every message?
   */
  void nameAndNumber(bool flag) {
    nameNum=flag;
  }

  void add(const my_string& str) {
    buff+=str; buff+="\n";
  }

  void add(const char* str) {
    buff+=str; buff+="\n";
  }

  void print(const char* name,int no,const my_string& str) {
    if (nameNum)
      os << "**" << name << "(" << no << "):\n" << buff << str << endl;
    else
      os << buff << str << endl;
    buff="";
  }

  void print(const char* name,int no,const char* str) {
    if (nameNum)
      os << "**" << name << "(" << no << "):\n" << buff << str << endl;
    else
      os << buff << str << endl;
    buff="";
  }

  my_string& getBuff() const {
    return (my_string&) buff;
  }
};

#endif
