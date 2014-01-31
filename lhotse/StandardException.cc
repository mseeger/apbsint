//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Definition class StandardException
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#include "lhotse/StandardException.h"
#if defined(HAVE_DEBUG) && !defined(HAVE_NO_BLAS)
#include "lhotse/DebugVars.h"
#endif

// Static members

string StandardException::name;

// Constructors

/*
 * Supporting the messing up LHOTSE exceptions debug facility (see
 * 'DebugVars'. If this is activated, we execute some code in the constructor
 * leading to a segmentation fault.
 */
StandardException::StandardException(const char* nam,const char* mess,
				     const char* file,int line) :
  message() {
  name=nam;
  if (mess==0 || strlen(mess)==0) {
    message=name+": unspecified";
  } else
    message=mess;
  if (file!=0) {
    message+="\nFile: "; message+=file;
    char buff[40];
    sprintf(buff," (line %d)",line);
    message+=buff;
  }
#if defined(HAVE_DEBUG) && !defined(HAVE_NO_BLAS)
  if (DebugVars::doWePrintExcEarly()) {
    cout << "DEBUG: Exception created and thrown. Message:" << endl
	 << message << endl;
  }
  if (DebugVars::doWeMessUpExceptions(name)) {
    // Causes a segmentation fault
    int* breakIt=0;
    *breakIt=1234;
  }
#endif
}
