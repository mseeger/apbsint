/* -------------------------------------------------------------------
 * EPTWRAP_GETPOTNAME
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_GETPOTNAME_H
#define EPTWRAP_GETPOTNAME_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_getpotname(int ain,int aout,int pid,char** name,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
