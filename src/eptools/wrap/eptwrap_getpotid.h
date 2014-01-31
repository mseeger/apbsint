/* -------------------------------------------------------------------
 * EPTWRAP_GETPOTID
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_GETPOTID_H
#define EPTWRAP_GETPOTID_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_getpotid(int ain,int aout,char* name,int* pid,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
