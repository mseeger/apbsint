/* -------------------------------------------------------------------
 * EPTWRAP_GETPOTAGROUP
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_GETPOTAGROUP_H
#define EPTWRAP_GETPOTAGROUP_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_getpotagroup(int ain,int aout,int pid,int* agid,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
