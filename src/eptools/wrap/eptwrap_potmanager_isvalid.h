/* -------------------------------------------------------------------
 * EPTWRAP_POTMANAGER_ISVALID
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_POTMANAGER_ISVALID_H
#define EPTWRAP_POTMANAGER_ISVALID_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_potmanager_isvalid(int ain,int aout,W_IARRAY(potids),
				  W_IARRAY(numpot),W_DARRAY(parvec),
				  W_IARRAY(parshrd),W_ARRAY(annobj,void*),
				  int posoff,W_IARRAY(tauind),char** retstr,
				  W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
