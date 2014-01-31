/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_SINGLE
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_EPUPDATE_SINGLE_H
#define EPTWRAP_EPUPDATE_SINGLE_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_epupdate_single1(int ain,int aout,int pid,W_DARRAY(pars),
				double cmu,double crho,int* rstat,
				double* alpha,double* nu,double* logz,
				W_ERRORARGS);

  void eptwrap_epupdate_single2(int ain,int aout,char* pname,W_DARRAY(pars),
				double cmu,double crho,int* rstat,
				double* alpha,double* nu,double* logz,
				W_ERRORARGS);

  void eptwrap_epupdate_single3(int ain,int aout,W_IARRAY(potids),
				W_IARRAY(numpot),W_DARRAY(parvec),
				W_IARRAY(parshrd),int pind,double cmu,
				double crho,int* rstat,double* alpha,
				double* nu,double* logz,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
