/* -------------------------------------------------------------------
 * EPTWRAP_EPUPDATE_SINGLE_BVPREC
 * -------------------------------------------------------------------
 * Declaration wrapper function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTWRAP_EPUPDATE_SINGLE_BVPREC_H
#define EPTWRAP_EPUPDATE_SINGLE_BVPREC_H

#include "src/eptools/wrap/eptools_helper_macros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void eptwrap_epupdate_single_bvprec1(int ain,int aout,int pid,W_DARRAY(pars),
				       void* annobj,double cmu,double crho,
				       double ca,double cc,int* rstat,
				       double* alpha,double* nu,double* hata,
				       double* hatc,double* logz,W_ERRORARGS);

  void eptwrap_epupdate_single_bvprec2(int ain,int aout,char* pname,
				       W_DARRAY(pars),void* annobj,double cmu,
				       double crho,double ca,double cc,
				       int* rstat,double* alpha,double* nu,
				       double* hata,double* hatc,double* logz,
				       W_ERRORARGS);

  void eptwrap_epupdate_single_bvprec3(int ain,int aout,W_IARRAY(potids),
				       W_IARRAY(numpot),W_DARRAY(parvec),
				       W_IARRAY(parshrd),W_ARRAY(annobj,void*),
				       int pind,double cmu,double crho,
				       double ca,double cc,int* rstat,
				       double* alpha,double* nu,double* hata,
				       double* hatc,double* logz,W_ERRORARGS);

#ifdef __cplusplus
}
#endif

#endif
