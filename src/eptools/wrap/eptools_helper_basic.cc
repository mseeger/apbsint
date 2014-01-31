/* -------------------------------------------------------------------
 * Basic helper functions for EPTOOLS wrapper functions
 * These can be compiled independent of LHOTSE
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/eptools/wrap/eptools_helper_basic.h"

void fillVec(double* vec,int n,double val)
{
  int i;
  for (i=0; i<n; i++) vec[i]=val;
}

void fillVecStep(double* vec,int n,int incv,double val)
{
  int i;
  for (i=0; i<n; i++,vec+=incv) *vec=val;
}
