/* -------------------------------------------------------------------
 * Basic helper functions for EPTOOLS wrapper functions
 * These can be compiled independent of LHOTSE
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_HELPER_BASIC_H
#define EPTOOLS_HELPER_BASIC_H

// Helper functions

void fillVec(double* vec,int n,double val); 

void fillVecStep(double* vec,int n,int incv,double val);

#endif
